package org.seqcode.projects.sem.mixturemodel;

import java.io.*;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import org.seqcode.deepseq.StrandedPair;
import org.seqcode.gseutils.Pair;
import org.seqcode.projects.sem.events.*;
import org.seqcode.projects.sem.framework.*;
import org.seqcode.projects.sem.utilities.Timer;
import org.seqcode.projects.sem.utilities.ConsoleProgressBar;

import org.seqcode.projects.sem.utilities.NucleosomePoissonBackgroundModel;
import org.seqcode.projects.sem.utilities.EMmode;
import org.seqcode.deepseq.experiments.*;
import org.seqcode.genome.*;
import org.seqcode.genome.location.*;
import org.seqcode.gsebricks.verbs.location.ChromosomeGenerator;

/**
 * BindingMixture: defines a mixture model over binding events.
 * @author Jianyu Yang
 *
 */

public class BindingMixture {
	
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected EventsConfig evconfig;
	protected SEMConfig config;
	protected ExperimentManager manager;
	protected BindingManager bindingManager;
	protected PotentialRegionFilter potRegFilter;
	protected List<Region> testRegions;
	protected HashMap<Region, List<List<BindingComponent>>> activeComponents; //Components active after a round of execute()
	protected HashMap<ExperimentCondition, NucleosomePoissonBackgroundModel> conditionBackgrounds = new HashMap<ExperimentCondition, NucleosomePoissonBackgroundModel>(); //Genomic Background models for each condition -- used to set alpha values in sparse prior
	protected List<BindingEvent> bindingEvents;
	protected Pair<String, Integer> plotDyad; //plot region which contains this dyad
	protected int trainingRound = 0;
	protected double LAP = 0;
	protected double lastLAP = -Double.MAX_VALUE;
	protected boolean converged = false;
	protected double noisePerBase[];		//Defines global noise
	protected double relativeCtrlNoise[];	//Defines global noise
	protected HashMap<Region, Double[]> noiseResp = new HashMap<Region, Double[]>(); //noise responsibilities after a round of execute(). Hashed by Region, indexed by condition
	protected HashMap<ExperimentCondition, Map<Integer, Double>> noiseFragSizeFreq;
	protected int numRegions = 0;										//Number of all regions
	protected AtomicInteger passRegion = new AtomicInteger();			//Number of regions have finished
	protected ConsoleProgressBar progressBar;
	
	public BindingMixture(GenomeConfig gcon, ExptConfig econ, EventsConfig evcon, SEMConfig semconfig, ExperimentManager eMan, BindingManager bMan, PotentialRegionFilter filter) {
		gconfig = gcon;
		econfig = econ;
		evconfig = evcon;
		config = semconfig;
		manager = eMan;
		bindingManager = bMan;
		potRegFilter = filter;
		testRegions = filter.getPotentialRegions();
		numRegions = testRegions.size();
//		plotDyad = bindingManager.getBindingModel(manager.getIndexedCondition(0)).get(0).getIntialDyadByIndex(0);
		plotDyad = semconfig.getPlotDyad();
		bindingEvents = new ArrayList<BindingEvent>();
		BindingEvent.setExperimentManager(manager);
		BindingEvent.setConfig(evconfig);
		
		activeComponents = new HashMap<Region, List<List<BindingComponent>>>();
		for(ExperimentCondition cond: manager.getConditions()) {
			//System.out.println(config.getGenome().getGenomeLength()-potRegFilter.getPotRegionLengthTotal());
			conditionBackgrounds.put(cond, new NucleosomePoissonBackgroundModel(-1, config.getSigLogConf(), cond.getTotalSignalPairCount()*(1-cond.getTotalSignalPairVsNoisePairFrac()), config.getGenome().getGenomeLength()-potRegFilter.getPotRegionLengthTotal(), econfig.getMappableGenomeProp(), bindingManager.getMaxInfluenceRange(cond), '.', 1, true));
			// ignore fixed alpha when determining threshold for each nucleosome
			double alf = (double)conditionBackgrounds.get(cond).calcCountThreshold(bindingManager.getMaxInfluenceRange(cond));
			if(config.getFixedAlpha()<0)
				System.err.println("DynamicAlpha "+cond.getName()+"\tRange="+bindingManager.getMaxInfluenceRange(cond)+"\t"+alf);
			else
				System.err.println("FixedAlpha "+cond.getName()+"\t"+config.getFixedAlpha());
		}
		
		noisePerBase = new double[manager.getNumConditions()];
		relativeCtrlNoise = new double[manager.getNumConditions()];
		noiseFragSizeFreq = filter.getNonPotRegFragSizeFreqSigChannel();
		
		initializeGlobalNoise();
	}
	
	// Is the model converged?
	public boolean ifConverged() {
		//monitor
		System.err.println("LAP this round: " + LAP);
		converged = Math.abs((LAP - lastLAP)/lastLAP)<config.EM_CONVERGENCE;
		lastLAP = LAP;
		LAP = 0;
		return converged;
	}
	
	/**
     * Initialize the global noise parameters. Inferred either from:
     *  - non-potential region read counts, or
     *  - noise proportion of reads from SES
     */
    protected void initializeGlobalNoise(){
    	for(int e=0; e<manager.getNumConditions(); e++){
    		ExperimentCondition cond = manager.getIndexedCondition(e);
    		
    		// Part that deals with read counts from non-potential regions... calculate values anyway whether using them or not
    		double potRegLengthTotal = potRegFilter.getPotRegionLengthTotal();
    		double nonPotRegLengthTotal = config.getGenome().getGenomeLength() - potRegLengthTotal;
    		//Combine control channel counts (avoiding duplication)
    		double potRegCountsSigChannel=potRegFilter.getPotRegCountsSigChannel(cond);
    		double nonPotRegCountsSigChannel=potRegFilter.getNonPotRegCountsSigChannel(cond); 
    		double potRegCountsCtrlChannel=potRegFilter.getPotRegCountsCtrlChannel(cond);
    		double nonPotRegCountsCtrlChannel=potRegFilter.getNonPotRegCountsCtrlChannel(cond); 
    		List<Sample> ctrls = cond.getControlSamples();
    		
    		//relativeCtrlNoise just tells us if there is a systemic over/under representation of reads in potential regions (in the control)
    		//NOTE: not used for anything right now. 
    		relativeCtrlNoise[e] = (potRegCountsCtrlChannel==0 && nonPotRegCountsCtrlChannel==0) ? 
    				1 : (potRegCountsCtrlChannel/potRegLengthTotal)/(nonPotRegCountsCtrlChannel/nonPotRegLengthTotal);

    		
    		noisePerBase[e] = nonPotRegCountsSigChannel/nonPotRegLengthTotal;  //Signal channel noise per base
    		
    		if(config.isVerbose()) {
    			System.out.println("Condition: "+cond.getName()+"\n"
    					+ "potRegCountsSignalChannel: "+potRegCountsSigChannel+"\n"
    					+ "nonPotRegCountsSignalChannel: "+nonPotRegCountsSigChannel+"\n"
    					+ "potRegCountsCtrlChannel: "+potRegCountsCtrlChannel + "\n"
    					+ "nonPotRegCountsCtrlChannel: " +nonPotRegCountsCtrlChannel+"\n"
    					+ "potRegLengthTotal: " + potRegLengthTotal + "\n"
    					+ "nonPotRegLengthTotal: " + nonPotRegLengthTotal + "\n");
    			System.err.println("Global noise per base initialization for "+cond.getName()+" = "+String.format("%.4f", noisePerBase[e]));
    		}
    	}
    }
    
    /**
     * Update the global noise parameters, using both non-potential region counts and assigned noise responsibilities
     */
    public void updateGlobalNoise(){
    	converged = true;
    	for(int e=0; e<manager.getNumConditions(); e++){
    		ExperimentCondition cond = manager.getIndexedCondition(e);
    		
    		//Don't need to examine noise reads in the update
    		double noiseReads=potRegFilter.getNonPotRegCountsSigChannel(cond); 
    		
    		//Don't forget to update noise fragment size frequency here!!!!!!!!!!!!!!!!
    		
    		for(Region r : noiseResp.keySet())
    			noiseReads+=noiseResp.get(r)[e];
    		
    		double newNoisePerBase = noiseReads/config.getGenome().getGenomeLength();
    		noisePerBase[e] = newNoisePerBase;  //Signal channel noise per base
    		
    		//monitor print updated noisePerBase
    		System.out.println("training round: "+trainingRound+"\tnoise per base: "+noisePerBase[e]);
    	}
    }
     
    /**
     * Update condition backgrounds for alpha
     */
    public void updateAlpha() {
    	for(ExperimentCondition cond: manager.getConditions()) {
    		int c = cond.getIndex();
    		conditionBackgrounds.put(cond, new NucleosomePoissonBackgroundModel(-1, config.getSigLogConf(), noisePerBase[c]*config.getGenome().getGenomeLength(), 
    				config.getGenome().getGenomeLength(), econfig.getMappableGenomeProp(), bindingManager.getMaxInfluenceRange(cond), '.', 1, true));
			double alf = (double)conditionBackgrounds.get(cond).calcCountThreshold(bindingManager.getMaxInfluenceRange(cond));
			if(config.getFixedAlpha()<0) 
				System.err.println("DynamicAlpha "+cond.getName()+"\tRange="+bindingManager.getMaxInfluenceRange(cond)+"\t"+alf);
			else
				System.err.println("FixedAlpha "+cond.getName()+"\t"+config.getFixedAlpha());
    	}
    }
    
    
	public void execute(boolean EM, boolean uniformBindingComponents, EMmode mode) {
		trainingRound++;
		progressBar = new ConsoleProgressBar(0, 100, 30, '#');
		passRegion.set(0);
		
		//Have to split the test regions up by chromosome in order to maintain compatibility with experiment file cache loading
		//There will be some performance hit here, as all threads have to finish in a given chromosome before moving on to the next one. 
		Iterator<Region> chroms = new ChromosomeGenerator().execute(gconfig.getGenome());
		while (chroms.hasNext()) {
			Region currChr = chroms.next();
			List<Region> currChrTestReg = new ArrayList<Region>();
			for(Region r : testRegions)
				if(currChr.overlaps(r))
					currChrTestReg.add(r);
			
			if(currChrTestReg.size()>0){
				int numThreads = config.getMaxThreads()>currChrTestReg.size() ?  currChrTestReg.size() : config.getMaxThreads(); 
				Thread[] threads = new Thread[numThreads];
		        ArrayList<Region> threadRegions[] = new ArrayList[numThreads];
		        int i = 0;
		        for (i = 0 ; i < threads.length; i++) {
		            threadRegions[i] = new ArrayList<Region>();
		        }i=0;
		        for(Region r : currChrTestReg){
		            threadRegions[(i++) % numThreads].add(r);
		        }
		
		        for (i = 0 ; i < threads.length; i++) {
		            Thread t = new Thread(new BindingMixtureThread(threadRegions[i], EM, uniformBindingComponents, mode));
		            t.start();
		            threads[i] = t;
		        }
		        boolean anyrunning = true;
		        while (anyrunning) {
		            anyrunning = false;
		            try {
		                Thread.sleep(5000);
		            } catch (InterruptedException e) { }
		            for (i = 0; i < threads.length; i++) {
		                if (threads[i].isAlive()) {
		                    anyrunning = true;
		                    break;
		                }
		            }
		        }
			}
		}
	}
	
    /**
     * Print all components active at the current time to a file.
     */
    public void printActiveComponentsToFile(EMmode mode){
    	try {
    		int totalActiveNuc = 0;
    		String filename = "";
    		if(mode.equals(EMmode.NORMAL))
    			filename = config.getOutputIntermediateDir()+File.separator+config.getOutBase()+"_t"+trainingRound+".components";
    		else if(mode.equals(EMmode.ALTERNATIVE))
    			filename = config.getOutputIntermediateDir()+File.separator+config.getOutBase()+"_alternative_t"+trainingRound+".components";
    		else if(mode.equals(EMmode.CONSENSUS))
    			filename = config.getOutputIntermediateDir()+File.separator+config.getOutBase()+"_consensus_t"+trainingRound+".components";
			FileWriter fout = new FileWriter(filename);
			fout.write("#region\tchromosome\tdyad\tpi\tsumResp\tfuzziness\ttau\tisPair\n");
			for(Region rr : activeComponents.keySet()){
	    		List<List<BindingComponent>> comps = activeComponents.get(rr);
	    		for(ExperimentCondition cond : manager.getConditions()){
	    			for(BindingComponent comp : comps.get(cond.getIndex())){
	    				totalActiveNuc += 1;
	    				fout.write(rr.getLocationString()+"\t"+cond.getName()+"\t"+comp.toString());
	    				if(config.getFixedAlpha()<0)
	    					fout.write("\t"+comp.getPValue());
	    				fout.write("\n");
	    			}
	    		}
	    	}
			fout.close();
			System.err.println("Total active nucleosome number: "+totalActiveNuc);
		} catch (IOException e) {
			e.printStackTrace();
		}
    }
    
    /**
     * Print nucleosome information per condition to a file
     * this will be called when EM ends
     */
    public void printNucleosomeInfoToFile() {
    	try {
    		for(ExperimentCondition cond: manager.getConditions()) {
    			String filename = config.getOutputParentDir() + File.separator + config.getOutBase() + "_" + cond.getName() + "_nucleosome_info.tsv";
    			BufferedWriter fout = new BufferedWriter(new FileWriter(filename));
    			//header
    			fout.write("chromosome\tdyad\toccupancy\tfuzziness\tsubtype\tProb_subtypes");
    			if(config.getFixedAlpha()<0) fout.write("\tP-value");
    			fout.write("\n");
    			for(Region rr: activeComponents.keySet()) {
    				List<List<BindingComponent>> comps = activeComponents.get(rr);
    				for(BindingComponent comp: comps.get(cond.getIndex())) {
    					fout.write(comp.standardInfo());
	    				if(config.getFixedAlpha()<0)
	    					fout.write("\t"+comp.getPValue());
	    				fout.write("\n");
    				}
    			}
    			fout.close();
    		}
    	} catch (IOException e) {
    		e.printStackTrace();
    	}
    }
    
    /**
     * Print nucleosome comparison information to a file (now only for two conditions)
     */
    public void printNucleosomeComparisonToFile() {
    	try {
    		for(ExperimentCondition cond: manager.getConditions()) {
	    		String filename = config.getOutputIntermediateDir() + File.separator + config.getOutBase() + "_" + cond.getName() + "_comparison_info.tsv";
	    		BufferedWriter fout = new BufferedWriter(new FileWriter(filename));
	    		//writer header information
	    		fout.write("#nucleosome comparison result info of "+cond.getName()+"\n");
	    		for(Region rr: activeComponents.keySet()) {
	    			List<List<BindingComponent>> comps = activeComponents.get(rr);
	    			for(BindingComponent comp: comps.get(cond.getIndex())) {
	    				fout.write(rr.getLocationString()+"\t"+manager.getIndexedCondition(cond.getIndex()).getName()+"\t"+comp.toString()+"\t");
	    				for(Pair<Integer, Integer> index: comp.getCompareRestulsConvert().keySet()) {
	    					fout.write(manager.getIndexedCondition(index.car()).getName()+"\t"+comps.get(index.car()).get(index.cdr()).toString());
	    					fout.write("\t"+Arrays.toString(comp.getCompareRestulsConvert().get(index).cdr()));
	    					fout.write("\t"+Arrays.toString(comp.getCompareRestulsConvert().get(index).car()));
	    				}
	    				fout.write("\n");
	    			}
	    		}
	    		fout.close();
    		}
    	} catch (IOException e) {
    		e.printStackTrace();
    	}
    }
	
	/**
	 * BindingMixtureThread: run binding mixtures over a subset of regions
	 * uniformComponents: if true, components are initialized at uniform prob and spacing.
	 * 					otherwise, components are initialized from activeComponents
	 */
	class BindingMixtureThread implements Runnable{
		private Collection<Region> regions;
		private int numBindingComponents = 1; 
		private boolean runEM = true;
		private boolean uniformBindingComponents = false;
		private EMmode mode = EMmode.NORMAL;	// function index: 0:EM, 1:alternative(30bp exclusion range), 2:consensus(147bp exclusion)
		
		public BindingMixtureThread(Collection<Region> regs, boolean EM, boolean uniformBindingComponents, EMmode mode) {
			regions = regs;
			this.uniformBindingComponents = uniformBindingComponents;
			runEM = EM;
			this.mode = mode;
		}
		
		//Run the binding mixture over each test region
		public void run() {
			//For each region in the test set
			for(Region rr: regions) {
				try {
					//Initialize array of binding component lists, indexed by condition
					List<List<BindingComponent>> currComps = new ArrayList<List<BindingComponent>>();
					for(int e=0; e<manager.getNumConditions(); e++)
						currComps.add(new ArrayList<BindingComponent>());
					
					ArrayList<Region> windows = new ArrayList<Region>();
					windows.add(rr);
					
					if(runEM) {
						//Run EM
						Double[] noiseRSums = new Double[manager.getNumConditions()];
						for(int e=0; e<manager.getNumConditions(); e++) { noiseRSums[e] = 0.0;}
						for(Region w: windows) {
							Pair<List<NoiseComponent>, List<List<BindingComponent>>> wComps = analyzeWindowEM(w);
							for(int e=0; e<manager.getNumConditions(); e++) {
								noiseRSums[e] += wComps.car().get(e).getSumResponsibility();
								currComps.get(e).addAll(wComps.cdr().get(e));
							}
							//update the finshied region number and progress bar
							synchronized(passRegion) {
								int numPass = passRegion.incrementAndGet();
								if(numPass % 1e2 == 0) {
									progressBar.show(100l*(long)numPass/(long)numRegions);
								} else if (numPass == numRegions) {
									progressBar.show(100l);
								}
							}
							//Free memory
							wComps.car().clear();
							wComps.cdr().clear();
							wComps = null;
						}
						
						//Only non-zero components are returned by analyzeWindow, so add them to the recorded active components
						synchronized(activeComponents) { activeComponents.put(rr, currComps);}
						
						//Add the sum of noise responsibilities to this region
						synchronized(noiseResp) { noiseResp.put(rr, noiseRSums);}
						
					} 
				} catch(Exception e) {
					System.err.println("ERROR: Exception when analyzing region"+rr.toString());
					e.printStackTrace(System.err);
					System.exit(-1);
				}
			}
		}
	
		private Pair<List<NoiseComponent>, List<List<BindingComponent>>> analyzeWindowEM(Region w) throws Exception {
//			System.err.println("Region: "+w.getChrom()+":"+w.getStart()+"-"+w.getEnd());
			Timer timer = new Timer();
			
			timer.extra_start();
			// Determine which BindingEM method will be used
			BindingEM_interface EM;
			EM = new BindingEM(config, manager, bindingManager, conditionBackgrounds, potRegFilter.getPotentialRegions().size());
				
			List<List<BindingComponent>> bindingComponents = null;
			List<NoiseComponent> noiseComponents = null;
			List<List<BindingComponent>> nonZeroComponents = new ArrayList<List<BindingComponent>>();
			for(int e=0; e<manager.getNumConditions(); e++)
				nonZeroComponents.add(new ArrayList<BindingComponent>());
				
			//check if plot this region
			boolean plotEM = false;
			if(w.getChrom().equals(plotDyad.car()) && 
					w.getStart() <= plotDyad.cdr() && 
					w.getEnd() >= plotDyad.cdr()) {
				plotEM = true;
			}
			
			//monitor: count time
			timer.start();
			
			//Load signal data
			List<List<StrandedPair>> signals = loadSignalData(w);
			if (signals==null)
				return new Pair<List<NoiseComponent>, List<List<BindingComponent>>>(noiseComponents, nonZeroComponents);
			
			//Load control data
			List<List<StrandedPair>> controls = loadControlData(w);
			
			//Initialize noise components
			noiseComponents = initializeNoiseComponents(w, signals, controls);
			
			//Initialize binding components
        	if(uniformBindingComponents)
        		bindingComponents = initializeBindingComponentsUniformly(w, noiseComponents);
        	else
        		bindingComponents = initializeBindingComponentsFromAllConditionActive(w, noiseComponents, mode);
            
            //ATAC-seq prior
            double[][] atacPrior = null;
            
			//monitor: count time
            timer.end("load");

            //EM learning: resulting binding components list will only contain non-zero components
            nonZeroComponents = EM.train(signals, w, noiseComponents, bindingComponents, numBindingComponents, atacPrior, trainingRound, mode, timer, plotEM);
            
            //Add the log likelihood to the whole model
            LAP += EM.getLAP();
            
            //Free memory
            signals.clear();
            signals = null;
            bindingComponents.clear();
            bindingComponents = null;
            EM = null;
            
            timer.extra_end();
            return new Pair<List<NoiseComponent>, List<List<BindingComponent>>>(noiseComponents, nonZeroComponents);
		}
		
		
		/**
		 * Load all signal read hits in a region by condition. 
		 * Note(jy): Change StrandedBaseCount to StrandedPair
		 * @param w
		 * @return List of List of StrandedPair, indexed by replicate index
		 * @author Jianyu Yang
		 */
		private List<List<StrandedPair>> loadSignalData(Region w){
			List<List<StrandedPair>> data = new ArrayList<List<StrandedPair>>();
			for(ExperimentCondition cond : manager.getConditions()){
				for(ControlledExperiment rep : cond.getReplicates())
					data.add(new ArrayList<StrandedPair>());
			}
			for(ExperimentCondition cond : manager.getConditions()){
				for(ControlledExperiment rep : cond.getReplicates()){
					data.get(rep.getIndex()).addAll(rep.getSignal().getPairsByMid(w));
				}
			}
			return data;
		}
		
		/**
		 * Load all control read hits in a region by condition. 
		 * Note(jy): Change StrandedBaseCount to StrandedPair
		 * @param w
		 * @return List of List of StrandedPair, indexed by replicate index
		 * @author Jianyu Yang
		 */
		private List<List<StrandedPair>> loadControlData(Region w){
			List<List<StrandedPair>> data = new ArrayList<List<StrandedPair>>();
			for(ExperimentCondition cond : manager.getConditions()){
				for(ControlledExperiment rep : cond.getReplicates())
					data.add(new ArrayList<StrandedPair>());
			}
			for(ExperimentCondition cond : manager.getConditions()){
				for(ControlledExperiment rep : cond.getReplicates()){
					if(rep.hasControl())
						data.get(rep.getIndex()).addAll(rep.getControl().getPairsByMid(w));
				}
			}
			return data;
		}	
	
	
		//Initialize the global noise parameters, inferred either from:
		// - non-potential region read counts, or
		// - noise proportion of reads from SES
		protected void initializeGlobalNoise() {
			for(int e=0; e<manager.getNumConditions(); e++) {
				ExperimentCondition cond = manager.getIndexedCondition(e);
			
				//Part that deals with read counts from non-potential regions... calculate values anyway whether using or not
				double potRegLengthTotal = potRegFilter.getPotRegionLengthTotal();
				double nonPotRegLengthTotal = config.getGenome().getGenomeLength() - potRegLengthTotal;
				//Combine control channel counts (avoiding duplication)
				double potRegCountsSigChannel = potRegFilter.getPotRegCountsSigChannel(cond);
				double nonPotRegCountsSigChannel=potRegFilter.getNonPotRegCountsSigChannel(cond); 
				double potRegCountsCtrlChannel=potRegFilter.getPotRegCountsCtrlChannel(cond);
				double nonPotRegCountsCtrlChannel=potRegFilter.getNonPotRegCountsCtrlChannel(cond); 
				List<Sample> ctrls = cond.getControlSamples();
			
				//relativeCtrlNoise just tells us if there is a systemic over/under representation of reads in potential regions (in the control)
				//NOTE: not used for anything right now. 
				relativeCtrlNoise[e] = (potRegCountsCtrlChannel==0 && nonPotRegCountsCtrlChannel==0) ? 
						1 : (potRegCountsCtrlChannel/potRegLengthTotal)/(nonPotRegCountsCtrlChannel/nonPotRegLengthTotal);
			
				noisePerBase[e] = nonPotRegCountsSigChannel/nonPotRegLengthTotal;  //Signal channel noise per base
				System.err.println("Global noise per base initialization for "+cond.getName()+" = "+String.format("%.4f", noisePerBase[e]));
			
			}
		}
		
		/**
		 * Set the conditions in which a given binding event is still active after EM training
		 * @param b
		 * @param currReg
		 */
		private void setFoundInConditions(BindingEvent b, Region currReg){
			for(int e=0; e<manager.getNumConditions(); e++)
        		for(BindingComponent comp : activeComponents.get(currReg).get(e)){
        			if(comp.getPosition() == b.getPoint().getLocation()){
        				b.setIsFoundInCondition(e, true);
        			}
        		}
		}
	
		/**	
		 * Initializes the noise components.
		 * Noise distributions are set per replicate, so control channel reads are not combined into a single list. 
		 *
		 * @param currReg
		 */
		private List<NoiseComponent> initializeNoiseComponents(Region currReg, List<List<StrandedPair>> sigHits, List<List<StrandedPair>> ctrlHits){
			List<NoiseComponent> noise = new ArrayList<NoiseComponent>();
			int numReps = manager.getReplicates().size();		
			double [] localSigRepCounts=new double [numReps];
			double [] localCtrlRepCounts=new double [numReps];
    	
			//Calculate expected noise distributions
			double [][] distribs=new double[numReps][];
			for(ExperimentCondition cond : manager.getConditions())
				for(ControlledExperiment rep : cond.getReplicates()){
					if(rep.hasControl() && ctrlHits.get(rep.getIndex()).size()>0){
						distribs[rep.getIndex()] = smoothNoiseDistribs(currReg, ctrlHits.get(rep.getIndex()));
						localCtrlRepCounts[rep.getIndex()]=0;
						for(StrandedPair b : ctrlHits.get(rep.getIndex()))
							localCtrlRepCounts[rep.getIndex()]+= b.getWeight();
					}else
						distribs[rep.getIndex()] = null;
				}
			
			//Calculate expected fragment size frequency distribution
			Map<ExperimentCondition, HashMap<Integer, Double>> fragSizeFreq = new HashMap<ExperimentCondition, HashMap<Integer, Double>>();
			for(ExperimentCondition cond: manager.getConditions()) {
				fragSizeFreq.put(cond, normalizeFragSizeFreq(noiseFragSizeFreq.get(cond)));
			}
			    	
			//Initialize the noise component
			for(int e=0; e<manager.getNumConditions(); e++){
				ExperimentCondition cond = manager.getIndexedCondition(e);
				double emission = 0;
    		
				//Sum signal reads & set local region experiment counts
				double sigCounts=0;
				for(ControlledExperiment rep : cond.getReplicates()){
					localSigRepCounts[rep.getIndex()]=0;
					for(StrandedPair b : sigHits.get(rep.getIndex())){
						sigCounts+=b.getWeight(); localSigRepCounts[rep.getIndex()]+=b.getWeight();
					}
				}
    		
				//Calculate a local noise factor to check for expected over-representation of noise reads, as specified in the control.
				//We assume that only noise generates the control channel. Therefore, the first two terms in the localNoiseFactor calculation
				//test for over-representation in the observed control read counts in the local window. 
				//The last term in the calculation is a local replicate weight, to account for the fact that some replicate signal channels are 
				//producing more of the signal reads (and of course, the noise we are actually accounting for is in the signal channel). 
				double localNoiseFactor = 0;
				for(ControlledExperiment rep : cond.getReplicates()){
					if(rep.hasControl() && ctrlHits.get(rep.getIndex()).size()>0){
						localNoiseFactor+=(localCtrlRepCounts[rep.getIndex()]/(double)currReg.getWidth()) /
										(rep.getControl().getHitCount()/(double)currReg.getWidth())     *
										(localSigRepCounts[rep.getIndex()]/sigCounts); //over-rep x weight
					}else{
						localNoiseFactor+=(localSigRepCounts[rep.getIndex()]/sigCounts); //1 x weight
					}
				}
    		
				//Calculate expected noise emission = number of expected noise reads over the total signal reads in this region
				// If local noise factor is above 1, account for it. Otherwise, meh.
				emission = (noisePerBase[e] * (double)currReg.getWidth()) / sigCounts;
				if(localNoiseFactor>1)
					emission*=localNoiseFactor;
				if(emission>config.NOISE_EMISSION_MAX)
					emission = config.NOISE_EMISSION_MAX;
				if(emission<config.NOISE_EMISSION_MIN)
					emission = config.NOISE_EMISSION_MIN;
				
//				System.err.println("localNoiseFactor: " + localNoiseFactor);
//				System.err.println("emission: " + emission);
				
				//Add the noise component
				NoiseComponent n = new NoiseComponent(emission, distribs, currReg, numReps, fragSizeFreq.get(cond));
				noise.add(n);
			}
			return noise;
		}//end of initializeNoiseComponents method
    
		/**
		 * Initializes the components uniformly: i.e. space them evenly along the region.
		 *
		 * @param currReg
		 */
		private List<List<BindingComponent>> initializeBindingComponentsUniformly(Region currReg, List<NoiseComponent> noise){
			List<List<BindingComponent>> components = new ArrayList<List<BindingComponent>>();
			for(int e=0; e<manager.getNumConditions(); e++)
				components.add(new ArrayList<BindingComponent>());
		
			//Place components along region
			int componentSpacing =  config.INIT_COMPONENT_SPACING;
			if(componentSpacing >= currReg.getWidth())
				System.err.println("Error:  region width less than component spacing in "+currReg.getLocationString());
			
			//TODO: I think here should be +1 because if there is an overhang, this region can hold +1 components.
			if(currReg.getWidth()%componentSpacing!=0) 
				numBindingComponents = currReg.getWidth()/componentSpacing + 1;
			else
				numBindingComponents = currReg.getWidth()/componentSpacing;

			//Set up the components
			for(int e=0; e<manager.getNumConditions(); e++){
				double numC=0;
				for(int i=0; i<numBindingComponents; i++){
					Point pos = new Point(config.getGenome(), currReg.getChrom(), currReg.getStart()+(i*componentSpacing));
					BindingComponent currComp = new BindingComponent(pos, manager.getReplicates().size());
					currComp.setIndex(i);
					numC++;
    	           components.get(e).add(currComp);
				}
				//Initialize normalized mixing probabilities (subtracting the noise emission probability)
				//Add initialization step for fuzziness and tau here @Jianyu Yang
				double emission = (1-noise.get(e).getPi())/numC;
				int numBindingSubtype = bindingManager.getBindingSubtypes(manager.getIndexedCondition(e)).size();
				double[] tauInit = new double[numBindingSubtype];
				for(int i=0; i<numBindingSubtype; i++) {
					tauInit[i] = (double)1/(double)numBindingSubtype;
				}
				for(BindingComponent b : components.get(e)){
					b.uniformInit(emission);
					b.setFuzziness(bindingManager.getBindingModel(manager.getIndexedCondition(e)).get(0).getIntialFuzziness());
					b.setTau(tauInit);
				}
			}
			return components; 
		}//end of initializeComponents method
		
		/**
         * Initializes the components from all active components in all conditions: 
         * 		Uses all active component locations from the last round of training in each condition,
         * 		i.e. not just the conditions in which those active components were active. Also adds in 
         * 		extra components flanking the active locations in case the binding distribution update
         * 		has made more joint events separable. If no components exist, a rescue component is added. 
         *
         * @param currReg
         */
        private List<List<BindingComponent>> initializeBindingComponentsFromAllConditionActive(Region currReg, List<NoiseComponent> noise, EMmode mode){
        	//Initialize component positions with active locations
        	List<Integer> componentPositions = new ArrayList<Integer>();
        	for(int e=0; e<manager.getNumConditions(); e++)
        		for(BindingComponent comp : activeComponents.get(currReg).get(e)){
        			if(!componentPositions.contains(comp.getPosition()) && comp.getPosition()>=currReg.getStart() && comp.getPosition()<currReg.getEnd())
        				componentPositions.add(comp.getPosition());
        		}
        	
        	//sort all component positions
        	Collections.sort(componentPositions);
        	
        	//If no components exist in region, add one to the center to allow rescues
        	numBindingComponents = componentPositions.size();
        	if(numBindingComponents==0){
        		componentPositions.add(currReg.getMidpoint().getLocation());
        		numBindingComponents++;
        	}
        	
        	//if the overhang region is longer than 1/2 maxIR, add one rescue nucleosome on the start or end of region
        	int maxIR = (int)Math.round(2 * 2.58 * Math.sqrt(config.INIT_FUZZINESS));
			int half_maxIR = maxIR / 2;
			if((currReg.getEnd()-componentPositions.get(componentPositions.size()-1))>half_maxIR) {
				componentPositions.add(currReg.getEnd());
				numBindingComponents++;
			}
        	if(componentPositions.get(0)-currReg.getStart()>half_maxIR) {
        		componentPositions.add(currReg.getStart());
        		numBindingComponents++;
        	}
        	
        	//if distance between two adjacent nucleosomes is > maxIR (defined by initial fuzziness), add one rescue nucleosome in the middle
        	//to rescue the fragments between them
        	int size = componentPositions.size();
        	for(int i=0; i<(size-1); i++) {
        		int currentPosition = componentPositions.get(i);
        		int nextPosition = componentPositions.get(i+1);
        		if((nextPosition-currentPosition)>maxIR) {
        			componentPositions.add((currentPosition + nextPosition)/2);
        			numBindingComponents++;
        		}
        	}

        	//Make new components with these locations
        	List<List<BindingComponent>> components = new ArrayList<List<BindingComponent>>();
        	for(int e=0; e<manager.getNumConditions(); e++)
        		components.add(new ArrayList<BindingComponent>());
        	
        	//Set up the components
        	for(int e=0; e<manager.getNumConditions(); e++){
        		double numC=(double)numBindingComponents; int index=0;
        		double emission = (1-noise.get(e).getPi())/numC;
	    		for(Integer i : componentPositions){
	    			Point pos = new Point(config.getGenome(), currReg.getChrom(), i);
	    			BindingComponent currComp = new BindingComponent(pos, manager.getReplicates().size());
	    			currComp.setIndex(index);
	    			index++;
	    			
    				//Initialize normalized mixing probabilities (subtracting the noise emission probability)
        			currComp.uniformInit(emission);
    				components.get(e).add(currComp);
    			}
	    		int numBindingSubtype = bindingManager.getBindingSubtypes(manager.getIndexedCondition(e)).size();
	    		double[] tauInit = new double[numBindingSubtype];
				for(int i=0; i<numBindingSubtype; i++) {
					tauInit[i] = (double)1/(double)numBindingSubtype;
				}
				for(BindingComponent b : components.get(e)){
					b.uniformInit(emission);
					b.setFuzziness(bindingManager.getBindingModel(manager.getIndexedCondition(e)).get(0).getIntialFuzziness());
					b.setTau(tauInit);
				}
    		}
        	return components; 
        }//end of initializeComponents method
		
		/**
         * Initializes the components from all active components in all conditions: 
         * 		Uses all active component locations from the last round of training in each condition,
         * 		i.e. not just the conditions in which those active components were active. Also adds in 
         * 		extra components flanking the active locations in case the binding distribution update
         * 		has made more joint events separable. If no components exist, a rescue component is added. 
         *
         * @param currReg
         */
        private List<List<BindingComponent>> initializeBindingComponentsWithExclusion(Region currReg, List<NoiseComponent> noise, EMmode mode){
        	//Sort active components by responsibilities
			Comparator<BindingComponent> bcCompare = new Comparator<BindingComponent>() {
				@Override
				public int compare(BindingComponent s1, BindingComponent s2) {
					if(s1==null && s2==null)
						return 0;
					if(s1==null)
						return -1;
					if(s2==null)
						return 1;
					if(s1.getResponsibility()<s2.getResponsibility())
						return -1;
					if(s1.getResponsibility()>s2.getResponsibility())
						return 1;
					return 0;
				}
			};
			for(int e=0; e<manager.getNumConditions(); e++)
				Collections.sort(activeComponents.get(currReg).get(e), bcCompare);
        	
        	//Initialize component positions with active locations
        	List<Integer> componentPositions = new ArrayList<Integer>();
        	for(int e=0; e<manager.getNumConditions(); e++)
        		for(BindingComponent comp : activeComponents.get(currReg).get(e)){
        			if(!componentPositions.contains(comp.getPosition()) && comp.getPosition()>=currReg.getStart() && comp.getPosition()<currReg.getEnd()) {
        				componentPositions.add(comp.getPosition());
        			}
        		}
        	
        	//If no components exist in region, add one to the center to allow rescues
        	if(numBindingComponents==0){
        		componentPositions.add(currReg.getMidpoint().getLocation());
        		numBindingComponents++;
        	}
        	
        	//TODO: Test to find consensus nucleosome location
        	//Use max influence range to exclude those overlapping nucleosome with low responsibility
        	boolean[][] componentActive = new boolean[manager.getNumConditions()][componentPositions.size()];
        	for(int e=0; e<manager.getNumConditions(); e++) {
        		boolean[] exclusionZone = new boolean[currReg.getWidth()];       		
        		for(BindingComponent comp: activeComponents.get(currReg).get(e)) {
        			//Use exclusion zone to exclude those nucleosomes who are too close to other nucleosomes
        			int exclusion = 0;
        			if(mode.equals(EMmode.ALTERNATIVE))
        				exclusion = config.getAlternativeExclusionZone();
        			else if(mode.equals(EMmode.CONSENSUS))
        				exclusion = config.getConsensusExclusionZone();
    				int start = Math.max(0 , comp.getPosition()-currReg.getStart()-exclusion/2);
    				int end = Math.min(currReg.getWidth()-1, comp.getPosition()-currReg.getStart()+exclusion/2);
    				boolean isOverlap = false;
    				for(int i=start; i<=end; i++) {
    					if(exclusionZone[i])
    						isOverlap = true;
    				}
        			if(!isOverlap) {
        				for(int z=start; z<end; z++) {
        					exclusionZone[z] = true;
        				}
        				componentActive[e][componentPositions.indexOf(comp.getPosition())] = true;
        			} else {
        				componentActive[e][componentPositions.indexOf(comp.getPosition())] = false;
        			}
        		}
        	}

        	numBindingComponents = componentPositions.size();
        	
        	//Make new components with these locations
        	List<List<BindingComponent>> components = new ArrayList<List<BindingComponent>>();
        	for(int e=0; e<manager.getNumConditions(); e++)
        		components.add(new ArrayList<BindingComponent>());
        	
        	//Set up the components
        	for(int e=0; e<manager.getNumConditions(); e++){
        		double numActiveC=0;
        		for(boolean isActive: componentActive[e])
        			if(isActive)
        				numActiveC++;
        		int index=0;
        		double emission = (1-noise.get(e).getPi())/numActiveC;
	    		int numBindingSubtype = bindingManager.getBindingSubtypes(manager.getIndexedCondition(e)).size();
	    		double[] tauInit = new double[numBindingSubtype];
				for(int i=0; i<numBindingSubtype; i++) {
					tauInit[i] = (double)1/(double)numBindingSubtype;
				}
	    		for(Integer i : componentPositions){
	    			Point pos = new Point(config.getGenome(), currReg.getChrom(), i);
	    			BindingComponent currComp = new BindingComponent(pos, manager.getReplicates().size());
	    			currComp.setIndex(index);
	    			index++;
	    			
    				//Initialize normalized mixing probabilities (subtracting the noise emission probability)
	    			if(componentActive[e][componentPositions.indexOf(i)])
	    				currComp.uniformInit(emission);
	    			else
	    				currComp.uniformInit(0);
        			currComp.setFuzziness(bindingManager.getBindingModel(manager.getIndexedCondition(e)).get(0).getIntialFuzziness());
        			currComp.setTau(tauInit);
    				components.get(e).add(currComp);
    			}
    		}
        	return components; 
        }//end of initializeComponents method
        
        /**
         * Initializes components from active components in a single condition: 
         * 		Uses active component locations from the last round of training in one condition,
         * 		No flanking components or resuce components added here, since resulting components will only be used
         * 		in ML assignment.  
         *
         * @param currReg
         */
        private List<BindingComponent> initializeBindingComponentsFromOneConditionActive(Region currReg, NoiseComponent noise, int conditionIndex){
        	//Initialize component positions with active locations
        	List<Integer> componentPositions = new ArrayList<Integer>();
        	for(BindingComponent comp : activeComponents.get(currReg).get(conditionIndex)){
        		if(!componentPositions.contains(comp.getPosition()) && comp.getPosition()>=currReg.getStart() && comp.getPosition()<currReg.getEnd())
        			componentPositions.add(comp.getPosition());
        	}

        	numBindingComponents = componentPositions.size();

        	//Make new components with these locations
        	List<BindingComponent> components = new ArrayList<BindingComponent>();
        	
        	//Set up the components
        	double numC=(double)numBindingComponents; int index=0;
    		double emission = (1-noise.getPi())/numC;
    		for(Integer i : componentPositions){
    			Point pos = new Point(config.getGenome(), currReg.getChrom(), i);
    			BindingComponent currComp = new BindingComponent(pos, manager.getReplicates().size());
    			currComp.setIndex(index);
    			index++;
    			//Initialize normalized mixing probabilities (subtracting the noise emission probability)
    			currComp.uniformInit(emission);
				components.add(currComp);
			}
    		
        	return components; 
        }//end of initializeComponents method
        
        /**
         * Smooth the distribution of the control reads over the window to make a probability distribution of noise over the region
         * @param currReg
         * @param ctrlHits
         * @return
         */
        private double[] smoothNoiseDistribs(Region currReg, List<StrandedPair> ctrlHits) {
        	double[] distrib = new double[currReg.getWidth()];
        	double[] counts = new double[currReg.getWidth()];
        	//Pseudocounts for distrib
        	for(int d=0; d<currReg.getWidth(); d++)
        		counts[d] = 1;
        	//Add int count weights
        	for(StrandedPair hit: ctrlHits) {
        		int index = hit.getMidpoint().getLocation() - currReg.getStart();
        		if(index>=0 && index<currReg.getWidth())
        			counts[index] += hit.getWeight();
        	}
        	
        	//Smooth
        	for(int i=0; i<currReg.getWidth(); i++) {
        		double sum=0, num=0;
        		for(int s=i-(config.NOISE_DISTRIB_SMOOTHING_WIN/2); s<i+(config.NOISE_DISTRIB_SMOOTHING_WIN/2); s++) {
        			if(s>=0 && s<currReg.getWidth()) {
        				num++; sum+=counts[s];
        			}
        		}
        		distrib[i] = sum/num;
        	}
        	
        	//Normalize
        	double total = 0;
        	for(int d=0; d<currReg.getWidth(); d++)
        		total += distrib[d];
        	for(int d=0; d<currReg.getWidth(); d++)
        		distrib[d] = distrib[d]/total;
        	
        	return distrib;
        }
        
        /**
         * Normalize the distribution of the fragment size
         */
        private HashMap<Integer, Double> normalizeFragSizeFreq(Map<Integer, Double> freq){
        	HashMap<Integer, Double> normalizedFreq = new HashMap<Integer, Double>();
        	//Add pseudo fragment size into frequency map
        	for(int i=1; i<=1000; i++) {
        		normalizedFreq.put(i, 1d);
        	}
        	//Add fragment size frequency in freq to normalizedFreq
        	for(int key: freq.keySet()) {
        		double count = normalizedFreq.containsKey(key)? normalizedFreq.get(key):1;
        		normalizedFreq.put(key, count+freq.get(key));
        	}
        	//Normalize        	
			double totalCount = 0;
			for(int key: normalizedFreq.keySet()) {
				totalCount += normalizedFreq.get(key); 
			}
			for(int key: normalizedFreq.keySet()) {
				double oldValue = normalizedFreq.get(key);
				normalizedFreq.put(key, oldValue/totalCount);
			}	
        	return normalizedFreq;
        }
        
        
    	/**
    	 * Consolidate and edit binding events.
    	 * Follows a strict definition of binding event quantification - 
    	 * if the event is not present in the condition, it gets a zero count assigned.
    	 * Also merges positional duplicate events. 
    	 * @param ev
    	 * @return
    	 */
    	private List<BindingEvent> consolidateBindingEvents(List<BindingEvent> ev){
    		List<BindingEvent> newEvents = new ArrayList<BindingEvent>();
    		HashMap<Point, Integer> eventMap = new HashMap<Point, Integer>();
    		int count=0;
    		for(BindingEvent be : ev){
    			if(!eventMap.containsKey(be.getPoint())){
    				eventMap.put(be.getPoint(), count);
    				count++;
    				
    				newEvents.add(be);
    				//First time an event is added, clear out inactive events
    				for(ExperimentCondition cond : manager.getConditions()){
    					if(!be.isFoundInCondition(cond)){
    						be.setCondSigHits(cond, 0.0);
    		        		for(ControlledExperiment rep : cond.getReplicates()){
    		        			be.setRepSigHits(rep, 0.0);
    		        		}if(evconfig.CALC_EVENTS_LL)
    			            	be.setLLd(cond, 0);
    		        }	}
    			}else{
    				int index = eventMap.get(be.getPoint());
    				//For events that are already in the list, just update the active condition read counts
    				for(ExperimentCondition cond : manager.getConditions()){
    					if(be.isFoundInCondition(cond)){
    						newEvents.get(index).setCondSigHits(cond, be.getCondSigHits(cond));
    						newEvents.get(index).setCondCtrlHits(cond, be.getCondCtrlHits(cond));
    		        		for(ControlledExperiment rep : cond.getReplicates()){
    		        			newEvents.get(index).setRepSigHits(rep, be.getRepSigHits(rep));
    							newEvents.get(index).setRepCtrlHits(rep, be.getRepCtrlHits(rep));
    		        		}if(evconfig.CALC_EVENTS_LL)
    		        			newEvents.get(index).setLLd(cond, be.getLLd(cond));
    		        }	}
    			}
    		}
    		return newEvents;
    	}
        
        /**
         * ComponentConfiguration: represents a configuration of binding components as an array of positions 
         * @author Jianyu Yang
         * 
         */
        protected class ComponentConfiguration{
        	int[] positions = null;
        	int parentCondIndex;
        	//Constructor
        	public ComponentConfiguration(List<BindingComponent> comps, int parentCondition) {
        		Collections.sort(comps);
        		positions = new int[comps.size()];
        		for(int p=0; p<comps.size(); p++)
        			positions[p] = comps.get(p).getPosition();
        		parentCondIndex = parentCondition;
        	}
        	//Return the index of the condition where this configuration comes from
        	public int getParentCondition() {return parentCondIndex;}
        	//Compare two configurations
        	public boolean isSameAs(ComponentConfiguration cc) {
        		if(positions.length != cc.positions.length)
        			return false;
        		boolean isEqual = true;
        		for(int p=0; p<positions.length; p++)
        			isEqual = isEqual && positions[p]==cc.positions[p];
        		return isEqual;
        	}
        }
	}
}
