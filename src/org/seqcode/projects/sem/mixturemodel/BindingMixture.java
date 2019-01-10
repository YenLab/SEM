package org.seqcode.projects.sem.mixturemodel;

import java.io.*;
import java.util.*;
import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.StrandedPair;
import org.seqcode.gseutils.Pair;
import org.seqcode.projects.sem.mixturemodel.BindingComponent;
import org.seqcode.projects.sem.mixturemodel.NoiseComponent;
import org.seqcode.projects.sem.events.*;
import org.seqcode.projects.sem.framework.*;
import org.seqcode.deepseq.experiments.*;
import org.seqcode.deepseq.stats.BackgroundCollection;
import org.seqcode.deepseq.stats.PoissonBackgroundModel;
import org.seqcode.genome.*;
import org.seqcode.genome.location.*;
import org.seqcode.gsebricks.verbs.location.ChromosomeGenerator;
import org.seqcode.deepseq.stats.BackgroundCollection;
import org.seqcode.deepseq.stats.PoissonBackgroundModel;


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
	protected HashMap<ExperimentCondition, BackgroundCollection> conditionBackgrounds = new HashMap<ExperimentCondition, BackgroundCollection>(); //Genomic Background models for each condition -- used to set alpha values in sparse prior
	protected List<BindingEvent> bindingEvents;
	protected List<Region> regionsToPlot;
	protected int trainingRound = 0;
	protected double noisePerBase[];		//Defines global noise
	protected double relativeCtrlNoise[];	//Defines global noise
	protected HashMap<Region, Double[]> noiseResp = new HashMap<Region, Double[]>(); //noise responsibilities after a round of execute(). Hashed by Region, indexed by condition
	
	public BindingMixture(GenomeConfig gcon, ExptConfig econ, EventsConfig evcon, SEMConfig semconfig, ExperimentManager eMan, BindingManager bMan, PotentialRegionFilter filter) {
		gconfig = gcon;
		econfig = econ;
		evconfig = evcon;
		config = semconfig;
		manager = eMan;
		bindingManager = bMan;
		potRegFilter = filter;
		testRegions = filter.getPotentialRegions();
		regionsToPlot = config.getRegionsToPlot();
		bindingEvents = new ArrayList<BindingEvent>();
		BindingEvent.setExperimentManager(manager);
		BindingEvent.setConfig(evconfig);
		
		activeComponents = new HashMap<Region, List<List<BindingComponent>>>();
		for(ExperimentCondition cond: manager.getConditions()) {
			conditionBackgrounds.put(cond, new BackgroundCollection());
			conditionBackgrounds.get(cond).addBackgroundModel(new PoissonBackgroundModel(-1, config.getSigLogConf(), cond.getTotalSignalCount()*(1-cond.getTotalSignalVsNoiseFrac()), config.getGenome().getGenomeLength(), econfig.getMappableGenomeProp(), bindingManager.getMaxInfluenceRange(cond), '.', 1, true));
			double alf = config.getFixedAlpha()>0 ? config.getFixedAlpha() : (double)conditionBackgrounds.get(cond).getMaxThreshold('.');
			System.err.println("Alpha "+cond.getName()+"\tRange="+bindingManager.getMaxInfluenceRange(cond)+"\t"+alf);
		}
		
		noisePerBase = new double[manager.getNumConditions()];
		relativeCtrlNoise = new double[manager.getNumConditions()];
		initializeGlobalNoise();
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
    		System.err.println("Global noise per base initialization for "+cond.getName()+" = "+String.format("%.4f", noisePerBase[e]));
    	}
    }
    
	
	public void execute(boolean EM, boolean uniformBindingComponents) {
		trainingRound++;
		
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
		            Thread t = new Thread(new BindingMixtureThread(threadRegions[i], EM, uniformBindingComponents));
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
	 * BindingMixtureThread: run binding mixtures over a subset of regions
	 * uniformComponents: if true, components are initialized at uniform prob and spacing.
	 * 					otherwise, components are initialized from activeComponents
	 */
	class BindingMixtureThread implements Runnable{
		private Collection<Region> regions;
		private int numBindingComponents = 1; 
		private boolean runEM = true;
		private boolean uniformBindingComponents = false;
		
		public BindingMixtureThread(Collection<Region> regs, boolean EM, boolean uniformBindingComponents) {
			regions = regs;
			this.uniformBindingComponents = uniformBindingComponents;
			runEM = EM;
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
						}
						
						//Only non-zero components are returned by analyzeWindow, so add them to the recorded active components
						synchronized(activeComponents) { activeComponents.put(rr, currComps);}
						
						//Add the sum of noise responsibilities to this region
						synchronized(noiseResp) { noiseResp.put(rr, noiseRSums);}
					} else {
						//Run ML assignment
						List<BindingEvent> windowBindingEvents = new ArrayList<BindingEvent>();
						for (Region w: windows) {
							windowBindingEvents.addAll(analyzeWindowML(w));
						}
						synchronized(bindingEvents) {bindingEvents.addAll(windowBindingEvents);}
					}
				} catch(Exception e) {
					System.err.println("ERROR: Exception when analyzing region"+rr.toString());
					e.printStackTrace(System.err);
					System.exit(-1);
				}
			}
		}
	
		private Pair<List<NoiseComponent>, List<List<BindingComponent>>> analyzeWindowEM(Region w) {
			BindingEM EM = new BindingEM(config, manager, bindingManager, conditionBackgrounds, potRegFilter.getPotentialRegions().size());
			List<List<BindingComponent>> bindingComponents = null;
			List<NoiseComponent> noiseComponents = null;
			List<List<BindingComponent>> nonZeroComponents = new ArrayList<List<BindingComponent>>();
			for(int e=0; e<manager.getNumConditions(); e++)
				nonZeroComponents.add(new ArrayList<BindingComponent>());
			
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
            	bindingComponents = initializeBindingComponentsFromAllConditionActive(w, noiseComponents, true);
            
            //ATAC-seq prior
            double[][] atacPrior;
            
            //EM learning: resulting binding components list will only contain non-zero components
            nonZeroComponents = EM.train(signals, w, noiseComponents, bindingComponents, numBindingComponents, atacPrior, trainingRound);
            
            return new Pair<List<NoiseComponent>, List<List<BindingComponent>>>(noiseComponents, nonZeroComponents);
		}
		
		/**
		 * Assign BindingComponents over a given window with ML solution
		 * 
		 * @param w
		 * @return Pair of component lists (noise components and binding components) indexed by condition
		 */
		private List<BindingEvent> analyzeWindowML(Region w){
			BindingMLAssignment ML = new BindingMLAssignment(econfig, evconfig, config, manager, bindingManager, conditionBackgrounds, potRegFilter.getPotentialRegions());
			List<BindingComponent> bindingComponents = null;
			List<NoiseComponent> noiseComponents = null;
			List<BindingEvent> currEvents = new ArrayList<BindingEvent>();
			
			//Load signal data
			List<List<StrandedPair>> signals = loadSignalData(w);
			if (signals==null)
				return currEvents;
			
			//Load control data
			List<List<StrandedPair>> controls = loadControlData(w);
			
			//Initialize noise components
			noiseComponents = initializeNoiseComponents(w, signals, controls);
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
					data.get(rep.getIndex()).addAll(rep.getSignal().getPairs(w));
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
						data.get(rep.getIndex()).addAll(rep.getControl().getPairs(w));
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
    		
				//Add the noise component
				NoiseComponent n = new NoiseComponent(emission, distribs, currReg, numReps);
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
				double emission = (1-noise.get(e).getPi())/numC;
				for(BindingComponent b : components.get(e)){
					b.uniformInit(emission);	                
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
        private List<List<BindingComponent>> initializeBindingComponentsFromAllConditionActive(Region currReg, List<NoiseComponent> noise, boolean addFlanking){
        	//Initialize component positions with active locations
        	List<Integer> componentPositions = new ArrayList<Integer>();
        	for(int e=0; e<manager.getNumConditions(); e++)
        		for(BindingComponent comp : activeComponents.get(currReg).get(e)){
        			if(!componentPositions.contains(comp.getPosition()) && comp.getPosition()>=currReg.getStart() && comp.getPosition()<currReg.getEnd())
        				componentPositions.add(comp.getPosition());
        			if(addFlanking){
        				if(!componentPositions.contains(comp.getPosition()-config.getAddFlankingComponentSpacing())
        						 && comp.getPosition()-config.getAddFlankingComponentSpacing()>=currReg.getStart())
        					componentPositions.add(comp.getPosition()-config.getAddFlankingComponentSpacing());
        				if(!componentPositions.contains(comp.getPosition()+config.getAddFlankingComponentSpacing())
        						&& comp.getPosition()+config.getAddFlankingComponentSpacing()<currReg.getEnd())
        					componentPositions.add(comp.getPosition()+config.getAddFlankingComponentSpacing());
        			}
        		}

        	numBindingComponents = componentPositions.size();
        	
        	//If no components exist in region, add one to the center to allow rescues
        	if(numBindingComponents==0 && addFlanking){
        		componentPositions.add(currReg.getMidpoint().getLocation());
        		numBindingComponents++;
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
        		counts[d] = 0;
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
	}
}
