package org.seqcode.projects.sem.framework;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.seqcode.deepseq.StrandedPair;
import org.seqcode.projects.sem.events.BindingManager;
import org.seqcode.projects.sem.events.BindingModel;
import org.seqcode.projects.sem.events.EventsConfig;
import org.seqcode.projects.sem.utilities.PotentialRegionPoissonBackgroundModel;

import umontreal.ssj.util.Systeme;

import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.deepseq.experiments.Sample;
import org.seqcode.deepseq.stats.BackgroundCollection;
import org.seqcode.deepseq.stats.PoissonBackgroundModel;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.gsebricks.verbs.location.ChromosomeGenerator;
import org.seqcode.gseutils.RealValuedHistogram;


/**
 * PotentialRegionFilter: Find a set of regions that are above a threshold in at least one condition. 
 * 		A region the size of the model span (i.e. 2x model range) potentially contains a binding site if 
 * 		it passes the Poisson threshold in at least one condition.
 * 		The Poisson thresholds are based on the model span size to keep consistent with the final used thresholds. 
 * Overall counts for reads in potential regions and outside potential regions are maintained to assist noise model initialization.  
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class PotentialRegionFilter {

	protected ExperimentManager manager; 
	protected BindingManager bindingManager;
	protected EventsConfig evconfig;
	protected SEMConfig config;
	protected ExptConfig econfig;
	protected Genome gen;
	protected float maxBinWidth=0, binStep, winExt;
	protected boolean loadControl=true; 
	protected boolean stranded=false;
	protected List<Region> potentialRegions = new ArrayList<Region>();
	protected double potRegionLengthTotal=0;
	protected HashMap<ExperimentCondition, BackgroundCollection> conditionBackgrounds=new HashMap<ExperimentCondition, BackgroundCollection>(); //Background models for each replicate
	protected HashMap<ExperimentCondition, Double> potRegCountsSigChannel = new HashMap<ExperimentCondition, Double>();
	protected HashMap<ExperimentCondition, Double> nonPotRegCountsSigChannel = new HashMap<ExperimentCondition, Double>();
	protected HashMap<ExperimentCondition, Double> potRegCountsCtrlChannel = new HashMap<ExperimentCondition, Double>();
	protected HashMap<ExperimentCondition, Double> nonPotRegCountsCtrlChannel = new HashMap<ExperimentCondition, Double>();	
	protected HashMap<ExperimentCondition, Map<Integer, Double>> nonPotRegFragSizeFreqSigChannel = new HashMap<ExperimentCondition, Map<Integer, Double>>(); //Fragment size frequency in non potential region
	protected HashMap<ControlledExperiment, Double> potRegCountsSigChannelByRep = new HashMap<ControlledExperiment, Double>();
	protected HashMap<ControlledExperiment, Double> nonPotRegCountsSigChannelByRep = new HashMap<ControlledExperiment, Double>();
	
	public PotentialRegionFilter(EventsConfig ec, SEMConfig c, ExptConfig econ, ExperimentManager eman, BindingManager bman){
		manager = eman;
		bindingManager = bman;
		evconfig = ec;
		config = c; 
		econfig = econ;
		gen = config.getGenome();
		//Initialize background models
		for(ExperimentCondition cond : manager.getConditions()){
			conditionBackgrounds.put(cond, new BackgroundCollection());
			int maxIR = bindingManager.getMaxInfluenceRange(cond); 
			boolean hasControls=true; 
			float binWidth = maxIR;
    		if(binWidth>maxBinWidth){maxBinWidth=binWidth;}
    			
    		//global threshold
    		conditionBackgrounds.get(cond).addBackgroundModel(new PotentialRegionPoissonBackgroundModel(-1, config.getPRLogConf(), cond.getTotalSignalPairCount(), config.getGenome().getGenomeLength(), econfig.getMappableGenomeProp(), binWidth, '.', 1, true));
    		//local windows won't work since we are testing per condition and we don't have a way to scale signal vs controls at the condition level (at least at this stage of execution)
    		
    		double thres = conditionBackgrounds.get(cond).getGenomicModelThreshold();
    		System.err.println("PotentialRegionFilter: genomic threshold for "+cond.getName()+" with bin width "+binWidth+" = "+thres);
    			
    		//Initialize counts
    		potRegCountsSigChannel.put(cond, 0.0);
    		nonPotRegCountsSigChannel.put(cond, 0.0);
    		potRegCountsCtrlChannel.put(cond, 0.0);
    		nonPotRegCountsCtrlChannel.put(cond, 0.0);
    		for(ControlledExperiment rep : cond.getReplicates()){
    			potRegCountsSigChannelByRep.put(rep, 0.0);
        		nonPotRegCountsSigChannelByRep.put(rep, 0.0);
    		}
    		nonPotRegFragSizeFreqSigChannel.put(cond, new HashMap<Integer, Double>());
    	}
		binStep = config.POTREG_BIN_STEP;
		if(binStep>maxBinWidth/2)
			binStep=maxBinWidth/2;
		winExt = maxBinWidth/2;
	}
	
	//Accessors for read counts
	public Double getPotRegCountsSigChannel(ExperimentCondition e){ return potRegCountsSigChannel.get(e);}
	public Double getNonPotRegCountsSigChannel(ExperimentCondition e){ return nonPotRegCountsSigChannel.get(e);}
	public Double getPotRegCountsCtrlChannel(ExperimentCondition e){ return potRegCountsCtrlChannel.get(e);}
	public Double getNonPotRegCountsCtrlChannel(ExperimentCondition e){ return nonPotRegCountsCtrlChannel.get(e);}
	public Double getPotRegCountsSigChannelByRep(ControlledExperiment e){ return potRegCountsSigChannelByRep.get(e);}
	public Double getNonPotRegCountsSigChannelByRep(ControlledExperiment e){ return nonPotRegCountsSigChannelByRep.get(e);}
	public HashMap<ExperimentCondition, Map<Integer, Double>> getNonPotRegFragSizeFreqSigChannel() { return nonPotRegFragSizeFreqSigChannel;}
	public List<Region> getPotentialRegions(){return potentialRegions;}
	public double getPotRegionLengthTotal(){return potRegionLengthTotal;}
	
	/**
	 * Find list of potentially enriched regions 
	 * (windows that contain the minimum number of reads needed to pass the Poisson backgrounds).
	 * @param testRegions
	 */
	public List<Region> execute(){
		//TODO: check config for defined subset of regions
		Iterator<Region> testRegions = new ChromosomeGenerator().execute(config.getGenome());
		
		//Threading divides analysis over entire chromosomes. This approach is not compatible with file caching. 
		int numThreads = econfig.getCacheAllData() ? config.getMaxThreads() : 1;
				
		Thread[] threads = new Thread[numThreads];
        ArrayList<Region> threadRegions[] = new ArrayList[numThreads];
        int i = 0;
        for (i = 0 ; i < threads.length; i++) {
            threadRegions[i] = new ArrayList<Region>();
        }i=0;
        while(testRegions.hasNext()){
        	Region r = testRegions.next(); 
            threadRegions[(i++) % numThreads].add(r);
        }

        for (i = 0 ; i < threads.length; i++) {
            Thread t = new Thread(new PotentialRegionFinderThread(threadRegions[i]));
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
        
        //Initialize signal & noise counts based on potential region calls
        for(ExperimentCondition cond : manager.getConditions()){
    		for(ControlledExperiment rep : cond.getReplicates()){
    			if(rep.getSignalVsNoiseFraction()==0) { //Only update if not already initialized
    				System.out.println(potRegCountsSigChannelByRep.get(rep));
    				System.out.println(potRegCountsSigChannelByRep.get(rep)+nonPotRegCountsSigChannelByRep.get(rep));
    				rep.setSignalVsNoiseFraction(potRegCountsSigChannelByRep.get(rep)/(potRegCountsSigChannelByRep.get(rep)+nonPotRegCountsSigChannelByRep.get(rep)));
    			}
    		}
        }
        
        for(Region r : potentialRegions)
        	potRegionLengthTotal+=(double)r.getWidth();
        
     	return potentialRegions;
	}
	
	/**
     * Print potential regions to a file.
     * TESTING ONLY 
     */
    public void printPotentialRegionsToFile(){
    	try {
    		String filename = config.getOutputIntermediateDir()+File.separator+config.getOutBase()+".potential.regions";
			FileWriter fout = new FileWriter(filename);
			for(Region r : potentialRegions){
	    		fout.write(r.getLocationString()+"\n");			
	    	}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
    }
	
    class PotentialRegionFinderThread implements Runnable {
        private Collection<Region> regions;
        private float[][] landscape=null;
        private float[][] starts=null;
        private List<Region> threadPotentials = new ArrayList<Region>();
        
        public PotentialRegionFinderThread(Collection<Region> r) {
            regions = r;
        }
        
        public void run() {
        	int expansion = (int)(winExt + maxBinWidth/2);
        	for (Region currentRegion : regions) {
            	Region lastPotential=null;
                List<List<StrandedPair>> ipHits = new ArrayList<List<StrandedPair>>();
                List<List<StrandedPair>> backHits = new ArrayList<List<StrandedPair>>();
                List<List<StrandedPair>> ipHitsByRep = new ArrayList<List<StrandedPair>>();
                //Split the job up into large chunks
                for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=config.MAXSECTION){
                    int y = (int) (x+config.MAXSECTION+(expansion)); //Leave a little overhang to handle enriched regions that may hit the border. Since lastPotential is defined above, a region on the boundary should get merged in.
                    if(y>currentRegion.getEnd()){y=currentRegion.getEnd();}
                    Region currSubRegion = new Region(gen, currentRegion.getChrom(), x, y);
                    
                    List<Region> currPotRegions = new ArrayList<Region>();
                    ipHits = new ArrayList<List<StrandedPair>>();
                    backHits = new ArrayList<List<StrandedPair>>();
                    ipHitsByRep = new ArrayList<List<StrandedPair>>();
                    
                    synchronized(manager){
	                    //Initialize the read lists
                    	for(ExperimentCondition cond : manager.getConditions()){
                    		ipHits.add(new ArrayList<StrandedPair>());
                			backHits.add(new ArrayList<StrandedPair>());
                    		for(ControlledExperiment rep : cond.getReplicates())
                    			ipHitsByRep.add(new ArrayList<StrandedPair>());
                    	}
                    	//Load signal reads by condition and by replicate, so that signal proportion estimates can be assigned to each replicate 
                    	for(ExperimentCondition cond : manager.getConditions()){
                    		for(ControlledExperiment rep : cond.getReplicates()){
                    			ipHits.get(cond.getIndex()).addAll(rep.getSignal().getPairsByMid(currSubRegion));
                    			ipHitsByRep.get(rep.getIndex()).addAll(rep.getSignal().getPairsByMid(currSubRegion));
                    		}for(Sample ctrl : cond.getControlSamples())
                    			backHits.get(cond.getIndex()).addAll(ctrl.getPairsByMid(currSubRegion));
                    		Collections.sort(ipHits.get(cond.getIndex()));
                    		Collections.sort(backHits.get(cond.getIndex()));
                    	}
                    }
                    
            		int numStrandIter = stranded ? 2:1;
                    for(int stranditer=1; stranditer<=numStrandIter; stranditer++){
                        //If stranded peak-finding, run over both strands separately
                        char str = !stranded ? '.' : (stranditer==1 ? '+' : '-');
					 
                        makeHitLandscape(ipHits, currSubRegion, maxBinWidth, binStep, str);
                        float ipHitCounts[][] = landscape.clone();
                        float ipBinnedStarts[][] = starts.clone();
                        float backBinnedStarts[][] = null;
                        if (loadControl) {
                            makeHitLandscape(backHits, currSubRegion, maxBinWidth, binStep, str);
                            backBinnedStarts = starts.clone();
                        }
					
                        //Scan regions
                        int currBin=0;
                        for(int i=currSubRegion.getStart(); i<currSubRegion.getEnd()-(int)maxBinWidth; i+=(int)binStep){
                        	boolean regionPasses=false;
                        	for(ExperimentCondition cond : manager.getConditions()){
                        		double ipWinHits=ipHitCounts[cond.getIndex()][currBin];
                        		//First Test: is the read count above the genome-wide thresholds?
                        		//If there is a fixed alpha, we should use that as the only threshold
//                        		if(config.getFixedAlpha()>0) {
//                        			if(ipWinHits>config.getFixedAlpha()){
//                        				regionPasses=true;
//                        				break;
//                        			}
//                        		}else if(conditionBackgrounds.get(cond).passesGenomicThreshold((int)ipWinHits, str)){
//                        			//Second Test: refresh all thresholds & test again
//                        			conditionBackgrounds.get(cond).updateModels(currSubRegion, i-x, ipBinnedStarts[cond.getIndex()], backBinnedStarts==null ? null : backBinnedStarts[cond.getIndex()], binStep);
//                        			if(conditionBackgrounds.get(cond).passesAllThresholds((int)ipWinHits, str)){
//                        				//If the region passes the thresholds for one condition, it's a potential
//                        				regionPasses=true;
//                        				break;
//		                            }
//		                        }
                        		// Test: is the read count above the alpha threshold
                        		// if fixedAlpha is set (>=0), potential region needs counts above fixedAlpha, else needs counts above 1
                        		if(config.getFixedAlpha()>=0) {
	                        		if(ipWinHits>config.getFixedAlpha()) {
	                        			regionPasses = true;
	                        			break;
	                        		}
                        		}
	                        	else {
	                        		if(ipWinHits>1) {	// if doesn't set fixedalpha, select region having at least 2 fragments as potential region
	                        			regionPasses = true;
	                        			break;
	                        		}
	                        	}
                        	}
                        	if(regionPasses){
                        		Region currPotential = new Region(gen, currentRegion.getChrom(), Math.max(i-expansion, 1), Math.min((int)(i-1+expansion), currentRegion.getEnd()));
                        		if(lastPotential!=null && currPotential.overlaps(lastPotential)){
                        			lastPotential = lastPotential.expand(0, currPotential.getEnd()-lastPotential.getEnd());
                        		}else{
                        			//Add the last recorded region to the list
                        			if(lastPotential!=null){
                        				if(lastPotential.getWidth()<=config.getBMAnalysisWindowMax()){
                        					currPotRegions.add(lastPotential);
                        					threadPotentials.add(lastPotential);
                        				}else{
                        					//Break up long windows
                        					List<Region> parts = breakWindow(lastPotential, ipHits, config.getBMAnalysisWindowMax(), str);
                        					for(Region p : parts){
                        						currPotRegions.add(p);
                        						threadPotentials.add(p);
                        					}
                        				}
                        			}lastPotential = currPotential;
                        		}
                        	}
                            currBin++;
                        }
					}
                    //Count all "signal" reads overlapping the regions in currPotRegions (including the lastPotential)
                    if(lastPotential!=null)
                    	currPotRegions.add(lastPotential);
                    currPotRegions = filterExcluded(currPotRegions);
                    countReadsInRegions(currPotRegions, ipHits, backHits, y==currentRegion.getEnd() ? y : y-expansion);
                    countReadsInRegionsByRep(currPotRegions, ipHitsByRep, y==currentRegion.getEnd() ? y : y-expansion);
                    //Note: it looks like currPotRegions and threadPotentials are redundant in the above, but they are not.
                    //currPotRegions is only used to count sig/noise reads in the current section. threadPotentials stores regions over the entire run.
                }
                //Add the final recorded region to the list
                //Warning: For SEM, it is possible that lastPotential region is quite large (e.g., millions of base pairs)
                //TODO: I add break step here to avoid a too large potential region
                if(lastPotential!=null) {
                	//break lastPotential
					List<Region> parts = breakWindow(lastPotential, ipHits, config.getBMAnalysisWindowMax(), '.');
					for(Region p: parts) {
						threadPotentials.add(p);
					}
                }
                threadPotentials = filterExcluded(threadPotentials);
            }
        	if(threadPotentials.size()>0){
        		synchronized(potentialRegions){
        			potentialRegions.addAll(threadPotentials);
        		}
        	}	
        }
        
        //Break up a long window into parts
        //For now, we just choose the break points as the bins with the lowest total signal read count around the desired length.
        //TODO: improve?
        protected List<Region> breakWindow(Region lastPotential, List<List<StrandedPair>> ipHits, int preferredWinLen, char str) {
			List<Region> parts = new ArrayList<Region>();
			makeHitLandscape(ipHits, lastPotential, maxBinWidth, binStep, str);
            float ipHitCounts[][] = landscape.clone();
            
            int currPartStart = lastPotential.getStart();
            double currPartTotalMin=Double.MAX_VALUE; int currPartTotalMinPos = -1;
            int currBin=0;
            for(int i=lastPotential.getStart(); i<lastPotential.getEnd()-(int)maxBinWidth; i+=(int)binStep){
            	if(lastPotential.getEnd()-currPartStart < (preferredWinLen*1.5))
            		break;
            	float currBinTotal=0;
            	for(ExperimentCondition cond : manager.getConditions())
                	currBinTotal+=ipHitCounts[cond.getIndex()][currBin];
            	
            	if(i>(currPartStart+preferredWinLen-1000) && i<(currPartStart+preferredWinLen+1000)){ 
            		if(currBinTotal<currPartTotalMin){
            			currPartTotalMin=currBinTotal;
            			currPartTotalMinPos=i;
            		}
            	}
            	//Add a new part
            	if(i>=(currPartStart+preferredWinLen+1000)){
            		parts.add(new Region(lastPotential.getGenome(), lastPotential.getChrom(), currPartStart, currPartTotalMinPos));
            		currPartStart = currPartTotalMinPos+1;
            		currPartTotalMin=Double.MAX_VALUE; currPartTotalMinPos = -1;
            	}
            	currBin++;
            }
            parts.add(new Region(lastPotential.getGenome(), lastPotential.getChrom(), currPartStart, lastPotential.getEnd()));
            
			return parts;
		}

		//Filter out pre-defined regions to ignore (e.g. tower regions)
        protected List<Region> filterExcluded(List<Region> testRegions) {
			List<Region> filtered = new ArrayList<Region>();
			if(config.getRegionsToIgnore().size()==0)
				return testRegions;
			
			for(Region t : testRegions){
				boolean ignore = false;
				for(Region i : config.getRegionsToIgnore()){
					if(t.overlaps(i)){
						ignore = true; break;
					}
				}
				if(!ignore)
					filtered.add(t);
			}
			return filtered;
		}

		//Makes integer arrays corresponding to the read landscape over the current region.
        //Reads are semi-extended out to bin width to account for the bin step
        //No needlefiltering here as that is taken care of during read loading (i.e. in Sample)
    	protected void makeHitLandscape(List<List<StrandedPair>> hits, Region currReg, float binWidth, float binStep, char strand){
    		int numBins = (int)(currReg.getWidth()/binStep);
    		landscape = new float[hits.size()][numBins+1];
    		starts = new float[hits.size()][numBins+1];
    		float halfWidth = binWidth/2;

    		for(ExperimentCondition cond : manager.getConditions()){
        		List<StrandedPair> currHits = hits.get(cond.getIndex());
    			for(int i=0; i<=numBins; i++){landscape[cond.getIndex()][i]=0; starts[cond.getIndex()][i]=0; }
	    		for(StrandedPair r : currHits){
	    			int offset=inBounds(r.getMidpoint().getLocation()-currReg.getStart(),0,currReg.getWidth());
	    			int binoff = inBounds((int)(offset/binStep), 0, numBins);
	    			starts[cond.getIndex()][binoff]+=r.getWeight();
	    			int binstart = inBounds((int)((double)(offset-halfWidth)/binStep), 0, numBins);
	    			int binend = inBounds((int)((double)(offset+halfWidth)/binStep), 0, numBins);
	    			for(int b=binstart; b<=binend; b++)
	    				landscape[cond.getIndex()][b]+=r.getWeight();
	    		}
           	}
    	}
    	protected final int inBounds(int x, int min, int max){
    		if(x<min){return min;}
    		if(x>max){return max;}
    		return x;
    	}
    	
    	/**
    	 * Count the total reads within potential regions via semi binary search (by condition).
    	 * Assumes both regs and ipHits are sorted.
    	 * We don't have to check chr String matches, as the hits were extracted from the chromosome
    	 * EndCoord accounts for the extra overhang added to some wide regions
    	 * We also ignore strandedness here -- the object is to count ALL reads that will be loaded for analysis later
    	 * (and that thus will not be accounted for by the global noise model)  
    	 * @param regs
    	 * @param ipHits
    	 * @param ctrlHits
    	 * @param endCoord
    	 */
    	protected void countReadsInRegions(List<Region> regs, List<List<StrandedPair>> ipHits, List<List<StrandedPair>> ctrlHits, int endCoord){
    		//Iterate through experiments
    		for(ExperimentCondition cond : manager.getConditions()){
    			double currPotWeightSig=0, currNonPotWeightSig=0, currPotWeightCtrl=0, currNonPotWeightCtrl=0;
    			//Iterate through signal hits
    			for(StrandedPair hit : ipHits.get(cond.getIndex())){
    				if(regs.size()==0)
    					currNonPotWeightSig+=hit.getWeight();
    				else{
    					//Binary search for closest region start
        				int hpoint = hit.getMidpoint().getLocation();
        				if(hpoint<endCoord){ //Throw this check in for the overhang
	        				int l = 0, r = regs.size()-1;
	        	            while (r - l > 1) {
	        	                int c = (l + r) / 2;
	        	                if (hpoint >= regs.get(c).getStart()) {
	        	                    l = c;
	        	                } else {
	        	                    r = c;
	        	                }
	        	            }
	        	            boolean inPot = false;
	        	            for(int x=l; x<=r; x++){
	        	            	if(hpoint >= regs.get(x).getStart() && hpoint <= regs.get(x).getEnd()){
	        	            		currPotWeightSig+=hit.getWeight(); 
	        	            		inPot=true;
	        	            		break;
	        	            	}
	        	            }
	        	            if(!inPot) {
	        	            	currNonPotWeightSig+=hit.getWeight();
	        	            	//Add hit to frequency channel
	        	            	try {
	        	            		double frequency = nonPotRegFragSizeFreqSigChannel.get(cond).containsKey(hit.getFragmentSize())? nonPotRegFragSizeFreqSigChannel.get(cond).get(hit.getFragmentSize()):0;
        	            		nonPotRegFragSizeFreqSigChannel.get(cond).put(hit.getFragmentSize(), frequency+hit.getWeight());
	        	            	} catch (Exception e) {
	        	            		e.printStackTrace();
	        	            		System.out.println(cond.getName());
	        	            		System.out.println(hit.getFragmentSize());
	        	            		System.out.println(hit.toString());
	        	            		System.exit(1);
	        	            	}
	        	            }
        				}
    				}
    			}
    			//Iterate through control hits
    			for(StrandedPair hit : ctrlHits.get(cond.getIndex())){
    				if(regs.size()==0)
    					currNonPotWeightCtrl+=hit.getWeight();
    				else{
        				//Binary search for closest region start
        				int hpoint = hit.getMidpoint().getLocation();
        				if(hpoint<endCoord){ //Throw this check in for the overhang
	        				int l = 0, r = regs.size()-1;
	        	            while (r - l > 1) {
	        	                int c = (l + r) / 2;
	        	                if (hpoint >= regs.get(c).getStart()) {
	        	                    l = c;
	        	                } else {
	        	                    r = c;
	        	                }
	        	            }
	        	            boolean inPot = false;
	        	            for(int x=l; x<=r; x++){
	        	            	if(hpoint >= regs.get(x).getStart() && hpoint <= regs.get(x).getEnd()){
	        	            		currPotWeightCtrl+=hit.getWeight(); inPot=true; break;
	        	            	}
	        	            }
	        	            if(!inPot) {
	        	            	currNonPotWeightCtrl+=hit.getWeight();
	        	            }
        				}
    				}
    			}
    			synchronized(potRegCountsSigChannel){
    				potRegCountsSigChannel.put(cond, potRegCountsSigChannel.get(cond)+currPotWeightSig);
    			}
    			synchronized(nonPotRegCountsSigChannel){
    				nonPotRegCountsSigChannel.put(cond, nonPotRegCountsSigChannel.get(cond)+currNonPotWeightSig);
    			}
    			synchronized(potRegCountsCtrlChannel){
    				potRegCountsCtrlChannel.put(cond, potRegCountsCtrlChannel.get(cond)+currPotWeightCtrl);
    			}
    			synchronized(nonPotRegCountsCtrlChannel){
    				nonPotRegCountsCtrlChannel.put(cond, nonPotRegCountsCtrlChannel.get(cond)+currNonPotWeightCtrl);
    			}
    		}
	    }
    }
    
    /**
	 * Count the total reads within potential regions via semi binary search (by replicate).
	 * Assumes both regs and ipHitsByRep are sorted.
	 * We don't have to check chr String matches, as the hits were extracted from the chromosome
	 * EndCoord accounts for the extra overhang added to some wide regions
	 * We also ignore strandedness here -- the object is to count ALL reads that will be loaded for analysis later
	 * (and that thus will not be accounted for by the global noise model)  
	 * @param regs
	 * @param ipHitsByRep
	 * @param ctrlHits
	 * @param endCoord
	 */
	protected void countReadsInRegionsByRep(List<Region> regs, List<List<StrandedPair>> ipHitsByRep, int endCoord){
		//Iterate through experiments
		for(ExperimentCondition cond : manager.getConditions()){
			for(ControlledExperiment rep : cond.getReplicates()){
				double currPotWeightSig=0, currNonPotWeightSig=0;
				//Iterate through signal hits
				for(StrandedPair hit : ipHitsByRep.get(rep.getIndex())){
					if(regs.size()==0)
						currNonPotWeightSig+=hit.getWeight();
					else{
						//Binary search for closest region start
	    				int hpoint = hit.getMidpoint().getLocation();
	    				if(hpoint<endCoord){ //Throw this check in for the overhang
	        				int l = 0, r = regs.size()-1;
	        	            while (r - l > 1) {
	        	                int c = (l + r) / 2;
	        	                if (hpoint >= regs.get(c).getStart()) {
	        	                    l = c;
	        	                } else {
	        	                    r = c;
	        	                }
	        	            }
	        	            boolean inPot = false;
	        	            for(int x=l; x<=r; x++){
	        	            	if(hpoint >= regs.get(x).getStart() && hpoint <= regs.get(x).getEnd()){
	        	            		currPotWeightSig+=hit.getWeight(); inPot=true; break;
	        	            	}
	        	            }
	        	            if(!inPot)
	        	            	currNonPotWeightSig+=hit.getWeight();
	    				}
					}
				}
				
				synchronized(potRegCountsSigChannelByRep){
					potRegCountsSigChannelByRep.put(rep, potRegCountsSigChannelByRep.get(rep)+currPotWeightSig);
				}
				synchronized(nonPotRegCountsSigChannelByRep){
					nonPotRegCountsSigChannelByRep.put(rep, nonPotRegCountsSigChannelByRep.get(rep)+currNonPotWeightSig);
				}
			}
		}
    }

}

