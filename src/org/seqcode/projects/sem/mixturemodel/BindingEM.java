package org.seqcode.projects.sem.mixturemodel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Comparator;

import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.StrandedPair;
import org.seqcode.deepseq.events.BindingModelPerBase;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.stats.BackgroundCollection;
import org.seqcode.genome.location.Region;
import org.seqcode.projects.sem.framework.SEMConfig;
import org.seqcode.projects.sem.events.BindingManager;
import org.seqcode.projects.sem.events.BindingModel;
import org.seqcode.projects.sem.events.BindingSubtype;
import org.seqcode.projects.sem.utilities.Timer;
import org.seqcode.gseutils.Pair;
import org.seqcode.math.stats.StatUtil;

public class BindingEM {
	protected ExperimentManager manager;
	protected BindingManager bindingManager;
	protected SEMConfig semconfig;
	protected List<List<BindingComponent>> components;
	protected List<NoiseComponent> noise;
	protected int numComponents;
	protected int numConditions;
	protected int numSubtypes;
	protected int trainingRound=0;
	protected HashMap<ExperimentCondition, BackgroundCollection> conditionBackgrounds;
//	protected HashMap<ExperimentCondition, Integer> numBindingSubtypes;
	protected Timer timer;
	
	// EM VARIABLES
	// H function and responsibility have to account for all reads in region now, as they will be updated
	// once the component positions change
	protected double[][]	hitCounts;			// Hit weights
	protected int[][]		hitPos;				// Hit positions
	protected int[][]		hitSize;			// &Hit fragment sizes 
	protected int[]			hitNum;				// Number of hits in each condition
	protected int[][]		repIndices;			// Index of replicate for the hit
	protected double[][][]	hAll;				// H function values for all positions in the current window (precomputed)
	protected double[][][]	h;					// H function (binding component probability per read)
	protected double[][]	n;					// N function (noise component probability per read)
	protected double[][][]	rBind;				// Binding component responsibilities
	protected double[][]	rNoise;				// Noise component responsibilities
	protected double[][]	sumResp;			// Sum of responsibility each bindingComponent
	protected double[][]	pi;					// pi: emission probabilities for binding components
	protected double[]		piNoise;			// piNoise: emission probabilities for noise components (fixed)
	protected int[][]		mu;					// mu: positions of the binding components (indexed by condition & binding component index)
	protected double[][]	fuzz;				// &fuzz: fuzziness of the binding components (indexed by condition & binding component index)
	protected double[][][]	tau;				// &tau: fragment size subtype probabilities (indexed by subtype index)
	protected double[]		alphaMax;			// Maximum alpha (indexed by condition)
	protected double[][]	atacPrior;			// &ATAC-seq prior (indexed by condition & base)
	protected double[][][]	lastRBind;			// Last responsibilities (monitor convergence)
	protected double[][]	lastPi;				// Last Pi (monitor convergence)
	protected int[][]		lastMu;				// Last positions (monitor convergence)
	protected double[][]	lastFuzz;			// &Last fuzziness (monitor convergence)
	protected double[][][]	lastTau;				// &Last tau (monitor convergence)
	protected double 		lastLAP, LAP;		// log-likelihood monitoring
	protected boolean 		plotEM=false;		// Plot the current region components
	protected Region		plotSubRegion=null;	// Sub region to plot
	protected double		numPotentialRegions;
	protected double		probAgivenB, probAgivenNOTB;
	protected int stateEquivCount = 0;
	
	public BindingEM(SEMConfig s, ExperimentManager eMan, BindingManager bMan, HashMap<ExperimentCondition, BackgroundCollection> condBacks, int numPotReg) {
		semconfig = s;
		manager = eMan;
		bindingManager = bMan;
		conditionBackgrounds = condBacks;
		numConditions = manager.getNumConditions();
		numPotentialRegions = (double)numPotReg;
		
		//Penalty probability applied on independent/shared conditions
		//TODO: how to adjust it when using it on nucleosomes?
		double N = numPotentialRegions;
		double S = N * semconfig.getProbSharedBinding();
		double L = (double)semconfig.getGenome().getGenomeLength();
		probAgivenB = Math.log(semconfig.getProbSharedBinding())/Math.log(2);
		probAgivenNOTB = Math.log((N-S)/(L-N))/Math.log(2);
	}
	
	public List<List<BindingComponent>> train(List<List<StrandedPair>> signals,
												Region w,
												List<NoiseComponent> noise,
												List<List<BindingComponent>> comps,
												int numComp,
												double[][] atacPrior,
												int trainingRound,
												Timer timer
												) throws Exception {
		components = comps;
		this.noise = noise;
		numComponents = numComp;
		this.atacPrior = atacPrior;
		this.trainingRound = trainingRound;
		this.plotSubRegion = plotSubRegion;
		this.timer = timer;
//		numBindingSubtypes = new HashMap<ExperimentCondition, Integer>();
		// Matrix initializations
		hitCounts = new double[numConditions][];		// Hit weights
		hitPos = new int[numConditions][];				// Hit positions
		hitSize = new int[numConditions][];				// Hit fragment sizes
		hitNum = new int[numConditions];				// Number of hits in each condition
		repIndices = new int[numConditions][];			// Index of replicate for the hit
		hAll = new double[numConditions][][];			// H functions for all positions in the current window
		h = new double[numConditions][][];				// H function (binding component probability per read)
		n = new double[numConditions][];				// N function (noise component probability per read)
		rBind = new double[numConditions][][];			// Binding component responsibilities
		rNoise = new double[numConditions][];			// Noise component responsibilities
		sumResp = new double[numConditions][];			// Sum of responsibilities each bindingComponent
		pi = new double[numConditions][numComponents];	// pi: emission probabilities for binding components
		piNoise = new double[numConditions];			// pi: emission probabilities for noise components (fixed)
		alphaMax = new double[numConditions];			// Maximum alpha
		mu = new int[numConditions][numComponents];		// mu: positions of the binding components
		fuzz = new double[numConditions][numComponents];	// &fuzz: fuzziness of the binding components
		tau = new double[numConditions][numComponents][];				// &tau: fragment size subtype probabilities (indexed by subtype index)
		plotEM = (plotSubRegion!=null && plotSubRegion.overlaps(w));
		// Monitor state convergence using the following last variables
		lastRBind = new double[numConditions][][];
		lastPi = new double[numConditions][numComponents];
		lastMu = new int[numConditions][numComponents];
		lastFuzz = new double[numConditions][numComponents];
		lastTau = new double[numConditions][numComponents][];
		
		// Initializing data structures
		for(ExperimentCondition cond: manager.getConditions()) {
			int c = cond.getIndex();
			
			// Set maximum alphas (TODO: Now we set 0 for SEM because we haven't found a good way to determine ahpha)
			alphaMax[c] = semconfig.getFixedAlpha()>0 ? semconfig.getFixedAlpha() :
					semconfig.getAlphaScalingFactor() * (double)conditionBackgrounds.get(cond).getMaxThreshold('.');
			
			//monitor: count time
			timer.start();
			
			// sort all read pairs by midpoint/size per replicate
			// TODO: move sort step to HitCache to reduce computation
			Comparator signalCompare = new Comparator<StrandedPair>() {
				@Override
				public int compare(StrandedPair s1, StrandedPair s2) {
					if(s1==null && s2==null)
						return 0;
					if(s1==null)
						return -1;
					if(s2==null)
						return 1;
					if(s1.getMidpoint().getLocation()<s2.getMidpoint().getLocation())
						return -1;
					if(s1.getMidpoint().getLocation()==s2.getMidpoint().getLocation()){
						if(s1.getFragmentSize()<s2.getFragmentSize())
							return -1;
						if(s1.getFragmentSize()==s2.getFragmentSize())
							return 0;
						return 1;
					}
					if(s1.getMidpoint().getLocation()>s2.getMidpoint().getLocation())
						return 1;
					return 0;
				}
			};	
			
			// Load Read pairs (merge from all replicates), sort it
			List<StrandedPair> pairs = new ArrayList<StrandedPair>();
			for(ControlledExperiment rep: cond.getReplicates()) {
				List<StrandedPair> repSignals = signals.get(rep.getIndex());
				Collections.sort(repSignals, signalCompare);
				pairs.addAll(repSignals);
//				pairs.addAll(signals.get(rep.getIndex()));
			}
			int numPairs = pairs.size();
			hitNum[c] = numPairs;

			// Load replicate index for each read
			
			repIndices[c] = new int[numPairs];
			int y=0, z=0;
			for(ControlledExperiment rep: cond.getReplicates()) {
				z=0;
				while(z<signals.get(rep.getIndex()).size()) {
					repIndices[c][y] = rep.getIndex();
					z++; y++;
				}
			}
			
			// Load read pair info
			// TODO: It looks like I can merge pairs with the same fragment size and midpoint
			double[] countc = new double[numPairs];
			int[] posc = new int[numPairs];
			int[] sizec = new int[numPairs];
			for(int i=0; i<numPairs; i++) {
				posc[i] = pairs.get(i).getMidpoint().getLocation();
				sizec[i] = pairs.get(i).getFragmentSize();
				countc[i] = pairs.get(i).getWeight();
			}
			hitPos[c] = posc;
			hitCounts[c] = countc;
			hitSize[c] = sizec;
			
			//test: collapse duplicate fragments when fragment density is higher than 10/bp
			if(numPairs/w.getWidth()>=10) {
				System.out.println("Collapsing high density region");
				int uniquePos=1;
				for(int k=0; k < hitPos[c].length-1; k++) {
					if(hitSize[c][k+1] != hitSize[c][k] ||
							hitPos[c][k+1] != hitPos[c][k] ||
							repIndices[c][k+1] != repIndices[c][k]) {
						uniquePos++;
					}
				}
				int[] tmphitPos = new int[uniquePos];
				int[] tmphitSize = new int[uniquePos];
				double[] tmphitCounts = new double[uniquePos];
				int[] tmprepIndices = new int[uniquePos];
				int x=0;
				tmphitPos[x] = hitPos[c][0];
				tmphitSize[x] = hitSize[c][0];
				tmphitCounts[c] += hitCounts[c][0];
				tmprepIndices[c] = repIndices[c][0];
				for(int k=1; k<hitPos[c].length; k++) {
					if(hitSize[c][k-1] != hitSize[c][k] ||
							hitPos[c][k-1] != hitPos[c][k] ||
							repIndices[c][k-1] != repIndices[c][k]) {
						x++;
					}
					tmphitPos[x] = hitPos[c][k];
					tmphitSize[x] = hitSize[c][k];
					tmphitCounts[x] += hitCounts[c][k];
					tmprepIndices[x] = repIndices[c][k];
				}
				hitPos[c] = tmphitPos;
				hitSize[c] = tmphitSize;
				hitCounts[c] = tmphitCounts;
				repIndices[c] = tmprepIndices;
				
				numPairs = hitPos[c].length;
				hitNum[c] = numPairs;
			}
			
			// Load pi for binding components
			for(int j=0;j<numComp; j++) {
				pi[c][j] = components.get(c).get(j).getPi();
			}
			
			// Load pi for noise components
			piNoise[c] = noise.get(c).getPi();
			
			// Load binding components positions, fuzziness, tau
//			numBindingSubtypes.put(cond, bindingManager.getBindingSubtypes(cond).size());
			for(int j=0; j<numComp; j++) {
				mu[c][j] = components.get(c).get(j).getPosition();
				fuzz[c][j] = components.get(c).get(j).getFuzziness();
				tau[c][j] = components.get(c).get(j).getTau();
//				lastTau[c][j] = new double[fragSizePDF.keySet().size()];
			}
			
			// Initialize H function for all positions in the current region
			// Due to different binding components have different fuzziness, we cannot precompute H function here for SEM
			
			// Initialize responsibility functions
			rBind[c] = new double[numComp][numPairs];
			rNoise[c] = new double[numPairs];
			lastRBind[c] = new double[numComp][numPairs];
		}
		// End of data structure initialization
		
		//monitor: count time
		timer.end("initialize");

		/////////
		// Run EM steps
		/////////
		EM_MAP(w);
		
		//////////
        // re-assign EM result back to component objects
        //////////
		List<List<BindingComponent>> activeComponents = new ArrayList<List<BindingComponent>>();
		for(ExperimentCondition cond: manager.getConditions()) {
			// Binding Components
			List<BindingComponent> currActiveComps = new ArrayList<BindingComponent>();
			int c = cond.getIndex();
			for(int j=0; j<numComp; j++) {
				BindingComponent comp = components.get(c).get(j);
				comp.setPi(pi[c][j]);
				comp.setPosition(mu[c][j]);
				comp.setFuzziness(fuzz[c][j]);
				comp.setTau(tau[c][j]);
				double sum_resp = 0.0;	
				for(int i=0;i<hitNum[c];i++){
                    sum_resp += hitCounts[c][i]*rBind[c][j][i];
                }
	            comp.setSumResponsibility(sum_resp);
				if(pi[c][j]>0.0) {
					currActiveComps.add(comp);
				}
			}
			activeComponents.add(currActiveComps);
			// Noise components
			double noise_resp = 0.0;
			for(int i=0; i<hitNum[c]; i++)
				noise_resp += hitCounts[c][i]*rNoise[c][i];
			noise.get(c).setSumResponsibility(noise_resp);
		}
		return activeComponents;
	}
	
	/**
     * Core EM iterations with sparse prior (component elimination) & multi-condition positional priors.
     * Assumes H function, pi, and responsibilities have all been initialized
     */
	private void EM_MAP(Region currRegion) throws Exception {
		
		int numComp = numComponents;
		double [][] totalResp = new double[numConditions][];
		int regStart = currRegion.getStart();
		
		// Variables for tracking mu maximization. Defined early to avoid memory assignment during main EM loop.
		double[][][] muSums = new double[numConditions][numComp][];		//Results of mu maximization summations for individual components across genome
        int[][] muSumStarts = new int[numConditions][numComp]; //Start positions of muSum arrays (start of maximization window).
        int[][] muSumWidths = new int[numConditions][numComp]; //Effective widths of muSum arrays (width of maximization window).
        int[][] muMax = new int[numConditions][numComp]; //Positions of maxima in mu maximization summations
        double[][] fuzzMax = new double[numConditions][numComp]; //Fuzziness of maxima in fuzziness maximization summations
        double[][] muSum = new double[numConditions][numComp]; //Sum of hitPos*rBind*hitCounts for each binding component
        double[][] fuzzSum = new double[numConditions][numComp]; //Sum of Variance*rBind*hitCounts for each binding component
        int[] muJoinClosestComps = new int[numConditions]; //Indices of nearest components in other conditions
        boolean[] muJoinSharedBetter = new boolean[numConditions]; //Indicator that sharing components across conditions is better than not
        int[][] newMu = new int[numConditions][numComponents];// mu update
        double[][] newFuzz = new double[numConditions][numComponents]; // fuzziness update 
        
        // Initialize responsibilities
        for(int c=0; c<numConditions; c++) {
        	int numPairs = hitNum[c];
        	totalResp[c] = new double[numPairs];
        	for(int i=0; i<numPairs; i++)
        		totalResp[c][i] = 0;
        }
        
        // Alpha is annealed in. Alpha=0 during ML steps
        double[] currAlpha = new double[numConditions];
        for(int c=0; c<numConditions; c++)
        	currAlpha[c] = 0;
        
    	//////////////////////////////////////////////////////////////////////////////////////////
        //Run EM while not converged
        // Note: iterations during which we eliminate a binding component don't count towards "t"
    	//////////////////////////////////////////////////////////////////////////////////////////
        int t=0, iter=0;
        
        //monitor code:
        long starttime = System.currentTimeMillis();     
        
        while(t<semconfig.MAX_EM_ITER) {
        	
        	//monitor
//        	System.out.println("In EM Region: "+currRegion.getChrom()+":"+currRegion.getStart()+"-"+currRegion.getEnd());
//        	System.out.println("In EM reads number: "+hitPos[0].length);
//        	System.out.println("\tE-step");
        	
    		////////
    		//E-step
    		////////
        	
        	//I want to get the range of fragment midpoint coordinates that will be influenced by each BindingComponent
        	//It's a trick to improve the performance because when encountering quite long region, calculate each fragment
        	//for each BindingComponent will waste a lot of time.
        	
        	//monitor
//        	System.out.println("\t\tMarking range");
        	
			//monitor: count time
        	timer.start();
        	
        	// Marking range for each bindingComponent
        	List<Map<Integer, Pair<Integer, Integer>>> pairIndexAroundMu = new ArrayList<Map<Integer, Pair<Integer, Integer>>>(); 
        	for(int c=0; c<numConditions; c++) {

        		pairIndexAroundMu.add(new HashMap<Integer, Pair<Integer, Integer>>());
        		int numPairs = hitNum[c];
        		for(int j=0; j<numComp; j++) {
        			int start=0; int end=numPairs-1;
        			boolean bigger=false; boolean smaller=false;
        			for(int i=0; i<numPairs; i++) {
        				if(!smaller) {
        					if((mu[c][j]-hitPos[c][i])<=BindingModel.getMaxInfluenceRange()/2) {
        						start = i;
        						smaller = true;
        					}
        				}
        				if(smaller && !bigger) {
        					if((hitPos[c][i]- mu[c][j])>BindingModel.getMaxInfluenceRange()/2) {
        						end = i-1;
        						bigger = true;
        					}
        				}
        			}
        			if(hitPos[c][start] <=(mu[c][j]+BindingModel.getMaxInfluenceRange()/2)) {
        				pairIndexAroundMu.get(c).put(j, new Pair<Integer, Integer>(start, end));
        			} else {
        				pairIndexAroundMu.get(c).put(j, new Pair<Integer, Integer>(-1, -1));
        			}
        		}
        	}
        	
			//monitor: count time
        	timer.end("mark");
        	timer.start();
        	
        	for(int c=0; c<numConditions; c++) {
        		int numPairs = hitNum[c];
        		ExperimentCondition cond = manager.getIndexedCondition(c);
        		
    			// Load binding fragment size PDF cache (indexed by type index)
    			Map<Integer, List<Double>> fragSizePDF = bindingManager.getCachePDF(cond);
    			
            	//monitor
//            	System.out.println("\t\tCompute H and N function");
    			
    			// Compute H and N function
    			double[][] hc = new double[numComp][numPairs];
    			double[] nc = new double[numPairs];
//            	//monitor
//            	System.out.println("\t\t\tIn loop");
//    			for(int i=0; i<numPairs; i++) {
//    				for(int j=0; j<numComp; j++) {
////    					System.out.println("\t\t\tcompute fuzzProb");
//    					double fuzzProb = BindingModel.probability(fuzz[c][j], mu[c][j]-hitPos[c][i]);
//    					double fragSizeProb = 0;
////    					System.out.println("\t\t\tcompute fragSizeProb");
//    					for(int index: fragSizePDF.keySet()) {
//    						fragSizeProb += tau[c][j][index] * fragSizePDF.get(index).get(hitSize[c][i]);
//    					}
//    					hc[j][i] = fuzzProb * fragSizeProb;
//    				}
//    				nc[i] = noise.get(c).score(hitPos[c][i], hitSize[c][i], repIndices[c][i]);
//    			}
    			
            	//test: reduce computation
    			for(int j=0; j<numComp; j++) {
    				if(pi[c][j]>0 && pairIndexAroundMu.get(c).get(j).car()!=-1) {
    					for(int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++) {
//    						System.out.println("\t\t\tcompute fuzzProb");
        					double fuzzProb = BindingModel.probability(fuzz[c][j], mu[c][j]-hitPos[c][i]);
        					double fragSizeProb = 0;
//        					System.out.println("\t\t\tcompute fragSizeProb");
        					for(int index: fragSizePDF.keySet()) {
        						fragSizeProb += tau[c][j][index] * fragSizePDF.get(index).get(hitSize[c][i]);
        					}
        					hc[j][i] = fuzzProb * fragSizeProb;
        					
        					if(Double.isNaN(hc[j][i])) {
        						System.out.println("Nan happens in H function: "+"\tfuzzProb: "+fuzzProb+"\tfragSizeProb: "+fragSizeProb+"fuzz: "+fuzz[c][j]+"\tsize: "+hitSize[c][i]+"\tdistance: "+(mu[c][j]-hitPos[c][i]));
        					}
    					}
    				}
    			}
    			for(int i=0; i<numPairs; i++) {
    				nc[i] = noise.get(c).score(hitPos[c][i], hitSize[c][i], repIndices[c][i]);
    			}
    			
    			h[c] = hc;
    			n[c] = nc;
    			
    			//monitor: count time
    			timer.end("HN");
    			timer.start();
        		
            	//monitor
//            	System.out.println("\t\tCompute responsibilities");
    			
        		// Compute responsibilities
//        		for(int i=0; i<numPairs; i++) {
//        			totalResp[c][i] = 0;
//        			for(int j=0; j<numComp; j++) {
//        				if(pi[c][j]>0) {
//        					rBind[c][j][i] = h[c][j][i]*pi[c][j];
//        					totalResp[c][i] += rBind[c][j][i];
//        				}
//        			}
//        			rNoise[c][i] = n[c][i] * piNoise[c];
//        			totalResp[c][i] += rNoise[c][i];
//        		}
            	
        		// test: reduce computation
    			totalResp[c] = new double[numPairs];
    			rBind[c] = new double[numComp][numPairs];
    			rNoise[c] = new double[numPairs];
        		for(int j=0; j<numComp; j++) {
        			if(pi[c][j]>0 && pairIndexAroundMu.get(c).get(j).car()!=-1) {
        				for(int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++) {
        					rBind[c][j][i] = h[c][j][i]*pi[c][j];
        					totalResp[c][i] += rBind[c][j][i];
        			}
        			}
        		}
        		for(int i=0; i<numPairs; i++) {
        			rNoise[c][i] = n[c][i] * piNoise[c];
        			totalResp[c][i] += rNoise[c][i];
        			//test: move normalize noise here
        			rNoise[c][i]/=totalResp[c][i];
        			
        			//monitor
        			if(Double.isNaN(rNoise[c][i])) {
        				System.out.println("NaN happens in normalize: "+"\tBefore normalize: "+n[c][i]*piNoise[c]+"\ttotalResp: "+totalResp[c][i]);
        			}
        		}
        		
    			//monitor: count time
        		timer.end("cResp");
        		timer.start();
        		
            	//monitor
//            	System.out.println("\t\tNormalize responsibilities");
        		
        		// Normalize responsibilities
//        		for(int i=0; i<numPairs; i++) {
//        			for(int j=0; j<numComp; j++) {
//        				if(pi[c][j]>0) {
//        					rBind[c][j][i] /= totalResp[c][i];
//        				}
//        			}
//        			rNoise[c][i]/=totalResp[c][i];
//        		}
        		
        		//test: reduce computation
        		for(int j=0; j<numComp; j++) {
        			if(pi[c][j]>0 && pairIndexAroundMu.get(c).get(j).car()!=-1) 
        				for(int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++) {
        						rBind[c][j][i] /= totalResp[c][i];
        				}
        		}
        	}
        	
			//monitor: count time
        	timer.end("nResp");
        	timer.start();
        	
        	//Compute sum of responsibility
        	for(int c=0; c<numConditions; c++) {
        		sumResp[c] = new double[numComp];
        		for(int j=0; j<numComp; j++) {
        			if(pi[c][j]>0 && pairIndexAroundMu.get(c).get(j).car()!=-1) {
        				for(int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++) {
        					sumResp[c][j] += rBind[c][j][i] * hitCounts[c][i];
        				}
        			}
        		}
        	}
        	
        	//monitor
//        	System.out.println("\tM-step");
        	
    		/////////////////////
    		//M-step: maximize mu (positions), fuzz (fuzziness), tau (fragment size subtype), pi (strength)
    		/////////////////////
        	if(numConditions>1 && t==semconfig.ALPHA_ANNEALING_ITER)
        		for(int c=0; c<numConditions; c++)
        			for(int j=0; j<numComp; j++)
        				if(pi[c][j]>0)
        					muSums[c][j] = new double[semconfig.EM_MU_UPDATE_WIN*2];
        	
        	//monitor
//        	System.out.println("\t\tMaximize mu");
        	
        	// Part1.1: Maximize mu (Because we use Gaussian distribution to describe tag distribution around dyad, we can get dyad directly by the proportional sum of tag position)
        	for(int c=0; c<numConditions; c++) {
        		int numPairs = hitNum[c];
        		for(int j=0; j<numComp; j++) {
        			if(pi[c][j]>0 && pairIndexAroundMu.get(c).get(j).car()!=-1) {
        				// Get the limit of new dyad locations
        				int start = Math.max(mu[c][j]-semconfig.EM_MU_UPDATE_WIN, regStart);
        				int end = Math.min(currRegion.getEnd(), mu[c][j]+semconfig.EM_MU_UPDATE_WIN);
        				// Assign variables for pairwise comparison
        				if(numConditions>1 && t>semconfig.ALPHA_ANNEALING_ITER) {
        					muSumStarts[c][j] = start; muSumWidths[c][j] = end-start;
        				}
        				// Score the current window
        				double currScore=0, maxScore=-Double.MAX_VALUE;
        				int maxPos = 0;
        				
        				//test here
//        				for(int i=0; i<numPairs; i++) {
        				for(int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++) {
        					muSum[c][j] += rBind[c][j][i] * hitCounts[c][i] * hitPos[c][i];
        				}

        				// Does maximize position exceed update window?
        				double mu_max = muSum[c][j]/sumResp[c][j];
        				if(mu_max > end) {
        					muMax[c][j] = end;
//        					mu[c][j] = end;
        				} else if(mu_max < start) {
        					muMax[c][j] = start;
//        					mu[c][j] = start;
        				} else {
        					muMax[c][j] = (int)Math.rint(mu_max);
//        					mu[c][j] = (int)Math.rint(sumTag/sumWeight);
        				}
        				
//        				for(int x=start; x<end; x++) {
//        					currScore=0;
//        					for(int i=0; i<numPairs; i++) {
//    							int dist = hitPos[c][i] - x;
//    							currScore += (rBind[c][j][i]*hitCounts[c][i]) * BindingModel.logProbability(fuzz[c][j], dist);
//        					}
//        					if(atacPrior!=null && semconfig.useAtacPrior())
//        						currScore += atacPrior[c][x-regStart];
//        					
//        					if(numConditions>1 && t>semconfig.ALPHA_ANNEALING_ITER) // Save the score
//        						muSums[c][j][x-start] = currScore;
//        					
//        					if(currScore>maxScore) {
//        						maxPos = x;
//        						maxScore = currScore;
//        					}
//        					
//        					if(j==2 && (x==216 || x==226)) {
//        						System.out.println(x + " prob: "+currScore);
//        					}
//        				}
//        				
//        				muSumMax[c][j] = maxPos;
//        				mu[c][j] = maxPos;
        			}
        		}

        	}
        	
			//monitor: count time
        	timer.end("mu");

        	timer.start();
        	
        	//monitor
//        	System.out.println("\t\tMaximize fuzziness");
        	
        	// Part2: Maximize fuzziness
        	if(t/semconfig.FUZZINESS_ANNEALING_ITER==0) {
        		for(int c=0; c<numConditions; c++) {
        			for(int j=0; j<numComp; j++) {if(pi[c][j]>0 && pairIndexAroundMu.get(c).get(j).car()!=-1) {
	        			for(int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++)
	        				fuzzSum[c][j] += (rBind[c][j][i] * hitCounts[c][i]) * Math.pow(hitPos[c][i]-mu[c][j], 2);
	        			fuzzMax[c][j] = fuzzSum[c][j]/sumResp[c][j];
//	        			fuzz[c][j] = fuzzSum[c][j]/sumResp[c][j];
        			}
        			}
        		}
        	}
        	
			//monitor: count time
        	timer.end("fuzz");
        	
        	//Insert pairwise bindingcomponents comparison here
        	//Due to I haven't found a way to compare fragment size distribution across different experiment condition
        	//I only consider position and fuzziness here.
        	
        	//evaluate whether joining nearby components across conditions is more favorable
        	
        	//First calculate independent log likelihood score for each binding component
        	double[][] indepScorePerBC = new double[numConditions][numComp];
        	for(int c=0; c<numConditions; c++) {
        		for(int j=0; j<numComp; j++) {
        			if(pi[c][j]>0 && pairIndexAroundMu.get(c).get(j).car()!=-1) {
        				for (int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++)
        					indepScorePerBC[c][j] += rBind[c][j][i] * hitCounts[c][i] * BindingModel.logProbability(fuzzMax[c][j], muMax[c][j]-hitPos[c][i]);
        			}
        		}
        	}
        	
        	for(int c=0; c<numConditions; c++) {
        		for(int j=0; j<numComp; j++) {
        			if(pi[c][j]>0 && pairIndexAroundMu.get(c).get(j).car()!=-1) {
        				if(numConditions>1 && t>semconfig.ALPHA_ANNEALING_ITER) {
        				//a: find the closest components to j in each condition
        				int closestComp=-1; int closestDist=Integer.MAX_VALUE;
        				for(int d=0; d<numConditions; d++) {
        					if(d!=c) {
        						closestComp=-1; closestDist=Integer.MAX_VALUE;
        						for(int k=0; k<numComp; k++) {
        							if(pi[d][k]>0 && pairIndexAroundMu.get(d).get(k).car()!=-1) {
        								int dist = Math.abs(mu[c][j]-mu[d][k]);
        								if(dist<closestDist && dist<semconfig.EM_MU_UPDATE_WIN)
        									closestDist = dist; closestComp = k;
        							}
        						}
        						muJoinClosestComps[d] = closestComp;
        					}
        				}
        				
        				//b: evaluate each pair of conditions, asking if a shared event involving j and its closest component would be better than independent events
        				int numSharedBetter = 0;
        				int maxMuStart = muSumStarts[c][j];
        				int minMuEnd = muSumStarts[c][j] + muSumWidths[c][j];
        				for(int d=0; d<numConditions; d++) {
        					if(d!=c) {
        						int k = muJoinClosestComps[d];
        						if(k==-1)
        							muJoinSharedBetter[d] = false;
        						else {
        							//Case 1: two independent components
        							double indepScore = indepScorePerBC[c][j] + probAgivenNOTB +
        												indepScorePerBC[d][k] + probAgivenNOTB;
        							//Case 2: single shared components
        							//Case 2.1: first get maximum shared position
        							double maxSharedScore = 0;
        							int maxSharedPos = (int)Math.rint((muSum[c][j] + muSum[d][k])/(sumResp[c][j] + sumResp[d][k]));
        							double maxSharedFuzz = 0;
        							maxSharedPos = maxSharedPos>Math.max(maxMuStart, muSumStarts[d][k]) ? maxSharedPos : Math.max(maxMuStart, muSumStarts[d][k]);
        							maxSharedPos = maxSharedPos<Math.min(minMuEnd, muSumStarts[d][k]+muSumWidths[d][k]) ? maxSharedPos : Math.min(minMuEnd, muSumStarts[d][k]+muSumWidths[d][k]);
        							//Case 2.2: then get maximum shared fuzziness
        							for (int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++) {
        								maxSharedFuzz += (rBind[c][j][i] * hitCounts[c][i]) * Math.pow(hitPos[c][i]-maxSharedPos, 2);
        							}
        							for (int i=pairIndexAroundMu.get(d).get(k).car(); i<=pairIndexAroundMu.get(d).get(k).cdr(); i++) {
        								maxSharedFuzz += (rBind[d][k][i] * hitCounts[d][i]) * Math.pow(hitPos[d][i]-maxSharedPos, 2);
        							}
        							maxSharedFuzz /= (sumResp[c][j] + sumResp[d][k]);
        							//Case 2.3: then get log likelihood score
        							for (int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++) {
        								maxSharedScore += rBind[c][j][i] * hitCounts[c][i] * BindingModel.logProbability(maxSharedFuzz, maxSharedPos-hitPos[c][i]);
        							}
        							for (int i=pairIndexAroundMu.get(d).get(k).car(); i<=pairIndexAroundMu.get(d).get(k).cdr(); i++) {
        								maxSharedScore += rBind[d][k][i] * hitCounts[d][i] * BindingModel.logProbability(maxSharedFuzz, maxSharedPos-hitPos[d][i]);
        							}
        							maxSharedScore += 2*probAgivenB;
        							//Is shared better? get overlapping region if true
        							muJoinSharedBetter[d] = maxSharedScore > indepScore ? true : false;
        							maxMuStart = muJoinSharedBetter[d] ? Math.max(maxMuStart, muSumStarts[d][k]) : maxMuStart; 
        							minMuEnd = muJoinSharedBetter[d] ? Math.min(minMuEnd, muSumStarts[d][k]+muSumWidths[d][k]) : minMuEnd;
        							if(muJoinSharedBetter[d])
        								numSharedBetter++;
        							
        							//update mu (Shortcut for numConditions==2)
        							if(numConditions==2) {
        								if(muJoinSharedBetter[d]) {
        									newMu[c][j] = maxSharedPos;
        									newFuzz[c][j] = maxSharedFuzz;
        								}
        								else {
        									newMu[c][j] = muMax[c][j];
        									newFuzz[c][j] = fuzzMax[c][j];
        								}
        							}
        						}
        					}
        				}
        				
        				//c: for all conditions that passed the pairwise test, evaluate if a single shared event is better than all independent
        				if(numConditions>2) { //Shortcut for 2 conditions above
        					//Case 1: sum of all independent components
        					double allIndepScore = indepScorePerBC[c][j] + probAgivenNOTB;
        					for(int d=0; d<numConditions; d++) {
        						if(d!=c) {
        							int k = muJoinClosestComps[d];
        							if(k!=-1) {
        								allIndepScore += indepScorePerBC[d][k] + probAgivenNOTB;
        							}
        						}
        					}
        					//Case 2: sum of shared component and non-shared
        					//Case 2.1: first get the maximum shared position for shared better components
        					int maxSomeSharedPos = 0;
        					double someSharedMuSum = muSum[c][j];
        					double someSharedSumResp = sumResp[c][j];
        					for(int d=0; d<numConditions; d++) {
        						if(d!=c) {
        							int k = muJoinClosestComps[d];
        							if(k!=-1 && muJoinSharedBetter[d]) {
        								someSharedMuSum += muSum[d][k];
        								someSharedSumResp += sumResp[d][k];
        							}
        						}
        					}
        					maxSomeSharedPos = (int)Math.rint(someSharedMuSum/someSharedSumResp);
        					maxSomeSharedPos = maxSomeSharedPos<maxMuStart ? maxMuStart : maxSomeSharedPos;
        					maxSomeSharedPos = maxSomeSharedPos>minMuEnd   ? minMuEnd   : maxSomeSharedPos;
        					//Case 2.2: then get maximum shared fuzziness
        					double maxSomeSharedFuzz = 0;
        					for (int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++) {
								maxSomeSharedFuzz += (rBind[c][j][i] * hitCounts[c][i]) * Math.pow(hitPos[c][i]-maxSomeSharedPos, 2);
							}
        					for(int d=0; d<numConditions; d++) {
        						if(d!=c) {
        							int k = muJoinClosestComps[d];
        							if(k!=-1 && muJoinSharedBetter[d]) {
        								for (int i=pairIndexAroundMu.get(d).get(k).car(); i<=pairIndexAroundMu.get(d).get(k).cdr(); i++) {
        									maxSomeSharedFuzz += (rBind[d][k][i] * hitCounts[d][i]) * Math.pow(hitPos[d][i]-maxSomeSharedPos, 2);
        								}
        							}
        						}
        					}
        					maxSomeSharedFuzz /= someSharedSumResp;							
        					//Case 2.2: then compute log likelihood
        					double maxSomeSharedScore = 0;
        					for(int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++) {
								maxSomeSharedScore += rBind[c][j][i] * hitCounts[c][i] * BindingModel.logProbability(maxSomeSharedFuzz, maxSomeSharedPos-hitPos[c][i]);
        					}
							maxSomeSharedScore += probAgivenB;
        					for(int d=0; d<numConditions; d++) {
        						if(d!=c) {
        							int k = muJoinClosestComps[d];
        							if(k!=-1) {
        								if(muJoinSharedBetter[d]) {
        									for(int i=pairIndexAroundMu.get(d).get(k).car(); i<=pairIndexAroundMu.get(d).get(k).cdr(); i++)
        										maxSomeSharedScore += rBind[d][k][i] * hitCounts[d][i] * BindingModel.logProbability(maxSomeSharedFuzz, maxSomeSharedPos-hitPos[d][i]);
        									maxSomeSharedScore += probAgivenB;
        								} else {
        									maxSomeSharedScore += indepScorePerBC[d][k] += probAgivenNOTB;
        								}
        							}
        						}
        					}
        					//Case 3: single shared position, regardless of what happened in the pairwise tests
        					//shortcut if all conditions are shared better
        					double maxAllSharedScore = 0; int maxAllSharedPos = 0;
        					if(numSharedBetter==numConditions-1) {
        						maxAllSharedScore = maxSomeSharedScore;
        						maxAllSharedPos = maxSomeSharedPos;
        					} else {
        						//Case 3.1: first get the maximum shared position for all nearby components
        						double allSharedMuSum = muSum[c][j];
            					double allSharedSumResp = sumResp[c][j];
            					for(int d=0; d<numConditions; d++) {
            						if(d!=c) {
            							int k = muJoinClosestComps[d];
            							if(k!=-1) {
            								allSharedMuSum += muSum[d][k];
            								allSharedSumResp += sumResp[d][k];
            							}
            						}
            					}
            					maxAllSharedPos = (int)Math.rint(allSharedMuSum/allSharedSumResp);
            					maxAllSharedPos = maxAllSharedPos<maxMuStart ? maxMuStart : maxAllSharedPos;
            					maxAllSharedPos = maxAllSharedPos>minMuEnd   ? minMuEnd   : maxAllSharedPos;
            					//Case 3.2: then get the maximum shared fuzziness 
            					double maxAllSharedFuzz = 0;
            					for (int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++) {
    								maxAllSharedFuzz += (rBind[c][j][i] * hitCounts[c][i]) * Math.pow(hitPos[c][i]-maxAllSharedPos, 2);
    							}
            					for(int d=0; d<numConditions; d++) {
            						if(d!=c) {
            							int k = muJoinClosestComps[d];
            							if(k!=-1) {
            								for (int i=pairIndexAroundMu.get(d).get(k).car(); i<=pairIndexAroundMu.get(d).get(k).cdr(); i++) {
            									maxSomeSharedFuzz += (rBind[d][k][i] * hitCounts[d][i]) * Math.pow(hitPos[d][i]-maxAllSharedPos, 2);
            								}
            							}
            						}
            					}
            					maxAllSharedFuzz /= allSharedSumResp;	
            					//Case 3.2: then log likelihood
            					for(int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++) {
    								maxAllSharedScore += rBind[c][j][i] * hitCounts[c][i] * BindingModel.logProbability(maxAllSharedFuzz, maxAllSharedPos-hitPos[c][i]);
            					}
    							maxAllSharedScore += probAgivenB;
            					for(int d=0; d<numConditions; d++) {
            						if(d!=c) {
            							int k = muJoinClosestComps[d];
            							if(k!=-1) {
            								for(int i=pairIndexAroundMu.get(d).get(k).car(); i<=pairIndexAroundMu.get(d).get(k).cdr(); i++)
            									maxAllSharedScore += rBind[d][k][i] * hitCounts[d][i] * BindingModel.logProbability(maxAllSharedFuzz, maxAllSharedPos-hitPos[d][i]);
            								maxAllSharedScore += probAgivenB;
            							}
            						}
            					}
            					
            					//e: update mu
    	    					if(maxAllSharedScore >=allIndepScore && maxAllSharedScore >=maxSomeSharedScore) {
    	    						newMu[c][j] = maxAllSharedPos;
    	    						newFuzz[c][j] = maxAllSharedFuzz;
    	    					} else if(maxSomeSharedScore >=allIndepScore) {
    	    						newMu[c][j] = maxSomeSharedPos;
    	    						newFuzz[c][j] = maxSomeSharedFuzz;
    	    					} else {
    	    						newMu[c][j] = muMax[c][j];
    	    						newFuzz[c][j] = fuzzMax[c][j];
    	    					}
        					}
        				}
        			}} else {
        				//Ignore other conditions in first phase of training
        				newMu[c][j] = muMax[c][j];
        				newFuzz[c][j] = fuzzMax[c][j];
        			}
        		}
        	}
        	
        	
        	timer.start();
        	
        	//monitor
//        	System.out.println("\t\tResolve duplicate");
        	
        	// Part1.2: Resolve duplicate binding components (share the same position) -> combine & delete one copy
        	for(int c=0; c<numConditions; c++) {
        		int numPairs = hitNum[c];
        		HashMap<Integer, Integer> pos2index = new HashMap<Integer, Integer>();	// Position to array index map
        		for(int j=0; j<numComp; j++) {
        			if(pi[c][j]>0) {
        				if(pos2index.containsKey(mu[c][j])) {
        					int orig = pos2index.get(mu[c][j]);
        					// Combine
        					pi[c][orig] += pi[c][j];
        					for(int i=0; i<numPairs; i++)
        						rBind[c][orig][i] += rBind[c][j][i];
        					// Delete
        					pi[c][j] = 0.0;
        					for(int i=0; i<numPairs; i++)
        						rBind[c][j][i] = 0;
        				} else {
        					pos2index.put(mu[c][j], j);
        				}
        			}
        		}
        	}

			//monitor: count time
        	timer.end("resolve");
        	timer.start();

        	//monitor
//        	System.out.println("\t\tMaximize tau");
        	
        	// Part3: Maximize tau
        	if(t/semconfig.TAU_ANNEALING_ITER==0) {
        	for(ExperimentCondition cond: manager.getConditions()) {
        		int c = cond.getIndex();
        		int numPairs = hitNum[c];
        		int numSubtypes = bindingManager.getBindingSubtypes(cond).size();
        		for(int j=0; j<numComp; j++) {if(pi[c][j]>0) {
        			tau[c][j] = new double[numSubtypes];
        			double probSum = 0;
        			// Compute tau probability of each fragment size subtype for each component
        			for(BindingSubtype b: bindingManager.getBindingSubtypes(cond)) {
//        				for(int i=0; i<numPairs; i++) {
        				if(pairIndexAroundMu.get(c).get(j).car()!=-1)
        					for(int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++) {
        						tau[c][j][b.getIndex()] += (rBind[c][j][i]* hitCounts[c][i]) * b.probability(hitSize[c][i]);
        					}
        				probSum += tau[c][j][b.getIndex()];
        			}
        			// Normalize tau for each component
        			for(BindingSubtype b: bindingManager.getBindingSubtypes(cond)) 
        				tau[c][j][b.getIndex()] /= probSum;
        			// Eliminate subtype with not enough ratio of tags (Assumption2: each nucleosome should not be associated with too many subtypes)
        			for(BindingSubtype b: bindingManager.getBindingSubtypes(cond)) {
        				if(tau[c][j][b.getIndex()] < semconfig.SPARSE_PRIOR_SUBTYPE)
        					tau[c][j][b.getIndex()] = 0;
        			}
        			// Re-normalize tau for each component after eliminating low ratio subtypes
        			probSum = 0;
        			for(BindingSubtype b: bindingManager.getBindingSubtypes(cond)) {
       					probSum += tau[c][j][b.getIndex()];

        			}
        			for(BindingSubtype b: bindingManager.getBindingSubtypes(cond)) 
        				tau[c][j][b.getIndex()] /= probSum;
        		}
        		}
        	}
        	}
        	
			//monitor: count time
        	timer.end("tau");
        	timer.start();

        	//monitor
//        	System.out.println("\t\tMaximize pi");
        	
        	//TODO: should we remove bindingComponents with fuzziness==0? because 0 fuzziness binding component won't "eat" other component's responsibility
        	//But it's still meaningful that a binding component has a fuzziness==0
        	// Part4: Maximize pi
        	boolean componentEliminated = false;
        	for(int c=0; c<numConditions; c++) {
        		int numPairs = hitNum[c];
        		// Maximize pi
        		double[] test = new double[numComponents];
        		double[] sumR = new double[numComponents];
        		for(int j=0; j<numComp; j++) {
        			if(pi[c][j]>0) {
            			if(pairIndexAroundMu.get(c).get(j).car()!=-1) {
            				for(int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++) {
            					sumR[j] += rBind[c][j][i] * hitCounts[c][i];
            				}
            			}
        			}
        		}
//        		for(int j=0; j<numComp; j++) {
//        			System.out.println("position:" + mu[c][j]);
//        			System.out.println("start:"+pairIndexAroundMu.get(c).get(j).car());
//        			System.out.println("end:"+pairIndexAroundMu.get(c).get(j).cdr());
//        			System.out.println("test: "+test[j]);
//        			System.out.println("sumR: "+sumR[j]);
//        		}
//				System.out.println("sum of responsibility: " +sumR[0]);
        		int minIndex=0; double minVal=Double.MAX_VALUE;
        		for(int j=0; j<numComp; j++) {
        			if(pi[c][j]>0) {
        				if(sumR[j]<minVal) {minVal=sumR[j]; minIndex=j;}
        			}
        		}
        		if(minVal>currAlpha[c]) {
        			// No component will be eliminated, update pi[j]
        			for(int j=0; j<numComp; j++) {
        				if(pi[c][j]>0) {
//        					pi[c][j] = Math.max(0, sumR[j]-currAlpha[c]);
        					pi[c][j] = Math.max(0, sumR[j]);
        				}
        			}
        		} else {
        			// Eliminate worst binding component Q: why only the worst one?
        			// Responsibilities will be redistributed in the E step
        			pi[c][minIndex]=0.0; sumR[minIndex]=0.0;
//        			for(int i=0; i<numPairs; i++)
        			if(pairIndexAroundMu.get(c).get(minIndex).car()!=-1)
        			for(int i=pairIndexAroundMu.get(c).get(minIndex).car(); i<=pairIndexAroundMu.get(c).get(minIndex).cdr(); i++)
        				rBind[c][minIndex][i] = 0;
        			// Re-estimate pi values for non-eliminated components using the current responsibility assignments
        			for(int j=0; j<numComp; j++) {
        				if(j!=minIndex)
        					pi[c][j] = Math.max(0, sumR[j]); // Q: I think here should be sumR[j]-currAlpha[c].
        			}
        			componentEliminated = true;
        		}
        		// Normalize pi (accounting for piNoise)
        		double totalPi = 0;
        		for(int j=0; j<numComp; j++) {
        			if(pi[c][j]>0)
        				totalPi += pi[c][j];
        		}
        		for(int j=0; j<numComp; j++) {
        			if(totalPi>0)
        				pi[c][j] = pi[c][j]/(totalPi/(1-piNoise[c]));
        		}
        		
        		/////////////
            	//Anneal alpha Q: Looks like annealling means reduce alpha each turn, I don't know why
            	//////////////
        		if (t >semconfig.EM_ML_ITER && t <= semconfig.ALPHA_ANNEALING_ITER)
        			currAlpha[c] = alphaMax[c] * (t-semconfig.EM_ML_ITER)/(semconfig.ALPHA_ANNEALING_ITER-semconfig.EM_ML_ITER);
        		else if(t > semconfig.ALPHA_ANNEALING_ITER)
        			currAlpha[c] = alphaMax[c];
        	}
        	
			//monitor: count time
        	timer.end("pi");
        	timer.start();
        	
        	// Non-zero components count
        	
        	int nonZeroComps = 0;
        	for(int c=0; c<numConditions; c++)
        		for(int j=0; j<numComp; j++)
        			if(pi[c][j]>0)
        				nonZeroComps++;
        	
        	//monitor
//        	System.out.println("\tCompute LL");
        	
        	////////////
        	//Compute LL
        	////////////
        	LAP=0;
        	if(semconfig.CALC_LL) {
        		// Log-likelihood calculation
        		double LL = 0;
        		for(int c=0; c<numConditions; c++) {
        			int numPairs = hitNum[c];
//        			for(int i=0; i<numPairs; i++) {
//        				// for each read pair, each component will give a conditional prob or bg prob
//        				double j_sum = 0;
//        				for(int j=0; j<numComp; j++) {
//        					if(pi[c][j]>0.0)
//        						j_sum += Math.log(rBind[c][j][i]/semconfig.LOG2); // Q: In this part, rBind has been normalized per read, I think it cannot be used to calculate log likelihood.
//        				}
//        				j_sum += Math.log(rNoise[c][i])/semconfig.LOG2;
//        			}
        			//test reduce conmputation
        			for(int j=0; j<numComp; j++) {
        				if(pairIndexAroundMu.get(c).get(j).car()!=-1)
                			for(int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++) {
                				LL += Math.log(rBind[c][j][i]/semconfig.LOG2) * hitCounts[c][i];
                			}
        			}
        			for(int i=0; i<numPairs; i++) {
        				LL += Math.log(rNoise[c][i])/semconfig.LOG2 * hitCounts[c][i];
        			}
        		}
        		// Log priors
        		double LP=0;
        		for(int c=0; c<numConditions; c++) {
        			// sum of log pi (Assumption1: Dirichlet prior on nucleosome occupancy)
        			double sum_log_pi = 0;
        			for(int j=0; j<numComp; j++) {
        				if(pi[c][j]>0.0)
        					sum_log_pi += Math.log(pi[c][j])/semconfig.LOG2;
        			}
        			// Position prior (Assumption3: Bernouli prior on nucleosome positions)
        			double sum_pos_prior = 0;
        			if(atacPrior!=null && semconfig.useAtacPrior())
        				for(int j=0; j<numComp; j++)
        					sum_pos_prior += atacPrior[c][mu[c][j]-currRegion.getStart()];
        		}
        		LAP = LL + LP;
        		
//        		System.out.println("EM: "+t+"\t"+LAP+"\t("+nonZeroComps+" non-zero components).");
        	}
        	
			//monitor: count time
        	timer.end("LL");
        	
        	//Tick the clock forward
    		if(!componentEliminated)
    			t++;
    		iter++;
    		
    		//monitor code
//    		long endtime = System.currentTimeMillis();
//    		if((endtime - starttime > 15000)) {
//    			System.out.println("Stack Warning:");
//    			System.out.println("Mu position: " + Arrays.toString(mu[0]));
//    			System.out.println("Pi: " + Arrays.toString(pi[0]));
//    		}
    		
    		 ////////////
          	//Check Stopping condition
          	////////////   		
            if (nonZeroComps>0 && ((numConditions>1 && t<=semconfig.POSPRIOR_ITER) || (numConditions==1 && t<=semconfig.ALPHA_ANNEALING_ITER) || (semconfig.CALC_LL && Math.abs(LAP-lastLAP)>semconfig.EM_CONVERGENCE) || stateEquivCount<semconfig.EM_STATE_EQUIV_ROUNDS)){
                copyStateToLast();
                lastLAP = LAP;
                continue;
            }else{
            	copyStateToLast();
            	lastLAP = LAP;
            	if(semconfig.isVerbose())
            		System.err.println("\tRegTrain:"+trainingRound+"\t"+currRegion.getLocationString()+"\t"+currRegion.getWidth()+"\t"+t+"\t"+iter+"\t"+nonZeroComps);
            	break;
            }
        } //LOOP: Run EM while not converged
        
        //monitor code: show binding component information after EM loop
//        for(int c=0; c<numConditions; c++) {
//        	for(int j=0; j<numComp; j++) {
//        		if(pi[c][j]>0) {
//        			System.out.println("\tBinding Component"+j);
//        			System.out.println("\t\tpi: "+pi[c][j]);
//        			System.out.println("\t\tmu: "+mu[c][j]);
//        			System.out.println("\t\tfuzziness: "+fuzz[c][j]);
//        			System.out.println("\t\ttau: "+Arrays.toString(tau[c][j]));
//        		}
//        	}
//        }
//        System.exit(1);
    } // end of EM_MAP method
	
	/**
     * Copy current variables to last variables (lastRBind, lastPi, lastMu).
     * Assumes visibility of both.
     */
    private void copyStateToLast(){
    	int numC = manager.getNumConditions();
    	for(int c=0; c<numC; c++){
    		for(int j=0; j<numComponents; j++){
    			lastPi[c][j] = pi[c][j];
    			lastMu[c][j] = mu[c][j];
    			lastFuzz[c][j] = fuzz[c][j];
    			lastTau[c][j] = tau[c][j];
    			for(int x=0; x<rBind[c][j].length; x++){
    				lastRBind[c][j][x] = rBind[c][j][x];
    			}
    		}
    	}
    }
    
    /**
     * Compare last variables to current (lastRBind, lastPi, lastMu).
     * Assumes visibility of both.
     * @return
     */
    private boolean lastEquivToCurr(){
    	int numC = manager.getNumConditions();
    	int currNZ=0, lastNZ=0;
    	for(int c=0; c<numConditions; c++)
    		for(int j=0;j<numComponents;j++){
    			if(pi[c][j]>0)
    				currNZ++;
    			if(lastPi[c][j]>0)
    				lastNZ++;
    		}
    	boolean numCompEqual = currNZ==lastNZ;
    	boolean compPosEqual=true;
    	if(numCompEqual){
    		for(int c=0; c<numC; c++)
    			for(int j=0; j<mu[c].length; j++){if(pi[c][j]>0){
    				compPosEqual = compPosEqual && (mu[c][j] == lastMu[c][j]);
    			}}
    	}else{
    		compPosEqual=false;
    	}
    	boolean piBindEquivalent=true;
    	for(int c=0; c<numC; c++)
			for(int j=0; j<pi[c].length; j++){if(pi[c][j]>0){
				piBindEquivalent = piBindEquivalent && (Math.abs(pi[c][j]-lastPi[c][j])<semconfig.EM_STATE_EQUIV_THRES);
			}}
    	boolean fuzzBindEquivalent=true;
    	for(int c=0; c<numC; c++)
    		for(int j=0; j<pi[c].length; j++) {
    			if(pi[c][j]>0)
    				fuzzBindEquivalent = fuzzBindEquivalent && (Math.abs(fuzz[c][j]-lastFuzz[c][j]) < semconfig.EM_STATE_EQUIV_THRES);
    		}
    	boolean tauBindEquivalent=true;
    	for(int c=0; c<numC; c++)
    		for(int j=0 ;j<pi[c].length; j++) {
    			if(pi[c][j]>0)
    				for(int t=0; t<tau[c][j].length; t++)
    					tauBindEquivalent = tauBindEquivalent && (Math.abs(tau[c][j][t] - lastTau[c][j][t]) < semconfig.EM_STATE_EQUIV_THRES);
    		}
    	boolean rBindEquivalent=true;
    	for(int c=0; c<numC; c++)
    		for(int j=0; j<pi[c].length; j++){if(pi[c][j]>0){
    			for(int x=0; x<rBind[c][j].length; x++){
    				rBindEquivalent = rBindEquivalent && (Math.abs(rBind[c][j][x]-lastRBind[c][j][x])<semconfig.EM_STATE_EQUIV_THRES);
    			}
			}}
		return numCompEqual && compPosEqual && piBindEquivalent && fuzzBindEquivalent && tauBindEquivalent && rBindEquivalent;
    }
}












































