package org.seqcode.projects.sem.mixturemodel;

/**
 * Statistic version of BindingEM
 * Nucleosome comparison procedure would be:
 * 1. F test: Check whether two nucleosomes share the same fuzziness (homogeneity of variance)
 * 2. T test or Welch's T test: To determine whether two nucleosomes share the same dyad location (mean)
 * 3. According to the above statistic results to determine two boolean
 * 		1) fuzzinessSharedBetter ?
 * 		2) positionSharedBetter ?
 * 4. In multi-condition (conditions>3) test
 * 		1) mean of locations of nucleosomes marked with positionSharedBetter will have a shared position
 * 		2) variance of nucleosomes marked with fuzzinessSharedBetter will have a shared fuzziness
 */

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Comparator;
import java.util.stream.Stream;
import java.util.stream.Collectors;
import java.math.BigDecimal;
import java.math.RoundingMode;

import org.seqcode.deepseq.StrandedPair;
import org.seqcode.deepseq.experiments.ControlledExperiment;

import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.genome.location.Region;
import org.seqcode.projects.sem.framework.SEMConfig;
import org.seqcode.projects.sem.events.BindingManager;
import org.seqcode.projects.sem.events.BindingModel;
import org.seqcode.projects.sem.events.BindingSubtype;
import org.seqcode.projects.sem.utilities.Timer;
import org.seqcode.projects.sem.utilities.EMStepPlotter;
import org.seqcode.projects.sem.utilities.Statistics;
import org.seqcode.projects.sem.utilities.GaleShapley;
import org.seqcode.projects.sem.utilities.NucleosomePoissonBackgroundModel;
import org.seqcode.projects.sem.utilities.EMmode;
import org.seqcode.gseutils.Pair;

public class BindingEM_Statistic implements BindingEM_interface {
	
	protected ExperimentManager manager;
	protected BindingManager bindingManager;
	protected SEMConfig semconfig;
	protected List<List<BindingComponent>> components;
	protected List<NoiseComponent> noise;
	protected int numComponents;
	protected int numConditions;
	protected int numSubtypes;
	protected int trainingRound=0;
	protected HashMap<ExperimentCondition, NucleosomePoissonBackgroundModel> conditionBackgrounds;
//	protected HashMap<ExperimentCondition, Integer> numBindingSubtypes;
	protected Timer timer;
	protected EMmode mode;
	
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
	protected double[][][]	resp;				// Responsibility each bindingComponent to each fragments
	protected double[][]	sumResp;			// Sum of responsibility each bindingComponent
	protected double[][]	pValue;				// p value of each components if alpha is determined by poisson background
	protected double[][]	pi;					// pi: emission probabilities for binding components
	protected double[]		piNoise;			// piNoise: emission probabilities for noise components (fixed)
	protected int[][]		mu;					// mu: positions of the binding components (indexed by condition & binding component index)
	protected double[][]	fuzz;				// &fuzz: fuzziness of the binding components (indexed by condition & binding component index)
	protected double[][][]	tau;				// &tau: fragment size subtype probabilities (indexed by subtype index)
	protected double[]		alphaMax;			// Maximum alpha (indexed by condition)
	protected double[][]	atacPrior;			// &ATAC-seq prior (indexed by condition & base)
	protected double[][][]	lastRBind;			// Last responsibilities (monitor convergence)
	protected double[][]	lastSumResp;		// Last sum of responsibilities (moniotr convergence)
	protected double[][]	lastPi;				// Last Pi (monitor convergence)
	protected int[][]		lastMu;				// Last positions (monitor convergence)
	protected double[][]	lastFuzz;			// &Last fuzziness (monitor convergence)
	protected double[][][]	lastTau;				// &Last tau (monitor convergence)
	protected int[][]		maxIR;				// Max influence range for each nucleosome
	protected double 		lastLAP, LAP;		// log-likelihood monitoring
	protected boolean 		plotEM=false;		// Plot the current region components
	protected Region		plotSubRegion=null;	// Sub region to plot
	protected double		numPotentialRegions;
	protected double		probAgivenB, probAgivenNOTB;
	protected int stateEquivCount = 0;
	protected Map<Pair<Integer, Integer>, Map<Pair<Integer, Integer>, Pair<double[], boolean[]>>> compareStore;
	
	public BindingEM_Statistic(SEMConfig s, ExperimentManager eMan, BindingManager bMan, HashMap<ExperimentCondition, NucleosomePoissonBackgroundModel> condBacks, int numPotReg) {
		semconfig = s;
		manager = eMan;
		bindingManager = bMan;
		conditionBackgrounds = condBacks;
		numConditions = manager.getNumConditions();
		numPotentialRegions = (double)numPotReg;		
		mode = EMmode.NORMAL;
	}
	
	public List<List<BindingComponent>> train(List<List<StrandedPair>> signals,
												Region w,
												List<NoiseComponent> noise,
												List<List<BindingComponent>> comps,
												int numComp,
												double[][] atacPrior,
												int trainingRound,
												EMmode mode,
												Timer timer,
												boolean plotEM
												) throws Exception {
		components = comps;
		this.noise = noise;
		numComponents = numComp;
		this.atacPrior = atacPrior;
		this.trainingRound = trainingRound;
		this.plotEM = plotEM;
		this.timer = timer;
		this.mode = mode;
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
		resp = new double[numConditions][][];
		sumResp = new double[numConditions][];			// Sum of responsibilities each bindingComponent
		pValue = new double[numConditions][];			// p value of each components if alpha is determined by poisson background
		pi = new double[numConditions][numComponents];	// pi: emission probabilities for binding components
		piNoise = new double[numConditions];			// pi: emission probabilities for noise components (fixed)
		alphaMax = new double[numConditions];			// Maximum alpha
		mu = new int[numConditions][numComponents];		// mu: positions of the binding components
		fuzz = new double[numConditions][numComponents];	// &fuzz: fuzziness of the binding components
		tau = new double[numConditions][numComponents][];				// &tau: fragment size subtype probabilities (indexed by subtype index)
		maxIR = new int[numConditions][numComponents];	// max influence range for each nucleosome
		// Monitor state convergence using the following last variables
		lastRBind = new double[numConditions][][];
		lastSumResp = new double[numConditions][numComponents];
		lastPi = new double[numConditions][numComponents];
		lastMu = new int[numConditions][numComponents];
		lastFuzz = new double[numConditions][numComponents];
		lastTau = new double[numConditions][numComponents][];
		
		// Initializing data structures
		for(ExperimentCondition cond: manager.getConditions()) {
			int c = cond.getIndex();
			
			// Set maximum alphas (TODO: Now we set 0 for SEM because we haven't found a good way to determine ahpha)
//			alphaMax[c] = semconfig.getFixedAlpha()>0 ? semconfig.getFixedAlpha() :
//					semconfig.getAlphaScalingFactor() * (double)conditionBackgrounds.get(cond).getMaxThreshold('.');
			
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
			
			// if there are no pairs in this condition, add a pseudo pair
			if(pairs.size()==0) {
				pairs.add(new StrandedPair(semconfig.getGenome(), semconfig.getGenome().getChromID(w.getChrom()), -11000, '+', 
						semconfig.getGenome().getChromID(w.getChrom()), -10000, '-', 1));
				repIndices[c] = new int[1];
				repIndices[c][0] = cond.getReplicates().get(0).getIndex();
				numPairs = 1;
				hitNum[c] = 1;
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
			resp[c] = new double[numComp][numPairs];
			sumResp[c] = new double[numComp];
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
		List<Map<Integer, Integer>> indexConverter = new ArrayList<Map<Integer, Integer>>();
		for(int c=0; c<numConditions; c++) {
			indexConverter.add(new HashMap<Integer, Integer>());
			// Binding Components
			List<BindingComponent> currActiveComps = new ArrayList<BindingComponent>();
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
	            if(semconfig.getFixedAlpha()<0)
	            	comp.setPValue(pValue[c][j]);
				if(pi[c][j]>0.0) {
					if(numConditions > 1)
						comp.setCompareResults(compareStore.get(new Pair<Integer, Integer>(c, j)));
					currActiveComps.add(comp);
					indexConverter.get(c).put(j, currActiveComps.size()-1);
				}
			}
			activeComponents.add(currActiveComps);
			// Noise components
			double noise_resp = 0.0;
			for(int i=0; i<hitNum[c]; i++)
				noise_resp += hitCounts[c][i]*rNoise[c][i];
			noise.get(c).setSumResponsibility(noise_resp);
		}
		// convert index of pairwise component to the index in activeComponents
		for(int c=0; c<numConditions; c++) {
			for(BindingComponent bc: activeComponents.get(c)) {
				bc.convertIndex(indexConverter);
				
				//monitor
//				if(bc.isPair()) {
//					System.out.println("original nucleosome info");
//					System.out.println("\tmu: "+bc.position);
//					System.out.println("\tfuzz: "+bc.fuzziness);
//					System.out.println("compare nucleosome info");
//					for(Pair<Integer, Integer> key: bc.getCompareRestulsConvert().keySet()) {
//						BindingComponent bcPair = activeComponents.get(key.car()).get(key.cdr());
//						System.out.println("\tmu: "+bcPair.position);
//						System.out.println("\tfuzz: "+bcPair.fuzziness);
//						System.out.println("\tcompareResults: sharePos " + bc.getCompareRestulsConvert().get(key).car() + "\tshareFuzz " + bc.getCompareRestulsConvert().get(key).cdr());
//					}
//				}
			}
		}
//		System.exit(1);
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
        double[][] muSum = new double[numConditions][]; //Sum of hitPos*rBind*hitCounts for each binding component
        double[][] fuzzSum = new double[numConditions][]; //Sum of Variance*rBind*hitCounts for each binding component
        int[] muJoinClosestComps = new int[numConditions]; //Indices of nearest components in other conditions
        boolean[] muSharedBetter = new boolean[numConditions]; //Indicator that sharing dyad locations across conditions is better than not
        boolean[] fuzzSharedBetter = new boolean[numConditions]; //Indicator that sharing fuzziness across conditions is better than not
        int[][] newMu = new int[numConditions][numComponents];// mu update
        double[][] newFuzz = new double[numConditions][numComponents]; // fuzziness update
                
        // Initialize responsibilities
        for(int c=0; c<numConditions; c++) {
        	int numPairs = hitNum[c];
        	totalResp[c] = new double[numPairs];
        	for(int i=0; i<numPairs; i++)
        		totalResp[c][i] = 0;
        }
        
        // Alpha and p value initialization
        double[][] currAlpha = new double[numConditions][numComponents];
        double alphaCoefficient = 0;
        for(int c=0; c<numConditions; c++)
        	for(int j=0; j<numComponents; j++) 
        		currAlpha[c][j] = 0;
        
        //record the initial parameters if plotting
        if(plotEM) {
        	EMStepPlotter plot = new EMStepPlotter(currRegion, semconfig, hitPos, hitCounts, hitSize, trainingRound);
        	EMStepPlotter.excute(mu, pi, fuzz, tau, 0, 0);
        }
        
    	//////////////////////////////////////////////////////////////////////////////////////////
        //Run EM while not converged
        // Note: iterations during which we eliminate a binding component don't count towards "t"
    	//////////////////////////////////////////////////////////////////////////////////////////
        int t=0, iter=0;   
        
        while(t<semconfig.MAX_EM_ITER) {
    		////////
    		//E-step
    		////////
        	
        	//I want to get the range of fragment midpoint coordinates that will be influenced by each BindingComponent
        	//It's a trick to improve the performance because when encountering quite long region, calculate each fragment
        	//for each BindingComponent will waste a lot of time.
        	
			//monitor: count time
        	timer.start();
        	
        	// Marking range for each bindingComponent
        	List<Map<Integer, Pair<Integer, Integer>>> pairIndexAroundMu = new ArrayList<Map<Integer, Pair<Integer, Integer>>>(); 
        	for(int c=0; c<numConditions; c++) {

        		pairIndexAroundMu.add(new HashMap<Integer, Pair<Integer, Integer>>());
        		for(int j=0; j<numComp; j++) {
        			// Get half 95% influence range for each nucleosome
        			maxIR[c][j] = (int)(Math.sqrt(fuzz[c][j]) * 1.96 * 2);
        			int half_maxIR = maxIR[c][j] / 2;
        			
        			// Increase speed by binary search
        			int left_bound = mu[c][j] - half_maxIR;
        			int right_bound = mu[c][j] + half_maxIR;
        			int start = Arrays.binarySearch(hitPos[c], left_bound);
        			int end = Arrays.binarySearch(hitPos[c], right_bound);
        			
        			if(start < 0) { start = -start - 1;}
        			if(end < 0) { end = -end - 1;}
        			
        			while(start > 0 && hitPos[c][start-1] >= left_bound) {
        				start--;
        			}
        			while(end < hitPos[c].length && hitPos[c][end] <= right_bound) {
        				end++;
        			}
        			
        			// if there is no fragment around dyad, mark start=-1
        			end -= 1;
        			if(start>end) {
        				start = -1;
        			}
        			pairIndexAroundMu.get(c).put(j, new Pair<Integer, Integer>(start, end));
        		}
        	}
        	
        	//compute current alpha for each component
	        for(int c=0; c<numConditions; c++) {
	        	ExperimentCondition cond = manager.getIndexedCondition(c);
	        	for(int j=0; j<numComponents; j++) 
	       			currAlpha[c][j] = semconfig.getFixedAlpha()<0 ? Math.max(conditionBackgrounds.get(cond).calcCountThreshold(maxIR[c][j]) * alphaCoefficient, 1):semconfig.getFixedAlpha();
	       	}
        	
			//monitor: count time
        	timer.end("mark");
        	
        	for(int c=0; c<numConditions; c++) {
        		int numPairs = hitNum[c];
        		ExperimentCondition cond = manager.getIndexedCondition(c);
        		
    			// Load binding fragment size PDF cache (indexed by type index)
    			double[][] fragSizePDF = bindingManager.getCachePDF(cond);

    			//monitor: count time
            	timer.start();
    			
    			// Compute H and N function
    			double[][] hc = new double[numComp][numPairs];
    			double[] nc = new double[numPairs];
    			
    			double fuzzProb;
    			double fragSizeProb;
    			for(int j=0; j<numComp; j++) {
    				if(pi[c][j]>0 && pairIndexAroundMu.get(c).get(j).car()!=-1) {
    					for(int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++) {
    						fuzzProb = BindingModel.probability(fuzz[c][j], mu[c][j]-hitPos[c][i]);
        					fragSizeProb = 0;
        					for(int index=0; index< fragSizePDF.length; index++) {
        						fragSizeProb += tau[c][j][index] * fragSizePDF[index][hitSize[c][i]];
        					}
        					hc[j][i] = fuzzProb * fragSizeProb;
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
    			
        		// Compute responsibilities
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
        			//normalize noise here
        			rNoise[c][i]/=totalResp[c][i];
        			
        			//monitor
        			if(Double.isNaN(rNoise[c][i])) {
        				System.out.println("NaN happens in normalize: "+"\tBefore normalize: "+n[c][i]*piNoise[c]+"\ttotalResp: "+totalResp[c][i]);
        			}
        		}
        		
    			//monitor: count time
        		timer.end("cResp");
        		timer.start();
        		
        		// Normalize responsibilities
        		for(int j=0; j<numComp; j++) {
        			if(pi[c][j]>0 && pairIndexAroundMu.get(c).get(j).car()!=-1) 
        				for(int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++) {
        					rBind[c][j][i] /= totalResp[c][i];
        					rBind[c][j][i] = new BigDecimal(String.valueOf(rBind[c][j][i])).setScale(2, RoundingMode.HALF_UP).doubleValue();
        				}
        		}
            	
    			//monitor: count time
            	timer.end("nResp");
        	}
        	
        	timer.start();
        	
        	//monitor
//        	System.out.println("\t\tResolve duplicate");
        	
        	// Resolve duplicate binding components (share the same position) -> combine & delete one copy
        	//TODO: I don't know how to combine binding components with different fuzziness at the same location, so I decide to combine binding components before computing fuzziness
        	for(int c=0; c<numConditions; c++) {
        		HashMap<Integer, Integer> pos2index = new HashMap<Integer, Integer>();	// Position to array index map
        		for(int j=0; j<numComp; j++) {
        			if(pi[c][j]>0 && pairIndexAroundMu.get(c).get(j).car()!=-1) {
        				if(pos2index.containsKey(mu[c][j])) {
        					int orig = pos2index.get(mu[c][j]);
        					// Combine
        					pi[c][orig] += pi[c][j];
        					for(int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++)
        						rBind[c][orig][i] += rBind[c][j][i];
        					// Delete
        					pi[c][j] = 0.0;
        					for(int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++)
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
        	
        	//Compute sum of responsibility and compute p value
        	for(int c=0; c<numConditions; c++) {
        		ExperimentCondition cond = manager.getIndexedCondition(c);
        		pValue[c] = new double[numComp];
        		sumResp[c] = new double[numComp];
        		resp[c] = new double[numComp][hitNum[c]];
        		for(int j=0; j<numComp; j++) {
        			if(pi[c][j]>0 && pairIndexAroundMu.get(c).get(j).car()!=-1) {
        				for(int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++) {
        					resp[c][j][i] = rBind[c][j][i] * hitCounts[c][i];
        					sumResp[c][j] += resp[c][j][i];
        				}
        				pValue[c][j] = conditionBackgrounds.get(cond).calcPValue(maxIR[c][j], sumResp[c][j]);
        			}
        		}
        	}
        	
        	timer.end("sResp");
        	timer.start();

    		/////////////////////
    		//M-step: maximize mu (positions), fuzz (fuzziness), tau (fragment size subtype), pi (strength)
    		/////////////////////
        	if(numConditions>1 && t==semconfig.ALPHA_ANNEALING_ITER)
        		for(int c=0; c<numConditions; c++)
        			for(int j=0; j<numComp; j++)
        				if(pi[c][j]>0)
        					muSums[c][j] = new double[semconfig.EM_MU_UPDATE_WIN*2];
        	
        	// Part1.1: Maximize mu (Because we use Gaussian distribution to describe tag distribution around dyad, we can get dyad directly by the proportional sum of tag position)
        	for(int c=0; c<numConditions; c++) {
        		muSum[c] = new double[numComp];
        		for(int j=0; j<numComp; j++) {
        			if(pi[c][j]>0 && pairIndexAroundMu.get(c).get(j).car()!=-1) {
        				// Get the limit of new dyad locations
        				int start = Math.max(mu[c][j]-semconfig.EM_MU_UPDATE_WIN, regStart);
        				int end = Math.min(currRegion.getEnd(), mu[c][j]+semconfig.EM_MU_UPDATE_WIN);
        				// Assign variables for pairwise comparison
        				if(numConditions>1 && t>semconfig.ALPHA_ANNEALING_ITER) {
        					muSumStarts[c][j] = start; muSumWidths[c][j] = end-start;
        				}
        				
        				for(int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++) {
        					muSum[c][j] += resp[c][j][i] * hitPos[c][i];
        				}

        				// Does maximize position exceed update window?
        				double mu_max = muSum[c][j]/sumResp[c][j];
        				muMax[c][j] = (int)Math.rint(mu_max);
        				muMax[c][j] = Math.min(muMax[c][j], end);
        				muMax[c][j] = Math.max(start, muMax[c][j]);
        			}
        		}
        	}
        	
			//monitor: count time
        	timer.end("mu");
        	timer.start();
        	
        	// Part2: Maximize fuzziness
        	if(t/semconfig.FUZZINESS_ANNEALING_ITER==0 || t>semconfig.ALPHA_ANNEALING_ITER) {
        		for(int c=0; c<numConditions; c++) {
        			fuzzSum[c] = new double[numComp];
        			for(int j=0; j<numComp; j++) {if(pi[c][j]>0 && pairIndexAroundMu.get(c).get(j).car()!=-1) {
        				double V1 = 0; double V2 = 0;
	        			for(int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++) {
	        				fuzzSum[c][j] += resp[c][j][i] * Math.pow(hitPos[c][i]-muMax[c][j], 2);
	        				V1 += resp[c][j][i];
	        				V2 += Math.pow(resp[c][j][i], 2);
	        			}
	        			//set fuzzMax = 0  if V1^2 == V2
	        			if((Math.pow(V1, 2) - V2) == 0) 
	        				fuzzMax[c][j] = 0;
	        			else
	        				fuzzMax[c][j] = (V1*fuzzSum[c][j])/(Math.pow(V1, 2) - V2);
	        			fuzzMax[c][j] = new BigDecimal(String.valueOf(fuzzMax[c][j])).setScale(2, RoundingMode.HALF_UP).doubleValue();
	        			if(fuzzMax[c][j]<0) {
	        				System.out.println("V1:" + V1);
	        				System.out.println("V1^2:" + Math.pow(V1, 2));
	        				System.out.println("V2:" + V2);
	        				System.out.println("V2/V1:" + (V2/V1));
	        				System.out.println("scale: "+(V1-(V2/V1)));
		        			for(int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++) 
		        				System.out.println(resp[c][j][i]);
	        			}
        			}
        			}
        		}
        	}
        	
			//monitor: count time
        	timer.end("fuzz");
        	
        	timer.start();
        	//Insert pairwise bindingcomponents comparison here
        	//Due to I haven't found a way to compare fragment size distribution across different experiment condition
        	//I only consider position and fuzziness here.
        	
        	//test code: employ Gale/Shapley algorithm to find pairwise nucleosomes
        	//TODO:Can we designate an anchor sample ?
        	//1. for each condition, get active components index
        	int ac = 0;		// anchor condition index
        	List<List<Integer>> activeIndex = new ArrayList<List<Integer>>();
        	for(int c=0; c<numConditions; c++) {
        		activeIndex.add(new ArrayList<Integer>());
        		for(int j=0; j<numComp; j++) {
        			if(pi[c][j]>0 && pairIndexAroundMu.get(c).get(j).car()!=-1)
        				activeIndex.get(c).add(j);
        		}
        	}
        	//2. for each nucleosome, get preference list of nucleosomes from other conditions
        	Map<Pair<Integer, Integer>, Map<Integer, List<Integer>>> preference = new HashMap<Pair<Integer, Integer>, Map<Integer, List<Integer>>>();
        	for(int c=0; c<numConditions; c++) {
        		for(int d=0; d<numConditions; d++) {
        			if(c!=d) {
        				Pair<Integer, Integer> key = new Pair<Integer, Integer>(c, d);
        				preference.put(key, new HashMap<Integer, List<Integer>>());
        			}
        		}
        	}
        	for(int c=0; c<numConditions; c++) {
	        	for(int j=0; j<numComp; j++) {if(pi[c][j]>0 && pairIndexAroundMu.get(c).get(j).car()!=-1 && sumResp[c][j]>Math.max(currAlpha[c][j], 1)) {
	        		for(int d=0; d<numConditions; d++) {
	        			if(c!=d) {
	        				Pair<Integer, Integer> key = new Pair<Integer, Integer>(c, d);
	        				preference.get(key).put(j, new ArrayList<Integer>());
	        				List<Integer> preferenceInter = new ArrayList<Integer>();
	        				List<Integer> distStore = new ArrayList<Integer>();
	        				for(int k=0; k<numComp; k++) {
	        					if(pi[d][k]>0 && pairIndexAroundMu.get(d).get(k).car()!=-1 && sumResp[d][k]>Math.max(currAlpha[d][k], 1)) {
	        						int dist = Math.abs(mu[c][j]-mu[d][k]);
	        						if(dist<semconfig.EM_MU_UPDATE_WIN) {
	        							preferenceInter.add(k);
	        							distStore.add(dist);
	        						}
	        					}
	        				}
	        				// sort preference list by distStore
	        				List<Integer> index = Stream.iterate(0, n->n+1).limit(distStore.size()).collect(Collectors.toList());
	        				Collections.sort(index, new Comparator<Integer>() {
	        					@Override public int compare(final Integer o1, final Integer o2) {
	        						return Float.compare(distStore.get(o1), distStore.get(o2));
	        					}        					
	        				});
	        				for(int i: index) {
	        					preference.get(key).get(j).add(preferenceInter.get(i));
	        				}
	        			}
	        		}
	        	}}
        	}
        	//3. put each pair of conditions into Gale-Shapley algorithm
        	Map<Pair<Integer, Integer>, GaleShapley<Integer>> gsMap = new HashMap<Pair<Integer, Integer>, GaleShapley<Integer>>();
        	for(int c=0; c<numConditions; c++) {
	        	for(int d=0; d<numConditions; d++) {
	        		if(d!=c) {
	        			Pair<Integer, Integer> menKey = new Pair<Integer, Integer>(c, d);
	        			Pair<Integer, Integer> womenKey	= new Pair<Integer, Integer>(d, c);
	        			GaleShapley<Integer> gs = new GaleShapley<Integer>(preference.get(menKey), preference.get(womenKey));
	        			gs.match();
	        			//save GaleShapley class
	        			gsMap.put(new Pair<Integer, Integer>(c, d), gs);
	        		}
	        	}
        	}
        	
        	// Store pairwise comparison results to reduce duplicated computation
        	Map<PairwiseKey, Pair<double[], boolean[]>> pairwiseComparisonResults = new HashMap<PairwiseKey, Pair<double[], boolean[]>>();
        	compareStore = new HashMap<Pair<Integer, Integer>, Map<Pair<Integer, Integer>, Pair<double[], boolean[]>>>();

        	// pairwise comparison, skip all nucleosomes with sumResp <= max(currAlpha, 1)
        	for(int c=0; c<numConditions; c++) {
        		for(int j=0; j<numComp; j++) {
        			if(pi[c][j]>0 && pairIndexAroundMu.get(c).get(j).car()!=-1) {
        				if(numConditions>1 && t>semconfig.ALPHA_ANNEALING_ITER && sumResp[c][j] > Math.max(currAlpha[c][j], 1)) {
        				
        				//a: get the paired nucleosome to j in each condition
        				for(int d=0; d<numConditions; d++) {
        					Pair<Integer, Integer> key = new Pair<Integer, Integer>(c, d);
        					if(d!=c) {
        						GaleShapley<Integer> gs = gsMap.get(key);
        						if(gs.isManSpecific(j)) {
        							muJoinClosestComps[d] = -1;
        						} else {
        							muJoinClosestComps[d] = gs.getWife(j);
        						}
        					}
        				}
        				
        				//b: evaluate each pair of conditions, asking if a shared event involving j and its closest component would be better than independent events
        				Pair<Integer, Integer> storeKey = new Pair<Integer, Integer>(c, j);
						compareStore.put(storeKey, new HashMap<Pair<Integer, Integer>, Pair<double[], boolean[]>>());
        				int numMuSharedBetter = 0;
        				int numFuzzSharedBetter = 0;
        				int maxMuStart = muSumStarts[c][j];
        				int minMuEnd = muSumStarts[c][j] + muSumWidths[c][j];
        				for(int d=0; d<numConditions; d++) {
        					if(d!=c) {
        						int k = muJoinClosestComps[d];
        						if(k==-1) {
        							muSharedBetter[d] = false;
        							fuzzSharedBetter[d] = false;
        						}
        						else {
        							PairwiseKey key;
        							if(c<d) {
        								key = new PairwiseKey(c, j, d, k);
        							} else {
        								key = new PairwiseKey(d, k, c, j);
        							}
        							int maxSharedPos; double maxSharedFuzz;
        							double[] confidence; boolean[] sharedBetter; 
        							if(pairwiseComparisonResults.containsKey(key)) {
        								//get shared position, fuzziness directly from stored map
        								confidence = pairwiseComparisonResults.get(key).car();
        								sharedBetter = pairwiseComparisonResults.get(key).cdr();
        							} else {
        								//Invoke comparison method in Statistics.java
        								//1. prepare array for comparison
        								List<Integer> hitPos1 = new ArrayList<Integer>();	List<Integer> hitPos2 = new ArrayList<Integer>();
        								List<Double> resp1 = new ArrayList<Double>();		List<Double> resp2 = new ArrayList<Double>();
        								for(int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++) {
        									for(int z=0; z<hitCounts[c][i]; z++) {
        										hitPos1.add(hitPos[c][i]);
        										resp1.add(resp[c][j][i]);
        									}
        								}
        								for(int i=pairIndexAroundMu.get(d).get(k).car(); i<=pairIndexAroundMu.get(d).get(k).cdr(); i++) {
        									for(int z=0; z<hitCounts[d][i]; z++) {
        										hitPos2.add(hitPos[d][i]);
        										resp2.add(resp[d][k][i]);
        									}
        								}

        								//2. get comparison results
        								Pair<double[], boolean[]> comparisonResults = Statistics.comparison(hitPos1, hitPos2, resp1, resp2);
        								confidence = comparisonResults.car();
        								sharedBetter = comparisonResults.cdr();
        		        				
        							}
        							
        							//get mu and fuzz share information
    								muSharedBetter[d] = sharedBetter[0];
    								fuzzSharedBetter[d] = sharedBetter[1];
        							
        							//save pairwise comparison results
        							compareStore.get(storeKey).put(new Pair<Integer, Integer>(d, k), new Pair<double[], boolean[]>(confidence, sharedBetter));
        							
        							maxMuStart = muSharedBetter[d] ? Math.max(maxMuStart, muSumStarts[d][k]) : maxMuStart; 
        							minMuEnd = muSharedBetter[d] ? Math.min(minMuEnd, muSumStarts[d][k]+muSumWidths[d][k]) : minMuEnd;
        							if(muSharedBetter[d])
        								numMuSharedBetter++;
        							if(fuzzSharedBetter[d])
        								numFuzzSharedBetter++;
        						}
        					}
        				}
        				
        				//get shared dyad location based on nucleosomes marked by muSharedBetter
        				int maxSharedMu;
        				double sharedMuSum = muSum[c][j];
        				double sharedMuSumResp = sumResp[c][j];
        				for(int d=0; d<numConditions; d++) {
        					if(d!=c) {
        						int k = muJoinClosestComps[d];
        						if(k!=-1 && muSharedBetter[d]) {
        							sharedMuSum += muSum[d][k];
        							sharedMuSumResp += sumResp[d][k];
        						}
        					}
        				}
        				maxSharedMu = (int)Math.rint(sharedMuSum/sharedMuSumResp);
        				maxSharedMu = maxSharedMu<maxMuStart ? maxMuStart : maxSharedMu;
        				maxSharedMu = maxSharedMu>minMuEnd   ? minMuEnd   : maxSharedMu;
        				//get shared fuzziness based nucleosomes marked by fuzzSharedBetter
        				double maxSharedFuzz;    
        				double sharedFuzzSum = 0;
//        				double sharedFuzzSumResp = sumResp[c][j];
        				double V1 = 0; double V2 = 0;
        				for (int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++) {
							sharedFuzzSum += resp[c][j][i] * Math.pow(hitPos[c][i]-maxSharedMu, 2);
							V1 += resp[c][j][i];
							V2 += Math.pow(resp[c][j][i], 2);
						}
    					for(int d=0; d<numConditions; d++) {
    						if(d!=c) {
    							int k = muJoinClosestComps[d];
    							if(k!=-1 && fuzzSharedBetter[d] && muSharedBetter[d]) {
    								for (int i=pairIndexAroundMu.get(d).get(k).car(); i<=pairIndexAroundMu.get(d).get(k).cdr(); i++) {
    									sharedFuzzSum += resp[d][k][i] * Math.pow(hitPos[d][i]-maxSharedMu, 2);
    									V1 += resp[d][k][i];
    									V2 += Math.pow(resp[d][k][i], 2);
    								}
    								//TODO: I remove this else branch because if two nucleosomes don't share position, they won't share fuzziness
//    								else {
//    									for (int i=pairIndexAroundMu.get(d).get(k).car(); i<=pairIndexAroundMu.get(d).get(k).cdr(); i++) {
//    										sharedFuzzSum += resp[d][k][i] * Math.pow(hitPos[d][i]-muMax[d][k], 2);
//    									}
//    								}
//    								sharedFuzzSumResp += sumResp[d][k];
    							}
    						}
    					}
	        			//set maxSharedFuzz = 0  if V1^2 == V2
	        			if((Math.pow(V1, 2) - V2) == 0) 
	        				maxSharedFuzz = 0;
	        			else
	        				maxSharedFuzz = (V1*sharedFuzzSum)/(Math.pow(V1, 2) - V2);
	        			maxSharedFuzz = new BigDecimal(String.valueOf(maxSharedFuzz)).setScale(2, RoundingMode.HALF_UP).doubleValue();
        				//update mu
    					newMu[c][j] = maxSharedMu;
    					newFuzz[c][j] = maxSharedFuzz;
        			}else {
        				//Ignore other conditions in first phase of training
        				newMu[c][j] = muMax[c][j];
        				newFuzz[c][j] = fuzzMax[c][j];
        			}
        			} 
        		}
        	}
        	
        	// update mu 
        	for(int c=0; c<numConditions; c++) {
        		for(int j=0; j<numComp; j++) {
        			if(pi[c][j]>0 && pairIndexAroundMu.get(c).get(j).car()!=-1) {
        				mu[c][j] = newMu[c][j];
        				fuzz[c][j] = newFuzz[c][j];
        			}
        		}
        	}
        	
        	timer.end("pair");
      
        	timer.start();

        	//monitor
//        	System.out.println("\t\tMaximize tau");
        	
        	// Part3: Maximize tau
        	if(t/semconfig.TAU_ANNEALING_ITER==0) {
        	for(ExperimentCondition cond: manager.getConditions()) {
        		int c = cond.getIndex();
        		int numSubtypes = bindingManager.getBindingSubtypes(cond).size();
				List<BindingSubtype> bindingSubtypes = bindingManager.getBindingSubtypes(cond);
        		for(int j=0; j<numComp; j++) {if(pi[c][j]>0 && pairIndexAroundMu.get(c).get(j).car()!=-1) {
        			tau[c][j] = new double[numSubtypes];
        			double probSum = 0;
        			
        			// Compute tau probability of each fragment size subtype for each component
        			for(int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++) {
        				double[] subProb = new double[bindingManager.getBindingSubtypes(cond).size()];
        				double subProbSum = 0;
        				
        				for(BindingSubtype b: bindingSubtypes) {
        					subProb[b.getIndex()] = b.probability(hitSize[c][i]);
        					subProbSum += subProb[b.getIndex()];
        				}
        				
        				for(BindingSubtype b: bindingSubtypes) {
        					subProb[b.getIndex()] /= subProbSum;
        					tau[c][j][b.getIndex()] += subProb[b.getIndex()] * resp[c][j][i];
        				}
        			}
        			
        			// Normalize tau for each component
        			for(BindingSubtype b: bindingSubtypes) 
        				tau[c][j][b.getIndex()] /= sumResp[c][j];
        			// Eliminate subtype with not enough ratio of tags (Assumption2: each nucleosome should not be associated with too many subtypes)
        			for(BindingSubtype b: bindingSubtypes) {
        				if(tau[c][j][b.getIndex()] < semconfig.SPARSE_PRIOR_SUBTYPE)
        					tau[c][j][b.getIndex()] = 0;
        			}
        			// Re-normalize tau for each component after eliminating low ratio subtypes
        			probSum = 0;
        			for(BindingSubtype b: bindingSubtypes) {
       					probSum += tau[c][j][b.getIndex()];

        			}
        			for(BindingSubtype b: bindingSubtypes) 
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
        		// Maximize pi
        		double[] sumR = new double[numComponents];
        		sumR = sumResp[c];
        		boolean[] isEliminated = new boolean[numComponents];
        		boolean ifEliminate = false;
        		for(int j=0; j<numComp; j++) {
        			if(pi[c][j]>0) {
        				if(sumR[j]<=currAlpha[c][j]) {
        					isEliminated[j] = true;
        					ifEliminate = true;
        				}
        			}
        		}
        		if(!ifEliminate) {
        			// No component will be eliminated, update pi[j]
        			for(int j=0; j<numComp; j++) {
        				if(pi[c][j]>0) {
//        					pi[c][j] = Math.max(0, sumR[j]-currAlpha[c]);
        					pi[c][j] = Math.max(0, sumR[j]);
        				}
        			}
        		} else {
        			// Select eliminated binding components following two criteria:
        			// 1. sumR < currAlpha		
        			// 2. no components eliminated in reassigned responsibility range (set as initialize fuzziness 95% interval)
        			boolean[] exclusionZone = new boolean[currRegion.getWidth()];
        			for(int j=0; j<numComp; j++	) {
        				if(isEliminated[j]) {
        					if(!exclusionZone[mu[c][j]-currRegion.getStart()]) {
        						int start = Math.max(0, mu[c][j]-currRegion.getStart()-maxIR[c][j]/2);
        						int end = Math.min(currRegion.getWidth()-1, mu[c][j]-currRegion.getStart()+maxIR[c][j]/2);
        						for(int z=start; z<=end; z++) {
        							exclusionZone[z] = true;
        						}
        					} else {
        						isEliminated[j] = false;
        					}
        				}
        			}
        			
        			// Eliminate selected eliminated binding components
        			// Responsibilities will be redistributed in the E step
        			for(int j=0; j<numComp; j++) {
        				if(isEliminated[j]) {
        					pi[c][j]=0.0; sumR[j]=0.0;
        					if(pairIndexAroundMu.get(c).get(j).car()!=-1)
        						for(int i=pairIndexAroundMu.get(c).get(j).car(); i<=pairIndexAroundMu.get(c).get(j).cdr(); i++)
        							rBind[c][j][i] = 0;
        				}
        			}
        			// Re-estimate pi values for non-eliminated components using the current responsibility assignments
        			for(int j=0; j<numComp; j++) {
        				if(!isEliminated[j])
        					pi[c][j] = Math.max(0, sumR[j]); // Q: I think here should be sumR[j]-currAlpha[c].
        			}
        			componentEliminated = true;
        		}
        		
        		// 	if there is no component eliminated because of resp < currAlpha, use exclusion zone to exclude nucleosome
        		// 	note: I don't want exclusion zone and alpha eliminate nucleosome at the same round because I want to make sure before
        		// any component is eliminated due to exclusion zone, all responsibilities have been assigned to components
        		if(t>=semconfig.ALPHA_ANNEALING_ITER && !ifEliminate && !mode.equals(EMmode.NORMAL)) {
            		boolean[] exclusionZone = new boolean[currRegion.getWidth()];
        			int exclusion = mode.equals(EMmode.ALTERNATIVE)? semconfig.getAlternativeExclusionZone() : semconfig.getConsensusExclusionZone();
        			//sort component according to their pi then get the excluded nucleosome with the lowest pi
        			ArrayIndexComparator comparator = new ArrayIndexComparator(pi[c]);
        			Integer[] indexes = comparator.createIndexArray();
        			Arrays.sort(indexes, comparator);
        			int minExIndex = -1;
        			for(int j : indexes) {
        				if(pi[c][j]>0) {
        					int start = Math.max(0, mu[c][j]-currRegion.getStart()-exclusion/2);
        					int end = Math.min(currRegion.getWidth()-1, mu[c][j]-currRegion.getStart()+exclusion/2);
        					boolean isOverlap = false;
        					for(int i=start; i<=end; i++) {
        						if(exclusionZone[i])
        							isOverlap = true;
        					}
        					if(!isOverlap) {
        						for(int z=start; z<end; z++)
        							exclusionZone[z] = true;
        					} else {
        						minExIndex = j;
        					}
        				}
        			}
        			// Eliminate the nucleosome in exclusion zone with the lowest pi
        			if(minExIndex!=-1) {
        				pi[c][minExIndex]=0.0; sumR[minExIndex]=0.0;
    					if(pairIndexAroundMu.get(c).get(minExIndex).car()!=-1)
    						for(int i=pairIndexAroundMu.get(c).get(minExIndex).car(); i<=pairIndexAroundMu.get(c).get(minExIndex).cdr(); i++)
    							rBind[c][minExIndex][i] = 0;
            			componentEliminated = true;
        			}
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
        	}
        	
    		/////////////
        	//Anneal alpha Q: Looks like annealling means reduce alpha each turn, I don't know why
        	//////////////
    		if (t >semconfig.EM_ML_ITER && t <= semconfig.ALPHA_ANNEALING_ITER)
    			alphaCoefficient = (double)(t-semconfig.EM_ML_ITER)/(double)(semconfig.ALPHA_ANNEALING_ITER-semconfig.EM_ML_ITER);
    		else if(t > semconfig.ALPHA_ANNEALING_ITER)
    			alphaCoefficient = 1;
        	
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
        			for(int j=0; j<numComp; j++) {
        				if(pi[c][j]>0 && pairIndexAroundMu.get(c).get(j).car()!=-1)
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
//        			double sum_log_pi = 0;
        			double sum_alpha = 0;
        			for(int j=0; j<numComp; j++) {
        				if(pi[c][j]>0.0) {
//        					sum_log_pi += Math.log(pi[c][j])/semconfig.LOG2;
        					sum_alpha += currAlpha[c][j] * Math.log(pi[c][j])/semconfig.LOG2;
        				}
        			}
        			// Position prior (Assumption3: Bernouli prior on nucleosome positions)
        			double sum_pos_prior = 0;
        			if(atacPrior!=null && semconfig.useAtacPrior())
        				for(int j=0; j<numComp; j++)
        					sum_pos_prior += atacPrior[c][mu[c][j]-currRegion.getStart()];
        			
        			LP += -sum_alpha + sum_pos_prior;
        		}
        		LAP = LL + LP;
        		
//        		System.out.println("EM: "+t+"\t"+LAP+"\t("+nonZeroComps+" non-zero components).");
        	}
        	
			//monitor: count time
        	timer.end("LL");
        	
        	//Is current state equivalent to the last?
            if(((numConditions>1 && t>semconfig.POSPRIOR_ITER) || (numConditions==1 && t>semconfig.ALPHA_ANNEALING_ITER)) && 
            		lastEquivToCurr())
            	stateEquivCount++;
            else
            	stateEquivCount=0;
        	
        	//Tick the clock forward
    		if(!componentEliminated)
    			t++;
    		iter++;
    		
    		//Save parameters if plotEM
    		if(plotEM) {
    			EMStepPlotter.excute(mu, pi, fuzz, tau, iter, t);
    		}
    		
    		//Are there any components sharing the same location from the same experiment condition?
    		boolean componentOverlapping = false;
    		for(int c=0; c<numConditions; c++) {
    			List<Integer> compLoc = new ArrayList<Integer>();
    			for(int j=0; j<numComp; j++) {
    				if(pi[c][j]>0) {
    				if(!compLoc.contains(mu[c][j])) {
    					compLoc.add(mu[c][j]);
    				} else {
    					componentOverlapping = true;
    					break;
    				}
    				}
    			}
    			if(componentOverlapping) {
    				break;
    			}
    		}
    		
    		 ////////////
          	//Check Stopping condition TODO: I don't know whether it is right to use Math.abs
          	////////////   		
            if (nonZeroComps>0 && (componentEliminated || componentOverlapping || (numConditions>1 && t<=semconfig.POSPRIOR_ITER) || (numConditions==1 && t<=semconfig.ALPHA_ANNEALING_ITER) 
            		|| (semconfig.CALC_LL && Math.abs(LAP-lastLAP)>Math.abs(semconfig.EM_CONVERGENCE*lastLAP)) || (stateEquivCount<semconfig.EM_STATE_EQUIV_ROUNDS))){
            	if(t==semconfig.MAX_EM_ITER) {
	            	System.out.println("\tcriteria 1: "+(numConditions>1 && t<=semconfig.POSPRIOR_ITER));
	            	System.out.println("\tcriteria 2: "+(numConditions==1 && t<=semconfig.ALPHA_ANNEALING_ITER));
	            	System.out.println("\tcriteria 3: "+(semconfig.CALC_LL && Math.abs(LAP-lastLAP)>Math.abs(semconfig.EM_CONVERGENCE*lastLAP)));
	            	System.out.println("\tcriteria 4: "+(stateEquivCount<semconfig.EM_STATE_EQUIV_ROUNDS));
	            	if(semconfig.CALC_LL && Math.abs(LAP-lastLAP)>Math.abs(semconfig.EM_CONVERGENCE*lastLAP)){
	            		System.out.println("\t\tlastLAP: "+lastLAP);
	            		System.out.println("\t\tLAP: "+LAP);
	            	}
            	}
                copyStateToLast();
                lastLAP = LAP;
                continue;
            }else{
            	copyStateToLast();
            	lastLAP = LAP;
            	if(semconfig.isVerbose())
            		System.err.println("\tRegTrain:"+trainingRound+"\t"+currRegion.getLocationString()+"\twidth:"+currRegion.getWidth()+"\tnoEliminated:"+t+"\toverall:"+iter+"\t"+nonZeroComps);
            	break;
            }
        } //LOOP: Run EM while not converged
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
    			lastSumResp[c][j] = sumResp[c][j];
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
    				compPosEqual = compPosEqual && Math.abs(mu[c][j] - lastMu[c][j])<3;
    			}}
    	}else{
    		compPosEqual=false;
    	}
    	boolean piBindEquivalent=true;
    	for(int c=0; c<numC; c++)
			for(int j=0; j<pi[c].length; j++){if(pi[c][j]>0){
				piBindEquivalent = piBindEquivalent && (Math.abs(pi[c][j]-lastPi[c][j])<semconfig.EM_STATE_EQUIV_PI_THRES);
			}}
    	boolean fuzzBindEquivalent=true;
    	for(int c=0; c<numC; c++)
    		for(int j=0; j<pi[c].length; j++) {
    			if(pi[c][j]>0)
    				fuzzBindEquivalent = fuzzBindEquivalent && (Math.abs(fuzz[c][j]-lastFuzz[c][j]) < semconfig.EM_STATE_EQUIV_FUZZ_THRES*lastFuzz[c][j]);
//    			if(Math.abs(fuzz[c][j]-lastFuzz[c][j]) >= semconfig.EM_STATE_EQUIV_FUZZ_THRES*lastFuzz[c][j]) {
//    				System.out.println("\t\t\tunequal fuzziness: "+fuzz[c][j]+"-"+lastFuzz[c][j]);
//    			}
    		}
    	boolean tauBindEquivalent=true;
    	for(int c=0; c<numC; c++)
    		for(int j=0 ;j<pi[c].length; j++) {
    			if(pi[c][j]>0)
    				for(int t=0; t<tau[c][j].length; t++)
    					tauBindEquivalent = tauBindEquivalent && (Math.abs(tau[c][j][t] - lastTau[c][j][t]) < semconfig.EM_STATE_EQUIV_TAU_THRES);
    		}
//    	System.out.println("\t\tnumCompEqual: "+numCompEqual);
//    	System.out.println("\t\tcompPosEqual: "+compPosEqual);
//    	System.out.println("\t\tpiBindEquivalent: "+piBindEquivalent);
//    	System.out.println("\t\tfuzzBindEquivalent: "+fuzzBindEquivalent);
//    	System.out.println("\t\ttauBindEquivalent: "+tauBindEquivalent);
		return numCompEqual && compPosEqual && piBindEquivalent && fuzzBindEquivalent && tauBindEquivalent;
    }
    
    private class PairwiseKey{
    	int cond1;
    	int index1;
    	int cond2;
    	int index2;
    	public PairwiseKey(int c1, int i1, int c2, int i2) {
    		cond1 = c1;
    		index1 = i1;
    		cond2 = c2;
    		index2 = i2;
    	}
    	
    	public int firstCond() {return cond1;}
    	public int firstIndex() {return index1;}
    	public int secondCond() {return cond2;}
    	public int secondIndex() {return index2;}
    	
    	public int hashCode() {
    		int code = 17;
    		code = code*31 + cond1;
    		code = code*31 + index1;
    		code = code*31 + cond2;
    		code = code*31 + index2;
    		return code;
    	}
    	
    	public boolean equals(Object o) {
    		if(!(o instanceof PairwiseKey)) {
    			return false;
    		}
    		PairwiseKey p = (PairwiseKey) o;
    		return p.firstCond()==cond1 && p.firstIndex()==index1 && p.secondCond()==cond2 && p.secondIndex()==index2;
    	}
    	
    	public String toString() {
    		return cond1 + " " + index1 + " " + cond2 + " " + index2;
    	}
    }
    
    private class ArrayIndexComparator implements Comparator<Integer>
    {
        private final double[] array;

        public ArrayIndexComparator(double[] array)
        {
            this.array = array;
        }

        public Integer[] createIndexArray()
        {
            Integer[] indexes = new Integer[array.length];
            for (int i = 0; i < array.length; i++)
            {
                indexes[i] = i; // Autoboxing
            }
            return indexes;
        }

        @Override
        public int compare(Integer index1, Integer index2)
        {
             // Autounbox from Integer to int to use as array indexes
            return -Double.compare(array[index1], array[index2]);
        }
    }
}














































