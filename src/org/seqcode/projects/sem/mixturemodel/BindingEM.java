package org.seqcode.projects.sem.mixturemodel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
	protected int stateEquivCount = 0;
	
	public BindingEM(SEMConfig s, ExperimentManager eMan, BindingManager bMan, HashMap<ExperimentCondition, BackgroundCollection> condBacks, int numPotReg) {
		semconfig = s;
		manager = eMan;
		bindingManager = bMan;
		conditionBackgrounds = condBacks;
		numConditions = manager.getNumConditions();
		numPotentialRegions = (double)numPotReg;

	}
	
	public List<List<BindingComponent>> train(List<List<StrandedPair>> signals,
												Region w,
												List<NoiseComponent> noise,
												List<List<BindingComponent>> comps,
												int numComp,
												double[][] atacPrior,
												int trainingRound
												) {
		components = comps;
		this.noise = noise;
		numComponents = numComp;
		this.atacPrior = atacPrior;
		this.trainingRound = trainingRound;
		this.plotSubRegion = plotSubRegion;
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
			
			// Set maximum alphas
			alphaMax[c] = semconfig.getFixedAlpha()>0 ? semconfig.getFixedAlpha() :
					semconfig.getAlphaScalingFactor() * (double)conditionBackgrounds.get(cond).getMaxThreshold('.');
					
			// Load Read pairs (merge from all replicates)
			List<StrandedPair> pairs = new ArrayList<StrandedPair>();
			for(ControlledExperiment rep: cond.getReplicates())
				pairs.addAll(signals.get(rep.getIndex()));
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
	private void EM_MAP(Region currRegion) {
		
		//monitor code: Mark EM_MAP start
		System.out.println("EM_MAP start.......................");
		
		int numComp = numComponents;
		double [][] totalResp = new double[numConditions][];
		int regStart = currRegion.getStart();
		
		// Variables for tracking mu maximization. Defined early to avoid memory assignment during main EM loop.
		double[][][] muSums = new double[numConditions][numComp][];		//Results of mu maximization summations for individual components across genome
        int[][] muSumStarts = new int[numConditions][numComp]; //Start positions of muSum arrays (start of maximization window).
        int[][] muSumWidths = new int[numConditions][numComp]; //Effective widths of muSum arrays (width of maximization window).
        int[][] muSumMaxPos = new int[numConditions][numComp]; //Positions of maxima in mu maximization summations
        int[] muJoinClosestComps = new int[numConditions]; //Indices of nearest components in other conditions
        boolean[] muJoinSharedBetter = new boolean[numConditions]; //Indicator that sharing components across conditions is better than not
        int[][] newMu = new int[numConditions][numComponents];// mu update
        
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
        while(t<semconfig.MAX_EM_ITER) {
        	
        	//monitor code: show EM round
        	System.err.println("EM round "+iter);
        	
    		////////
    		//E-step
    		////////
        	for(int c=0; c<numConditions; c++) {
        		int numPairs = hitNum[c];
        		ExperimentCondition cond = manager.getIndexedCondition(c);
        		
    			// Load binding fragment size PDF cache (indexed by type index)
    			Map<Integer, List<Double>> fragSizePDF = bindingManager.getCachePDF(cond);
    			
    			// Compute H and N function
    			double[][] hc = new double[numComp][numPairs];
    			double[] nc = new double[numPairs];
    			for(int i=0; i<numPairs; i++) {
    				for(int j=0; j<numComp; j++) {
    					double fuzzProb = BindingModel.probability(fuzz[c][j], mu[c][j]-hitPos[c][i]);
    					double fragSizeProb = 0;
    					for(int index: fragSizePDF.keySet()) {
    						fragSizeProb += tau[c][j][index] * fragSizePDF.get(index).get(hitSize[c][i]);
    					}
    					hc[j][i] = fuzzProb * fragSizeProb;
    				}
    				nc[i] = noise.get(c).score(hitPos[c][i], hitSize[c][i], repIndices[c][i]);
    			}
    			h[c] = hc;
    			n[c] = nc;
        		
        		// Compute responsibilities
        		for(int i=0; i<numPairs; i++) {
        			totalResp[c][i] = 0;
        			for(int j=0; j<numComp; j++) {
        				if(pi[c][j]>0) {
        					rBind[c][j][i] = h[c][j][i]*pi[c][j];
        					totalResp[c][i] += rBind[c][j][i];
        				}
        			}
        			rNoise[c][i] = n[c][i] * piNoise[c];
        			totalResp[c][i] += rNoise[c][i];
        		}
        		// Normalize responsibilities
        		for(int i=0; i<numPairs; i++) {
        			for(int j=0; j<numComp; j++) {
        				if(pi[c][j]>0) {
        					rBind[c][j][i] /= totalResp[c][i];
        				}
        			}
        			rNoise[c][i]/=totalResp[c][i];
        		}
        	}

    		/////////////////////
    		//M-step: maximize mu (positions), fuzz (fuzziness), tau (fragment size subtype), pi (strength)
    		/////////////////////
        	if(numConditions>1 && t==semconfig.ALPHA_ANNEALING_ITER)
        		for(int c=0; c<numConditions; c++)
        			for(int j=0; j<numComp; j++)
        				if(pi[c][j]>0)
        					muSums[c][j] = new double[semconfig.EM_MU_UPDATE_WIN*2];
        	// Part1.1: Maximize mu
        	for(int c=0; c<numConditions; c++) {
        		int numPairs = hitNum[c];
        		for(int j=0; j<numComp; j++) {
        			if(pi[c][j]>0) {
        				int start = Math.max(mu[c][j]-semconfig.EM_MU_UPDATE_WIN, regStart);
        				int end = Math.min(currRegion.getEnd(), mu[c][j]+semconfig.EM_MU_UPDATE_WIN);
        				// Assign special variables
        				if(numConditions>1 && t>semconfig.ALPHA_ANNEALING_ITER) { // I don't know why needs this if-
        					muSumStarts[c][j] = start; muSumWidths[c][j] = end - start;
        				}
        				// Score the current window
        				double currScore=0, maxScore=-Double.MAX_VALUE;
        				int maxPos = 0;
        				for(int x=start; x<end; x++) {
        					currScore=0;
        					for(int i=0; i<numPairs; i++) {
    							int dist = hitPos[c][i] - x;
    							currScore += (rBind[c][j][i]*hitCounts[c][i]) * BindingModel.logProbability(fuzz[c][j], dist);
        					}
        					if(atacPrior!=null && semconfig.useAtacPrior())
        						currScore += atacPrior[c][x-regStart];
        					
        					if(numConditions>1 && t>semconfig.ALPHA_ANNEALING_ITER) // Save the score
        						muSums[c][j][x-start] = currScore;
        					
        					if(currScore>maxScore) {
        						maxPos = x;
        						maxScore = currScore;
        					}
        				}
        				muSumMaxPos[c][j] = maxPos;
        				mu[c][j] = maxPos;
        			}
        		}
        	}
        	
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
        	

        	// Part2: Maximize fuzziness
        	for(int c=0; c<numConditions; c++) {
        		int numPairs = hitNum[c];
        		for(int j=0; j<numComp; j++) {if(pi[c][j]>0) {
        			double posVar = 0;
        			double sumWeight = 0;
        			for(int i=0; i<numPairs; i++) 
        				sumWeight += rBind[c][j][i] * hitCounts[c][i];
        			for(int i=0; i<numPairs; i++)
        				posVar += (rBind[c][j][i] * hitCounts[c][i]) * Math.pow(hitPos[c][i]-mu[c][j], 2);
        			posVar /= sumWeight;
        			fuzz[c][j] = Math.sqrt(posVar);
        		}
        		}
        	}

        	// Part3: Maximize tau
        	for(ExperimentCondition cond: manager.getConditions()) {
        		int c = cond.getIndex();
        		int numPairs = hitNum[c];
        		int numSubtypes = bindingManager.getBindingSubtypes(cond).size();
        		for(int j=0; j<numComp; j++) {if(pi[c][j]>0) {
        			tau[c][j] = new double[numSubtypes];
        			double probSum = 0;
        			// Compute tau probability of each fragment size subtype for each component
        			for(BindingSubtype b: bindingManager.getBindingSubtypes(cond)) {
        				for(int i=0; i<numPairs; i++) {
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
        	
        	
        	// Part4: Maximize pi
        	boolean componentEliminated = false;
        	for(int c=0; c<numConditions; c++) {
        		int numPairs = hitNum[c];
        		// Maximize pi
        		double[] sumR = new double[numComponents];
        		for(int j=0; j<numComp; j++) {
        			if(pi[c][j]>0) {
        				for(int i=0; i<numPairs; i++)
        					sumR[j] += rBind[c][j][i] * hitCounts[c][i];
        			}
        		}
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
        			for(int i=0; i<numPairs; i++)
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
        	
        	// Non-zero components count
        	int nonZeroComps = 0;
        	for(int c=0; c<numConditions; c++)
        		for(int j=0; j<numComp; j++)
        			if(pi[c][j]>0)
        				nonZeroComps++;
        	
        	////////////
        	//Compute LL
        	////////////
        	LAP=0;
        	if(semconfig.CALC_LL) {
        		// Log-likelihood calculation
        		double LL = 0;
        		for(int c=0; c<numConditions; c++) {
        			int numPairs = hitNum[c];
        			for(int i=0; i<numPairs; i++) {
        				// for each read pair, each component will give a conditional prob or bg prob
        				double j_sum = 0;
        				for(int j=0; j<numComp; j++) {
        					if(pi[c][j]>0.0)
        						j_sum += Math.log(rBind[c][j][i]/semconfig.LOG2); // Q: In this part, rBind has been normalized per read, I think it cannot be used to calculate log likelihood.
        				}
        				j_sum += Math.log(rNoise[c][i])/semconfig.LOG2;
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
        		
        		System.out.println("EM: "+t+"\t"+LAP+"\t("+nonZeroComps+" non-zero components).");
        	}
        	
        	//Tick the clock forward
    		if(!componentEliminated)
    			t++;
    		iter++;
    		
    		 ////////////
          	//Check Stopping condition
          	////////////   		
//            if (nonZeroComps>0 && ((numConditions>1 && t<=semconfig.POSPRIOR_ITER) || (numConditions==1 && t<=semconfig.ALPHA_ANNEALING_ITER) || (semconfig.CALC_LL && Math.abs(LAP-lastLAP)>semconfig.EM_CONVERGENCE) || stateEquivCount<semconfig.EM_STATE_EQUIV_ROUNDS)){
//                copyStateToLast();
//                lastLAP = LAP;
//                continue;
//            }else{
//            	copyStateToLast();
//            	lastLAP = LAP;
//            	//if(config.isVerbose())
//            		//System.err.println("\tRegTrain:"+trainingRound+"\t"+currRegion.getLocationString()+"\t"+currRegion.getWidth()+"\t"+t+"\t"+iter+"\t"+nonZeroComps);
//            	break;
//            }
        } //LOOP: Run EM while not converged
        
        //monitor code: show binding component information after EM loop
        for(int c=0; c<numConditions; c++) {
        	for(int j=0; j<numComp; j++) {
        		if(pi[c][j]>0) {
        			System.out.println("\tBinding Component"+j);
        			System.out.println("\t\tpi: "+pi[c][j]);
        			System.out.println("\t\tmu: "+mu[c][j]);
        			System.out.println("\t\tfuzziness: "+fuzz[c][j]);
        			System.out.println("\t\ttau: "+Arrays.toString(tau[c][j]));
        		}
        	}
        }
        
    	//monitor code: exit here
    	System.exit(1);
    	
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












































