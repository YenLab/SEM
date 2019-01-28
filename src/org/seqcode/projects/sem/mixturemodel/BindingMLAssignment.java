package org.seqcode.projects.sem.mixturemodel;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.seqcode.deepseq.StrandedPair;
import org.seqcode.deepseq.events.BindingModelPerBase;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.deepseq.experiments.Sample;
import org.seqcode.deepseq.stats.BackgroundCollection;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.sequence.SequenceUtils;
import org.seqcode.projects.sem.framework.SEMConfig;
import org.seqcode.projects.sem.events.BindingManager;
import org.seqcode.projects.sem.events.BindingModel;
import org.seqcode.projects.sem.events.EventsConfig;
import org.seqcode.projects.sem.events.BindingEvent;

/**
 * BindingMLAssignment: Maximum likelihood assignment of reads to a configuration of binding components.
 * @author Jianyu Yang
 */
public class BindingMLAssignment {

	protected ExperimentManager manager;
	protected BindingManager bindingManager;
	protected ExptConfig econfig;
	protected EventsConfig evconfig;
	protected SEMConfig config;
	protected List<BindingComponent> components;
	protected List<NoiseComponent> noise;
	protected int numComponents;	//Assumes the same number of active+inactive components in each condition
	protected int numConditions;
	protected HashMap<ExperimentCondition, BackgroundCollection> conditionBackgrounds; 
	
	// EM VARIABLES
	protected double[][]	sigHitCounts;		// Hit weights
	protected int[][]		sigHitPos;			// Hit positions
	protected int[][]		sigHitSize;			// Hit fragment size
	protected int[]			sigHitNum;			// Number of hits in each condition
	protected int[][]		sigRepIndices;		// Index of replicate for each hit
	protected double[][]	ctrlHitCounts;		// Hit weights
	protected int[][]		ctrlHitPos;			// Hit positions
	protected int[][]		ctrlHitSize;		// Hit fragment size
	protected int[]			ctrlHitNum;			// Number of hits in each condition
	protected int[][]		ctrlRepIndices;		// Index of replicate for each hit
	protected double[][][]	h;					// H function (binding component probability per read)
	protected double[][]	n;					// N function (noise component probability per read)
	protected double[][][]	rBindSig;			// Binding component responsibilities (signal reads)
	protected double[][]	rNoiseSig;			// Noise component responsibilities (control reads)
	protected double[][][]	rBindCtrl;			// Binding component responsibilities (control reads)
	protected double[][]	rNoiseCtrl;			// Noise component responsibilities (control reads)
	protected double[][]	pi;					// pi: emission probabilities for binding components
	protected double[]		piNoise;			// piNoise: emission probabilities for noise components (fixed)
	protected int[][]		mu;					// mu: positions of the binding components
	protected double[][]	fuzz;				// &fuzz: fuzziness of the binding components (indexed by condition & binding component index)
	protected double[][][]	tau;				// &tau: fragment size subtype probabilities (indexed by subtype index)
	protected double[][]	compLL;				// Log-likelihood for each component in each condition
	protected double[][][]	lastRBind;			// Last responsibilities (monitor convergence)
	protected double[][]	lastPi;				// Last Pi (monitor convergence)
	protected int[][]		lastMu;				// Last positions (monitor convergence)
	protected double[][]	lastFuzz;			// &Last fuzziness (monitor convergence)
	protected double[][][]	lastTau;				// &Last tau (monitor convergence)
	protected double[][]	tmp_pi;				// pi used in ML calc
	protected double[][][]	tmp_h;				// h used in ML calc
	protected double[][][]	tmp_rBindSig;		// rBindSig used in ML calc
	protected double[][]	tmp_rNoiseSig;		// rNoiseSig used in ML calc
	protected double[]		tmp_piNoise;		// piNoise used in ML calc
	
	protected double[]		sigRepHitCountTotals;	// Hit count totals counted by replicate
	protected double[]		uniformRepHitCountTotals;	// Hit count totals by replicate if signal read counts were distributed uniformly (used only if there is no control) 
	protected double 		numPotentialRegions;	// Q: Why double?
	
	/**
	 * Constructor
	 */
	public BindingMLAssignment(ExptConfig econ, EventsConfig evcon, SEMConfig s, ExperimentManager eMan, BindingManager bMan, HashMap<ExperimentCondition, BackgroundCollection> condBacks, int numPotReg) {
		config = s;
		evconfig = evcon;
		econfig = econ;
		manager = eMan;
		bindingManager = bMan;
		conditionBackgrounds = condBacks;
		numConditions = manager.getNumConditions();
		numPotentialRegions = (double)numPotReg;
	}
	
	/**
	 * ML assignment
	 * 
	 * Takes as input a SINGLE list of binding components (i.e. a single configuration)
	 * Returns lists of binding events indexed by condition.
	 * Pi values are calculated from the signal hits and then those same components are directly applied to the control hits.
	 */
	public List<BindingEvent> assign(List<List<StrandedPair>> signals,
									 List<List<StrandedPair>> controls,
									 Region w,
									 List<NoiseComponent> noise,
									 List<BindingComponent> comps,
									 int numComp) {
		components = comps;
		this.noise = noise;
		numComponents = numComp;
		// Matrix initializations
		sigHitCounts = new double[numConditions][];		//Hit weights
		sigHitPos = new int[numConditions][];			//Hit positions
		sigHitSize = new int[numConditions][];			//Hit fragment size
		sigHitNum = new int[numConditions];				//Number of hits in each condition
		sigRepIndices = new int[numConditions][];		//Index of replicate for the hit
		ctrlHitCounts = new double[numConditions][];	//Hit weights
		ctrlHitPos = new int[numConditions][];			//Hit positions
		ctrlHitSize = new int[numConditions][];			//Hit fragment size
		ctrlHitNum = new int[numConditions];			//Number of hits in each condition
		ctrlRepIndices = new int[numConditions][];		//Index of replicate for the hit
		h = new double[numConditions][][];				//H function (binding component probability per read)
		n = new double[numConditions][];				//N funciton (noise component probability per read)
		rBindSig = new double[numConditions][][];		//Binding component responsibilities (signal reads)
		rNoiseSig = new double[numConditions][];		//Noise component responsibilities (signal reads)
		rBindCtrl = new double[numConditions][][];		//Binding component responsibilities (control reads)
		rNoiseCtrl = new double[numConditions][];		//Noise component responsibilities (control read)
		pi = new double[numConditions][numComponents];	//pi: emission probabilities for binding components
		piNoise = new double[numConditions];			//piNoise: emission probabilities for noise components (fixed)
		mu = new int[numConditions][numComponents];		//mu: positions of the binding components
		fuzz = new double[numConditions][numComponents];	// &fuzz: fuzziness of the binding components
		tau = new double[numConditions][numComponents][];				// &tau: fragment size subtype probabilities (indexed by subtype index)
		compLL = new double[numConditions][numComponents];	//Log-likelihood for each component in each condition
		sigRepHitCountTotals = new double[manager.getReplicates().size()];	//Hit count totals counted by replicate (for convenience)
		uniformRepHitCountTotals = new double[manager.getReplicates().size()];	//Hit count totals by replicate if reads were distributed uniformly
		//Monitor state convergence using the following last variables
		lastRBind = new double[numConditions][][];
		lastPi = new double[numConditions][numComponents];
		lastMu = new int[numConditions][numComponents];
		lastFuzz = new double[numConditions][numComponents];
		lastTau = new double[numConditions][numComponents][];
		//Temporary variables
		tmp_pi = new double[numConditions][numComponents];	// pi used in ML calc
		tmp_rBindSig = new double[numConditions][][];		// rBindSig used in ML calc
		tmp_rNoiseSig = new double[numConditions][];		// rNoiseSig used in ML calc
		tmp_piNoise = new double[numConditions];			// piNoise used in ML calc
		tmp_h = new double[numConditions][][];				// H function used in ML calc
		
		//Initializing data structures
		for(ExperimentCondition cond: manager.getConditions()) {
			int c = cond.getIndex();
			
			// Load binding fragment size PDF cache
			Map<Integer, List<Double>> fragSizePDF = bindingManager.getCachePDF(cond);
			
			//Load reads (merge from all replicates) & Indexed by condition
			List<StrandedPair> sigPairs = new ArrayList<StrandedPair>();
			List<StrandedPair> ctrlPairs = new ArrayList<StrandedPair>();
			for(ControlledExperiment rep: cond.getReplicates()) {
				sigPairs.addAll(signals.get(rep.getIndex()));
				if(controls.get(rep.getIndex())!=null)
					ctrlPairs.addAll(controls.get(rep.getIndex()));
			}
			sigHitNum[c] = sigPairs.size();
			ctrlHitNum[c] = ctrlPairs.size();
			
			//Count total weights each replicate for convenience & Indexed by replicate
			for(ControlledExperiment rep: cond.getReplicates()) {
				sigRepHitCountTotals[rep.getIndex()] = 0;
				for(StrandedPair s: signals.get(rep.getIndex()))
					sigRepHitCountTotals[rep.getIndex()] += s.getWeight();
				uniformRepHitCountTotals[rep.getIndex()] = (((rep.getSignal().getHitCount()*(1-rep.getSignalVsNoiseFraction())) /
						econfig.getMappableGenomeLength())*(double)w.getWidth())/rep.getControlScaling();
			}
			
			//Load replicate index for each read
			sigRepIndices[c] = new int[sigHitNum[c]];
			ctrlRepIndices[c] = ctrlHitNum[c]==0 ? null : new int[ctrlHitNum[c]];
			int ys=0, yc=0, z=0;
			for(ControlledExperiment rep: cond.getReplicates()) {
				z=0;
				while(z<signals.get(rep.getIndex()).size()) {
					sigRepIndices[c][ys] = rep.getIndex();
					z++; ys++;
				}
				z=0;
				while(z<controls.get(rep.getIndex()).size()) {
					ctrlRepIndices[c][yc] = rep.getIndex();
					z++; yc++;
				}
			}
			
			//Load signal read info
			sigHitCounts[c] = new double[sigHitNum[c]];
			sigHitPos[c] = new int[sigHitNum[c]];
			sigHitSize[c] = new int[sigHitNum[c]];
			for(int i=0; i<sigHitNum[c]; i++) {
				sigHitPos[c][i] = sigPairs.get(i).getMidpoint().getLocation();
				sigHitSize[c][i] = sigPairs.get(i).getFragmentSize();
				sigHitCounts[c][i] = 1;
			}
			
			//Load control read info
			if(ctrlHitNum[c]>0) {
				ctrlHitCounts[c] = new double[ctrlHitNum[c]];
				ctrlHitPos[c] = new int[ctrlHitNum[c]];
				ctrlHitSize[c] = new int[ctrlHitNum[c]];
				for(int i=0; i<ctrlHitNum[c]; i++) {
					ctrlHitPos[c][i] = ctrlPairs.get(i).getMidpoint().getLocation();
					ctrlHitSize[c][i] = sigPairs.get(i).getFragmentSize();
					ctrlHitCounts[c][i] = 1;
				}
			}
			
			//Load pi for binding component
			for(int j=0; j<numComp; j++) {
				BindingComponent comp = components.get(j);
				pi[c][j] = comp.getPi();
			}
			//Load pi for noise components
			piNoise[c] = noise.get(c).getPi();
			
			//Load binding component positions, fuzziness, tau
			for(int j=0; j<numComp; j++) {
				mu[c][j] = components.get(j).getPosition();
				fuzz[c][j] = components.get(j).getFuzziness();
				tau[c][j] = components.get(j).getTau();
				lastTau[c][j] = new double[fragSizePDF.keySet().size()];
			}
			
			//Initialize responsibility functions
			double[][] hc = new double[numComp][sigHitNum[c]];
			double[][] thc = new double[numComp][sigHitNum[c]];
			double[] nc = new double[sigHitNum[c]];
			for(int i=0; i<sigHitNum[c]; i++) {
				for(int j=0; j<numComp; j++) {
					double fuzzProb = BindingModel.probability(fuzz[c][j], mu[c][j]-sigHitPos[c][i]);
					double fragSizeProb = 0;
					for(int index: fragSizePDF.keySet()) {
						fragSizeProb += tau[c][j][index] * fragSizePDF.get(index).get(sigHitSize[c][i]);
					}
					hc[j][i] = fuzzProb * fragSizeProb;
					thc[j][i] = fuzzProb * fragSizeProb;
				}
				nc[i] = noise.get(c).scorePosition(sigHitPos[c][i], sigRepIndices[c][i]);
			}
			h[c] = hc;
			n[c] = nc;
			tmp_h[c] = thc;
			
			rBindSig[c] = new double[numComp][sigHitNum[c]];
			rNoiseSig[c] = new double[sigHitNum[c]];
			lastRBind[c] = new double[numComp][sigHitNum[c]];
			tmp_rBindSig[c] = new double[numComp][sigHitNum[c]];
			tmp_rNoiseSig[c] = new double[sigHitNum[c]];
		}//End of data structure initialization
		
        //////////
        // Run ML steps
        //////////
		ML(w);
		
        //////////
        // Assign ML result to BindingEvents
        //////////
		List<BindingEvent> events = new ArrayList<BindingEvent>();
		for(int j=0; j<numComp; j++) {
			BindingEvent event = new BindingEvent(components.get(j).getCoord(), w);
			
			for(ExperimentCondition cond: manager.getConditions()) {
				ArrayList<Sample> controlsSeen = new ArrayList<Sample>();
				boolean uniformBackAdded=false;
				double condSigResp = 0.0, condCtrlResp = 0.0;
				int c = cond.getIndex();
				for(ControlledExperiment rep: cond.getReplicates()) {
					int r = rep.getIndex();
					double repSigResp = 0.0, repCtrlResp = 0.0;
					if(pi[c][j]>0) {
						double scount = 0;
						for(int i=0; i<sigHitNum[c]; i++)
							if(sigRepIndices[c][i]==r)
								scount += sigHitCounts[c][i] * rBindSig[c][j][i];
						repSigResp += scount;
						condSigResp += scount;
						
						double ccount = 0;
						if(rep.hasControl()) {
							for(int i=0; i<ctrlHitNum[c];i++)
								if(ctrlRepIndices[c][i]==r)
									ccount += ctrlHitCounts[c][i] * rBindCtrl[c][j][i];
							repCtrlResp += ccount;
							if(!controlsSeen.contains(rep.getControl()))
								condCtrlResp+=ccount;
							controlsSeen.add(rep.getControl());
						}else {  //If there is no control channel, assign pseudo-control counts as if the noise reads in the IP channel were distributed perfectly uniformly
							repCtrlResp = uniformRepHitCountTotals[rep.getIndex()]*pi[c][j];
							if(!uniformBackAdded)
								condCtrlResp += repCtrlResp;
							uniformBackAdded=true;
						}
					}
					event.setRepSigHits(rep, repSigResp);
					event.setRepCtrlHits(rep, repCtrlResp);
				}
				event.setCondSigHits(cond, condSigResp);
				event.setCondCtrlHits(cond, condCtrlResp);
				
			}
		}
		return events;
	}//end of EMTrain method
	
	/**
	 * Core EM iterations with sparse prior (component elimination) & multi-condition positional priors.
	 * Assumes H function, pi, and responsibilities have all been initialized
	 */
	private void ML(Region currRegion) {
		int numComp = numComponents;
		double[][] totalRespSig = new double[numConditions][];
		double[][] totalRespCtrl = new double[numConditions][];
		
		//Initialize responsibilities
		for(int c=0; c<numConditions; c++) {
			totalRespSig[c] = new double[sigHitNum[c]];
			for(int i=0; i<sigHitNum[c]; i++)
				totalRespSig[c][i] = 0;
			if(ctrlHitNum[c]>0) {
				rBindCtrl[c] = new double[numComp][ctrlHitNum[c]];
				rNoiseCtrl[c] = new double[ctrlHitNum[c]];
				totalRespCtrl[c] = new double[ctrlHitNum[c]];
				for(int i=0; i<ctrlHitNum[c]; i++)
					totalRespCtrl[c][i] = 0;
			}
		}
		
    	////////////////////////////
        //Run ML -- this should only need one or two rounds
    	////////////////////////////
		for(int t=0; t<config.EM_ML_ITER; t++) {
			////////
			//E-step
			////////
			for(int c=0; c<numConditions; c++) {
				int numPairs = sigHitNum[c];
				// Load binding fragment size PDF cache
				Map<Integer, List<Double>> fragSizePDF = bindingManager.getCachePDF(manager.getIndexedCondition(c));
				//Recompute h function, given binding component positions (n function is constant because noise model doesn't move)
				for(int i=0; i<numPairs; i++)
					for(int j=0; j<numComp; j++) {
						if(pi[c][j]>0) {
							double fuzzProb = BindingModel.probability(fuzz[c][j], mu[c][j]-sigHitPos[c][i]);
							double fragSizeProb = 0;
							for(int index: fragSizePDF.keySet()) {
								fragSizeProb += tau[c][j][index] * fragSizePDF.get(index).get(sigHitSize[c][i]);
							}
							h[c][j][i] = fuzzProb * fragSizeProb;
						}
					}
				
				//Compute responsibilities
				for(int i=0; i<numPairs; i++)
					totalRespSig[c][i] = 0;
				for(int i=0; i<numPairs; i++) {
					for(int j=0; j<numComp; j++) {
						if(pi[c][j]>0) {
							rBindSig[c][j][i] = h[c][j][i]*pi[c][j];
							totalRespSig[c][i] += rBindSig[c][j][i];
						}
					}
					rNoiseSig[c][i] = n[c][i] * piNoise[c];
					totalRespSig[c][i] += rNoiseSig[c][i];
				}
				//Normalize responsibilities
				for(int i=0; i<numPairs; i++) {
					for(int j=0; j<numComp; j++) {
						if(pi[c][j]>0) {
							rBindSig[c][j][i]/=totalRespSig[c][i];
						}
					}
					rNoiseSig[c][i]/=totalRespSig[c][i];
				}
			}
			
    		/////////////////////
    		//M-step: maximize pi
    		/////////////////////
			for(int c=0; c<numConditions; c++) {
				int numPairs = sigHitNum[c];
				//Maximize pi
				double[] sumR = new double[numComponents]; //Consider counts, should equal to totalRespSig in SEM
				for(int j=0; j<numComp; j++) {
					if(pi[c][j]>0) {
						for(int i=0; i<numPairs; i++)
							sumR[j] += rBindSig[c][j][i] * sigHitCounts[c][i];
					}
				}
				
				//No components to be eliminated in ML, just update pi(j) (Q should we update fuzziness, tau in SEM?)
				for(int j=0; j<numComp; j++) {
					pi[c][j] = Math.max(0, sumR[j]);
				}
				
				//Normalize pi (accounting for piNoise)
				double totalPi = 0;
				for(int j=0; j<numComp; j++) {
					if(pi[c][j]>0) {
						totalPi += pi[c][j];
					}
				}
				for(int j=0; j<numComp; j++) {
					if(pi[c][j]>0) {
						pi[c][j]=pi[c][j]/(totalPi/(1-piNoise[c]));
					}
				}
			}
			
			//Non-zero components count
        	int nonZeroComps=0;
        	for(int c=0; c<numConditions; c++)
        		for(int j=0;j<numComp;j++)
        			if(pi[c][j]>0.0)
        				nonZeroComps++;
        	
        	////////////
          	//Check Stopping condition
          	////////////	
              if (nonZeroComps>0 && (t==0 || !lastEquivToCurr())){
              	copyStateToLast();
                  continue;
              }else{
              	copyStateToLast();
              	break;
              }
		}//LOOP: Run ML while not converged
		
		//Base log-likelihood calculation
		double[] baseLL = new double[numConditions];
		for(int c=0; c<numConditions; c++) {
			baseLL[c] = 0;
			int numPairs = sigHitNum[c];
			for(int i=0; i<numPairs; i++) {
				// for each read pair, each event will give a conditional prob or bg prob
				double j_sum = 0;
				for(int j=0; j<numComp; j++) {
					if(pi[c][j]>0.0) {
						j_sum += Math.log(rBindSig[c][j][i])/config.LOG2;
					}
				}
				j_sum += Math.log(rNoiseSig[c][i])/config.LOG2;
				baseLL[c] += j_sum * sigHitCounts[c][i];
			}
		}
		
		//ML assignment of signal reads to components is finished
		//Assign control reads with converged pi values here
		for(int c=0; c<numConditions; c++) {
			int numPairs = ctrlHitNum[c];
			double[][] hCtrl = new double[numComp][numPairs];
			double[] nCtrl = new double[numPairs];
			
			// Load binding fragment size PDF cache
			Map<Integer, List<Double>> fragSizePDF = bindingManager.getCachePDF(manager.getIndexedCondition(c));

			
			//Recompute h & n functions for control reads, given binding component positions
			for(int i=0; i<numPairs; i++) {
				for(int j=0; j<numComp; j++) {
					if(pi[c][j]>0) {
						double fuzzProb = BindingModel.probability(fuzz[c][j], mu[c][j]-ctrlHitPos[c][i]);
						double fragSizeProb = 0;
						for(int index: fragSizePDF.keySet()) {
							fragSizeProb += tau[c][j][index] * fragSizePDF.get(index).get(ctrlHitSize[c][i]);
						}
						hCtrl[j][i] = fuzzProb * fragSizeProb;
					}
				}
				nCtrl[i] = noise.get(c).scorePosition(ctrlHitPos[c][i], ctrlRepIndices[c][i]);
			}
			
			//Compute responsibilities
			for(int i=0;i<numPairs;i++)
	            totalRespCtrl[c][i] = 0;
			for(int i=0;i<numPairs;i++){
				for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
					rBindCtrl[c][j][i] = hCtrl[j][i]*pi[c][j];
					totalRespCtrl[c][i] +=rBindCtrl[c][j][i]; 
				}}
				rNoiseCtrl[c][i] = nCtrl[i] * piNoise[c];
				totalRespCtrl[c][i] +=rNoiseCtrl[c][i];
			}
			//Normalize responsibilities
			for(int i=0;i<numPairs;i++){
				for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
					rBindCtrl[c][j][i]/=totalRespCtrl[c][i];
				}}
				rNoiseCtrl[c][i]/=totalRespCtrl[c][i];
			}
		}
	}
	
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
    			for(int x=0; x<rBindSig[c][j].length; x++){
    				lastRBind[c][j][x] = rBindSig[c][j][x];
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
				piBindEquivalent = piBindEquivalent && (Math.abs(pi[c][j]-lastPi[c][j])<config.EM_STATE_EQUIV_THRES);
			}}
    	boolean fuzzBindEquivalent=true;
    	for(int c=0; c<numC; c++)
    		for(int j=0; j<pi[c].length; j++) {
    			if(pi[c][j]>0)
    				fuzzBindEquivalent = fuzzBindEquivalent && (Math.abs(fuzz[c][j]-lastFuzz[c][j]) < config.EM_STATE_EQUIV_THRES);
    		}
    	boolean tauBindEquivalent=true;
    	for(int c=0; c<numC; c++)
    		for(int j=0 ;j<pi[c].length; j++) {
    			if(pi[c][j]>0)
    				for(int t=0; t<tau[c][j].length; t++)
    					tauBindEquivalent = tauBindEquivalent && (Math.abs(tau[c][j][t] - lastTau[c][j][t]) < config.EM_STATE_EQUIV_THRES);
    		}
    	boolean rBindEquivalent=true;
    	for(int c=0; c<numC; c++)
    		for(int j=0; j<pi[c].length; j++){if(pi[c][j]>0){
    			for(int x=0; x<rBindSig[c][j].length; x++){
    				rBindEquivalent = rBindEquivalent && (Math.abs(rBindSig[c][j][x]-lastRBind[c][j][x])<config.EM_STATE_EQUIV_THRES);
    			}
			}}
		return numCompEqual && compPosEqual && piBindEquivalent && fuzzBindEquivalent && tauBindEquivalent && rBindEquivalent;
    }
}



































