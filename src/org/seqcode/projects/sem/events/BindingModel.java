package org.seqcode.projects.sem.events;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

import org.seqcode.gseutils.Pair;
import org.seqcode.math.stats.StatUtil;

import org.seqcode.deepseq.StrandedPair;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;

/**
 * BindingModel defines a (probabilistic) model of read occurrences around a binding event.
 * In SEM, we assume fragment midpoints follow a Gaussian distribution around nucleosome dyad location
 * with different variance (which we call "fuzziness")
 * @author Jianyu Yang
 *
 */
public class BindingModel {
	protected ExperimentManager manager;
	protected GenomeConfig gconfig;
	protected ExperimentCondition cond;
	
	protected double initialFuzziness;
	protected List<Pair<String, Integer>> initialDyad;
	protected Map<Integer, Double> pairFreqAroundInitialDyad;
	
	protected static boolean isCache = false;
	protected static Double[][] probCache;	// Store probability for each distance under each fuzziness to reduce computation
		
	protected static final double LOG2 = Math.log(2);
	protected static final double ROOT2PI = Math.sqrt(2*Math.PI);
	protected static final double bgProb=0;
	protected static final double logBgProb=-Double.MAX_VALUE;
	protected final int maxIR;	// 95% influence range computed by initialized fuzziness
	
	// Constructor: Read in dyad location of nucleosome to initialize fuzziness
	public BindingModel(String dyadFile, ExperimentManager eman, ExperimentCondition ec, GenomeConfig gc) {
		manager = eman;	
		gconfig = gc;
		cond = ec;
		initialDyad = new ArrayList<Pair<String, Integer>>();
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(dyadFile));
			String line = "";
			while( (line=br.readLine()) != null) {
				if(line.charAt(0) == '#') {continue;}
				// use \t as delimiter
				String[] info = line.split("\t");
				initialDyad.add(new Pair<String, Integer>(info[0].replaceFirst("^chromosome", "").
																	replaceFirst("^chrom", "").replaceFirst("^chr", ""), Integer.parseInt(info[1])));
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		initializeFuzziness();		
		
		maxIR = (int)Math.rint(Math.sqrt(initialFuzziness) * 1.96) * 2;
		
		//monitor
		System.out.println("Initialize Fuzziness: "+initialFuzziness);
		System.out.println("Initialize maxIR: "+maxIR);
	};
	
	//Accessors
	public double getIntialFuzziness() {return initialFuzziness;}
	public int getMaxInfluenceRange() {return maxIR;}
	public Pair<String, Integer> getIntialDyadByIndex(int index) {return initialDyad.get(index);}
	
	//Look up the probability corresponding to a distance
	//Distance should be defined as (Read position - Peak position)
	public static double probability(double variance, int distance) throws Exception {
		double prob;
//		if(Math.abs(distance) > max) {
//			prob = 0;
//		} else 
		if (variance > 0) {
			prob = 1/(Math.sqrt(variance)*ROOT2PI) * Math.exp(-Math.pow(distance, 2)/(2*variance));
		} else if (variance == 0 && distance == 0) {
			prob = 1;
		} else if (variance == 0 && distance != 0) {
			prob = bgProb;
		} else {
			throw new Exception("Variance must >= 0!");
		}
		
		if(prob >= 0) {
			return prob;
		} else {
			throw new Exception("Detected prob < 0!");
		}
	}
	
	public static double logProbability(double variance, int distance) throws Exception {	
		double prob = probability(variance, distance);
		if(prob > 0) {
			return Math.log(prob);
		} else {
			return logBgProb;
		}
	}
	
	private void initializeFuzziness() {
		// Initialize frequency map
		pairFreqAroundInitialDyad = new HashMap<Integer, Double>();
		for(int i=-75; i<=75; i++) {
			pairFreqAroundInitialDyad.put(i, 0.);
		}
		
		for(Pair<String, Integer> p: initialDyad) {
			//Get strandedPair around dyad +/- 75bp
			int start = (p.cdr()-75)>0 ? (p.cdr()-75) : 0;
			int end = p.cdr() + 75;
			Region r = new Region(gconfig.getGenome(), p.car(), start, end);
			List<StrandedPair> pairs = new ArrayList<StrandedPair>();
			for(ControlledExperiment rep: cond.getReplicates()) {
				pairs.addAll(rep.getSignal().getPairsByMid(r));
			}
			
			//Add strandedPair into frequency map
			for (StrandedPair pair: pairs) {
				int distance = p.cdr() - pair.getMidpoint().getLocation();
				pairFreqAroundInitialDyad.put(distance, pairFreqAroundInitialDyad.get(distance) + pair.getWeight());
			}
		}
		
		//Compute initial fuzziness for all nucleosomes
		initialFuzziness = 0;
		double sumWeight = 0;
		for(int dis: pairFreqAroundInitialDyad.keySet()) {
			initialFuzziness += Math.pow(dis, 2) * pairFreqAroundInitialDyad.get(dis);
			sumWeight += pairFreqAroundInitialDyad.get(dis);
		}
		initialFuzziness /= sumWeight;	
		if(Double.isNaN(initialFuzziness)) {
			System.err.println("NaN initialized fuzziness detected, set 100 as initial fuzziness");
			initialFuzziness = 100;
		}
	}
	
	public static void main(String[] args) {
		try {
			System.out.println(Math.log(1/(Math.sqrt(1000)*ROOT2PI) * Math.exp(-Math.pow(2000, 2)/(2*1000))));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

	
	
