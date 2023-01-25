package org.seqcode.projects.sem.events;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;
import org.apache.commons.math3.distribution.NormalDistribution;

import org.seqcode.gseutils.Pair;
import org.seqcode.deepseq.StrandedPair;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.projects.sem.framework.SEMConfig;

import org.seqcode.genome.location.Point;

/**
 * BindingModel defines a (probabilistic) model of read occurrences around a binding event.
 * In SEM, we assume fragment midpoints follow a Gaussian distribution around nucleosome dyad location
 * with different variance (which we call "fuzziness")
 * @author Jianyu Yang
 *
 */
public class BindingModel {
	protected SEMConfig	semconfig;
	protected ExperimentManager manager;
	protected GenomeConfig gconfig;
	protected ExperimentCondition cond;
	
	protected double initialFuzziness;
	protected final int maxIR;	// 99% influence range computed by initialized fuzziness
	protected Map<Integer, Double> pairFreqAroundInitialDyad;
	
	protected static boolean isCache = false;
	protected static Double[][] probCache;	// Store probability for each distance under each fuzziness to reduce computation
		
	protected static final double LOG2 = Math.log(2);
	protected static final double ROOT2PI = Math.sqrt(2*Math.PI);
	protected static final double bgProb=1e-5;
	protected static final double logBgProb=-5;
	
	// Constructor: Read in dyad location of nucleosome to initialize fuzziness
	public BindingModel(SEMConfig config, ExperimentManager eman, ExperimentCondition ec, GenomeConfig gc) {
		semconfig = config;
		manager = eman;	
		gconfig = gc;
		cond = ec;
		
		
		initialFuzziness = semconfig.INIT_FUZZINESS;
		maxIR = (int)Math.rint(Math.sqrt(initialFuzziness) * 2.58) * 2;
		
		//monitor
		if (semconfig.isVerbose()){
			System.out.println("Initialize Fuzziness: "+initialFuzziness);
			System.out.println("Initialize maxIR: "+maxIR);
		}		
	};
	
	//Accessors
	public double getIntialFuzziness() {return initialFuzziness;}
	public int getMaxInfluenceRange() {return maxIR;}
	
	//Look up the probability corresponding to a distance
	//Distance should be defined as (Read position - Peak position)
	public static double probability(double variance, int distance) throws Exception {
		double prob;
		if (variance > 0) {
			if(Math.abs(distance) > (2.58 * Math.sqrt(variance))) {
				prob = bgProb;
			}  else {
				prob = 1/(Math.sqrt(variance)*ROOT2PI) * Math.exp(-Math.pow(distance, 2)/(2*variance));
			}
		} else if (variance == 0 && distance == 0) {
			prob = 1;
		} else if (variance == 0 && distance != 0) {
			prob = bgProb;
		} else {
			throw new Exception("Detected variance " + variance + "Variance must >= 0!");
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
	
	public static void main(String[] args) {
		double a;
		double start; double end;
		start = System.currentTimeMillis();
		for (int i=0; i<10000000; i++)
			a = 1/(Math.sqrt(100)*ROOT2PI) * Math.exp(-Math.pow(20, 2)/(2*100));
		end = System.currentTimeMillis();
		System.out.println("direct: " + (end-start));
		NormalDistribution normal = new NormalDistribution(10, 10);
		start = System.currentTimeMillis();
		for (int i=0; i<10000000; i++)
			a = normal.density(30);
		end = System.currentTimeMillis();
		System.out.println("commons: " + (end-start));
	}
}

	
	
