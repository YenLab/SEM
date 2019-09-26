package org.seqcode.projects.sem.events;

import java.util.*;

import org.apache.commons.math3.linear.RealVector;

import org.seqcode.gseutils.Pair;
import org.seqcode.deepseq.experiments.ExperimentCondition;


/**
 * BindingSubtype class represent binding subtypes that are assocaited with fragment size distribution
 * 
 * @author Jianyu Yang
 *
 */

public class BindingSubtype {
	protected static final double ROOT2PI = Math.sqrt(2*Math.PI);
	protected ExperimentCondition condition;
	protected RealVector fragSizeDisPara;
	protected double mean;
	protected double var;
	protected double weight;
	protected double[] probStore = new double[1001];
	protected int index; 
	protected double min=0, max=1000;

	public BindingSubtype(double mean, double variance, double weight) {
		this.mean = mean;
		this.var = variance;
		this.weight = weight;
		
		// Pre-compute the probability of each fragment size to reduce time
		for(int i=0; i<=1000; i++) {
			probStore[i] = computeProbability(i);
		}
	}
	
	public BindingSubtype(ExperimentCondition cond, RealVector fragSizeDisPara, int index) {
		condition = cond;
		this.fragSizeDisPara = fragSizeDisPara;
		mean = fragSizeDisPara.getEntry(0);
		var = fragSizeDisPara.getEntry(1);
		weight = fragSizeDisPara.getEntry(2);
		this.index = index;
		
		// Initialize 95 confidence interval
		min = mean - 1.96 * Math.sqrt(var);
		max = mean + 1.96 * Math.sqrt(var);
		
		// Pre-compute the probability of each fragment size to reduce time
		for(int i=0; i<=1000; i++) {
			probStore[i] = computeProbability(i);
		}
	}
	
	//Accessors
	public double getMean() {return mean;}
	public double getVar() {return var;}
	public double getWeight() {return weight;}
	public int getIndex() {return index;}
	
	public double probability(int fragSize) {
		if(fragSize <= 1000)
			return probStore[fragSize];
		else
			return computeProbability(fragSize);
	}
	
	public double logProbability(int fragSize) {
		return Math.log(probability(fragSize));
	}
	
	public double computeProbability(int fragSize) {
		
		double prob;
		if(fragSize<min || fragSize>max) {
			prob = 1e-20;
		} else {
			prob = 1/(Math.sqrt(var)*ROOT2PI) * Math.exp(-Math.pow(fragSize-mean, 2)/(2*var));
			prob *= weight;
		}
		return prob;
	}
	
	@Override
	public String toString() {
		return new String("index: " + index + "\t" + "mean: "+Double.toString(mean)+"\tvar: "+Double.toString(var)+"\tweight: "+Double.toString(weight));
	}
	
	public static void main(String[] args) {
		BindingSubtype bs = new BindingSubtype(200, 400, 1);
		double start = System.currentTimeMillis();
		for(int i=0; i<1000000; i++) {
			int length = i%1000;
			bs.probability(length);
		}
		double end = System.currentTimeMillis();
		System.out.println("time 1: " + (end-start));
		
		start = System.currentTimeMillis();
		for(int i=0; i<1000000; i++) {
			int length = i%1000;
			bs.computeProbability(length);
		}
		end = System.currentTimeMillis();
		System.out.println("time 2: "+(end-start));
	}
	
}
