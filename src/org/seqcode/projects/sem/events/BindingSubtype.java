package org.seqcode.projects.sem.events;

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
	protected int index; 
	protected double min=0, max=1000;
	
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
	}
	
	//Accessors
	public double getMean() {return mean;}
	public double getVar() {return var;}
	public double getWeight() {return weight;}
	public int getIndex() {return index;}
	
	public double probability(int fragSize) {
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
		return new String("mean: "+Double.toString(mean)+"\nvar: "+Double.toString(var)+"\nweight: "+Double.toString(weight));
	}
}
