package org.seqcode.projects.sem.events;

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
	protected Pair<Double, Double> fragSizeDisPara;
	protected double mean;
	protected double var;
	protected int index; 
	protected int min=0, max=1000;
	
	public BindingSubtype(ExperimentCondition cond, Pair<Double, Double> fragSizeDisPara, int index) {
		condition = cond;
		this.fragSizeDisPara = fragSizeDisPara;
		mean = fragSizeDisPara.car();
		var = fragSizeDisPara.cdr();
		this.index = index;
	}
	
	//Accessors
	public double getMean() {return mean;}
	public double getVar() {return var;}
	public int getIndex() {return index;}
	
	public double probability(int fragSize) {
		if(fragSize<min || fragSize>max) {
			return 0;
		} else {
			double prob = 1/(var*ROOT2PI) * Math.pow(Math.E, -Math.pow(fragSize-mean, 2)/(2*var*var));
			return prob;
		}
	}
}
