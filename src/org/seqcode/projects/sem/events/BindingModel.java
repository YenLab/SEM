package org.seqcode.projects.sem.events;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import org.seqcode.gseutils.Pair;
import org.seqcode.math.stats.StatUtil;

/**
 * BindingModel defines a (probabilistic) model of read occurrences around a binding event.
 * In SEM, we assume fragment midpoints follow a Gaussian distribution around nucleosome dyad location
 * with different variance (which we call "fuzziness")
 * @author Jianyu Yang
 *
 */
public class BindingModel {
	protected static final double LOG2 = Math.log(2);
	protected static final double ROOT2PI = Math.sqrt(2*Math.PI);
	protected static int min = -200;			// the start position
	protected static int max = 200;				// the end position
	protected static double bgProb, logBgProb;
	protected int influenceRange;	// 95% probability range 
	
	// Constructor
	protected BindingModel() {};
	
	//Accessors
	public int getMin() {return min;}
	public int getMax() {return max;}
	public int getInfluenceRange() {return influenceRange;}
	
	//Look up the probability corresponding to a distance
	//Distance should be defined as (Read position - Peak position)
	public static double probability(double variance, int distance) {
		if(distance<min || distance>max) {
			return(bgProb);
		} else {
			double prob = 1/(variance*ROOT2PI) * Math.pow(Math.E, -Math.pow(distance, 2)/(2*variance*variance));
			return prob;
		}
	}
	
	public static double logProbability(double variance, int distance) {
		if(distance<min || distance>max){
		  return(logBgProb);
		}else{
		  return(Math.log(probability(variance, distance)));
		}	  
	}
	
}

	
	
