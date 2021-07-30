package org.seqcode.projects.sem.utilities;

import org.seqcode.deepseq.stats.PoissonBackgroundModel;

import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;

/**
 * Due to in SEM, different nucleosome will have different influence range, this class
 * has a `setRange` method to apply for nucleosomes with different fuzziness
 * @author Jianyu Yang
 *
 */
public class NucleosomePoissonBackgroundModel extends PoissonBackgroundModel {
	
	protected double lambda;
	
	public NucleosomePoissonBackgroundModel(int modelType, double logconfidence, double totalReads, double genomeLength, double mappableGenome, double binWidth, char strand, double scaling, boolean useThisExpt) {
		super(modelType, logconfidence, totalReads, genomeLength, mappableGenome, binWidth, strand, scaling, useThisExpt);
	}
	
	//Setters
	protected void setRange(int r) {binWidth = r;}
	
	//Get Poison thresholds according to the given influence range
	public int calcCountThreshold(int range) {
		setRange(range);
		return calcCountThreshold();
	}
	
	//Calculate the p value given influence range and counts
	public double calcPValue(int range, double counts) {
		setRange(range);
		int icounts = (int)Math.floor(counts);
		DRand re = new DRand();
		Poisson P = new Poisson(0, re);
		lambda = (totalReads*binWidth)/(regionLength*mappableRegion); 
		P.setMean(lambda);
		return 1 - P.cdf(icounts);
	}
	
	public static void main(String[] args) {
		Poisson p = new Poisson(0, new DRand());
		p.setMean(50);
		System.out.println(1-p.cdf(60));
	}
}
