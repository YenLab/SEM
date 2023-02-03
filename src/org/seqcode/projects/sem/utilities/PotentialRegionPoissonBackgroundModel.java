package org.seqcode.projects.sem.utilities;

import org.seqcode.deepseq.stats.PoissonBackgroundModel;

import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;

/**
 * PotentialRegionPoissonBackgroundModel: A background model using the Poisson distribution specified for SEM
 * Due to nucleosomes are all over the entire genome, I assume if the reads in the bin are significantly
 * lower than average reads number, this bin will not be considered as a potential region
 * @author Jianyu Yang
 *
 */
public class PotentialRegionPoissonBackgroundModel extends PoissonBackgroundModel {

	protected double lambda;
	
	public PotentialRegionPoissonBackgroundModel(int modelType, double logconfidence, double totalReads, double genomeLength, double mappableGenome, double binWidth, char strand, double scaling, boolean useThisExpt) {
		super(modelType, logconfidence, totalReads, genomeLength, mappableGenome, binWidth, strand, scaling, useThisExpt);
	}
	
	//Set Poisson thresholds (lower bound for nucleosome)
	@Override
	protected int calcCountThreshold() {
		int countThres = 0;
		DRand re = new DRand();
		Poisson p = new Poisson(0, re);
		lambda = (totalReads*binWidth)/(regionLength*mappableRegion);
		
		p.setMean(lambda);
		double l = 0;
		for(int b=1; l<confThreshold; b++) {
			l = p.cdf(b);
			countThres = b;
		}
		return(Math.max(1, countThres));
	}
}
