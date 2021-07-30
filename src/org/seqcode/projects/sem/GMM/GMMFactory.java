package org.seqcode.projects.sem.GMM;

import java.util.HashMap;
import java.util.List;

import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.projects.sem.framework.SEMConfig;

public class GMMFactory {
	
	/**
	 * Return Infinite GMM class if numMix==-1, else return Finite GMM class
	 * @param cond
	 * @param s
	 * @param frequency
	 * @param numMix
	 * @return
	 */
	public static AbstractCluster getGMMClass(ExperimentCondition cond, SEMConfig s, List<HashMap<Integer, Integer>> frequency, int numClusters) {
		if(numClusters == -1) {
			return new InfiniteGaussianMixture(cond, s, frequency);
		} else {
			return new FiniteGaussianMixture(cond, s, frequency, numClusters);
		}
	}
}
