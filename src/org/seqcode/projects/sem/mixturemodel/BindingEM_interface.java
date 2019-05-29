package org.seqcode.projects.sem.mixturemodel;

import java.util.List;

import org.seqcode.deepseq.StrandedPair;

import org.seqcode.genome.location.Region;
import org.seqcode.projects.sem.utilities.Timer;
import org.seqcode.projects.sem.utilities.EMmode;

/**
 * Interface for group of BindingEM class
 * @author Jianyu Yang
 *
 */

public interface BindingEM_interface {
	public abstract List<List<BindingComponent>> train(List<List<StrandedPair>> signals,
			Region w,
			List<NoiseComponent> noise,
			List<List<BindingComponent>> comps,
			int numComp,
			double[][] atacPrior,
			int trainingRound,
			EMmode mode,
			Timer timer,
			boolean plotEM) throws Exception;
}
