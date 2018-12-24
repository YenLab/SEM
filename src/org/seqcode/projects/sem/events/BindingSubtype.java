package org.seqcode.projects.sem.events;

/**
 * BindingSubtype class represent binding subtypes that are assocaited with fragment size distribution
 * 
 * @author Jianyu Yang
 *
 */

public class BindingSubtype {
	ExperimentCondition condition;
	FragSizeProbabilityDensity[] repBindingModel£»
	
	public BindingSubtype(ExperimentCondition cond, List<StrandedPair> modelRefs, int bindingModelWidth) {
		condition = cond;
		
	}
	
}
