package org.seqcode.projects.sem;

import java.util.*;

import org.seqcode.math.diff.Normalization;

import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.projects.sem.events.EventsConfig;
import org.seqcode.projects.sem.events.BindingManager;
import org.seqcode.projects.sem.events.BindingModel;
import org.seqcode.projects.sem.events.BindingSubtype;
import org.seqcode.projects.sem.mixturemodel.BindingMixture;
import org.seqcode.projects.sem.GMM.GaussianMixture;
import org.seqcode.projects.sem.framework.SEMConfig;
import org.seqcode.projects.sem.framework.OutputFormatter;
import org.seqcode.projects.sem.framework.PotentialRegionFilter;

public class SEM {
	
	protected ExperimentManager manager;
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected EventsConfig evconfig;
	protected SEMConfig semconfig;
	protected BindingManager bindingManager;
	protected BindingMixture mixtureModel;
	protected PotentialRegionFilter potentialFilter;
	protected OutputFormatter outFormatter;
	protected Normalization normalizer;
	protected GaussianMixture gmm;
	protected Map<ControlledExperiment, List<BindingModel>> repBindingModels;
	protected Map<ExperimentCondition, List<BindingSubtype>> condModels;
	protected Map<ExperimentCondition, List<HashMap<Integer, Integer>>> condFragSizeFrequency;
	
	public SEM(GenomeConfig gcon, ExptConfig econ, EventsConfig evcon, SEMConfig c, ExperimentManager eMan) {
		gconfig = gcon;
		econfig = econ;
		evconfig = evcon;
		manager = eMan;
		semconfig = c;
		semconfig.makeSEMOutputDirs(true);
		bindingManager = new BindingManager(evconfig, manager);
		
		//Initialize binding models & binding model record
		repBindingModels = new HashMap<ControlledExperiment, List<BindingModel>>();
		for(ControlledExperiment rep: manager.getReplicates()) {
			if(evconfig.getDefaultBindingModel()!=null)
				bindingManager.setBindingModel(rep, evconfig.getDefaultBindingModel());
			repBindingModels.put(rep, new ArrayList<BindingModel>());
			repBindingModels.get(rep).add(bindingManager.getBindingModel(rep));
		}
		for(ExperimentCondition cond: manager.getConditions())
			bindingManager.updateMaxInfluenceRange(cond);
		
		//Insert fragment size distribution initializing here
		System.err.println("Doing GMM on fragment size");
		condModels = new HashMap<ExperimentCondition, List<BindingSubtype>>();
		//Read fragment size distribution from each replicate, indexed by experiment condition
		condFragSizeFrequency = new HashMap<ExperimentCondition, List<HashMap<Integer, Integer>>>();
		for(ExperimentCondition cond: manager.getConditions()) {
			condFragSizeFrequency.put(cond, new ArrayList<HashMap<Integer, Integer>>());
			for(ControlledExperiment rep: cond.getReplicates()) {
				condFragSizeFrequency.get(cond).add(rep.getSignal().getFragSizeFrequency());
			}
		}
		//Employ GMM on each experiment condition's fragment size distribution
		for(ExperimentCondition cond: manager.getConditions())
			gmm = new GaussianMixture(cond, condFragSizeFrequency.get(cond));
			List<BindingSubtype> subtypes = gmm.getBindingSubtypes();
			condModels.put(cond, subtypes);
		
		//Find potential binding regions
		System.err.println("Finding potential binding regions.");
		potentialFilter = new PotentialRegionFilter(evconfig, semconfig, econfig, manager, bindingManager);
		List<Region> potentials = potentialFilter.execute();
		System.err.println(potentials.size()+" potential regions found.");
		if(potentials.size()==0) {
			System.err.println("No potential regions - exiting.");
			System.exit(1);
		}
		potentialFilter.printPotentialRegionsToFile();
	}
	
	//Run the mixture model to find binding events
	public void runMixtureModel() {
		Double[] kl;
		System.err.println("Initialize mixture model");
		mixtureModel = new BindingMixture(gconfig, econfig, evconfig, semconfig, manager, bindingManager, potentialFilter);
		
	}
	
}
