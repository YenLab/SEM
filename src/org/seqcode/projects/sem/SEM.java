package org.seqcode.projects.sem;

import java.util.*;

import org.seqcode.math.diff.Normalization;

import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.deepseq.events.EnrichmentSignificance;
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
import org.seqcode.gseutils.Pair;

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
		for(ExperimentCondition cond: manager.getConditions()) {
			gmm = new GaussianMixture(cond, semconfig, condFragSizeFrequency.get(cond));
			gmm.excute();
			List<BindingSubtype> fragSizeSubtypes = new ArrayList<BindingSubtype>();
			int index=0;
			for (Pair<Double, Double> para: gmm.getParameters()) {
				fragSizeSubtypes.add(new BindingSubtype(cond, para, index));
				index++;
			}
			bindingManager.setBindingSubtypes(cond, fragSizeSubtypes);
		}
		
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
		
		System.err.println("Initialize mixture model");
		mixtureModel = new BindingMixture(gconfig, econfig, evconfig, semconfig, manager, bindingManager, potentialFilter);
		
		int round = 0;
		boolean converged = false;
		while (!converged) {
		
		System.err.println("\n============================== Round "+round+" =====================");
		
		//Execute the SEM mixture model
		if(round==0)
			mixtureModel.execute(true, true); //EM
		else
			mixtureModel.execute(true, false); //EM
		
		//Update binding models in multiGPS, I think SEM doesn't need it
		
		//Update motifs in multiGPS, I think SEM doesn't need it
		
		//Update noise models
		mixtureModel.updateGlobalNoise();
		
		//Print current components
		mixtureModel.printActiveComponentsToFile();
		
		round++;
		
		//Check for convergence
		if(round>semconfig.getMaxModelUpdateRounds()) {
			converged=true;
		}else {
			converged=true;
		}
		}
		
	}
	
}

























