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
//		mixtureModel.updateGlobalNoise();
		
		//Print current components
//		mixtureModel.printActiveComponentsToFile();
		
		round++;
		
		//Check for convergence
		if(round>semconfig.getMaxModelUpdateRounds()) {
			converged=true;
		}else {
			converged=true;
		}
		}
		
	}
	
	/**
	 * Main driver method for SEM
	 */
	public static void main(String[] args) {
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		EventsConfig evconfig = new EventsConfig(gcon, args);
		SEMConfig config = new SEMConfig(gcon, args);
		if(config.helpWanted()){
			System.err.println(SEM.getSEMArgsList());
		}else{
			ExperimentManager manager = new ExperimentManager(econ);
			SEM sem = new SEM(gcon, econ, evconfig, config, manager);
			sem.runMixtureModel();
		
			manager.close();
		}
	}
	
	/**
	 * returns a string describing the arguments for the public version of MultiGPS. 
	 * @return String
	 */
	public static String getSEMArgsList(){
		return(new String("" +
				"Copyright (C) Shaun Mahony 2012-2016\n" +
				"<http://mahonylab.org/software/multigps>\n" +
				"\n" +
				"MultiGPS comes with ABSOLUTELY NO WARRANTY.  This is free software, and you\n"+
				"are welcome to redistribute it under certain conditions.  See the MIT license \n"+
				"for details.\n"+
				"\n OPTIONS:\n" +
				" General:\n"+
				"\t--out <output file prefix>\n" +
				"\t--threads <number of threads to use (default=1)>\n" +
				"\t--verbose [flag to print intermediate files and extra output]\n" +
				"\t--config <config file: all options here can be specified in a name<space>value text file, over-ridden by command-line args>\n" +
				" Genome:\n" +
				"\t--geninfo <genome info file> AND --seq <fasta seq directory reqd if using motif prior>\n" +
				" Loading Data:\n" +
				"\t--expt <file name> AND --format <SAM/BED/SCIDX>\n" +
				"\t--ctrl <file name (optional argument. must be same format as expt files)>\n" +
				"\t--design <experiment design file name to use instead of --expt and --ctrl; see website for format>\n"+
				"\t--fixedpb <fixed per base limit (default: estimated from background model)>\n" +
				"\t--poissongausspb <filter per base using a Poisson threshold parameterized by a local Gaussian sliding window>\n" +
				"\t--nonunique [flag to use non-unique reads]\n" +
				"\t--mappability <fraction of the genome that is mappable for these experiments (default=0.8)>\n" +
				"\t--nocache [flag to turn off caching of the entire set of experiments (i.e. run slower with less memory)]\n" +
				"Scaling control vs signal counts:\n" +
				"\t--noscaling [flag to turn off auto estimation of signal vs control scaling factor]\n" +
				"\t--medianscale [flag to use scaling by median ratio (default = scaling by NCIS)]\n" +
				"\t--regressionscale [flag to use scaling by regression (default = scaling by NCIS)]\n" +
				"\t--sesscale [flag to use scaling by SES (default = scaling by NCIS)]\n" +
				"\t--fixedscaling <multiply control counts by total tag count ratio and then by this factor (default: NCIS)>\n" +
				"\t--scalewin <window size for scaling procedure (default=10000)>\n" +
				"\t--plotscaling [flag to plot diagnostic information for the chosen scaling method]\n" +
				" Running MultiGPS:\n" +
				"\t--d <binding event read distribution file>\n" +
				"\t--r <max. model update rounds, default=3>\n" +
				"\t--nomodelupdate [flag to turn off binding model updates]\n" +
				"\t--minmodelupdateevents <minimum number of events to support an update (default=500)>\n" +
				"\t--nomodelsmoothing [flag to turn off binding model smoothing]\n" +
				"\t--splinesmoothparam <spline smoothing parameter (default=30)>\n" +
				"\t--gaussmodelsmoothing [flag to turn on Gaussian model smoothing (default = cubic spline)]\n" +
				"\t--gausssmoothparam <Gaussian smoothing std dev (default=3)>\n" +
				"\t--jointinmodel [flag to allow joint events in model updates (default=do not)]\n" +
				"\t--fixedmodelrange [flag to keep binding model range fixed to inital size (default: vary automatically)]\n" +
				"\t--prlogconf <Poisson log threshold for potential region scanning(default=-6)>\n" +
				"\t--alphascale <alpha scaling factor(default=1.0>\n" +
				"\t--fixedalpha <impose this alpha (default: set automatically)>\n" +
				"\t--mlconfignotshared [flag to not share component configs in the ML step]\n" +
				"\t--exclude <file of regions to ignore>\n" +
				" MultiGPS priors:\n"+
				"\t--noposprior [flag to turn off inter-experiment positional prior (default=on)]\n" +
				"\t--probshared <probability that events are shared across conditions (default=0.9)>\n" +
				"\t--nomotifs [flag to turn off motif-finding & motif priors]\n" +
				"\t--nomotifprior [flag to turn off motif priors only]\n" +
				"\t--memepath <path to the meme bin dir (default: meme is in $PATH)>\n" +
				"\t--memenmotifs <number of motifs MEME should find for each condition (default=3)>\n" +
				"\t--mememinw <minw arg for MEME (default=6)>\n"+
				"\t--mememaxw <maxw arg for MEME (default=18)>\n"+
				"\t--memeargs <additional args for MEME (default=  -dna -mod zoops -revcomp -nostatus)>\n"+
				"\t--meme1proc [flag to enforce non-parallel version of MEME]\n"+
				" Reporting binding events:\n" +
				"\t--q <Q-value minimum (default=0.001)>\n" +
				"\t--minfold <minimum event fold-change vs scaled control (default=1.5)>\n" +
				"\t--nodifftests [flag to turn off differential enrichment tests]\n" +
				"\t--rpath <path to the R bin dir (default: R is in $PATH). Note that you need to install edgeR separately>\n" +
				"\t--edgerod <EdgeR overdispersion parameter (default=0.15)>\n" +
				"\t--diffp <minimum p-value for reporting differential enrichment (default=0.01)>\n" +
				"\t--eventsaretxt [add .txt to events file extension]\n" +
				""));
	}
}

























