package org.seqcode.projects.sem;

import java.util.*;

import org.seqcode.math.diff.Normalization;
import org.apache.commons.math3.linear.RealVector;

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
import org.seqcode.projects.sem.GMM.AbstractCluster;
import org.seqcode.projects.sem.GMM.GMMFactory;
import org.seqcode.projects.sem.framework.SEMConfig;
import org.seqcode.projects.sem.framework.OutputFormatter;
import org.seqcode.projects.sem.framework.PotentialRegionFilter;
import org.seqcode.projects.sem.utilities.EMmode;
import org.seqcode.projects.sem.utilities.Timer;

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
	protected AbstractCluster gmm;
	protected Map<ExperimentCondition, List<BindingModel>> condBindingModels;
	protected Map<ExperimentCondition, List<BindingSubtype>> condModels;
	protected Map<ExperimentCondition, List<HashMap<Integer, Integer>>> condFragSizeFrequency;
	
	public static String logo = 				
			"\r\n"
			+ "███████╗███████╗███╗   ███╗\r\n"
			+ "██╔════╝██╔════╝████╗ ████║\r\n"
			+ "███████╗█████╗  ██╔████╔██║\r\n"
			+ "╚════██║██╔══╝  ██║╚██╔╝██║\r\n"
			+ "███████║███████╗██║ ╚═╝ ██║\r\n"
			+ "╚══════╝╚══════╝╚═╝     ╚═╝\r\n"
			+ "                           \r\n"
			+ "\n" +
			"Copyright (C) Jianyu Yang 2019-2022\n"
//			"<http://mahonylab.org/software/multigps>\n"
			;
	
	public SEM(GenomeConfig gcon, ExptConfig econ, EventsConfig evcon, SEMConfig c, ExperimentManager eMan) {
		gconfig = gcon;
		econfig = econ;
		evconfig = evcon;
		manager = eMan;
		semconfig = c;
		semconfig.makeSEMOutputDirs(true);
		bindingManager = new BindingManager(semconfig, evconfig, manager);
		
		//Insert fragment size distribution initializing here
		System.err.println("\nGMM on fragment size frequency...");
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
			int numClusters = semconfig.getNumClusters();
			// Use DPMM to determine the number of clusters first if numClusters not specified
			if(numClusters<=0) {
				gmm = GMMFactory.getGMMClass(cond, semconfig, condFragSizeFrequency.get(cond), -1);
				gmm.excute();
				numClusters = gmm.getNumClusters();
			}
			gmm = GMMFactory.getGMMClass(cond, semconfig, condFragSizeFrequency.get(cond), numClusters);
			gmm.excute();
			List<BindingSubtype> fragSizeSubtypes = new ArrayList<BindingSubtype>();
			int index=0;
			for (RealVector para: gmm.getParameters()) {
				fragSizeSubtypes.add(new BindingSubtype(cond, para, index));
				index++;
			}
			bindingManager.setBindingSubtypes(cond, fragSizeSubtypes);
			bindingManager.cache();
			bindingManager.updateNumBindingTypes();
		}
		
		//Insert bindingModel initialization here
		condBindingModels = new HashMap<ExperimentCondition, List<BindingModel>>();
		for(ExperimentCondition cond:manager.getConditions()) {
			condBindingModels.put(cond, new ArrayList<BindingModel>());
			condBindingModels.get(cond).add(new BindingModel(semconfig.getInitialDyad(),semconfig, manager, cond, gconfig));
			bindingManager.setMaxInfluenceRange(cond, condBindingModels.get(cond).get(0).getMaxInfluenceRange());
		}
		bindingManager.setBindingModels(condBindingModels);
		
		//Find potential binding regions
		//Check if user already provided the potential regions
		System.err.println("Finding potential binding regions...");
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
		
		mixtureModel = new BindingMixture(gconfig, econfig, evconfig, semconfig, manager, bindingManager, potentialFilter);
		int round = 0;
		boolean converged = false;
		while (!converged) {
		
			System.err.println("\n============================== EM Round "+(round+1)+" =====================");
			long start = System.currentTimeMillis();
			//Execute the SEM mixture model, only EM
			if(round==0)
				mixtureModel.execute(true, true, EMmode.CONSENSUS); //EM
			else
				mixtureModel.execute(true, false, EMmode.CONSENSUS); //EM
		
			//Update binding models in multiGPS
			//TODO: add statistical test for fuzziness distribution for SEM?
		
			//Update noise models
			mixtureModel.updateGlobalNoise();
			 
			mixtureModel.updateAlpha();
		
			//Print current components
			mixtureModel.printActiveComponentsToFile(EMmode.CONSENSUS);
			
			long end = System.currentTimeMillis();
			System.err.println("Round "+round+"\tOverall time: "+(end-start)/60000+"min");
			round++;
		
			//Check for convergence
			if(round>=semconfig.getMaxModelUpdateRounds() || (mixtureModel.ifConverged() && round>=semconfig.getMinModelUpdateRounds())) {
				converged=true;
			}else {
				converged=false;
			}
			
			//monitor: count time
			if(semconfig.isVerbose()) {
				Timer.summary();
				Timer.showTime();
				Timer.reset();
			}
		}
		//Final output
		// print subtypes info
		bindingManager.printSubtypes();
		// print nucleosome information per condition
		mixtureModel.printNucleosomeInfoToFile();
		// print nucleosome comparison results to file
		mixtureModel.printNucleosomeComparisonToFile();
		
		
		// find alternative and consensus nucleosomes after EM mode has converged
		// alternative nucleosome calling
//		for(int i=0; i<2; i++) {
//			System.err.format("============================== Finding alternative nucleosome (%d/2) =====================\n", i+1);
//			mixtureModel.execute(true, false, EMmode.ALTERNATIVE);
//			// print consensus components
//			mixtureModel.printActiveComponentsToFile(EMmode.ALTERNATIVE);
//			Timer.summary();
//			Timer.showTime();
//			Timer.reset();
//		}
		// find consensus nucleosome after model has converged
//		for(int i=0; i<2; i++) {
//			System.err.format("============================== Finding consensus nucleosome (%d/2) =====================\n", i+1);
//			mixtureModel.execute(true, false, EMmode.CONSENSUS);
//			// print consensus components
//			mixtureModel.printActiveComponentsToFile(EMmode.CONSENSUS);
//			mixtureModel.printNucleosomeComparisonToFile();
//			Timer.summary();
//			Timer.showTime();
//			
//		}
	}
	
	/**
	 * Main driver method for SEM
	 */
	public static void main(String[] args) {
		// Logo
		System.out.println(logo);
		
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args); econ.setLoadPairs(true); econ.setSortMid(true);			// SEM needs midpoint sorted read pairs
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
				logo +
				"\n" +
				"SEM comes with ABSOLUTELY NO WARRANTY.  This is free software, and you\n"+
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
				"\t--providedPotenialRegions <bed file to restrict nucleosome detection regions>" +
				" Loading Data:\n" +
				"\t--expt <file name> AND --format <SAM/BED/SCIDX>\n" +
				"\t--ctrl <file name (optional argument. must be same format as expt files)>\n" +
				"\t--design <experiment design file name to use instead of --expt and --ctrl; see website for format>\n"+
				"\t--nonunique [flag to use non-unique reads]\n" +
				"\t--mappability <fraction of the genome that is mappable for these experiments (default=0.8)>\n" +
				"\t--nocache [flag to turn off caching of the entire set of experiments (i.e. run slower with less memory)]\n" +
				" Detecting nucleosome type:\n" +
				"\t--numClusters <number of nucleosome types> \n\t\tnumber of clusters for finite GMM on fragment size distribution, if set -1, GMM with Dirichlet prior will be used to determine number of types\n" +
				" Running SEM:\n" +
				"\t--r <max. model update rounds, default=3>\n" +
				"\t--alphascale <alpha scaling factor(default=1.0>\n" +
				"\t--fixedalpha <minimum number of fragments a nucleosome should have (default=1, must >= 1)>\n" +
				"\t--exclude <file of regions to ignore>\n" +
				"\t--alternativeExclusion <alternative exclusion zone>\n" +
				"\t--consensusExclusion <consensus exclusion zone>\n" +
				""));
	}
}

























