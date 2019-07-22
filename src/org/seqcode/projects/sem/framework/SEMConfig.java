package org.seqcode.projects.sem.framework;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.*;

import org.seqcode.data.io.BackgroundModelIO;
import org.seqcode.data.io.IOUtil;
import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.data.motifdb.CountsBackgroundModel;
import org.seqcode.data.motifdb.MarkovBackgroundModel;
import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.gseutils.*;
import org.seqcode.motifs.FreqMatrixImport;

/**
 * SEMConfig:
 * 		Maintains all constants needed by ChExMix.
 * 
 * @author Jianyu Yang
 * @version
 */
public class SEMConfig {
	
	public static String version = "0.1";
	public boolean isGPS = true; //???
	protected GenomeConfig gconfig;
	protected Genome gen = null;
	protected String outName = "sem", outBase = "sem";
	protected File outDir = null, interDir = null, imagesDir = null;
	protected boolean printHelp = false;
	protected double sigLogConf = -5; //???
	protected double prLogConf = -10; //???
	protected int minModelUpdateRounds = 1; //Minimum number of EM training rounds
	protected int maxModelUpdateRounds = 15; //Maximum number of EM training rounds (May increase @ Jianyu Yang)
	protected int posPriorScaling = 10; //???
	protected int maxThreads = 1;
	protected double alphaScalingFactor = 1.0; //Scaling the condition-specific alpha value by this factor
	protected double fixedAlpha = -1.0; //Fixed alpha value if >= 0 else automatic mode (use 1 to find potential regions then automatically adjust alpha each EM round)
	protected double betaScalingFactor = 0.05; //Scale the condition and component-specfic beta value by this factor (May change @ Jianyu Yang)
	protected double extendWindow = 500; //Range extension around gff points
	protected double prob_shared_binding = 0.9; //Prior probability that binding sites are shared between conditions
	protected int bmAnalysisWindowMax = 2000; //???
	protected int minComponentsForBMUpdate = 50;
	protected int minRefsForBMUpdate = 25;
	protected double minSubtypeFraction = 0.05; // A subtype needs to be associated with at least this fraction of binding events to be supported 
	protected double minComponentReadFactorForBM = 3;// Components must have (this factor times the condition alpha) number of reads assigned before being included in BM update
	protected boolean updateBM = true; //Set to false to turn off binding model update
	protected double gauss_smooth = 1; //Variance for Gaussian smoothing
	protected int addFlankingComponentSpacing = 20; //In non-first rounds of EM, the components are initialized using the positions from the last round with additional flanking components added at this spacing(???)
	protected boolean addFlankingComponent = true;
	protected List<Region> regionsToPlot = new ArrayList<Region>(); // List of regions that will be printed during EM training (for debugging/demonstration)
	protected List<Region> regionsToIgnore = new ArrayList<Region>(); //List of regions that will be ignored during EM training (i.e. known towers, etc)
	protected List<Point> initialPos = null; //List of points loaded from peak file and used to place binding components.
	protected boolean doReadFilter = false; // Turn on per base read filter in case of highly duplicated experiment
	protected int initComponentSpacing = 30; //Initial component spacing
	protected int modelRange = 100; // Window size to extract tag counts
	protected boolean verbose = false; //Print extra output
	protected boolean shareSubtypes=true; //Share subtypes across experiments
	protected boolean useAtacPrior=true;
	protected boolean MLSharedComponentConfiguration = true; //For ML assignment: use a component configuration shared across all conditions or have condition-specific configs.
	protected int alternativeExclusion = 30; // & Exclusion zone used to determine alternative nucleosome
	protected int consensusExclusion = 127; // & Exclusion zone used to determine consensus nucleosome
	protected int numClusters = -1; // & Number of clusters for GMM (if numCluster==-1 will use InfiniteGMM class to determine cluster number automatically)
	protected String initialDyad = "";	// & File containing the dyad locations for fuzziness initialization (format:chr	coordinate)
	protected int test = 0; // Determine whether to use BindingEM_test instead of BindingEM_Statistic
	
		
	//Constants
	public final double LOG2 = Math.log(2);
	public final int POTREG_BIN_STEP = 100; //Sliding window step in potential region scanner(?)
	public final int MAXSECTION = 5000000;
	public final double NOISE_EMISSION_MIN = 0.01; //Arbitrary floor on the emission probability of noise (must be non-zero to mop up noise reads)
    public final double NOISE_EMISSION_MAX = 0.95; //Arbitrary ceiling on the emission probability of noise
    public final int NOISE_DISTRIB_SMOOTHING_WIN = 50; //Smoothing window for the noise distribution used in the BindingMixture
	public final int EM_MU_UPDATE_WIN = 50; // &
	public final int MAX_EM_ITER = 100;
	public final int EM_ML_ITER = 5; // &
	public final int ALPHA_ANNEALING_ITER = 10; // &
	public final int POSPRIOR_ITER = 15; // &
	public final int FUZZINESS_ANNEALING_ITER = 2; // & Update fuzziness every ? turns
	public final int TAU_ANNEALING_ITER = 2; // & Update tau every ? turns
	public final double SPARSE_PRIOR_SUBTYPE = 0.05; // &
	public final boolean CALC_LL = true; // &
	public final double EM_CONVERGENCE = 0.01; // &
	public final int EM_STATE_EQUIV_ROUNDS = 0; // &
	public final double EM_STATE_EQUIV_PI_THRES = 0.01; // &
	public final double EM_STATE_EQUIV_FUZZ_THRES = 0.1; // &
	public final double EM_STATE_EQUIV_TAU_THRES = 0.05; // &
	public final int INIT_COMPONENT_SPACING = 100; // &
	
	protected String[] args;
	public String getArgs() {
		String a = "";
		for(int i = 0; i<args.length; i++)
			a = a + " " + args[i];
		return a;
	}
	
	public SEMConfig(GenomeConfig gcon, String[] arguments) {this(gcon, arguments, true);}
	public SEMConfig(GenomeConfig gcon, String[] arguments, boolean isGPS) {
		System.setProperty("java.awt.headless", "true"); //??
		gconfig = gcon;
		gen = gconfig.getGenome();
		this.args = arguments;
		this.isGPS = isGPS;
		ArgParser ap = new ArgParser(args);
		if(args.length==0 || ap.hasKey("h")) {
			printHelp=true;
		} else {
			try {
				//Test for a config file... if there is concatenate the contents into the args
				if(ap.hasKey("config")) {
					ArrayList<String> confArgs = new ArrayList<String>();
					String confName = ap.getKeyValue("config");
					File confFile = new File(confName);
					if(!confFile.isFile())
						System.err.println("\nCannot find configuration file: "+confName);
					BufferedReader reader = new BufferedReader(new FileReader(confFile));
				    String line;
			        while ((line = reader.readLine()) != null) {
			        	line = line.trim();
			        	String[] words = line.split("\\s+");
			        	if(!words[0].startsWith("--"))
			        		words[0] = new String("--"+words[0]);
			        	confArgs.add(words[0]); 
			        	if(words.length>1){
				        	String rest=words[1];
				        	for(int w=2; w<words.length; w++)
				        		rest = rest+" "+words[w];
				        	confArgs.add(rest);
			        	}
			        }
				
				String [] confArgsArr = confArgs.toArray(new String[confArgs.size()]);
		        String [] newargs =new String[args.length + confArgsArr.length];
		        System.arraycopy(args, 0, newargs, 0, args.length);
		        System.arraycopy(confArgsArr, 0, newargs, args.length, confArgsArr.length);
		        args = newargs;
		        ap = new ArgParser(args);
				}
		        
		        /****Miscellaneous arguments****/
				//Maximum number of model update rounds
				maxModelUpdateRounds = Args.parseInteger(args,"round", maxModelUpdateRounds);
				//Turn off binding model updates
				updateBM = Args.parseFlags(args).contains("nomodelupdate") ? false : true;
				//Minimum number of components to support a binding model update		
				minComponentsForBMUpdate = Args.parseInteger(args,"minmodelupdateevents",minComponentsForBMUpdate);
				//Minimum number of motif references  to support a binding model update		
				minRefsForBMUpdate = Args.parseInteger(args,"minmodelupdaterefs",minRefsForBMUpdate);
				//Parameter for Gaussian smoothing (std. dev.)
				gauss_smooth = Args.parseDouble(args,"gausssmoothparam",gauss_smooth);
				//Background model parameters		
				sigLogConf = Args.parseDouble(args,"highlogconf",sigLogConf);		
				prLogConf = Args.parseDouble(args,"prlogconf",prLogConf);
				//Threads
				maxThreads = Args.parseInteger(args,"threads",maxThreads);
				maxThreads = Math.min(maxThreads, java.lang.Runtime.getRuntime().availableProcessors());
				//Alpha scaling factor
				alphaScalingFactor = Args.parseDouble(args,"alphascale",alphaScalingFactor);
				//Fixed alpha value
				fixedAlpha = Args.parseDouble(args,"fixedalpha",fixedAlpha);
				//Beta scaling factor
				betaScalingFactor = Args.parseDouble(args,"betascale",betaScalingFactor);
				//Number of base pair to extend around gff
				extendWindow = Args.parseDouble(args, "extwin", extendWindow);
				//Initial component spacing
				initComponentSpacing = Args.parseInteger(args,"compspacing",initComponentSpacing);
				//Window size for extracting tag counts
				modelRange = Args.parseInteger(args,"mrange",modelRange);
				//Number of clusters to divide fragment size frequency distribution
				numClusters = Args.parseInteger(args, "numClusters", -1);
				//Initial dyad location file for fuzziness initialization
				initialDyad = Args.parseString(args, "initialDyad", "");
				//Run SEM using BindingEM_test.java instead of BindingEM_Statistic ?
				test = Args.parseInteger(args, "test", 0);
				if(test>0)
					System.err.println("SEM is in test mode" + test + ".............");
				//Fixed exclusion zone
				alternativeExclusion = Args.parseInteger(args, "alternativeExclusion", 30);
				consensusExclusion = Args.parseInteger(args, "consensusExclusion", 147);
				//Output path
				DateFormat df = new SimpleDateFormat("yyyy-MM-dd-hh-mm-ss");  
			    df.setTimeZone(TimeZone.getTimeZone("EST"));
				outName = Args.parseString(args, "out", outName+"_"+df.format(new Date()));
				outDir =  new File(outName); //Output directory
				outBase = outDir.getName(); //Last part of name
				
				if(ap.hasKey("plotregions"))
					regionsToPlot = RegionFileUtilities.loadRegionsFromFile(Args.parseString(args, "plotregions", null), gen, -1);
				//Regions to ignore during EM training
				if(ap.hasKey("exclude"))
					regionsToIgnore = RegionFileUtilities.loadRegionsFromFile(Args.parseString(args, "exclude", null), gen, -1);
				else if (ap.hasKey("excludebed"))
					regionsToIgnore = RegionFileUtilities.loadRegionsFromBEDFile(gen, Args.parseString(args, "excludebed", null), -1);
				//Initial peak file
				if (ap.hasKey("peakf"))
					initialPos = RegionFileUtilities.loadPeaksFromPeakFile(gen, Args.parseString(args, "peakf", null));
				
				//Turn off adding franking components
				addFlankingComponent = Args.parseFlags(args).contains("noflanking") ? false : true; 
				
				// Positional prior weights
				posPriorScaling = Args.parseInteger(args,"pospriorscale",posPriorScaling);
				// Turn on per base read filtering
				doReadFilter = Args.parseFlags(args).contains("readfilter") ? true : false;	
							
				//Extra output
				verbose = Args.parseFlags(args).contains("verbose") ? true : false;
				
				//Not share subtype motifs across experiments
				shareSubtypes = Args.parseFlags(args).contains("subtypenotshared") ? false : true;
				
				
			} catch (FileNotFoundException e){
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	//Accessors
	public Genome getGenome(){return gen;}
	public boolean helpWanted(){return printHelp;}
	public double getSigLogConf(){return sigLogConf;}
	public double getPRLogConf(){return prLogConf;}
	public int getMaxThreads(){return maxThreads;}
	public double getAlphaScalingFactor(){return alphaScalingFactor;}
	public double getBetaScalingFactor(){return betaScalingFactor;}
	public double getFixedAlpha(){return fixedAlpha;}
	public double getProbSharedBinding() {return prob_shared_binding;}
	public double getWindowExtension(){return extendWindow;}
	public int getBMAnalysisWindowMax(){return bmAnalysisWindowMax;}
	public int getAddFlankingComponentSpacing(){return addFlankingComponentSpacing;}
	public boolean getAddFlankingComponents(){return addFlankingComponent;}
	public List<Region> getRegionsToPlot(){return regionsToPlot;}
	public List<Region> getRegionsToIgnore(){return regionsToIgnore;}
	public List<Point> getInitialPos(){return initialPos;}
	public boolean doBMUpdate(){return updateBM;}
	public int getMinComponentsForBMUpdate(){return minComponentsForBMUpdate;}
	public int getMinRefsForBMUpdate(){return minRefsForBMUpdate;}
	public double getMinSubtypeFraction(){return minSubtypeFraction;}
	public double getMinComponentReadFactorForBM(){return minComponentReadFactorForBM;}
	public double getGaussSmoothParam(){return gauss_smooth;}
	public int getMinModelUpdateRounds(){return minModelUpdateRounds;}
	public int getMaxModelUpdateRounds(){return maxModelUpdateRounds;}
	public boolean useReadFilter(){return doReadFilter;}
	public double getPosPriorScaling(){return posPriorScaling;}
	public boolean isVerbose(){return verbose;}
	public int getInitialCompSpacing(){return initComponentSpacing;}
	public int getModelRange(){return modelRange;}
	public boolean getShareSubtypes(){return shareSubtypes;}
	public boolean useAtacPrior() {return useAtacPrior;}
	public boolean getMLSharedComponentConfiguration(){return MLSharedComponentConfiguration;}
	public int getNumClusters() {return numClusters;}
	public int getAlternativeExclusionZone() {return alternativeExclusion;}
	public int getConsensusExclusionZone() {return consensusExclusion;}
	public String getInitialDyad() {return initialDyad;}
	public int getTestMode() {return test;}
	
	/**
	 * Make output directories used by SEM
	 */
	public void makeSEMOutputDirs(boolean makeInterAndImageDirs) {
		//Test if output directory already exists. If it does,  recursively delete contents
		outDir =  new File(outName);
//		if(outDir.exists())
//			deleteDirectory(outDir);
		outBase = outDir.getName();
		//(re)make the output directory
		outDir.mkdirs();
		if(makeInterAndImageDirs){
			//Make the gps intermediate results output directory
			interDir = new File(outDir.getAbsolutePath()+File.separator+"intermediate-results");
			interDir.mkdirs();
			//Make the image results output directory
			imagesDir = new File(outDir.getAbsolutePath()+File.separator+"images");
			imagesDir.mkdirs();
		}
	}
	public String getOutName(){return outName;}
	public String getOutBase(){return outBase;}
	public File getOutputParentDir(){return outDir;}
	public File getOutputIntermediateDir(){return interDir;}
	public File getOutputImagesDir(){return imagesDir;}
	
	/**
	 * Delete a directory
	 */
	public boolean deleteDirectory(File path) {
	    if( path.exists() ) {
	      File[] files = path.listFiles();
	      for(int i=0; i<files.length; i++) {
	         if(files[i].isDirectory()) {
	           deleteDirectory(files[i]);
	         }
	         else {
	           files[i].delete();
	         }
	      }
	    }
	    return( path.delete() );
	}
	
	/**
	 * returns a string describing the arguments handled by this parser. 
	 * @return String
	 */
	public String getArgsList(){
		return(new String("" +
				"Genome:" +
				"\t--species <Species;Genome>\n" +
				"\tOR\n" +
				"\t--geninfo <genome info file> AND --seq <fasta seq directory>\n" +
				"General:\n" +
				"\t--round <max. model update rounds (default="+maxModelUpdateRounds+">\n" +
				"\t--out <out name (default="+outBase+">\n" +
				"\t--d <read distribution model file>\n" +
				"\t--nonunique [flag to use non-unique reads]\n" +
				"\t--threads <number of threads to use (default="+maxThreads+")>\n" +
				"Experiment Design File:\n" +
				"\t--design <file name>\n" +
				"SEM Model:" +
				"\t--model <filename>\n" +
				"Miscellaneous:\n" +
				"\t--initialDyad <File containing the dyad locations for fuzziness intialization>\n"+
				"\t--numClusters <Number of clusters will be used to do GMM on fragment size frequency, do infinite GMM if -1(default=-1)>\n"+
				"\t--prlogconf <Poisson log threshold for potential region scanning(default="+prLogConf+")>\n" +
				"\t--alphascale <alpha scaling factor(default="+alphaScalingFactor+")>\n" +
				"\t--fixedalpha <impose this alpha (default: set automatically)>\n" +
				"\t--extwin <number of bp expansion centered around gff points (default: 500)]\n" +
				"\t--nomodelupdate [flag to turn off binding model updates]\n" +
				"\t--gausssmoothparam <Gaussian smoothing std dev (default="+gauss_smooth+">\n" +
				"\t--exclude <file of regions to ignore>\n" +
				"\t--plotregions <regions to print during EM training>\n" +
				"\t--peakf <peak file used for component initialization>\n" +
				"\t--motifregions <regions to print component distribution histogram>\n" +
				"\t--eventbasecomp [flag to record event base compositions]\n"+
				"\t--nomotifs [flag to turn off motif-finding & motif priors]\n" +
				"\t--nomotifprior [flag to turn off motif priors only]\n" +
				"\t--memepath <path to the meme bin dir (default: meme is in $PATH)>\n" +
				"\t--memenmotifs <number of motifs MEME should find for each condition (default=3)>\n" +
				"\t--back <Markov background model>\n"+
				"\t--verbose [flag to print intermediate files and extra output]\n" +
				"\t--config <config file: all options can be specified in a name<space>value text file, over-ridden by command-line args>\n" +
				""));
	}
}

