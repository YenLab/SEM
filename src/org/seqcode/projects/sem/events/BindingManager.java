package org.seqcode.projects.sem.events;

import java.io.*;
import java.util.*;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealVector;

import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.projects.sem.framework.SEMConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.projects.sem.GMM.AbstractCluster;
import org.seqcode.projects.sem.GMM.GMMFactory;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.xy.XYDataset;


/**
 * BindingManager stores lists of binding events and motifs associated with experiment conditions,
 * and binding models associated with replicates.
 * @author Jianyu Yang
 *
 */

public class BindingManager {
	protected SEMConfig semconfig;
	protected EventsConfig config;
	protected GenomeConfig gconfig;
	protected ExperimentManager manager;
	protected Map<ExperimentCondition, List<BindingModel>> condBindingModels;
	protected Map<ExperimentCondition, List<BindingSubtype>> bindingSubtypes;
	protected Map<ExperimentCondition, List<HashMap<Integer, Integer>>> condFragSizeFrequency;
	protected AbstractCluster gmm;
	protected Map<ExperimentCondition, Integer> numBindingType;
	protected Map<ExperimentCondition, Double> alpha;
	protected Map<ExperimentCondition, List<BindingSubtype>> potentialBindingSubtypes;
	protected Map<ExperimentCondition, Integer> maxInfluenceRange;
	protected Map<ExperimentCondition, double[][]> cachePDF; // &Indexed by ExperimentCondition:BindingSubtype index:fragment size
	
	public BindingManager(SEMConfig sconfig, EventsConfig con, GenomeConfig gcon, ExperimentManager exptman) {
		semconfig = sconfig;
		config = con;
		manager = exptman;
		bindingSubtypes = new HashMap<ExperimentCondition, List<BindingSubtype>>();
		alpha = new HashMap<ExperimentCondition, Double>();
		numBindingType = new HashMap<ExperimentCondition, Integer>();
		potentialBindingSubtypes = new HashMap<ExperimentCondition, List<BindingSubtype>>();
		cachePDF = new HashMap<ExperimentCondition, double[][]>();
		maxInfluenceRange = new HashMap<ExperimentCondition, Integer>();
		for(ExperimentCondition cond: manager.getConditions()) {
			bindingSubtypes.put(cond, new ArrayList<BindingSubtype>());
			numBindingType.put(cond, 1);
			alpha.put(cond, 0.0);
			potentialBindingSubtypes.put(cond, new ArrayList<BindingSubtype>());
		}
	}
	
	//Accessors
	public List<BindingSubtype> getBindingSubtypes(ExperimentCondition ec){return bindingSubtypes.get(ec);}
	public Integer getNumBindingType(ExperimentCondition ec){return numBindingType.get(ec);}
	public Double getAlpha(ExperimentCondition ec){return alpha.get(ec);}
//	public List<List<StrandedPoint>> getAlignedEventPoints(ExperimentCondition ec){return alignedEventPoints.get(ec);}
	public List<BindingSubtype> getPotentialBindingSubtypes(ExperimentCondition ec){return potentialBindingSubtypes.get(ec);}
	public Integer getMaxInfluenceRange(ExperimentCondition ec) {return maxInfluenceRange.get(ec);}
	public double[][] getCachePDF(ExperimentCondition ec) {return cachePDF.get(ec);}
	public List<BindingModel> getBindingModel(ExperimentCondition ec) {return condBindingModels.get(ec);}
	
	//Setters
	public void setAlpha(ExperimentCondition ec, Double a){alpha.put(ec,a);}
//	public void setAlignedEventPoints(ExperimentCondition ec, List<List<StrandedPoint>> points){alignedEventPoints.put(ec, points);}
	public void addPotentialBindingSubtypes(ExperimentCondition ec, List<BindingSubtype> subtypes){potentialBindingSubtypes.get(ec).addAll(subtypes);}
	public void clearPotentialBindingSubtypes(ExperimentCondition ec){ potentialBindingSubtypes.put(ec, new ArrayList<BindingSubtype>());}
	public void setBindingModels(Map<ExperimentCondition, List<BindingModel>> condBMs) {condBindingModels = condBMs;}
	public void setMaxInfluenceRange(ExperimentCondition ec, int maxIR) {maxInfluenceRange.put(ec, maxIR);}
	public void setBindingSubtypes(ExperimentCondition ec, List<BindingSubtype> sub){
		bindingSubtypes.put(ec, sub); 
		numBindingType.put(ec, sub.size());
		}
	
	/**
	 *  Initialize binding subtypes using GMM
	 */
	public void initializeBindingSubtypes() {
		//Read fragment size distribution from each replicate, indexed by experiment condition
		condFragSizeFrequency = new HashMap<ExperimentCondition, List<HashMap<Integer, Integer>>>();
		for(ExperimentCondition cond: manager.getConditions()) {
			condFragSizeFrequency.put(cond, new ArrayList<HashMap<Integer, Integer>>());
			for(ControlledExperiment rep: cond.getReplicates()) {
				HashMap<Integer, Integer> freq = rep.getSignal().getFragSizeFrequency();
				freq.keySet().removeIf(size -> size>semconfig.getMaxFragmentLen()); 				// remove fragment sizes over the set max length
				condFragSizeFrequency.get(cond).add(freq);
			}
		}
		
		if (semconfig.getUserBindingSubtypes().equals("")) {
			System.out.println("GMM on fragment size distribution...");
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
			    setBindingSubtypes(cond, fragSizeSubtypes);
			}
		} else {
			System.out.println("User provided binding subtypes, loading...");
			//Load user provided Binding subtypes info
			for(ExperimentCondition cond: manager.getConditions()) {
				try {
					BufferedReader br = new BufferedReader(new FileReader(semconfig.getUserBindingSubtypes()));
					String line;
					List<BindingSubtype> fragSizeSubtypes = new ArrayList<BindingSubtype>();
					int index=0;
					while((line = br.readLine()) != null) {
						if(semconfig.isVerbose()) System.out.println(line);
						// Delimiter: Tab
						String[] entry = line.split("\t");
						fragSizeSubtypes.add(new BindingSubtype(cond, MatrixUtils.createRealVector(new double[] {
							Double.parseDouble(entry[0]), Double.parseDouble(entry[1]), Double.parseDouble(entry[2])	
						}), index));
						index++;
					}
					setBindingSubtypes(cond, fragSizeSubtypes);
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		cache();
		updateNumBindingTypes();
	}
	
	/**
	 * Initialize binding model (fuzziness)
	 */
	public void initializeBindingModels() {
		//Insert bindingModel initialization here
		condBindingModels = new HashMap<ExperimentCondition, List<BindingModel>>();
		for(ExperimentCondition cond:manager.getConditions()) {
			condBindingModels.put(cond, new ArrayList<BindingModel>());
			condBindingModels.get(cond).add(new BindingModel(semconfig, manager, cond, gconfig));
			setMaxInfluenceRange(cond, condBindingModels.get(cond).get(0).getMaxInfluenceRange());
		}
		setBindingModels(condBindingModels);
		
	}
	
	/**
	 * Cache PDF value for fragment size ranging from min to max
	 * Reduce computational cost
	 * @param min
	 * @param max
	 * @author Jianyu Yang
	 */
	public void cache() {
		for(ExperimentCondition cond: manager.getConditions()) {
			double[][] interCache = new double[numBindingType.get(cond)][semconfig.getMaxFragmentLen()+1];
			for(BindingSubtype b: bindingSubtypes.get(cond)) {
				for(int j=0; j<=semconfig.getMaxFragmentLen(); j++) {
					interCache[b.getIndex()][j]=b.probability(j);
				}
			}
			cachePDF.put(cond, interCache);
		}
	}
	
	public void updateNumBindingTypes() {
		int[] numBindingTypes = new int[manager.getNumConditions()];
		for (ExperimentCondition cond: manager.getConditions()) {
			numBindingTypes[cond.getIndex()] = getNumBindingType(cond);
		BindingEvent.setNumBindingTypes(numBindingTypes);
		}
	}
	
	/**
	 * Save fragment length probability distribution of each nucleosome type into a JPEG file
	 * @param outFile
	 * @param delimiter
	 * @throws IOException 
	 */
	public void printSubTypesToPNG() {
		for(ExperimentCondition cond: manager.getConditions()) {
			String filename = semconfig.getOutputIntermediateDir() + File.separator + semconfig.getOutBase() + "_" + cond.getName() + "_fragLenDist.jpeg";
			
			// Create a XYSeries dataset for each binding subtype
			final XYSeriesCollection dataset = new XYSeriesCollection();
			for(BindingSubtype b: bindingSubtypes.get(cond)) {
				final XYSeries fragLenDist = new XYSeries(b.toString());
				for(int i=0; i<=500; i++) {
					fragLenDist.add(i, b.probability(i));
				}
				dataset.addSeries(fragLenDist);
			}
			
			// Plot Line Chart
			JFreeChart xylineChart = ChartFactory.createXYLineChart(
					"Fragment Length Distribution on " + cond.getName(), 
					"Fragment Length", 
					"Frequency", 
					dataset, 
					PlotOrientation.VERTICAL, 
					true, true, false);
			
			// Save into jpeg
			int width = 640;
			int height = 480;
			File XYChart = new File(filename);
			try {
				ChartUtilities.saveChartAsJPEG(XYChart, xylineChart, width, height);
			} catch (IOException e) {
				e.printStackTrace();
			}
 		}
	}
	
	
	//Print the subtypes information into a file
	public void printSubtypesToFile() {
		String filename = semconfig.getOutputParentDir()+File.separator+semconfig.getOutBase() + "_subtypes.info";
		try {
			FileWriter fout = new FileWriter(filename);
			for(ExperimentCondition cond: manager.getConditions()) {
				fout.write(cond.getName() + "\n");
				for(BindingSubtype b: bindingSubtypes.get(cond)) {
					fout.write("\t" + b + "\n");
				}
			}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
