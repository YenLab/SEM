package org.seqcode.projects.sem.GMM;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Random;
import java.util.TreeMap;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.projects.sem.framework.SEMConfig;

import com.datumbox.framework.common.Configuration;
import com.datumbox.framework.common.dataobjects.AssociativeArray;
import com.datumbox.framework.core.common.dataobjects.Dataframe;
import com.datumbox.framework.core.common.dataobjects.Record;
import com.datumbox.framework.core.machinelearning.MLBuilder;
import com.datumbox.framework.core.machinelearning.clustering.GaussianDPMM;

public class InfiniteGaussianMixture extends AbstractCluster{
	protected ExperimentCondition cond;
	protected SEMConfig semconfig;
	protected Map<Integer, Integer> mergeFragSizeFrequency;
	protected List<HashMap<Integer, Integer>> fragSizeFrequency;
	
	protected Dataframe trainingData;

	protected static Configuration configuration = Configuration.getConfiguration();
	protected static final int dimension = 1;
	protected static double alpha = 0.4;
	protected static int kappa0 = 0;
	protected static int nu0 = 0;
	protected static RealVector mu0 = MatrixUtils.createRealVector(new double[1]);
	protected static RealMatrix psi0 = MatrixUtils.createRealIdentityMatrix(1);
	
	// Store number of clusters each round
	protected int round = 10;
	protected Map<Integer, Integer> numClusters = new HashMap<Integer, Integer>();
	
	/**
	 * Standard constructor
	 * @param cond
	 * @param s
	 * @param frequency
	 * @author Jianyu Yang
	 */
	public InfiniteGaussianMixture(ExperimentCondition cond, SEMConfig s, List<HashMap<Integer, Integer>> frequency) {
		this.cond = cond;
		this.semconfig = s;
		fragSizeFrequency = frequency;
		mergeFragSizeFrequency = new HashMap<Integer, Integer>();
		
		//Merge all fragment size frequency to a single frequency
		for(HashMap<Integer, Integer> d: fragSizeFrequency)
			d.forEach((k,v) -> mergeFragSizeFrequency.merge(k, v, (a,b)->a+b));
	}
	
	/**
	 * Test Constructor
	 */
	public InfiniteGaussianMixture(Dataframe trainingData) {
		this.trainingData = trainingData;
	}
	
	//Setters
	public void setAlpha(double alpha) {this.alpha = alpha;}
	public void setKappa(int kappa) {kappa0 = kappa;}
	public void setNu(int nu) {nu0 = nu;}
	public void setMu(RealVector mu) {mu0 = mu;}
	public void setPis(RealMatrix psi) {psi0 = psi;}
	
	//Accessors
	/**
	 * Return the number of clusters with the highest frequency
	 * @return
	 */
	public int getNumClusters() {
		Map.Entry<Integer, Integer> maxEntry = null;

		for (Map.Entry<Integer, Integer> entry : numClusters.entrySet())
		{
		    if (maxEntry == null || entry.getValue().compareTo(maxEntry.getValue()) > 0)
		    {
		        maxEntry = entry;
		    }
		}
		
		if (semconfig.isVerbose())
			System.out.println("Number of clusters reported by DPMM: " + maxEntry.getKey().toString());
		
		return maxEntry.getKey();
	}
	
	/**
	 * Fit Gaussian DPMM on trainingData
	 */
	@Override
	public void excute() {	
		for(int r=0; r<round; r++) {
			//Randomly select fragment size according to the weight then read data into Dataframe
			RandomCollection<Integer> rc = new RandomCollection<>();
			for(int fz: mergeFragSizeFrequency.keySet())
				rc.add(mergeFragSizeFrequency.get(fz), fz);
			trainingData = new Dataframe(configuration);
			for(int i=0; i<1e4; i++)
				trainingData.add(newDataVector(new Object[] {rc.next()}, null));
			
			long startTime = System.nanoTime();
			
			//Construct GaussianDPMM class, configurations store information about memory storage, multithread and other configurations
			//GaussianDPMM.TrainingParameters store model parameters (alpha and NIW distribution parameters)
			
	        String storageName = "temporal";
	        
	        GaussianDPMM.TrainingParameters param = new GaussianDPMM.TrainingParameters();
	        param.setAlpha(alpha);
	        param.setMaxIterations(50);
	        param.setInitializationMethod(GaussianDPMM.TrainingParameters.Initialization.RANDOM_ASSIGNMENT);
	        param.setKappa0(kappa0);
	        param.setNu0(nu0);
	        param.setMu0(mu0);
	        param.setPsi0(psi0);
	        
	        //Fit data
	        GaussianDPMM instance = MLBuilder.create(param, configuration);
	        
	        instance.fit(trainingData);
	        instance.save(storageName);                
	        
	        //Initialize
	    	Map<Integer, List<RealVector>> clusterVector = new HashMap<Integer, List<RealVector>>();
	    	Map<Integer, RealVector> clusterMean = new HashMap<Integer, RealVector>();
	    	Map<Integer, Integer> clusterKappa = new HashMap<Integer, Integer>();
	    	Map<Integer, Integer> clusterNu = new HashMap<Integer, Integer>();
	    	Map<Integer, RealMatrix> clusterPsi = new HashMap<Integer, RealMatrix>();
	    	Map<Integer, Double> clusterC = new HashMap<Integer, Double>();
	    	clusterSigma = new LinkedHashMap<Integer, RealMatrix>();
	    	clusterMu = new LinkedHashMap<Integer, RealVector>();
	    	clusterWeight = new LinkedHashMap<Integer, Double>();
	    	
	    	// Store each cluster's X data as RealVector
	        for(Object o: trainingData.toArray()) {
	        	Record a = (Record) o;
	        	double[] d = new double[a.getX().size()];
	        	for(int i=0; i<dimension; i++) {
	        		d[i] = a.getX().getDouble(i);
	        	}
	        	if (!clusterVector.containsKey(a.getYPredicted())) {
	        		clusterVector.put((Integer)a.getYPredicted(), new ArrayList<RealVector>());
	        		clusterVector.get((Integer)a.getYPredicted()).add(MatrixUtils.createRealVector(d));
	        	} else {
	        		clusterVector.get((Integer)a.getYPredicted()).add(MatrixUtils.createRealVector(d));
	        	}
	        }
	        	
	        // Calculate each cluster's X data mean
	        for(Integer key: clusterVector.keySet()) {
	        	RealVector sum = MatrixUtils.createRealVector(new double[clusterVector.get(key).get(0).getDimension()]);
	        	for(RealVector rv : clusterVector.get(key)) {
	        		sum = sum.add(rv);
	        	}
	        	RealVector mean = sum.mapDivide(clusterVector.get(key).size());
	        	clusterMean.put(key, mean);
	        }
	        
	        // Calculate each cluster's weight
	        for(Integer key: clusterVector.keySet()) {
	        	clusterWeight.put(key, (double)clusterVector.get(key).size()/(double)trainingData.size());
	        }
	        
	        // Calculate each clusters posterior distribution parameters
	        for(Integer key: clusterVector.keySet()) {
	        	int size  = clusterVector.get(key).size();
	        	
	        	clusterMu.put(key, mu0.mapMultiply(kappa0).add(clusterMean.get(key).mapMultiply(size)).mapDivide(kappa0+size));
	        	clusterKappa.put(key, kappa0 + size);
	        	clusterNu.put(key, nu0 + size);
	        	// C
	        	double sum = 0;
	        	for(RealVector rv : clusterVector.get(key)) {
	        		RealVector minusVector = rv.subtract(clusterMean.get(key));
	        		RealMatrix minusMatrix = MatrixUtils.createRealMatrix(new double[][] {
	        			minusVector.toArray()
	        		});
	        		sum += minusMatrix.multiply(minusMatrix.transpose()).getEntry(0, 0);
	        	}
	        	clusterC.put(key, sum);
	        	// PsiN
	        	RealMatrix intermediate = MatrixUtils.createRealMatrix(new double[][] {
	        		clusterMean.get(key).subtract(mu0).toArray()
	        	});        	
	        	double constant = intermediate.multiply(intermediate.transpose()).getEntry(0, 0)*(kappa0*size)/(kappa0+size);
	        	clusterPsi.put(key, MatrixUtils.createRealIdentityMatrix(dimension).scalarMultiply(clusterC.get(key)+constant));
	        	// SigmaN
	        	double coefficient = Double.valueOf(clusterKappa.get(key)+1)/(Double.valueOf(clusterKappa.get(key))*Double.valueOf(clusterNu.get(key)-dimension+1));
	        	clusterSigma.put(key, clusterPsi.get(key).scalarMultiply(coefficient));        	
	        }
	               
	        //Print cluster parameters
//	        System.err.println(cond.getName());
//	        for(Integer o: clusterMean.keySet()) {
//	        	System.err.println("\tCluster " + o);
//	        	System.err.println("\tCluster Weight = " + clusterWeight.get(o));
//	        	System.err.println("\tCluster Mean = " + clusterMu.get(o));
//	        	System.err.println("\tCluster Sigma = " + clusterSigma.get(o));
//	        }
	        
	        //Store cluster number
	        int nc = clusterWeight.size();
	        if(numClusters.containsKey(nc))
	        	numClusters.put(nc, numClusters.get(nc)+1);
	        else
	        	numClusters.put(nc, 1);
	        
	        //Report run time and close the instance
	        long endTime = System.nanoTime();
	        System.err.println("Run time: " + (endTime-startTime)/1000000000 + " seconds");
	        
	        trainingData.close();
	        instance.delete();
		}
	}
	
	/**
	 * Input a T[] array which contains x variables and an Object y which represents cluster (null in DPMM model), returns a Record Object for Datumbox package
	 * Each Record represents a single data point
	 * @param xArray
	 * @param y
	 * @return
	 * @author Jianyu Yang
	 */
	private static <T> Record newDataVector(T[] xArray, Object y) {
        AssociativeArray x = new AssociativeArray();
        for(int i=0;i<xArray.length;++i) {
            x.put(i, xArray[i]);
        }
        return new Record(x, y);
    }
	
	//Save cluster parameters to csv file
	public void save(String outFile, String csvSplitBy) {
		
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(outFile))) {
        	
        	for(Integer o: clusterMu.keySet()) {
        		bw.write(Double.toString(clusterMu.get(o).getEntry(0)));
        		bw.write(csvSplitBy);
        		bw.write(Double.toString(clusterSigma.get(o).getEntry(0, 0)));
        		bw.newLine();
        	}
        	
        	
        } catch (IOException e) {
        	e.printStackTrace();
        }
	}
	
	public class RandomCollection<E> {
	    private final NavigableMap<Double, E> map = new TreeMap<Double, E>();
	    private final Random random;
	    private double total = 0;

	    public RandomCollection() {
	        this(new Random());
	    }

	    public RandomCollection(Random random) {
	        this.random = random;
	    }

	    public RandomCollection<E> add(double weight, E result) {
	        if (weight <= 0) return this;
	        total += weight;
	        map.put(total, result);
	        return this;
	    }

	    public E next() {
	        double value = random.nextDouble() * total;
	        return map.higherEntry(value).getValue();
	    }
	}
	
	//DPMM on low MNase data from MPE-seq
	public static void main(String[] args) {
//		String csvSplitBy = "\t";
//		String folderName = "C:\\Users\\Administrator\\Dropbox\\Code\\GaussianDPMM_fragSize_sample\\LM_H3_rep2\\";
//		for(int index=0; index<10; index++) {
//			String csvFile = folderName + "fragSizeFrequencySample" + index;
//			String outFile = folderName + "ClusterParameters" + index;
//		
//			GaussianMixture g = new GaussianMixture(csvFile, csvSplitBy);
//			g.excute();
//			g.save(outFile, csvSplitBy);
//		}
		
		Dataframe trainingData = new Dataframe(configuration);
		NormalDistribution n1 = new NormalDistribution(300, 50);
		NormalDistribution n2 = new NormalDistribution(110, 30);
		NormalDistribution n3 = new NormalDistribution(200, 30);
        int observationsPerCluster = 1000;
        for(int i=0;i<observationsPerCluster;++i) {
            trainingData.add(newDataVector(new Object[] {n1.sample()}, "c1"));
        }
        
        for(int i=0;i<observationsPerCluster;++i) {
            trainingData.add(newDataVector(new Object[] {n2.sample()}, "c2"));
        }
        
        for(int i=0;i<observationsPerCluster;++i) {
            trainingData.add(newDataVector(new Object[] {n3.sample()}, "c3"));
        }
        
        InfiniteGaussianMixture g = new InfiniteGaussianMixture(trainingData);
        g.setAlpha(0.4);
        g.excute();
	}
}
