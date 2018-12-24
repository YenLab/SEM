package org.seqcode.projects.sem.GMM;

import java.util.*;
import java.io.*;

import org.seqcode.projects.sem.events.BindingSubtype;
import org.seqcode.projects.sem.framework.SEMConfig;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ControlledExperiment;

import com.datumbox.framework.common.Configuration;
import com.datumbox.framework.core.machinelearning.MLBuilder;
import com.datumbox.framework.core.machinelearning.clustering.GaussianDPMM;
import com.datumbox.framework.common.dataobjects.AssociativeArray;
import com.datumbox.framework.core.common.dataobjects.Dataframe;
import com.datumbox.framework.core.common.dataobjects.Record;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.distribution.NormalDistribution;

public class GaussianMixture {
	protected ExperimentCondition cond;
	protected SEMConfig semconfig;
	protected Map<Integer, Integer> mergeFragSizeFrequency;
	protected List<HashMap<Integer, Integer>> fragSizeFrequency;
	protected List<BindingSubtype> subtypes;
	
	protected Dataframe trainingData;
	protected Map<Object, RealVector> clusterMu;
	protected Map<Object, RealMatrix> clusterSigma;
	
	protected static Configuration configuration = Configuration.getConfiguration();
	protected static final int dimension = 1;
	protected static double alpha = 0.3;
	protected static int kappa0 = 0;
	protected static int nu0 = 0;
	protected static RealVector mu0 = MatrixUtils.createRealVector(new double[1]);
	protected static RealMatrix psi0 = MatrixUtils.createRealIdentityMatrix(1);
	
	
	/**
	 * Standard constructor
	 * @param cond
	 * @param s
	 * @param frequency
	 * @author Jianyu Yang
	 */
	public GaussianMixture(ExperimentCondition cond, SEMConfig s, List<HashMap<Integer, Integer>> frequency) {
		this.cond = cond;
		this.semconfig = s;
		fragSizeFrequency = frequency;
		mergeFragSizeFrequency = new HashMap<Integer, Integer>();
		
		//Merge all fragment size frequency to a single frequency
		for(HashMap<Integer, Integer> d: fragSizeFrequency)
			d.forEach((k,v) -> mergeFragSizeFrequency.merge(k, v, (a,b)->a+b));
		
		//Read data into Dataframe
		trainingData = new Dataframe(configuration);
		for(int fz: mergeFragSizeFrequency.keySet()) {
			for(int count=0; count<mergeFragSizeFrequency.get(fz); count++) {
				trainingData.add(newDataVector(new Object[] {fz}, null));
			}
		}
		
		subtypes = new ArrayList<BindingSubtype>();
	}
	
	/**
	 * Test constructor
	 * @param csvFile represents file containing fragment size frequency information
	 */
	public GaussianMixture(String csvFile, String csvSplitBy) {
		//Merge all fragment size frequency to a single frequency
		
        String line = "";
        Map<Integer, Integer> frequency = new HashMap<Integer, Integer>();
        
        try (BufferedReader br = new BufferedReader(new FileReader(csvFile))) {

            while ((line = br.readLine()) != null) {

                // use tab as separator
                String[] info = line.split(csvSplitBy);
                frequency.put(Integer.parseInt(info[0]), Integer.parseInt(info[1]));
            }
            System.out.println(frequency);

        } catch (IOException e) {
            e.printStackTrace();
        }
        
        mergeFragSizeFrequency  = frequency;
        
        //Read data into Dataframe
		trainingData = new Dataframe(configuration);
		for(int fz: mergeFragSizeFrequency.keySet()) {
			for(int count=0; count<mergeFragSizeFrequency.get(fz); count++) {
				trainingData.add(newDataVector(new Object[] {fz}, null));
			}
		}
        
	}
	
	/**
	 * Constructor
	 */
	public GaussianMixture(Dataframe trainingData) {
		this.trainingData = trainingData;
	}
	
	//Setters
	public void setAlpha(double alpha) {this.alpha = alpha;}
	public void setKappa(int kappa) {kappa0 = kappa;}
	public void setNu(int nu) {nu0 = nu;}
	public void setMu(RealVector mu) {mu0 = mu;}
	public void setPis(RealMatrix psi) {psi0 = psi;}
	
	//Accessors
	public Map<Object, RealVector> getMu() {return clusterMu;}
	public Map<Object, RealMatrix> getSigma() {return clusterSigma;}
	
	/**
	 * Fit Gaussian DPMM on trainingData
	 */
	public void excute() {	
		long startTime = System.nanoTime();
		
		//Construct GaussianDPMM class, configurations store information about memory storage, multithread and other configurations
		//GaussianDPMM.TrainingParameters store model parameters (alpha and NIW distribution parameters)
		
        String storageName = "temporal";
        
        GaussianDPMM.TrainingParameters param = new GaussianDPMM.TrainingParameters();
        param.setAlpha(alpha);
        param.setMaxIterations(200);
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
    	Map<Object, List<RealVector>> clusterVector = new HashMap<Object, List<RealVector>>();
    	Map<Object, RealVector> clusterMean = new HashMap<Object, RealVector>();
    	Map<Object, Integer> clusterKappa = new HashMap<Object, Integer>();
    	Map<Object, Integer> clusterNu = new HashMap<Object, Integer>();
    	Map<Object, RealMatrix> clusterPsi = new HashMap<Object, RealMatrix>();
    	Map<Object, Double> clusterC = new HashMap<Object, Double>();
    	clusterSigma = new HashMap<Object, RealMatrix>();
    	clusterMu = new HashMap<Object, RealVector>();
    	
    	
    	// Store each cluster's X data as RealVector
        for(Object o: trainingData.toArray()) {
        	Record a = (Record) o;
        	double[] d = new double[a.getX().size()];
        	for(int i=0; i<dimension; i++) {
        		d[i] = a.getX().getDouble(i);
        	}
        	if (!clusterVector.containsKey(a.getYPredicted())) {
        		clusterVector.put(a.getYPredicted(), new ArrayList<RealVector>());
        		clusterVector.get(a.getYPredicted()).add(MatrixUtils.createRealVector(d));
        	} else {
        		clusterVector.get(a.getYPredicted()).add(MatrixUtils.createRealVector(d));
        	}
        }
        	
        // Calculate each cluster's X data mean
        for(Object key: clusterVector.keySet()) {
        	RealVector sum = MatrixUtils.createRealVector(new double[clusterVector.get(key).get(0).getDimension()]);
        	for(RealVector rv : clusterVector.get(key)) {
        		sum = sum.add(rv);
        	}
        	RealVector mean = sum.mapDivide(clusterVector.get(key).size());
        	clusterMean.put(key, mean);
        }
        
        
        // Calculate each clusters posterior distribution parameters
        for(Object key: clusterVector.keySet()) {
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
        int count = 1;
        for(Object o: clusterMean.keySet()) {
        	System.out.println("Cluster " + count);
        	System.out.println("Number of data point = " + clusterVector.get(o).size());
        	System.out.println("Cluster Mean = " + clusterMu.get(o));
        	System.out.println("Cluster Sigma = " + clusterSigma.get(o));
        	count++;
        }
        
        //Report run time and close the instance
        long endTime = System.nanoTime();
        System.out.println("Run time: " + (endTime-startTime)/1000000000 + " seconds");
        
        trainingData.close();
        instance.delete();
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
        	
        	for(Object o: clusterMu.keySet()) {
        		bw.write(Double.toString(clusterMu.get(o).getEntry(0)));
        		bw.write(csvSplitBy);
        		bw.write(Double.toString(clusterSigma.get(o).getEntry(0, 0)));
        		bw.newLine();
        	}
        	
        	
        } catch (IOException e) {
        	e.printStackTrace();
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
        
        GaussianMixture g = new GaussianMixture(trainingData);
        g.setAlpha(0.4);
        g.excute();
	}
}
