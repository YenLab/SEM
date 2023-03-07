package org.seqcode.projects.sem.GMM;

import java.util.*;
import java.io.*;

import org.seqcode.projects.sem.framework.SEMConfig;
import org.seqcode.deepseq.experiments.ExperimentCondition;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class FiniteGaussianMixture extends AbstractCluster{
	protected ExperimentCondition cond;
	protected SEMConfig semconfig;
	protected Map<Integer, Integer> mergeFragSizeFrequency;
	protected List<HashMap<Integer, Integer>> fragSizeFrequency;
	
	protected int[] trainingData;
	protected int[] freq;
	protected int freqSum = 0;
	protected int size = 0;
	protected int mixNum = 3;
	protected double[] weights;
	protected double[] m_means;
	protected double[] m_vars;
	protected double[] m_minVars;
	
	protected static final double ROOT2PI = Math.sqrt(2*Math.PI);
	protected static final int dimension = 1;
	protected static final double MIN_VAR = 1E-10;
	
	/**
	 * Standard constructor
	 * @param cond
	 * @param s
	 * @param frequency
	 * @author Jianyu Yang
	 */
	public FiniteGaussianMixture(ExperimentCondition cond, SEMConfig s, List<HashMap<Integer, Integer>> frequency, int mixNum) {
		this.cond = cond;
		this.semconfig = s;
		this.mixNum = mixNum;
		fragSizeFrequency = frequency;
		mergeFragSizeFrequency = new HashMap<Integer, Integer>();
		
		//Merge all fragment size frequency to a single frequency
		for(HashMap<Integer, Integer> d: fragSizeFrequency) {
			d.forEach((k,v) -> mergeFragSizeFrequency.merge(k, v, (a,b)->a+b));
		}
		
		//Generate sorted dataset
		size = mergeFragSizeFrequency.keySet().size();
		trainingData = new int[size];
		freq = new int[size];
		int index=0;
		int[] keys = new ArrayList<Integer>(mergeFragSizeFrequency.keySet()).stream().mapToInt(Integer::valueOf).toArray();
		Arrays.sort(keys);
		for(int fz: keys) {
			trainingData[index] = fz;
			freq[index] = mergeFragSizeFrequency.get(fz);
			freqSum += mergeFragSizeFrequency.get(fz);
			index++;
		}
		
		//Initialize matrix
		weights = new double[mixNum];
		m_means = new double[mixNum];
		m_vars = new double[mixNum];
		m_minVars = new double[mixNum];
	}
	
	/**
	 * Test constructor
	 * @param csvFile represents file containing fragment size frequency information
	 */
	public FiniteGaussianMixture(String csvFile, String csvSplitBy) {
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
        
		//Generate dataset
        size = mergeFragSizeFrequency.keySet().size();
        trainingData = new int[size];
        freq = new int[size];
		int index=0;
		for(int fz: mergeFragSizeFrequency.keySet()) {
			trainingData[index] = fz;
			freq[index] = mergeFragSizeFrequency.get(fz);
			freqSum += mergeFragSizeFrequency.get(fz);
			index++;
		}
		
		//Initialize matrix
		weights = new double[mixNum];
		m_means = new double[mixNum];
		m_vars = new double[mixNum];
		m_minVars = new double[mixNum];
	}
	
	/**
	 * Fit finite Gaussian Mixture on trainingData
	 */
	@Override
	public void excute() {	
		int m_maxIterNum = 500;
		double err = 0.00001;
		
		boolean loop = true;
		double iterNum = 0;
		double lastL = 0;
		double currL = 0;
		int unchanged = 0;
		
		initParameters(trainingData, freq);
		
		double[] next_means = new double[mixNum];
		double[] next_weights = new double[mixNum];
		double[] next_vars = new double[mixNum];
		double[][] resp = new double[mixNum][size];
		List<DataNode> cList = new ArrayList<DataNode>();
		
		while(loop) {
			Arrays.fill(next_weights, 0);
			cList.clear();
			for(int i=0; i<mixNum; i++) {
				Arrays.fill(next_means, 0);
				Arrays.fill(next_vars, 0);
			}
			
			lastL = currL;
			currL = 0;
			for(int k=0; k<size; k++) {
				double p = getProbability(trainingData[k]);	// Sum of all probability density function on this data point
				DataNode dn = new DataNode(trainingData[k], freq[k]);
				dn.index = k;
				cList.add(dn);
				double maxResp = Double.MIN_VALUE;
				for(int j=0; j<mixNum; j++) {					
					resp[j][k] = (getProbability(trainingData[k], j) * weights[j] / p ) * freq[k];	// Proportion of a specific PDF on this data point		
					next_weights[j] += resp[j][k];
					if(resp[j][k]>maxResp) {
						maxResp = resp[j][k];
						dn.cindex = j;
					}
				}
				currL += ((p > 1E-20) ? Math.log(p) : -20) * freq[k];
			}
			currL /= freqSum;
			
			// Re-estimation: generate new weights, means and variances.
			for (int j=0; j<mixNum; j++) {
				weights[j] = next_weights[j] / freqSum;
			}
			
			// means
			for (DataNode dn: cList) {
				for(int j=0; j<mixNum; j++) {
				if(weights[j]>0) {
					next_means[j] += resp[j][dn.index] * trainingData[dn.index] / next_weights[j];
				}
				}
			}
			
			// variances
			for (DataNode dn: cList) {
				for(int j=0; j<mixNum; j++) {
				if(weights[j]>0) {
					next_vars[j] += resp[j][dn.index] * Math.pow((trainingData[dn.index]-next_means[j]), 2) / next_weights[j];
				}
				}
			}
			
			// Check termination			
			m_means = next_means.clone();
			m_vars = next_vars.clone();
			iterNum++;
			if(Math.abs(currL - lastL) < err * Math.abs(lastL)) {
				unchanged++;
			}
			if(iterNum >= m_maxIterNum || unchanged >= 5 ) {
				loop = false;
			}
		}
		
		// Print result
		System.out.println("======================Nucleosome subtype parameters===================");
		for(int j=0; j<mixNum; j++) {
				System.out.println("["+j+"]");
				System.out.println("means: "+m_means[j]);
				System.out.println("vars: "+m_vars[j]);
				System.out.println("weight: " + weights[j]);
				System.out.println();
		}
		
//		// Get cluster assignment
//		for(int i=0; i<size; i++) {
//			System.out.println("data[" + i + "]=" + trainingData[i] + " cindex : " + cList.get(i).cindex);
//		}
		
		// Put data into map
		clusterMu = new LinkedHashMap<Integer, RealVector>();
		clusterSigma = new LinkedHashMap<Integer, RealMatrix>();
		clusterWeight = new LinkedHashMap<Integer, Double>();
		
		for(int j=0; j<mixNum; j++) {
			clusterMu.put(j, MatrixUtils.createRealVector(new double[] {m_means[j]}));
			clusterSigma.put(j, MatrixUtils.createRealMatrix(new double[][] { new double[] {m_vars[j]}}));
			clusterWeight.put(j, weights[j]);
		}
	}

	/**
	 * @param data
	 */
	private void initParameters(int[] data, int[] freq) {
		// Initialize cluster according to cluster number (assign mean to percentile equally)
		List<DataNode> cList = new ArrayList<DataNode>();	
		int[] m_means_index = new int[mixNum];
		for(int i=0; i<mixNum; i++) {
			m_means_index[i] = (int)Math.ceil((double)(i+1)/(double)(mixNum+1) * freqSum);
		}
		int index = 0; int freq_sum = 0;
		for(int k=0; k<size; k++) {
			freq_sum += freq[k];
			if(freq_sum >= m_means_index[index]) {
				m_means[index++] = data[k];
				if(index>=m_means_index.length)
					break;
			}
		}
		
		// Assign data points to the closest cluster
		double[] next_means = new double[mixNum];
		double[] counts = new double[mixNum];
		int unchanged = 0;
		int loopNum = 0;
		while(true) {
			counts = new double[mixNum];
			next_means = new double[mixNum];
			cList.clear();
			for(int k=0; k<size; k++) {
				DataNode dn = new DataNode(data[k], freq[k]);
				dn.index = k;
				double min = Double.MAX_VALUE;
				for(int i=0; i<mixNum; i++) {
					double distance=0;
					distance = Math.abs(data[k]-m_means[i]);
					if(distance<min) {
						min=distance;
						dn.cindex = i;
					}
				}
				counts[dn.cindex]+=freq[k];
				cList.add(dn);
			}
			
			// compute new means
			for(DataNode dn: cList) {
				next_means[dn.cindex] += dn.freq * dn.value;
			}
			
			for(int i=0; i<mixNum; i++)
				next_means[i] /= counts[i];
						
			// check termination
			double diff = 0;
			for(int i=0; i<mixNum; i++)
				diff += Math.abs(m_means[i]-next_means[i]);
			if(diff<5)
				unchanged++;
			loopNum++;
			
			m_means = next_means.clone();
			
			if(unchanged>=5 || loopNum>=100)
				break;
		}
		
		// Compute weights
		for(int i=0; i<mixNum; i++) {
			weights[i] = counts[i] / freqSum;
		}
		
		// Compute variance
		for(DataNode dn: cList) {
			// Count each Gaussian
			m_vars[dn.cindex] += Math.pow((dn.value - m_means[dn.cindex]), 2) * dn.freq;
		}
		
		// Initialize each Gaussian.
		for(int i=0; i<mixNum; i++) {
			if(weights[i]>0) {
				m_vars[i] = m_vars[i] / counts[i];
			}
		}
		
		if (semconfig.isVerbose()) {
		System.out.println("=================GMM Initialization=================");
			for(int i=0; i<mixNum; i++) {
					System.out.println("[" + i + "]: ");
					System.out.println("means : " + m_means[i]);
					System.out.println("var : " + m_vars[i]);
					System.out.println("weights: "+weights[i]);
			}
		}
	}
	
	public double getProbability(double x) {
		double p = 0;
		for(int i=0; i<mixNum; i++) {
			p += weights[i] * getProbability(x, i);
		}
		return p;
	}
	
	public double getProbability(double x, int j) {
		double p = 1;
		double z = (x-m_means[j])/m_vars[j];
		if( z < -1.96 || z > 1.96 ) {
			p = 1e-20;
		} else {
			p *= 1/(Math.sqrt(m_vars[j])*ROOT2PI) * Math.exp(-Math.pow(x-m_means[j], 2)/(2*m_vars[j]));
		}
		return p;
	}
	
	//Save cluster parameters to csv file
	public void saveCSV(String outFile) {
		
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(outFile))) {
        	
        	for(Object o: clusterMu.keySet()) {
        		bw.write(Double.toString(clusterMu.get(o).getEntry(0)));
        		bw.write(',');
        		bw.write(Double.toString(clusterSigma.get(o).getEntry(0, 0)));
        		bw.newLine();
        	}
        	
        	
        } catch (IOException e) {
        	e.printStackTrace();
        }
	}
	
	public class DataNode {
		public int cindex; // cluster
		public int index;
		public int freq;
		public double value;
		
		public DataNode(double v, int f) {
			this.value = v;
			this.freq = f;
			cindex = -1;
			index = -1;
		}
	}
	
	//GMM on low MNase data from MPE-seq
	public static void main(String[] args) {
		String csvSplitBy = "\t";
		String folderName = "D:\\Dropbox\\Code\\GaussianDPMM_fragSize_sample\\LM_H3_rep2\\";
		for(int index=0; index<1; index++) {
			String csvFile = folderName + "fragSizeFrequencySample" + index;
			String outFile = folderName + "ClusterParameters" + index;
		
			FiniteGaussianMixture g = new FiniteGaussianMixture(csvFile, csvSplitBy);
			g.excute();
			g.saveCSV(outFile);
		}
	}
}
