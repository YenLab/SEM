package org.seqcode.projects.sem.GMM;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.LinkedHashMap;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.MatrixUtils;

public abstract class AbstractCluster {
	protected LinkedHashMap<Integer, RealVector> clusterMu;	// Means of clusters
	protected LinkedHashMap<Integer, RealMatrix> clusterSigma;	// Variances of clusters
	protected LinkedHashMap<Integer, Double> clusterWeight;	// Weights of clusters
	
	//Accessors
	public Map<Integer, RealVector> getMu() {return clusterMu;}
	public Map<Integer, RealMatrix> getSigma() {return clusterSigma;}
	public Map<Integer, Double> getWeight() {return clusterWeight;}
	public int getNumClusters() {return clusterMu.size();}
	
	public abstract void excute();
	
	public List<RealVector> getParameters() {
		List<RealVector> fragSizeParaList = new ArrayList<RealVector>();
		for (int i: clusterMu.keySet()) {
			fragSizeParaList.add(MatrixUtils.createRealVector(new double[] {
					clusterMu.get(i).getEntry(0), clusterSigma.get(i).getEntry(0, 0), clusterWeight.get(i)
			}));
		}
		return fragSizeParaList;
	}
}
