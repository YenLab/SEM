package org.seqcode.projects.sem.GMM;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.seqcode.gseutils.Pair;

public abstract class AbstractCluster {
	protected Map<Integer, RealVector> clusterMu;
	protected Map<Integer, RealMatrix> clusterSigma;
	
	public void excute() {}
	
	public List<Pair<Double, Double>> getParameters() {
		List<Pair<Double, Double>> fragSizeParaList = new ArrayList<Pair<Double, Double>>();
		for (Object o: clusterMu.keySet()) {
			fragSizeParaList.add(new Pair(clusterMu.get(o).getEntry(0), clusterSigma.get(0).getEntry(0, 0)));
		}
		return fragSizeParaList;
	}
}
