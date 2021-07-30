package org.seqcode.projects.sem.utilities;

/**
 * Statistical test for nucleosome comparison
 * @author Jianyu Yang *
 */

import java.util.*;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.distribution.FDistribution;

import org.seqcode.projects.sem.framework.SEMConfig;
import org.seqcode.gseutils.Pair;;

public class Statistics {
	
	protected static double fThreshold;
	protected static double tThreshold;
	
	public Statistics(SEMConfig semconfig) {
		fThreshold = semconfig.getFThreshold();
		tThreshold = semconfig.getTThreshold();
	}
	
	public static Pair<double[], boolean[]> comparison(List<Integer> hitPos1, List<Integer> hitPos2, List<Double> resp1, List<Double> resp2) throws Exception {
		//convert arraylist to array
		int[] ix = hitPos1.stream().mapToInt(i->i).toArray();
		int[] iy = hitPos2.stream().mapToInt(i->i).toArray();
		double[] wx = resp1.stream().mapToDouble(i->i).toArray();
		double[] wy = resp2.stream().mapToDouble(i->i).toArray();
		
		double[] x = copyFromIntArray(ix);
		double[] y = copyFromIntArray(iy);
		
		checkSamples(x, wx);
		checkSamples(y, wy);
		
		// comparison of two nucleosomes
		double fuzzConfidence = 0;
		double muConfidence = 0;
		boolean fuzzShared = false;
		boolean muShared = false;
		// step 1: f-test to check if it is better to share fuzziness
		fuzzConfidence = fTest(x, y, wx, wy);
		fuzzShared = fuzzConfidence > fThreshold;
		// step 2: (welch's) t-test to check if it is better to share dyad location
		if(fuzzShared) {
			// t-test if share fuzziness (equal variance)
			muConfidence = tTest(x, y, wx, wy, true);
		} else {
			// welch's t-test if don't share fuzziness (unequal variance)
			muConfidence = tTest(x, y, wx, wy, false);
		}
		muShared = muConfidence > tThreshold;
		
		return new Pair<double[], boolean[]>(new double[] {muConfidence, fuzzConfidence}, new boolean[] {muShared, fuzzShared});
	}
	
	// basic f-test
	public static double fTest(double f, int df1, int df2) {
		FDistribution dist = new FDistribution(df1, df2);
		return 1-dist.cumulativeProbability(f);
	}
	
	// f-test on weighted data
	public static double fTest(double[] x, double[] y, double[] wx, double[] wy) {
		if(x.length != wx.length || y.length != wy.length) {
			System.err.println("Unequal length of array detected!");
			System.exit(1);
		}
		int n1 = x.length;
		int n2 = y.length;
		
		double var1 = weightedVar(x, wx);
		double var2 = weightedVar(y, wy);
		
		double f;
		if(var1 > var2) {
			f = var1/var2;
			return fTest(f, n1-1, n2-1);
		} else {
			f = var2/var1;
			return fTest(f, n2-1, n1-1);
		}
	}
	
	// basic t-test/welch's t-test
	public static double tTest(long n1, long n2, double mu1, double mu2, double var1, double var2, boolean equalVariance) {
		if (equalVariance) {
            long df = n1 + n2 - 2;

            double svar = ((n1 - 1) * var1 + (n2 - 1) * var2) / df;
            double t = (mu1 - mu2) / Math.sqrt(svar * (1.0 / n1 + 1.0 / n2));
            TDistribution dist = new TDistribution(df);
            double p = 2.0 * dist.cumulativeProbability(-Math.abs(t));
            
            return p;
            
        } else {
            double df = Math.pow(var1/n1 + var2/n2, 2)/(var1*var1/(n1*n1*(n1-1)) + var2*var2/(n2*n2*(n2-1)));
            
            //monitor
            if(df<0) {
            	System.out.println(n1+"\t"+n2+"\t"+var1+"\t"+var2);
            }
            
            double t = (mu1 - mu2) / Math.sqrt(var1 / n1 + var2 / n2);
            TDistribution dist = new TDistribution(df);
            double p = 2.0 * dist.cumulativeProbability(-Math.abs(t));
            
            return p;
        }
	}
	
	// t-test/welch's t-test on weighted data
	public static double tTest(double[] x, double[] y, double[] wx, double[] wy, boolean equalVariance) {
		long n1 = x.length;
		long n2 = y.length;
		
		double mu1 = weightedMean(x, wx);
		double mu2 = weightedMean(y, wy);
		
		double var1 = weightedVar(x, wx);
		double var2 = weightedVar(y, wy);
		
		return tTest(n1, n2, mu1, mu2, var1, var2, equalVariance);
	}
	
	public static double sum(double[] arr) {
		double sum = 0;
		for(double d: arr) {
			sum += d;
		}
		return sum;
	}
	
	public static double mean(double[] arr) {
		return sum(arr)/arr.length;
	}
	
	public static double var(double[] arr) {
		double mean = mean(arr);
		double sumVar = 0;
		for(double d: arr) {
			sumVar += Math.pow(d-mean, 2);
		}
		return sumVar/(arr.length-1);
	}
	
	public static double weightedMean(double[] arr, double[] weight) {
		double sumWeight = sum(weight);
		double sum = 0;
		for(int i=0; i<arr.length; i++) {
			sum += arr[i] * weight[i];
		}
		return sum/sumWeight;
	}
	
	public static double weightedVar(double[] arr, double[] weight) {
		double weightedMean = weightedMean(arr, weight);
		double sum = 0;
		double V1 = 0; double V2 = 0;
		for(int i=0; i<arr.length; i++) {
			sum += weight[i] * Math.pow(arr[i]-weightedMean, 2);
			V1 += weight[i];
			V2 += Math.pow(weight[i], 2);
		}
		//return 0  if V1^2 == V2
		if((Math.pow(V1, 2) - V2) == 0) 
			return 0;
		else
			return (V1*sum)/(Math.pow(V1, 2) - V2);
	}
	
	public static double[] copyFromIntArray(int[] source) {
	    double[] dest = new double[source.length];
	    for(int i=0; i<source.length; i++) {
	        dest[i] = source[i];
	    }
	    return dest;
	}
	
	public static void checkSamples(double[] x, double[] wx) throws Exception {
		if(x.length != wx.length)
			throw new Exception("Length of x and wx are not equal!");
		if(x.length < 2 ) { 
			System.out.println(Arrays.toString(x));
			System.out.println(Arrays.toString(wx));
			throw new Exception("Length of x is too small!");
		}
	}
	
    public static void main(String[] args) {
        double x[] = {3.0, 4.0, 1.0, 2.1};
        double y[] = {490.2, 340.0, 433.9};
        double wx[] = {1, 1, 1, 1};
        double wy[] = {1, 1, 1};
        double p = tTest(x, y, wx, wy, false);
        System.out.println("p = " + p);
        System.out.println(fTest(2, 50, 50));
    }
}


