package org.seqcode.projects.sem.utilities;

import org.seqcode.projects.sem.framework.SEMConfig;

import org.seqcode.genome.location.Region;

import java.io.*;
import java.util.Arrays;

/**
 * EMStepPlotter: Used to plot a step in EM over a given region
 * @author Jianyu Yang
 *
 */
public class EMStepPlotter {
	protected int[][]		hitPos;
	protected double[][]	hitCount;
	protected int[][]		hitSize;
	
	protected static Region plotRegion;
	protected static SEMConfig semconfig;
	
	protected static File dir;
	protected static int trainingRound;
	
	public EMStepPlotter(Region w, SEMConfig config, int[][] hitPos, double[][] hitCount, int[][] hitSize, int tr) {
		plotRegion = w;
		semconfig = config;
		
		this.hitPos = hitPos;
		this.hitCount = hitCount;
		this.hitSize = hitSize;
		
		trainingRound = tr;
		
		dir = semconfig.getOutputImagesDir();
		
		//save hits information 
		try {
			String regStr = w.getLocationString().replaceAll(":", "-");
			String fileName = dir.getAbsolutePath()+File.separator+"EM_"+regStr+"_hitsInfo.txt";
			
			FileWriter fout = new FileWriter(fileName);
			fout.write("#"+w.getLocationString()+"\n");
			fout.write("#hitPos\thitSize\thitCount\n");
			for(int c=0; c<hitPos.length; c++) {
				for(int i=0; i<hitPos[c].length; i++) {
					fout.write(hitPos[c][i]+"\t"+hitSize[c][i]+"\t"+hitCount[c][i]+"\n");
				}
			}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Save 
	 * @param hitPos
	 * @param hitCount
	 * @param hitSize
	 * @param mu
	 * @param pi
	 * @param fuzz
	 * @param tau
	 */
	public static void excute(int[][] mu, double[][] resp, double[][] fuzz, double[][][] tau, int r, int t) {
		try {
			String regStr = plotRegion.getLocationString().replaceAll(":", "-");
			String fileName = dir.getAbsolutePath()+File.separator+"EM_"+regStr+"_trainingRound"+trainingRound+"_dyadInfo.txt";
			BufferedWriter bw = new BufferedWriter(new FileWriter(fileName, true));
			bw.write("#"+plotRegion.getLocationString()+"\n");
			bw.write("#Iter\tt\tcondition\tindex\tmu\tpi\tfuzz\n");
			for(int c=0; c<mu.length; c++) {	
				for(int j=0; j<mu[c].length; j++) {
					if(resp[c][j]>0) {
						bw.write(r+"\t"+t+"\t"+c+"\t"+j+"\t"+mu[c][j]+"\t"+resp[c][j]+"\t"+
								fuzz[c][j]+"\t"+Arrays.toString(tau[c][j])+"\n");
					}
				}
			}
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
	}
}
