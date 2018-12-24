package org.seqcode.projects.sem.events;

import org.seqcode.projects.sem.events.EventsConfig;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.projects.sem.framework.SEMConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.genome.location.Region;


/**
 * 
 * @author Jianyu Yang
 *
 */

public class BindingEvent implements Comparable<BindingEvent>{
	
	protected static ExperimentManager experiments = null;
	protected static EventsConfig config = null;
	protected static SEMConfig semconfig = null;
	protected static ExperimentCondition sortingCondition = null;
	protected static ExperimentCondition sortingConditionB = null;
	protected static final int numSingleCondCols = 4;
	protected static final int numInterCondCols = 3;
	protected static int[] numBindingTypes = null;
	protected Point point;
	protected StrandedPoint[][][] typePoints;
	protected double[][][] typeProbs;
	protected Region containingReg;
	protected boolean[] foundInCond; //The binding event was discovered in these conditions [Indexed by condition]
	protected double[] condSigHits; //Signal counts by condition (not scaled) [Indexed by condition]
	protected double[] condCtrlHits; //Control counts by condition (not scales) [Indexed by condition]
	protected double [] condSigVCtrlFold;  //Signal vs Control fold by condition      [indexed by condition]
	protected double [] condSigVCtrlP;  //Signal vs Control P-value by condition      [indexed by condition]
	protected double [] repSigHits;  //Signal counts by replicate  (not scaled)     [indexed by replicate]
	protected double [] repCtrlHits;  //Control counts by replicate  (not scaled)   [indexed by replicate]
	protected double [] repSigVCtrlFold;  //Signal vs Control fold by replicate     [indexed by replicate]
	protected double [] repSigVCtrlP;  //Signal vs Control P-value by replicate     [indexed by replicate]
	protected double [][] interCondScMean;   //Signal vs signal scaled mean (logCPM from EdgeR), inter-condition  [indexed by condition & condition]
	protected double [][] interCondFold;   //Signal vs signal fold difference  (logFold from EdgeR), inter-condition  [indexed by condition & condition]
	protected double [][] interCondP;   //Signal vs signal P, inter-condition         [indexed by condition & condition]
	protected double [][] interRepP;   //Signal vs signal P, inter-replicate        [indexed by replicate & replicate]
	protected double []   LLd; 			//Log-likelihood loss test statistic resulting from eliminating component [indexed by condition]
	protected double []   LLp;			//P-value for LL [indexed by condition]
	
	public BindingEvent(Point p, Region potentialReg) {
		point = p;
		containingReg = potentialReg;
		
		int numC = experiments.getConditions().size();
		int numR = experiments.getReplicates().size();
		typePoints = new StrandedPoint[numC][][];
		typeProbs = new double[numC][][];
		foundInCond = new boolean[numC];
		condSigHits = new double [numC];
		condCtrlHits = new double [numC];
		condSigVCtrlFold = new double [numC];
		condSigVCtrlP = new double [numC];
		interCondScMean = new double [numC][numC];
		interCondFold = new double [numC][numC];
		interCondP = new double [numC][numC];
		repSigHits = new double [numR];
		repCtrlHits = new double [numR];
		repSigVCtrlFold = new double [numR];
		repSigVCtrlP = new double [numR];
		interRepP = new double [numR][numR];
		LLd = new double [numC];
		LLp = new double [numC];
		
		for(int i=0; i<numC; i++) {
			foundInCond[i]=false; 
		}
	}
	
	//Accessors
	public StrandedPoint[][] getTypePoints(ExperimentCondition c){return typePoints[experiments.getConditionIndex(c)];}
	public double[][] getTypeProbs(ExperimentCondition c){return typeProbs[experiments.getConditionIndex(c)];}
 	public boolean isFoundInCondition(ExperimentCondition c){return foundInCond[experiments.getConditionIndex(c)];}
 	public boolean isFoundInCondition(int conditionIndex){return foundInCond[conditionIndex];}
 	public double getCondSigHits(ExperimentCondition c){return(condSigHits[experiments.getConditionIndex(c)]);}
 	public double getCondCtrlHits(ExperimentCondition c){return(condCtrlHits[experiments.getConditionIndex(c)]);}
 	public double getCondSigVCtrlP(ExperimentCondition c){return(condSigVCtrlP[experiments.getConditionIndex(c)]);}
 	public double getCondSigVCtrlFold(ExperimentCondition c){return(condSigVCtrlFold[experiments.getConditionIndex(c)]);}
 	public double getRepSigHits(ControlledExperiment r){return(repSigHits[r.getIndex()]);}
 	public double getRepCtrlHits(ControlledExperiment r){return(repCtrlHits[r.getIndex()]);}
 	public double getRepSigVCtrlP(ControlledExperiment r){return(repSigVCtrlP[r.getIndex()]);}
 	public double getRepSigVCtrlFold(ControlledExperiment r){return(repSigVCtrlFold[r.getIndex()]);}
 	public double getInterCondScMean(ExperimentCondition c1, ExperimentCondition c2){return(interCondScMean[experiments.getConditionIndex(c1)][experiments.getConditionIndex(c2)]);}
 	public double getInterCondFold(ExperimentCondition c1, ExperimentCondition c2){return(interCondFold[experiments.getConditionIndex(c1)][experiments.getConditionIndex(c2)]);}
 	public double getInterCondP(ExperimentCondition c1, ExperimentCondition c2){return(interCondP[experiments.getConditionIndex(c1)][experiments.getConditionIndex(c2)]);}
	public double getInterRepP(ControlledExperiment r1, ControlledExperiment r2){return interRepP[r1.getIndex()][r2.getIndex()];}
	public static int getNumSingleCondCols(){return numSingleCondCols;}
	public static int getNumInterCondCols(){return numInterCondCols;}
	public double getLLd(ExperimentCondition c1){return(LLd[experiments.getConditionIndex(c1)]);}
	public double getLLp(ExperimentCondition c1){return(LLp[experiments.getConditionIndex(c1)]);}
	public Point getPoint() {return point;}
	public Region getContainingRegion() {return containingReg;}
 	
 	//Setters
	public static void setExperimentManager(ExperimentManager e){experiments = e; sortingCondition = experiments.getConditions().get(0);}
	public static void setConfig(EventsConfig c){config = c;}
	
 	public void setTypePoints(ExperimentCondition c, StrandedPoint[][] points) {typePoints[experiments.getConditionIndex(c)] = points;}
 	public void setTypeProbs(ExperimentCondition c, double[][] probs) { typeProbs[experiments.getConditionIndex(c)]=probs;}
	public void setIsFoundInCondition(ExperimentCondition c, boolean found){foundInCond[experiments.getConditionIndex(c)] = found;}
	public void setIsFoundInCondition(int c, boolean found){foundInCond[c] = found;}
	public void setCondSigHits(ExperimentCondition c, double x){condSigHits[experiments.getConditionIndex(c)]=x;}
	public void setCondCtrlHits(ExperimentCondition c, double x){condCtrlHits[experiments.getConditionIndex(c)]=x;}
	public void setCondSigVCtrlFold(ExperimentCondition c, double x){condSigVCtrlFold[experiments.getConditionIndex(c)]=x;}
	public void setCondSigVCtrlP(ExperimentCondition c, double x){condSigVCtrlP[experiments.getConditionIndex(c)]=x;}
	public void setLLd(ExperimentCondition c, double l){LLd[experiments.getConditionIndex(c)]=l;}
	public void setLLp(ExperimentCondition c, double p){LLp[experiments.getConditionIndex(c)]=p;}
	public void setRepSigHits(ControlledExperiment r, double x){repSigHits[r.getIndex()]=x;}
	public void setRepCtrlHits(ControlledExperiment r, double x){repCtrlHits[r.getIndex()]=x;}
	public void setRepSigVCtrlFold(ControlledExperiment r, double x){repSigVCtrlFold[r.getIndex()]=x;}
	public void setRepSigVCtrlP(ControlledExperiment r, double x){repSigVCtrlP[r.getIndex()]=x;}
	public void setInterCondScMean(ExperimentCondition c1, ExperimentCondition c2, double x){interCondScMean[experiments.getConditionIndex(c1)][experiments.getConditionIndex(c2)]=x;}
	public void setInterCondFold(ExperimentCondition c1, ExperimentCondition c2, double x){interCondFold[experiments.getConditionIndex(c1)][experiments.getConditionIndex(c2)]=x;}
	public void setInterCondP(ExperimentCondition c1, ExperimentCondition c2, double x){interCondP[experiments.getConditionIndex(c1)][experiments.getConditionIndex(c2)]=x;}
	public void setInterRepP(ControlledExperiment r1, ControlledExperiment r2, double x){interRepP[r1.getIndex()][r2.getIndex()]=x;}
	
	//Rank according to location
	public int compareByLocation(BindingEvent f) {
		return getPoint().compareTo(f.getPoint());
	}
	
	//Comparable default method
	public int compareTo(BindingEvent f) {
		return compareByLocation(f);
	}
	
	//Record number of binding types per condition
	public static void setNumBindingTypes(int[] bt) {numBindingTypes = bt;}
}
