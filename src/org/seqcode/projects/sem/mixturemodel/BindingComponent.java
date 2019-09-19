package org.seqcode.projects.sem.mixturemodel;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.HashMap;

import org.seqcode.genome.location.Point;
import org.seqcode.gseutils.Pair;

/**
 * BindingComponents are used in mixture models to represent potential binding events
 * 
 * @author Jianyu Yang
 *
 */

public class BindingComponent implements Comparable<BindingComponent>{
	
	protected Point coord;			// Event coordinate (nucleosome dyad location)
	protected int position;			// Position without the chromosome name
	protected double fuzziness;		// Fuzziness score (Variance of Gaussian Distribution) per replicate
	protected int maxIR;			// Max influence range of this binding component (95 percent interval determined by fuzziness)
	protected double[] tau;			// Fragment size subtype probabilities (indexed by subtype index) per replicate
	protected double pi;			// Emission probability
	protected double sumResp;		// Sum of read responsibilities
	protected double pValue;		// p value of bindingComponent responsibility with respect to noise (only active when fixedAlpha<0)
	protected double[][] readProfile;	// Read responsibility for each read (indexed by replicate & read index) 
	protected int index = 0;
	protected int subtypeIndex = 0;	// Index of subtype with highest probability
	protected boolean isSubtype = false;
	protected boolean isPair = false;
	protected Map<Pair<Integer, Integer>, Pair<double[], boolean[]>> compareResults;			//paired nucleosome index in EM (muSharedBetter, fuzzSharedBetter)
	protected Map<Pair<Integer, Integer>, Pair<double[], boolean[]>> compareResultsConvert;	//paired nucleosome index in activeComponents
	
	public BindingComponent(Point pos, int numReps) {
		coord = pos;
		position  = coord.getLocation();
		pi = 1;
		
		compareResults = new HashMap<Pair<Integer, Integer>, Pair<double[], boolean[]>>();
		compareResults = new HashMap<Pair<Integer, Integer>, Pair<double[], boolean[]>>();
	}
	
	//Accessors
	public double getPi() {return pi;}
	public Point getCoord() {return coord;}
	public int getPosition() {return position;}
	public double getFuzziness() {return fuzziness;}
	public int getMaxIR() {return maxIR;}
	public int getIndex() {return index;}
	public double[] getTau() {return tau;}
	public boolean isSubtype() {return isSubtype;}
	public int getMaxSubtype() {return subtypeIndex; }
	public boolean isPair() {return isPair;}
	public boolean isNonZero() {return pi>0;}
	public double getResponsibility() {return sumResp;}
	public double getPValue() {return pValue;}
	public Map<Pair<Integer, Integer>, Pair<double[], boolean[]>> getCompareResults() {return compareResults;}
	public Map<Pair<Integer, Integer>, Pair<double[], boolean[]>> getCompareRestulsConvert() {return compareResultsConvert;}
	
	//Setters
	public void setPi(double p) {pi = p;}
	public void setPosition(int p) {position = p; updateCoordFromLocation();}
	public void setCoord(Point p) {coord = p; position = p.getLocation();}
	public void updateCoordFromLocation(){Point newCoord = new Point(coord.getGenome(), coord.getChrom(), position); coord=newCoord;}
	public void setIndex(int i){index=i;}
	public void setFuzziness(double f) {
		fuzziness=f;
		maxIR = (int)(Math.sqrt(fuzziness * 1.96 * 2));
	}
	public void setTau(double[] t) {
		tau=new double[t.length];
		for(int i=0; i<t.length; i++)
			tau[i] = t[i];
		isSubtype=true;
		
		//Set subtype with highest probability as the subtype of this component
		int maxIndex = -1; double maxProb = -Double.MAX_VALUE;
		for(int i=0; i<tau.length; i++) {
			if(tau[i] > maxProb) {
				maxIndex = i;
				maxProb = tau[i];
			}
		}
		subtypeIndex = maxIndex;
	}
	public void setSumResponsibility(double sumResp) { this.sumResp = sumResp;}
	public void setPValue(double p) {pValue = p;}
	public void setCompareResults(Map<Pair<Integer, Integer>, Pair<double[], boolean[]>> cr) {
		compareResults = cr;
		if(compareResults.size()>0) {
			isPair = true;
		}
	}
	
	// Convert index of paired nucleosome to new index used in activeComponents to make sure we can find paired nucleosome in activeComponents 
	public void convertIndex(List<Map<Integer, Integer>> indexConverter) {
		compareResultsConvert = new HashMap<Pair<Integer, Integer>, Pair<double[], boolean[]>>();
		if(isPair) {
			for(Pair<Integer, Integer> index: compareResults.keySet()) {
				try {
					int newIndex = indexConverter.get(index.car()).get(index.cdr());
					compareResultsConvert.put(new Pair<Integer, Integer>(index.car(), newIndex), compareResults.get(index));
				} catch (Exception e) {
					e.printStackTrace();
					System.out.println(index);
					System.exit(1);
				}
			}
		}
	}
	
	public void uniformInit(double initValue){
		pi=initValue;
	}
	
	//Comparable default method
	public int compareTo(BindingComponent m) {
		return getCoord().compareTo(m.getCoord());
	}
	
	//Compare by responsibility
	public int compareByResp(BindingComponent m) {
			return Double.compare(sumResp, m.sumResp);
	}
	
	public String toString(){
//		return "B\tcoor: "+coord.getLocationString()+"\tpi: "+String.format("%.3f",pi)+"\tsumResp: "+String.format("%.3f", sumResp)+
//				"\tfuzziness: "+fuzziness+"\ttau: "+Arrays.toString(tau)+"\tindex: "+index;
		return "chr"+coord.getChrom()+"\t"+position+"\t"+pi+"\t"+sumResp+"\t"+fuzziness+"\t"+subtypeIndex+"\t"+isPair;
	}
}
