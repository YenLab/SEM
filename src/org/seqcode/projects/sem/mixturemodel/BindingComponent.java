package org.seqcode.projects.sem.mixturemodel;

import org.seqcode.genome.location.Point;

/**
 * BindingComponents are used in mixture models to represent potential binding events
 * 
 * @author Jianyu Yang
 *
 */

public class BindingComponent implements Comparable<BindingComponent>{
	
	protected Point coord;			// Event coordinate (nucleosome dyad location)
	protected int position;			// Position without the chromosome name
	protected double fuzziness;		// Fuzziness score (Variance of Gaussian Distribution)
	protected double[] tau;			// Fragment size subtype probabilities (indexed by subtype index)
	protected double pi;			// Emission probability
	protected double sumResp = 0;	// Sum of read responsibilities
	protected double[][]	readProfile;	// Read responsibility for each read (indexed by replicate & read index) 
	protected int index = 0;
	protected boolean isSubtype = false;
	
	public BindingComponent(Point pos, int numReps) {
		coord = pos;
		position  = coord.getLocation();
		pi = 1;		
	}
	
	//Accessors
	public double getPi() {return pi;}
	public Point getCoord() {return coord;}
	public int getPosition() {return position;}
	public double getFuzziness() {return fuzziness;}
	public int getIndex() {return index;}
	public double[] getTau() {return tau;}
	public boolean isSubtype() {return isSubtype;}
	public boolean isNonZero() {return pi>0;}
	
	public double getResponsibility() {return sumResp;}
	
	//Setters
	public void setPi(double p) {pi = p;}
	public void setPosition(int p) {position = p; updateCoordFromLocation();}
	public void setCoord(Point p) {coord = p; position = p.getLocation();}
	public void updateCoordFromLocation(){Point newCoord = new Point(coord.getGenome(), coord.getChrom(), position); coord=newCoord;}
	public void setIndex(int i){index=i;}
	public void setFuzziness(double f) {fuzziness=f;}
	public void setTau(double[] t){tau=t; isSubtype=true;}
	public void setSumResponsibility(double sumResp) { this.sumResp = sumResp;}
	
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
}
