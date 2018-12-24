package org.seqcode.projects.sem.mixturemodel;

import org.seqcode.genome.location.Point;
import org.seqcode.projects.chexmix.mixturemodel.BindingSubComponents;

/**
 * BindingComponents are used in mixture models to represent potential binding events
 * 
 * @author Jianyu Yang
 *
 */

public class BindingComponent implements Comparable<BindingComponent>{
	
	protected Point coord;			//Event coordinate
	protected Point[][] coords; 	//Event coordinates
	protected int position;			//Position without the chromosome name
	protected int[][] positions;	//Position without the chromosome name
	protected double[][] tau;		//Subtype probabilities
	protected double pi;		//Emission probability
	protected double sumResp = 0;	//Sum of read responsibilities
	protected int index = 0;
	protected boolean isSubtype = false;
	
	public BindingComponent(Point pos, int numReps) {
		coord = pos;
		position  = coord.getLocation();
		positions = null;
		pi = 1;		
	}
	
	//Accessors
	public double getPi() {return pi;}
	public Point getCoord() {return coord;}
	public int getPosition() {return position;}
	public int getIndex() {return index;}
	public double[][] getTau() {return tau;}
	public boolean isSubtype() {return isSubtype;}
	
	public boolean isNonZero() {return pi>0;}
	
	public double getResponsibility() {return sumResp;}
	
	//Setters
	public void setPi(double p) {pi = p;}
	public void setPosition(int p) {position = p; updateCoordFromLocation();}
	public void setCoord(Point p) {coord = p; position = p.getLocation();}
	public void updateCoordFromLocation(){Point newCoord = new Point(coord.getGenome(), coord.getChrom(), position); coord=newCoord;}
	public void setIndex(int i){index=i;}
	public void setTau(double[][] t){tau=t;isSubtype=true;}
	public void setSumResponsibilities(double sumResp) { this.sumResp = sumResp;}
	
	public void setPositions(int[][] p){
		positions = p;
		coords = new Point[p.length][2];
		for (int bt=0; bt< p.length;bt++){
			for (int s=0; s< p[bt].length; s++){
				Point newCoord = new Point(coord.getGenome(), coord.getChrom(), positions[bt][s]); 
				coords[bt][s]= newCoord;
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
}
