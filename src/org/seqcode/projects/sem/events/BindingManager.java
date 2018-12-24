package org.seqcode.projects.sem.events;

import java.io.*;
import java.util.*;

import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.gseutils.Pair;
import org.seqcode.projects.sem.events.BindingEvent;
import org.seqcode.projects.sem.events.BindingModel;
import org.seqcode.projects.sem.events.BindingSubtype;

/**
 * BindingManager stores lists of binding events and motifs associated with experiment conditions,
 * and binding models associated with replicates.
 * @author Jianyu Yang
 *
 */

public class BindingManager {
	
	protected EventsConfig config;
	protected ExperimentManager manager;
	protected List<BindingEvent> events;
	protected Map<ExperimentCondition, List<BindingEvent>> conditionEvents;
	protected Map<ExperimentCondition, List<BindingSubtype>> bindingSubtypes;
	protected Map<ExperimentCondition, Integer> numBindingType;
	protected Map<ControlledExperiment, BindingModel> unstrandedModel;
	protected Map<ExperimentCondition, Double> alpha;
	protected Map<ExperimentCondition, List<BindingSubtype>> potentialBindingSubtypes;
	protected Map<ExperimentCondition, Integer> maxInfluenceRange;
	
	public BindingManager(EventsConfig con, ExperimentManager exptman) {
		config = con;
		manager = exptman;
		events = new ArrayList<BindingEvent>();
		conditionEvents = new HashMap<ExperimentCondition, List<BindingEvent>>();
		bindingSubtypes = new HashMap<ExperimentCondition, List<BindingSubtype>>();
		alpha = new HashMap<ExperimentCondition, Double>();
		numBindingType = new HashMap<ExperimentCondition, Integer>();
		unstrandedModel = new HashMap<ControlledExperiment, BindingModel>();
		potentialBindingSubtypes = new HashMap<ExperimentCondition, List<BindingSubtype>>();
		for(ExperimentCondition cond: manager.getConditions()) {
			conditionEvents.put(cond,  new ArrayList<BindingEvent>());
			bindingSubtypes.put(cond, new ArrayList<BindingSubtype>());
			numBindingType.put(cond, 1);
			alpha.put(cond, 0.0);
			potentialBindingSubtypes.put(cond, new ArrayList<BindingSubtype>());
		}
	}
	
	//Accessors
	public List<BindingEvent> getBindingEvents(){return events;}
	public List<BindingEvent> getConditionBindingEvents(ExperimentCondition ec){return conditionEvents.get(ec);}
	public List<BindingSubtype> getBindingSubtype(ExperimentCondition ec){return bindingSubtypes.get(ec);}
	public Integer getNumBindingType(ExperimentCondition ec){return numBindingType.get(ec);}
	public BindingModel getBindingModel(ControlledExperiment ce){return unstrandedModel.get(ce);}
	public Double getAlpha(ExperimentCondition ec){return alpha.get(ec);}
//	public List<List<StrandedPoint>> getAlignedEventPoints(ExperimentCondition ec){return alignedEventPoints.get(ec);}
	public List<BindingSubtype> getPotentialBindingSubtypes(ExperimentCondition ec){return potentialBindingSubtypes.get(ec);}
	public Integer getMaxInfluenceRange(ExperimentCondition ec) {return maxInfluenceRange.get(ec);}
	
	//Setters
	public void setBindingEvents(List<BindingEvent> e){events =e;}
	public void setConditionBindingEvents(ExperimentCondition ec, List<BindingEvent> e){conditionEvents.put(ec, e);}
	public void setBindingSubtypes(ExperimentCondition ec, List<BindingSubtype> sub){bindingSubtypes.put(ec, sub); numBindingType.put(ec, sub.size());}
	public void setAlpha(ExperimentCondition ec, Double a){alpha.put(ec,a);}
	public void setUnstrandedBindingModel(ControlledExperiment ce, BindingModel mod){unstrandedModel.put(ce, mod);}
//	public void setAlignedEventPoints(ExperimentCondition ec, List<List<StrandedPoint>> points){alignedEventPoints.put(ec, points);}
	public void addPotentialBindingSubtypes(ExperimentCondition ec, List<BindingSubtype> subtypes){potentialBindingSubtypes.get(ec).addAll(subtypes);}
	public void clearPotentialBindingSubtypes(ExperimentCondition ec){ potentialBindingSubtypes.put(ec, new ArrayList<BindingSubtype>());}
	public void setBindingModel(ControlledExperiment ce, BindingModel mod) {unstrandedModel.put(ce, mod);}
	
	public void updateNumBindingTypes() {
		int[] numBindingTypes = new int[manager.getNumConditions()];
		for (ExperimentCondition cond: manager.getConditions()) {
			numBindingTypes[cond.getIndex()] = getNumBindingType(cond);
		BindingEvent.setNumBindingTypes(numBindingTypes);
		}
	}
	
	public void updateMaxInfluenceRange(ExperimentCondition ec) {
		int max=0;
		for(ControlledExperiment rep: ec.getReplicates()) {
			if(getBindingModel(rep).getInfluenceRange()>max)
				max=getBindingModel(rep).getInfluenceRange();
		}
		maxInfluenceRange.put(ec, max);
	}
	
	//For each controlled experiment, simply calculate the proportion of reads in the provided list of
	//binding events to everything else.
	public void estimateSignalVsNoiseFractions(List<BindingEvent> signalEvents) {
		for(ExperimentCondition cond: manager.getConditions()) {
			for(ControlledExperiment r: cond.getReplicates()) {
				double repSigCount = 0, repNoiseCount = 0;
				for(BindingEvent event: signalEvents) {
					if(event.isFoundInCondition(cond)) {
						repSigCount += event.getRepSigHits(r);
					}
				}
				repNoiseCount = r.getSignal().getHitCount() - repSigCount;
				r.setSignalVsNoiseFraction(repSigCount/repNoiseCount); //??? I think the original code is wrong
				System.err.println(r.getName()+"\t"+r.getIndex()+"\tsignal-noise ratio:\t" + String.format("%.4f", r.getSignalVsNoiseFraction()));
			}
		}
	}
	
	//Count the binding events present in a given condition
	public int countSubtypeEventsInCondition(ExperimentCondition cond, double qMinThres) {
		int count = 0;
		for(BindingEvent e: events) {
			if(e.isFoundInCondition(cond) && e.getCondSigVCtrlP(cond)<=qMinThres) {
				count++;
			}
		}
		return count;
	}
	
	//Count the differential binding events present in a given pair of conditions
	public int countDiffEventsBetweenConditions(ExperimentCondition cond, ExperimentCondition othercond, double qMinThres, double diffPMinThres) {
		int count = 0;
		for(BindingEvent e: events) {
			if(e.isFoundInCondition(cond) && e.getCondSigVCtrlP(cond)<=qMinThres) 
				if(e.getInterCondP(cond, othercond)<=diffPMinThres && e.getInterCondFold(cond, othercond)>0)
					count++;
		}
		return count;
	}
}
