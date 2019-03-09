package org.seqcode.projects.sem.utilities;
import java.util.*;

/**
 * Count time for each step in BindingEM
 * @author Jianyu
 *
 */
public class Timer {
	protected static double load_pairs = 0;
	protected static double initialize_data = 0;
	protected static double mark_range = 0;
	protected static double compute_HN = 0;
	protected static double compute_resp = 0;
	protected static double normalize_resp = 0;
	protected static double maximize_mu = 0;
	protected static double resolve_duplicate = 0;
	protected static double maximize_fuzz = 0;
	protected static double maximize_tau = 0;
	protected static double maximize_pi = 0;
	protected static double compute_LL = 0;
	
	protected static double EM_MAP = 0;
	protected static double E_step = 0;
	protected static double M_step = 0;
	
	protected double start_time = 0;
	protected double end_time = 0;
	protected double duration = 0;
	
	public void start() {
		start_time = System.currentTimeMillis();
	}
	
	public void end(String step) {
		end_time = System.currentTimeMillis();
		duration = end_time - start_time;
		if(step=="load") {
			load_pairs += duration;
		} else if(step == "initialize") {
			initialize_data += duration;
		} else if(step == "mark") {
			mark_range += duration;
		} else if(step == "HN") {
			compute_HN += duration;
		} else if(step == "cResp") {
			compute_resp += duration;
		} else if(step == "nResp") {
			normalize_resp += duration;
		} else if(step == "mu") {
			maximize_mu += duration;
		} else if(step == "resolve") {
			resolve_duplicate += duration;
		} else if(step == "fuzz") {
			maximize_fuzz += duration;
		} else if(step == "tau") {
			maximize_tau += duration;
		} else if(step == "pi") {
			maximize_pi += duration;
		} else if(step == "LL") {
			compute_LL += duration;
		} else {
			System.out.println("Please input right step name!");
		}
	}
	
	public static void summary() {
		E_step = mark_range + compute_HN + compute_resp + normalize_resp;
		M_step = maximize_mu + resolve_duplicate + maximize_fuzz + maximize_tau + maximize_pi + compute_LL;
	}
	
	public String toString() {
		return "Load pairs: "+load_pairs/1000+"s"+"\n"+
				"Initalize data: "+initialize_data/1000+"s"+"\n" +
				"E step: "+E_step/1000+"s"+
				"\tMark range: "+mark_range/1000+"s"+"\n"+
				"\tCompute H N function: "+compute_HN/1000+"s"+"\n"+
				"\tCompute responsibility: "+compute_resp/1000+"s"+"\n"+
				"\tNormalize responsibility: "+normalize_resp/1000+"s"+"\n"+
				"M step: "+M_step/1000+"s"+
				"\tMaximize mu: "+maximize_mu/1000+"s"+"\n"+
				"\tResolve duplicate: "+resolve_duplicate/1000+"s"+"\n"+
				"\tMaximize fuzziness: "+maximize_fuzz/1000+"s"+"\n"+
				"\tMaximize tau: "+maximize_tau/1000+"s"+"\n"+
				"\tMaximize pi: "+maximize_pi/1000+"s"+"\n"+
				"\tCompute LL: "+compute_LL/1000+"s"+"\n";
	}
}
