package org.seqcode.projects.sem.utilities;

/**
 * Count time for each step in BindingMixture
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
	protected static double sum_resp = 0;
	protected static double maximize_fuzz = 0;
	protected static double pairwise = 0;
	protected static double maximize_tau = 0;
	protected static double maximize_pi = 0;
	protected static double compute_LL = 0;
	
	protected static double extra = 0;
	
	protected static double EM_MAP = 0;
	protected static double E_step = 0;
	protected static double M_step = 0;
	
	protected double start_time = 0;
	protected double end_time = 0;
	protected double duration = 0;
	
	protected double extra_start = 0;
	protected double extra_end = 0;
	
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
		} else if(step == "sResp") {
			sum_resp += duration;
		} else if(step == "fuzz") {
			maximize_fuzz += duration;
		} else if(step == "tau") {
			maximize_tau += duration;
		} else if(step == "pi") {
			maximize_pi += duration;
		} else if(step == "LL") {
			compute_LL += duration;
		} else if(step == "pair") {
			pairwise += duration;
		} else {
			System.out.println("Please input right step name!");
		}
	}
	
	public void extra_start() {
		extra_start = System.currentTimeMillis();
	}
	
	public void extra_end() {
		extra_end = System.currentTimeMillis();
		duration = extra_end - extra_start;
		extra += duration;
	}
	
	public static void summary() {
		E_step = mark_range + compute_HN + compute_resp + normalize_resp;
		M_step = maximize_mu + resolve_duplicate + maximize_fuzz + 
				maximize_tau + maximize_pi + compute_LL + pairwise;
	}
	
	public static void reset() {
		load_pairs = 0;
		initialize_data = 0;
		mark_range = 0;
		compute_HN = 0;
		compute_resp = 0;
		normalize_resp = 0;
		maximize_mu = 0;
		resolve_duplicate = 0;
		sum_resp = 0;
		maximize_fuzz = 0;
		pairwise = 0;
		maximize_tau = 0;
		maximize_pi = 0;
		compute_LL = 0;
		
		EM_MAP = 0;
		E_step = 0;
		M_step = 0;
		
		extra = 0;
	}
	
	public static void showTime() {
		System.out.println( "Load pairs: "+load_pairs/1000+"s"+"\n"+
				"Initalize data: "+initialize_data/1000+"s"+"\n" +
				"E step: "+E_step/1000+"s"+"\n"+
				"\tMark range: "+mark_range/1000+"s"+"\n"+
				"\tCompute H N function: "+compute_HN/1000+"s"+"\n"+
				"\tCompute responsibility: "+compute_resp/1000+"s"+"\n"+
				"\tNormalize responsibility: "+normalize_resp/1000+"s"+"\n"+
				"\tSum responsibility: " + sum_resp/1000+"s"+"\n"+
				"M step: "+M_step/1000+"s"+"\n"+
				"\tMaximize mu: "+maximize_mu/1000+"s"+"\n"+
				"\tResolve duplicate: "+resolve_duplicate/1000+"s"+"\n"+
				"\tMaximize fuzziness: "+maximize_fuzz/1000+"s"+"\n"+
				"\tPairwise comparison: "+pairwise/1000+"s"+"\n"+
				"\tMaximize tau: "+maximize_tau/1000+"s"+"\n"+
				"\tMaximize pi: "+maximize_pi/1000+"s"+"\n"+
				"\tCompute LL: "+compute_LL/1000+"s"+"\n"+
				"\textra: "+extra/1000+"s"+"\n");
	}
}
