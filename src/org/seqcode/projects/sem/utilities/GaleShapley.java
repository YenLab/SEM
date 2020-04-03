package org.seqcode.projects.sem.utilities;

import java.util.Map;
import java.util.List;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Queue;
import java.io.*;

import org.seqcode.gseutils.Pair;

/**
 * Gale-Shapley algorithm used to find pairwise nucleosomes across different conditions
 * @author Jianyu Yang
 *
 */

public class GaleShapley<T> {
	protected Map<T, List<T>> M;
	protected Map<T, List<T>> W;
	protected Map<T, T> wives;		// get husband according to wife
	protected Map<T, T> husbands;	// get wife according to husband
	protected List<Pair<T, T>> pairs;	// car: man, cdr: woman
	protected List<T> mspecific;
	protected List<T> wspecific;
	
	protected Map<T, Map<T, Integer>> mrank;
	protected Map<T, Map<T, Integer>> wrank;
	
	protected boolean isMatch = false;
	
	public GaleShapley(Map<T, List<T>> men, Map<T, List<T>> women) {
		this.M = men;
		this.W = women;
		
		husbands = new HashMap<T, T>();
		wives = new HashMap<T, T>();
		pairs = new ArrayList<Pair<T, T>>();
		mspecific = new ArrayList<T>();
		wspecific = new ArrayList<T>();
		
		// index spousal preferences at initialization to avoid expensive lookups when matching
		mrank = new HashMap<T, Map<T, Integer>>();
		wrank = new HashMap<T, Map<T, Integer>>();
		
		if(M.size()==0 || W.size()==0) {
			if(M.size()==0)
				wspecific.addAll(W.keySet());
			else
				mspecific.addAll(M.keySet());
			isMatch = true;
		} else {
			for(T m: M.keySet()) {
				if(M.get(m).size()>0) {
					mrank.put(m, new HashMap<T, Integer>());
					for(int i=0; i<M.get(m).size(); i++) {
						mrank.get(m).put(M.get(m).get(i), i);
					}
				} else {
					mspecific.add(m);
				}
			}
			for(T w: W.keySet()) {
				if(W.get(w).size()>0) {
					wrank.put(w, new HashMap<T, Integer>());
					for(int i=0; i<W.get(w).size(); i++) {
						wrank.get(w).put(W.get(w).get(i), i);
					}
				} else {
					wspecific.add(w);
				}
			}
		}
	}
	
	//Accessors
	public List<Pair<T, T>> getPairs() {return pairs;}
	public List<T> getMenSpecific() {return mspecific;}
	public List<T> getWoMenSpecific() {return wspecific;}
	public boolean isManSpecific(T query) {return mspecific.contains(query);}
	public boolean isWomanSpecific(T query) {return wspecific.contains(query);}
	public T getWife(T query) {return husbands.get(query);}
	
	protected boolean prefers(T w, T m, T h) {
		// Test whether w prefers m over h.
		return wrank.get(w).get(m) < wrank.get(w).get(h);
	}
	
	protected T after(T m, T w) {
		// Return the woman favored by m after w.
		int i = mrank.get(m).get(w) + 1;
		if(i<M.get(m).size()) {
			return M.get(m).get(i);
		} else {
			return null;
		}
	}
	
	public void match() {
		if(!isMatch) {
			Map<T, T> next= new HashMap<T, T>();
			for(T m: mrank.keySet()) {
				next.put(m, M.get(m).get(0));
			}
			Queue<T> unpairedMen = new LinkedList<T>();
			unpairedMen.addAll(mrank.keySet());
			match(unpairedMen, next, wives);
		}
	}
	
	protected void match(Queue<T> unpairedMen, Map<T, T> next, Map<T, T> wives) {
		// Try to match all men with their next preferred spouse
		if(unpairedMen.size()==0) {
			for(T w: wives.keySet()) {
				pairs.add(new Pair<T, T>(wives.get(w), w));
				husbands.put(wives.get(w), w);
			}
			// add unpaired women to the `wspecific`
			for(T w: W.keySet()) {
				if(!wives.containsKey(w)) {
					wspecific.add(w);
				}
			}
			return;
		}
		T m = unpairedMen.poll();
		T w = next.get(m);
		if(w!=null) {
			next.put(m, after(m, w));
			if(wives.containsKey(w)) {
				T h = wives.get(w);
				if(prefers(w, m, h)) {
					unpairedMen.offer(h);
					wives.put(w, m);
				} else {
					unpairedMen.offer(m);
				}
			} else {
				wives.put(w, m);
			}
		} else {
			mspecific.add(m);
		}
		match(unpairedMen, next, wives);
	}
	
	public boolean is_stable() {
		for(T w: wives.keySet()) {
			T m = wives.get(w);
			int i = M.get(m).indexOf(w);
			for(int index=0; index<i; index++) {
				T p = M.get(m).get(index);
				T h = wives.get(p);
				if(W.get(p).indexOf(m) < W.get(p).indexOf(h)) {
					System.out.println("Unstable pair results!");
					return false;
				}
			}
		}
		return true;
	}
	
	public void print_results() {
		System.out.println("Couples:");
		for(Pair<T, T> p: pairs) {
			System.out.println("\t"+p);
		}
		System.out.println("lonely men:");
		for(T m: mspecific) {
			System.out.println("\t"+m);
		}
		System.out.println("lonely women:");
		for(T w: wspecific) {
			System.out.println("\t"+w);
		}
	}
	
	public static void main(String[] args) {
		String menFile = "D:\\Dropbox\\Yen lab\\2017\\Analysis\\Code\\python\\practice\\men.txt";
		String womenFile = "D:\\Dropbox\\Yen lab\\2017\\Analysis\\Code\\python\\practice\\women.txt";
		Map<String, List<String>> men = new HashMap<String, List<String>>();
		Map<String, List<String>> women = new HashMap<String, List<String>>();
		try {
			BufferedReader mbr = new BufferedReader(new FileReader(menFile));
			BufferedReader wbr = new BufferedReader(new FileReader(womenFile));
			String line;
			while((line=mbr.readLine()) != null) {
				if(line.strip().startsWith("#"))
					continue;
				String[] split = line.strip().split(":");
				String key = split[0];
				List<String> values = new ArrayList<String>();
				for(String value: split[1].split(",")) {
					values.add(value.strip());
				}
				men.put(key, values);
			}
			while((line=wbr.readLine()) != null) {
				if(line.strip().startsWith("#"))
					continue;
				String[] split = line.strip().split(":");
				String key = split[0];
				List<String> values = new ArrayList<String>();
				for(String value: split[1].split(",")) {
					values.add(value.strip());
				}
				women.put(key, values);
			}
		} catch(Exception e) {
			e.printStackTrace();
		}
		
		GaleShapley gs = new GaleShapley(men, women);
		gs.match();
		System.out.println(gs.is_stable());
		gs.print_results();
	}
}
