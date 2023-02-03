package org.seqcode.projects.sem.events;

import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.gseutils.ArgParser;

public class EventsConfig {
	protected GenomeConfig gconfig;
	protected Genome gen = null;
	protected BindingModel defaultModel=null;
	protected String[] args;
	
	public final boolean CALC_EVENTS_LL=false; //Calculate component-wise log-likelihoods during ML
	
	public EventsConfig(GenomeConfig gcon, String[] arguments) {
		gconfig = gcon;
		gen = gconfig.getGenome();
		this.args = arguments;
		ArgParser ap = new ArgParser(args);
	}
	
	//Accessors
	public Genome getGenome() {return gen;}
	public BindingModel getDefaultBindingModel() {return defaultModel;}
	
	/**
	 * returns a string describing the arguments handled by this paper
	 */
	public String getArgsList() {
		return(new String("" +
					"BindingModels:\n"
				));
	}
}
