package org.seqcode.projects.sem.events;

import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;

public class EventsConfig {
	protected GenomeConfig gconfig;
	protected Genome gen = null;
	protected BindingModel defaultModel=null;
	
	//Accessors
	public BindingModel getDefaultBindingModel() {return defaultModel;}
}
