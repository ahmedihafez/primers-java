package org.primer3.libprimer3;

/* This struct captures informatin about the similarity of
   an oligo (primer) to elements in a mispriming (repeat)
   library (which is read in from a fasta file). */
public class RepSim {
	/** Name of the sequence format with maximum
       similarity to the oligo. */
	String name;       

	/** 
     * The minimum score in slot 'score' (below).
     * (Used when the objective function involves
     * minimization of mispriming possibilities.)
     */
	short min;        
	/** The index of the maximum score in slot 'score' (below). */
	short max;        
	
	/** 
     * Array of similarity (i.e. false-priming) scores,
     * one for each entry in the 'repeat_lib' slot
     * of the primargs struct.  In libprimer3.c,
     * score is set to NULL to indicate that
     * the rep_sim structure is uninitialized.
     */
	double[] score = null;    
}