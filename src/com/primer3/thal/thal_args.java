package com.primer3.thal;


public class thal_args {
	int debug; /* if non zero, print debugging info to stderr */
	public thal_alignment_type type; /* one of the
		      1 THAL_ANY, (by default)
		      2 THAL_END1,
		      3 THAL_END2,
		      4 THAL_HAIRPIN */
	public int maxLoop;  /* maximum size of loop to consider; longer than 30 bp are not allowed */
	public double mv; /* concentration of monovalent cations */
	public double dv; /* concentration of divalent cations */
	public double dntp; /* concentration of dNTP-s */
	public double dna_conc; /* concentration of oligonucleotides */
	public double temp; /* temperature from which hairpin structures will be calculated */
	public int temponly; /* if non zero, print only temperature to stderr */
	public int dimer; /* if non zero, dimer structure is calculated */
	
	
	public void set_thal_default_args() {
		this.debug = 0;
		this.type = thal_alignment_type.thal_any; /* thal_alignment_type THAL_ANY */
		this.maxLoop = thal.MAX_LOOP;
		this.mv = 50; /* mM */
		this.dv = 0.0; /* mM */
		this.dntp = 0.8; /* mM */
		this.dna_conc = 50; /* nM */
		this.temp = thal.TEMP_KELVIN; /* Kelvin */
		this.temponly = 1; /* return only melting temperature of predicted structure */
		this.dimer = 1; /* by default dimer structure is calculated */
	}
	
	/**
	 * Set default args for oligo
	 */
	public void set_thal_oligo_default_args(){
		this.debug = 0;
		this.type = thal_alignment_type.thal_any; /* thal_alignment_type THAL_ANY */
		this.maxLoop = thal.MAX_LOOP;
		this.mv = 50; /* mM */
		this.dv = 0.0; /* mM */
		this.dntp = 0.0; /* mM */
		this.dna_conc = 50; /* nM */
		this.temp = thal.TEMP_KELVIN; /* Kelvin */
		this.temponly = 1; /* return only melting temperature of predicted structure */
		this.dimer = 1; /* by default dimer structure is calculated */
	}
	
	
}
