package com.primer3.libprimer3;

import com.primer3.p3_seq_lib.seq_lib;

public class args_for_one_oligo_or_primer {
	  
	public seq_lib       repeat_lib = null;
	public oligo_weights weights = new oligo_weights();
	double opt_tm;
	double min_tm;
	double max_tm;
	double opt_gc_content;
	double max_gc;
  double min_gc;

  /* Warning: also used for product Tm (TO DO: factor this out) */
  public double salt_conc;

  public double divalent_conc;
  /*
    DIVALENT_CONC and DNTP_CONC are both needed for enabling use of
    divalent cations for calculation of melting temperature of short
    and long oligos.  The formula for converting the divalent cations
    to monovalent cations is in the paper [Ahsen von N, Wittwer CT,
    Schutz E (2001) "Oligonucleotide Melting Temperatures under PCR
    Conditions: Nearest-Neighbor Corrections for Mg^2+,
    Deoxynucleotide Triphosphate, and Dimethyl Sulfoxide
    Concentrations with Comparision to Alternative Empirical
    Formulas", Clinical Chemistry 47:1956-61
    http://www.clinchem.org/cgi/content/full/47/11/1956] The default
    is 0.0.  (New in v. 1.1.0, added by Maido Remm and Triinu
    Koressaar.)
  */

  public double dntp_conc;

  public double dna_conc;
  int    num_ns_accepted;
  int    opt_size;
  int    min_size;
  int    max_size;
  int    max_poly_x;      /* 
                       * Maximum length of mononucleotide sequence in an
                       * oligo.
                       */

  int    min_end_quality;
  int    min_quality;       /* Minimum quality permitted for oligo sequence.*/

  double max_self_any;  
  double max_self_end;
  double max_self_any_th;
  double max_self_end_th;
  double max_hairpin_th;
  double max_repeat_compl;   /* 
                          * Acceptable complementarity with repeat
                          * sequences.
                          */

  	double max_template_mispriming;
	double max_template_mispriming_th;

	// orginal char *
	char[] must_match_five_prime;
	// orginal char *
	char[] must_match_three_prime;
	/* Primers and Oligos must match this 5 prime and 3 prime sequences.
	     This allows to select a set of primer pair with an identical 3' end
	     to avoid primer dimers. On the 5 prime end bases quenching
	     flourochromes can be avoided.
	     The sequence must be 5 nucletides long and can contain the following letters:

	     N Any nucleotide

	     A Adenine
	     G Guanine
	     C Cytosine
	     T Thymine

	     R Purine (A or G)
	     Y Pyrimidine (C or T)
	     W Weak (A or T)
	     S Strong (G or C)
	     M Amino (A or C)
	     K Keto (G or T)
	     B Not A (G or C or T)
	     H Not G (A or C or T)
	     D Not C (A or G or T)
	     V Not T (A or G or C)
	  */
	
	/**
	 * PRIMER_MIN_QUALITY
	 * @param datum
	 */
	public void set_min_quality(String datum) {
		this.min_quality = Integer.parseInt(datum);		
	}
	
	/**
	 * PRIMER_MAX_SELF_ANY
	 * @param datum
	 */
	public void set_max_self_any(String datum) {
		this.max_self_any = Double.parseDouble(datum);		
	}
	
	/**
	 * PRIMER_MAX_NS_ACCEPTED
	 * @param datum
	 */
	public void set_num_ns_accepted(String datum) {
		this.num_ns_accepted = Integer.parseInt(datum);
	}
	
	/**
	 * PRIMER_DNA_CONC
	 * @param datum
	 */
	public void set_dna_conc(String datum) {
		this.dna_conc = Double.parseDouble(datum);		

	}
	
	/**
	 * PRIMER_DNTP_CONC
	 * @param datum
	 */
	public void set_dntp_conc(String datum) {
		this.dntp_conc = Double.parseDouble(datum);		
	}
	
	/**
	 * PRIMER_SALT_DIVALENT
	 * @param datum
	 */
	public void set_divalent_conc(String datum) {
		this.divalent_conc = Double.parseDouble(datum);
	}
	
	/**
	 * PRIMER_SALT_MONOVALENT
	 * @param datum
	 */
	public void set_salt_conc(String datum) {
		this.salt_conc = Double.parseDouble(datum);		
		
	}

	/**
	 * PRIMER_MAX_GC
	 * @param datum
	 */
	public void set_max_gc(String datum) {
		this.max_gc = Double.parseDouble(datum);		
		
	}
	
	/**
	 * PRIMER_MIN_GC
	 * @param datum
	 */
	public void set_min_gc(String datum) {
		this.min_gc = Double.parseDouble(datum);		
	}
		
	/**
	 * PRIMER_MAX_TM
	 * @param datum
	 */
	public void set_max_tm(String datum) {
		this.max_tm = Double.parseDouble(datum);		
	}
	
	/**
	 * PRIMER_MIN_TM
	 * @param datum
	 */
	public void set_min_tm(String datum) {
		this.min_tm = Double.parseDouble(datum);
		
	}
	
	/**
	 * PRIMER_OPT_GC_PERCENT
	 * @param datum
	 */
	public void set_opt_gc_content(String datum) {
		this.opt_gc_content = Double.parseDouble(datum);		
	}
	
	/**
	 * PRIMER_OPT_TM
	 * @param datum
	 */
	public void set_opt_tm(String datum) {
		this.opt_tm = Double.parseDouble(datum);		
	}
	
	/**
	 * PRIMER_MAX_POLY_X
	 * @param datum
	 */
	public void set_max_poly_x(String datum) {
		this.max_poly_x  = Integer.parseInt(datum);		
		
	}
	
	/**
	 * PRIMER_MAX_SIZE
	 * @param datum
	 */
	public void set_max_size(String datum) {
		this.max_size  = Integer.parseInt(datum);		
		
	}
	
	/**
	 * PRIMER_MIN_SIZE
	 * @param datum
	 */
	public void set_min_size(String datum) {
		this.min_size  = Integer.parseInt(datum);		
		
	}
	
	/**
	 * PRIMER_OPT_SIZE
	 * @param datum
	 */
	public void set_opt_size(String datum) {
		this.opt_size  = Integer.parseInt(datum);		
	}
	
	/**
	 * PRIMER_MIN_END_QUALITY
	 * @param datum
	 */
	public void set_min_end_quality(String datum) {
		this.min_end_quality =Integer.parseInt(datum);		
	}
	
	/**
	 * PRIMER_MAX_SELF_ANY_TH
	 * @param datum
	 */
	public void set_max_self_any_th(String datum) {
		this.max_self_any_th = Double.parseDouble(datum);
	}
	
	/**
	 * PRIMER_MAX_SELF_END_TH
	 * @param datum
	 */
	public void set_max_self_end_th(String datum) {
		this.max_self_end_th = Double.parseDouble(datum);
	}
	
	/**
	 * PRIMER_MAX_HAIRPIN_TH
	 * @param datum
	 */
	public void set_max_hairpin_th(String datum) {
		this.max_hairpin_th = Double.parseDouble(datum);
	}
	
	/**
	 * 
	 * @param datum
	 */
	public void set_max_repeat_compl(String datum) {
		this.max_repeat_compl = Double.parseDouble(datum);
	}
	
	public void set_max_template_mispriming(String datum) {
		this.max_template_mispriming = Double.parseDouble(datum);
	}
	
	public void set_max_template_mispriming_th(String datum) {
		this.max_template_mispriming_th = Double.parseDouble(datum);
	}
	
	public void set_must_match_five_prime(String datum) {
		this.must_match_five_prime = datum.toCharArray();
	}
	
	public void set_must_match_three_prime(String datum) {
		this.must_match_three_prime = datum.toCharArray();		
	}

	/**
	 * PRIMER_MAX_SELF_END/PRIMER_INTERNAL_MAX_SELF_END
	 * @param datum
	 */
	public void set_max_self_end(String datum) {
		this.max_self_end = Double.parseDouble(datum);		
	}
}