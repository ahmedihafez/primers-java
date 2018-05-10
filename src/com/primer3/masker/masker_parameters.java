package com.primer3.masker;

import java.util.List;


public class masker_parameters {
	/* strand to mask */
	public masking_direction mdir;
	
	/* primer failure rate cutoff used in primer design, 
	 * potential locations in a sequence for primers with PCR 
	 * failure rate over the given cutoff are masked
	 *(see function calculate_scores() from masker.c) */
	public double failure_rate;
	
	/* absolute value cutoff, this can be used for masking all the k-mers in a sequence
	 * that have the frequency over abs_cutoff in a k-mer list */
	int abs_cutoff;
	
	/* number of nucleotides masked in 5' and 3' direction with respect
	 * to the 3' end of a primer */
	public int nucl_masked_in_5p_direction;
	public int nucl_masked_in_3p_direction;
	
	/* If masker is used as a separate application then always print_sequence=1, 
	 * i.e the output is sent to stdout.
	 * If print_sequence=0 the output is written in a string variable and can be forwarded
	 * to the next function */
	public boolean print_sequence;
	
	/* if do_soft_masking=1, masked nucleotides and converted to lower-case, else 
	 * masked nucleotide are converted to masking_char ('N' by default) */
	public boolean do_soft_masking;
	char masking_char;
	
	/* size of the masking window */
	public int window_size;
	
	/* number of k-mer lists used in the masking formula */
	public int nlists;
	/* k-mer lists and all their parameters which are used in the masking formula */
	public String list_prefix;
	public List<formula_parameters> fp;
	public double formula_intercept;
	
	public void set_nucl_masked_in_3p_direction(String datum) {
		this.nucl_masked_in_3p_direction = Integer.parseInt(datum);
	}
	
	public void set_list_prefix(String datum) {
		this.list_prefix = datum;
	}

	public void set_nucl_masked_in_5p_direction(String datum) {
		this.nucl_masked_in_5p_direction = Integer.parseInt(datum);
	}
	
	
	
	public void set_failure_rate(String datum) {
		this.failure_rate = Double.parseDouble(datum);		
	} 
	
	
	



}
