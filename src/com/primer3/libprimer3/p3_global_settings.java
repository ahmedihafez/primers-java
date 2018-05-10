package com.primer3.libprimer3;

import java.util.regex.Pattern;

import com.primer3.masker.masker;
import com.primer3.masker.masker_parameters;
import com.primer3.masker.masking_direction;
import com.primer3.oligotm.salt_correction_type;
import com.primer3.oligotm.tm_method_type;

public class p3_global_settings {

	/* ================================================== */
	/* Arguments that control behavior of choose_primers() */

	/** 
	 * 2 if left primer only, 3 if right primer
	 * only, 4 if internal oligo only.  
	 */
	public task   primer_task; 

	public boolean    pick_left_primer;
	public boolean    pick_right_primer;
	public boolean    pick_internal_oligo;

	/* TODO, See if this can be factored out. */
	int    file_flag;   

	/** 
	 * The index of the first base in the input
	 * sequence.  This parameter is ignored within
	 * pr_choice; pr_choice's caller must assure that
	 * all indexes are 0-based.  However, this
	 * parameter should used by output routines to
	 * adjust base indexes.
	 */
	public int    first_base_index;  

	/** 
	 * If non-0 then turn characters other than
	 * [ATGCNatgcn] into N.
	 */
	boolean    liberal_base;   

	/** The number of best primer pairs to return. */
	public int    num_return; 

	/** Pick even if input primer or oligos
  violate constraints. */
	public boolean    pick_anyway;    

	/** If non-0, treat ambiguity codes in a mispriming/mishyb
	 * library as representing a consensus.  So, for example,
	 * S would match C or G.  N would match any nucleotide.
	 * It turns out that this _not_ what one normally wants,
	 * since many libraries contain strings of N, which then
	 * match every oligo (very bad).
	 */
	boolean    lib_ambiguity_codes_consensus;


	int    quality_range_min;
	int    quality_range_max;

	/* ================================================== */
	/* Arguments for individual oligos and/or primers */

	public args_for_one_oligo_or_primer p_args = new args_for_one_oligo_or_primer();
	public args_for_one_oligo_or_primer o_args = new args_for_one_oligo_or_primer();

	/** 
	 * tm_santalucia added by T.Koressaar for updated table
	 * thermodynamics.  Specifies details of melting temperature
	 * calculation.  (New in v. 1.1.0, added by Maido Remm and Triinu
	 * Koressaar.) 
	 * A value of 1 (recommended) directs primer3 to use the table of
	 * thermodynamic values and the method for melting temperature
	 * calculation suggested in the paper [SantaLucia JR (1998) "A
	 * unified view of polymer, dumbbell and oligonucleotide DNA
	 * nearest-neighbor thermodynamics", Proc Natl Acad Sci 95:1460-65
	 * http://dx.doi.org/10.1073/pnas.95.4.1460].
	 *  
	 * A value of 0 directs primer3 to a backward compatible calculation
	 * (in other words, the only calculation availble in previous
	 * version of primer3).
	 *  
	 * This backward compatible calculation uses the table of
	 * thermodynamic parameters in the paper [Breslauer KJ, Frank R,
	 * Bloecker H and Marky LA (1986) "Predicting DNA duplex stability
	 * from the base sequence" Proc Natl Acad Sci 83:4746-50
	 * http://dx.doi.org/10.1073/pnas.83.11.3746],
	 * and the method in the paper [Rychlik W, Spencer WJ and Rhoads
	 * RE (1990) "Optimization of the annealing temperature for DNA
	 * amplification in vitro", Nucleic Acids Res 18:6409-12
	 * http://www.pubmedcentral.nih.gov/articlerender.fcgi?tool=pubmed&pubmedid=2243783].
	 *   
	 * The default value is 0 only for backward compatibility.
	 */
	public tm_method_type tm_santalucia;  

	/** 
     salt_corrections added by T.Koressaar for salt correction for Tm
     calculation.  A value of 1 (recommended) directs primer3 to use
     the salt correction formula in the paper [SantaLucia JR (1998) "A
     unified view of polymer, dumbbell and oligonucleotide DNA
     nearest-neighbor thermodynamics", Proc Natl Acad Sci 95:1460-65
     http://dx.doi.org/10.1073/pnas.95.4.1460]

     A value of 0 directs primer3 to use the the salt correction
     formula in the paper [Schildkraut, C, and Lifson, S (1965)
     "Dependence of the melting temperature of DNA on salt
     concentration", Biopolymers 3:195-208 (not available on-line)].
     This was the formula used in previous version of primer3.

     A value of 2 directs primer3 to use the salt correction formula
     in the paper [Owczarzy R, You Y, Moreira BG, Manthey JA, Huang L,
     Behlke MA and Walder JA (2004) "Effects of sodium ions on DNA
     duplex oligomers: Improved predictions of melting temperatures",
     Biochemistry 43:3537-54 http://dx.doi.org/10.1021/bi034621r].

     The default is 0 only for backward compatibility.
	 */
	public salt_correction_type salt_corrections; 


	/* ================================================== */
	/* Start of arguments applicable to primers but not
     oligos */

	double max_end_stability;
	/* The maximum value allowed for the delta
	 * G of disruption for the 5 3' bases of
	 * a primer.
	 */

	int    max_end_gc;
	/* The maximum number of allowed G or C
	 * for the 5 3' bases of a primer.
	 */

	int    gc_clamp;  /* Required number of GCs at 3' end. */

	/* ================================================== */
	/* Start of arguments related to primer and/or oligo
     location in the template. */

	public boolean lowercase_masking; 
	/* 
     lowercase_masking added by T.Koressaar. Enables design of primers
     from lowercase masked template.  A value of 1 directs primer3 to
     reject primers overlapping lowercase a base exactly at the 3'
     end.

     This property relies on the assumption that masked features
     (e.g. repeats) can partly overlap primer, but they cannot overlap the
     3'-end of the primer.  In other words, lowercase bases at other
     positions in the primer are accepted, assuming that the masked
     features do not influence the primer performance if they do not
     overlap the 3'-end of primer.
	 */

	public sequencing_parameters sequencing = new sequencing_parameters();  /* Used to calculate the position
					of sequencing primers */

	double outside_penalty; /* Multiply this value times the number of NTs
	 * from the 3' end to the the (unique) target to
	 * get the 'position penalty'.
	 * Meaningless if there are multiple targets
	 * or if the primer cannot be part of a pair
	 * that spans the target.
	 */

	double inside_penalty;  /* Multiply this value times the number of NT
	 * positions by which the primer overlaps
	 * the (unique) target to the 'position penalty'.
	 * Meaningless if there are multiple targets
	 * or if the primer cannot be part of a pair
	 * that spans the target.
	 */

	/* ================================================== */
	/* Arguments for primer pairs and products. */

	/* Warning: Use p3_empty_gs_product_size_range and
     p3_add_to_gs_product_size_range to set these next
     three slots. */
	/** Minimum product sizes. */
	int[]    pr_min = new int[LibPrimer3.PR_MAX_INTERVAL_ARRAY];
	/** Maximum product sizes. */
	int[]    pr_max = new int[LibPrimer3.PR_MAX_INTERVAL_ARRAY]; 
	/** 
	 * Number of product size intervals
	 * (i.e. number of elements in pr_min and
	 * pr_max)
	 */
	int    num_intervals;        

	int    product_opt_size;
	double product_max_tm;
	double product_min_tm;
	double product_opt_tm;
	double pair_max_template_mispriming;
	double pair_max_template_mispriming_th;
	double pair_repeat_compl;
	double pair_compl_any;
	double pair_compl_any_th;
	double pair_compl_end;
	double pair_compl_end_th;




	// refactor to true / false
	public static final int USE_THERMODYNAMICS_ALIGNMENT = 1;
	public static final int DONOT_USE_THERMODYNAMICS_ALIGNMENT = 0;

	/** 
     Enables to use approach of thermodynamical alignment for dimer
     and hairpin calculations.

     <br>
     <code>DONOT_USE_THERMODYNAMICS_ALIGNMENT</code> 0 = Use alignment *not* based on thermodynamics - the only dimer
     calculation approach until primer3 2.2.0.  No hairpins are
     calculated.
	<br>
     <code>USE_THERMODYNAMICS_ALIGNMENT</code> 1 = Use alignment based on thermodynamics. Hairpins are calculated
     TODO :: PR_ASSERT ( thermodynamic_oligo_alignment == 0 || == 1 ) 
	 */

	public int thermodynamic_oligo_alignment;


	/** 
     Enables to use approach of thermodynamical alignment for template 
     mispriming calculations.
     0 = Use alignment *not* based on thermodynamics.
     1 = Use alignment based on thermodynamics.
	 */
	public int thermodynamic_template_alignment;


	/** Maximum allowed difference between temperature of primer and
     temperature of product.  Cannot be calculated until product is
     known. */
	double max_diff_tm; 


	public pair_weights  pr_pair_weights = new pair_weights();

	int    min_left_three_prime_distance;
	/** Minimum number of base pairs between the 3' ends of any two left
     or any two right primers when returning num_return primer pairs.
     The objective is get 'truly different' primer pairs.

     Please see the user documentation (primer3_manual.htm) for
     PRIMER_{LEFT,RIGHT}_MIN_THREE_PRIME_DISTANCE.
	 */
	int    min_right_three_prime_distance;

	/** The number of basepairs
  the primer has to
  overlap an overlap
  junction. 
	 */
	int    min_5_prime_overlap_of_junction;   
	int    min_3_prime_overlap_of_junction;

	public boolean    mask_template;

	/** Turn on masking of the trimmed_orig_seq (added by M. Lepamets)*/
	public boolean    masking_parameters_changed;

	/**
    a struct containing all masking parameters
	 */
	public masker_parameters mp = new masker_parameters();

	/** dump fields for global settings and seq args if dump == 1 */
	private boolean dump;




	public static p3_global_settings p3_create_global_settings(){
		p3_global_settings r = new p3_global_settings();
		r.pr_set_default_global_args_2();
		return r;
	}

	public static p3_global_settings p3_create_global_settings_default_version_1(){
		p3_global_settings r = new p3_global_settings();
		r.pr_set_default_global_args_1();
		return r;
	}

	/*
	 * Write the default values for default_values=2
	 */
	private void pr_set_default_global_args_2() {
		this.pr_set_default_global_args_1();
		this.tm_santalucia                    = tm_method_type.santalucia_auto;
		this.salt_corrections                 = salt_correction_type.santalucia;
		this.thermodynamic_oligo_alignment    = 1;
		this.thermodynamic_template_alignment = 0;
		this.p_args.divalent_conc             = 1.5;
		this.p_args.dntp_conc                 = 0.6;
		this.lib_ambiguity_codes_consensus    = false;
	}



	/**
	 * Write the default values for default_values=1
	 */
	private void pr_set_default_global_args_1() {
		/* Arguments for primers ================================= */
		this.p_args.opt_size          = 20;
		this.p_args.min_size          = 18;
		this.p_args.max_size          = 27;

		this.p_args.opt_tm            = 60;
		this.p_args.min_tm            = 57;
		this.p_args.max_tm            = 63;

		this.p_args.min_gc            = 20.0;
		this.p_args.opt_gc_content    = LibPrimer3.DEFAULT_OPT_GC_PERCENT;
		this.p_args.max_gc            = 80.0;
		this.p_args.salt_conc         = 50.0;
		this.p_args.divalent_conc     = 0.0;
		this.p_args.dntp_conc         = 0.0;

		this.p_args.dna_conc          = 50.0;
		this.p_args.num_ns_accepted   = 0;
		this.p_args.max_self_any      = 8.0;
		this.p_args.max_self_end      = 3.0;
		this.p_args.max_self_any_th   = 47.0;
		this.p_args.max_self_end_th   = 47.0;
		this.p_args.max_hairpin_th    = 47.0;
		this.p_args.max_poly_x        = 5;
		this.p_args.max_repeat_compl  = 12.0;
		this.p_args.min_quality       = 0;
		this.p_args.min_end_quality   = 0;
		this.p_args.max_template_mispriming = LibPrimer3.PR_UNDEFINED_ALIGN_OPT;
		this.p_args.max_template_mispriming_th = LibPrimer3.PR_UNDEFINED_ALIGN_OPT;
		/* The following apply only to primers (and not to internal
		     oligos). */
		this.gc_clamp                 = 0;
		this.max_end_gc               = 5;

		/* Weights for objective functions for oligos and pairs. */
		this.p_args.weights.compl_any     = 0;
		this.p_args.weights.compl_any_th  = 0;
		this.p_args.weights.compl_end     = 0;
		this.p_args.weights.compl_end_th  = 0;
		this.p_args.weights.end_quality   = 0;
		this.p_args.weights.end_stability = 0;
		this.p_args.weights.gc_content_gt = 0;
		this.p_args.weights.gc_content_lt = 0;
		this.p_args.weights.hairpin_th    = 0;
		this.p_args.weights.length_gt     = 1;
		this.p_args.weights.length_lt     = 1;
		this.p_args.weights.num_ns        = 0;
		this.p_args.weights.pos_penalty   = 1;
		this.p_args.weights.repeat_sim    = 0;
		this.p_args.weights.seq_quality   = 0;
		this.p_args.weights.temp_cutoff   = 5;
		this.p_args.weights.temp_gt       = 1;
		this.p_args.weights.temp_lt       = 1;
		this.p_args.weights.template_mispriming = 0.0;
		this.p_args.weights.template_mispriming_th = 0.0;
		this.p_args.must_match_five_prime  = null;
		this.p_args.must_match_three_prime = null;
		/* End of weights for objective functions for oligos and pairs. */

		/* End of arguments for primers =========================== */

		this.max_diff_tm         = 100.0;
		this.tm_santalucia       = tm_method_type.breslauer_auto;
		this.salt_corrections    = salt_correction_type.schildkraut;
		this.pair_compl_any      = 8.0;
		this.pair_compl_end      = 3.0;
		this.pair_compl_any_th   = 47.0;
		this.pair_compl_end_th   = 47.0;
		this.thermodynamic_oligo_alignment = 0;
		this.thermodynamic_template_alignment = 0;
		this.liberal_base        = false;
		this.primer_task         = task.generic;
		this.pick_left_primer    = true;
		this.pick_right_primer   = true;
		this.pick_internal_oligo = false;
		this.first_base_index    = 0;
		this.num_return          = 5;
		this.pr_min[0]           = 100;
		this.pr_max[0]           = 300;
		this.num_intervals       = 1;
		this.pair_repeat_compl   = 24.0;
		this.quality_range_min   = 0;
		this.quality_range_max   = 100;
		this.outside_penalty     = LibPrimer3.PR_DEFAULT_OUTSIDE_PENALTY;
		this.inside_penalty      = LibPrimer3.PR_DEFAULT_INSIDE_PENALTY;
		this.max_end_stability   = 100.0;
		this.lowercase_masking   = false;
		this.product_max_tm      = LibPrimer3.PR_DEFAULT_PRODUCT_MAX_TM;
		this.product_min_tm      = LibPrimer3.PR_DEFAULT_PRODUCT_MIN_TM;
		this.product_opt_tm      = LibPrimer3.PR_UNDEFINED_DBL_OPT;
		this.product_opt_size    = LibPrimer3.PR_UNDEFINED_INT_OPT;
		this.pair_max_template_mispriming = LibPrimer3.PR_UNDEFINED_ALIGN_OPT;
		this.pair_max_template_mispriming_th = LibPrimer3.PR_UNDEFINED_ALIGN_OPT;
		this.o_args.opt_size        = 20;
		this.o_args.min_size        = 18;
		this.o_args.max_size        = 27;
		this.o_args.opt_tm          = 60.0;
		this.o_args.min_tm          = 57.0;
		this.o_args.max_tm          = 63.0;
		this.o_args.min_gc          = 20.0;
		this.o_args.max_gc          = 80.0;
		this.o_args.opt_gc_content  = LibPrimer3.DEFAULT_OPT_GC_PERCENT;
		this.o_args.max_poly_x      = 5;
		this.o_args.salt_conc       = 50.0;
		this.o_args.divalent_conc   = 0.0;
		this.o_args.dntp_conc       = 0.0;
		this.o_args.dna_conc        = 50.0;
		this.o_args.num_ns_accepted = 0;
		this.o_args.max_self_any    = 12.0;
		this.o_args.max_self_end    = 12.0;
		this.o_args.max_self_any_th = 47.0;
		this.o_args.max_self_end_th = 47.0;
		this.o_args.max_hairpin_th  = 47.0;
		this.o_args.max_repeat_compl= 12.0;

		this.o_args.min_quality           = 0;
		this.o_args.min_end_quality       = 0;
		this.o_args.max_template_mispriming = LibPrimer3.PR_UNDEFINED_ALIGN_OPT;
		this.o_args.max_template_mispriming_th = LibPrimer3.PR_UNDEFINED_ALIGN_OPT;
		this.o_args.weights.temp_gt       = 1;
		this.o_args.weights.temp_lt       = 1;
		this.o_args.weights.length_gt     = 1;
		this.o_args.weights.length_lt     = 1;
		this.o_args.weights.gc_content_gt = 0;
		this.o_args.weights.gc_content_lt = 0;
		this.o_args.weights.compl_any     = 0;
		this.o_args.weights.compl_end     = 0;
		this.o_args.weights.compl_any_th  = 0;
		this.o_args.weights.compl_end_th  = 0;
		this.o_args.weights.hairpin_th    = 0;
		this.o_args.weights.num_ns        = 0;
		this.o_args.weights.repeat_sim    = 0;
		this.o_args.weights.seq_quality   = 0;
		this.o_args.weights.end_quality   = 0;
		this.o_args.must_match_five_prime  = null;
		this.o_args.must_match_three_prime = null;

		this.pr_pair_weights.primer_quality  = 1;
		this.pr_pair_weights.io_quality      = 0;
		this.pr_pair_weights.diff_tm         = 0;
		this.pr_pair_weights.compl_any       = 0;
		this.pr_pair_weights.compl_end       = 0;
		this.pr_pair_weights.compl_any_th    = 0;
		this.pr_pair_weights.compl_end_th    = 0;
		this.pr_pair_weights.temp_cutoff     = 5;
		this.pr_pair_weights.repeat_sim      = 0;

		this.pr_pair_weights.product_tm_lt   = 0;
		this.pr_pair_weights.product_tm_gt   = 0;
		this.pr_pair_weights.product_size_lt = 0;
		this.pr_pair_weights.product_size_gt = 0;

		this.lib_ambiguity_codes_consensus   = true;
		/*  Set to 1 for backward compatibility. This _NOT_ what
		      one normally wants, since many libraries contain
		      strings of N, which then match every oligo (very bad).
		 */

		this.min_left_three_prime_distance   = -1;
		this.min_right_three_prime_distance  = -1;

		this.sequencing.lead                 = 50;
		this.sequencing.spacing              = 500;
		this.sequencing.interval             = 250;
		this.sequencing.accuracy             = 20;

		this.min_5_prime_overlap_of_junction = 7;
		this.min_3_prime_overlap_of_junction = 4;

		this.mask_template                   = false;
		this.masking_parameters_changed      = false;
		this.mp.mdir                         = masking_direction.both_separately;
		this.mp.failure_rate                 = 0.1;
		this.mp.nucl_masked_in_5p_direction  = 1;
		this.mp.nucl_masked_in_3p_direction  = 0;
		this.mp.print_sequence               = false;
		this.mp.do_soft_masking              = true;
		this.mp.nlists                       = masker.DEFAULT_NLISTS;
		this.mp.list_prefix                  = masker.DEFAULT_LIST_FILE_PREFIX;
		this.mp.fp                           = null;
		this.mp.formula_intercept            = masker.DEFAULT_INTERCEPT;  		
	}





	/**
	 * @return the dump
	 */
	public boolean isDump() {
		return dump;
	}

	/**
	 * @param dump the dump to set
	 */
	public void setDump(boolean dump) {
		this.dump = dump;
	}


	public boolean PR_DEFAULT_POSITION_PENALTIES() {
		return (LibPrimer3.PR_DEFAULT_INSIDE_PENALTY == this.inside_penalty 
				&& LibPrimer3.PR_DEFAULT_OUTSIDE_PENALTY == this.outside_penalty);
	}

	public boolean _PR_DEFAULT_POSITION_PENALTIES() {
		return (LibPrimer3.PR_DEFAULT_INSIDE_PENALTY == this.inside_penalty  && LibPrimer3.PR_DEFAULT_OUTSIDE_PENALTY == this.outside_penalty);
	}

	public static p3_global_settings p3_create_global_settings(
			int default_version) {
		if (default_version == 1)
			return p3_create_global_settings_default_version_1();
		else if (default_version == 2) 
			return p3_create_global_settings();
		return null;
	}


	/**
	 * PRIMER_PAIR_MAX_DIFF_TM
	 * @param datum
	 */
	public void p3_set_pa_max_diff_tm(String datum) {
		this.max_diff_tm = Double.parseDouble(datum);		
	}  

	/**
	 * PRIMER_TM_FORMULA
	 */
	public void p3_set_pa_tm_method_type(String datum) {
		int typeInt = Integer.parseInt(datum);
		if(typeInt == tm_method_type.santalucia_auto.getValue() )
		{
			this.tm_santalucia = tm_method_type.santalucia_auto;
		}
		else if(typeInt == tm_method_type.breslauer_auto.getValue() )
		{
			this.tm_santalucia = tm_method_type.breslauer_auto;
		}

	}


	/**
	 * PRIMER_SALT_CORRECTIONS
	 * @param datum
	 */
	public void p3_set_pa_salt_corrections(String datum) {
		int typeInt = Integer.parseInt(datum);

		if(typeInt == salt_correction_type.owczarzy.getValue() )
		{
			this.salt_corrections = salt_correction_type.owczarzy;
		} 
		else if(typeInt == salt_correction_type.santalucia.getValue() ) 
		{
			this.salt_corrections = salt_correction_type.santalucia;
		}
		else if(typeInt == salt_correction_type.schildkraut.getValue() ) 
		{
			this.salt_corrections = salt_correction_type.schildkraut;
		}
		else {
			// FIXME :: what to do in this case
		}
	}


	/**
	 * PRIMER_PRODUCT_SIZE_RANGE
	 * @param datum
	 */
	public void p3_set_pa_product_size(String datum)
	{
		// TODO :: catch ex and clear values then throws it again
		String intervalSep = " ";
		String numSep = "-";
		String[] intervalStrs = datum.split(Pattern.quote(intervalSep));
		this.num_intervals=0;
		for(String intervalStr : intervalStrs){
			if(this.num_intervals >=  LibPrimer3.PR_MAX_INTERVAL_ARRAY)
				throw new IndexOutOfBoundsException("Too many elements for tag ");
			String[] intervalNums = intervalStr.split(Pattern.quote(numSep));
			this.pr_min[this.num_intervals] = Integer.parseInt(intervalNums[0]);
			this.pr_max[this.num_intervals] = Integer.parseInt(intervalNums[1]);
			this.num_intervals++;
		}		
	}

	/**
	 * PRIMER_PRODUCT_OPT_SIZE
	 * @param datum
	 */
	public void p3_set_pa_product_opt_size(String datum) {
		this.product_opt_size = Integer.parseInt(datum);
	}

	/**
	 * PRIMER_PAIR_MAX_COMPL_ANY
	 * @param datum
	 */
	public void p3_set_pa_pair_compl_any(String datum) {
		this.pair_compl_any = Double.parseDouble(datum);
	}

	/**
	 * PRIMER_PAIR_MAX_COMPL_END
	 * @param datum
	 */
	public void p3_set_pa_pair_compl_end(String datum) {
		this.pair_compl_end = Double.parseDouble(datum);
	}

	/**
	 * PRIMER_PAIR_MAX_COMPL_ANY_TH
	 * @param datum
	 */
	public void p3_set_pa_pair_compl_any_th(String datum) {
		this.pair_compl_any_th = Double.parseDouble(datum);
	}

	/**
	 * PRIMER_PAIR_MAX_COMPL_END_TH
	 * @param datum
	 */
	public void p3_set_pa_pair_compl_end_th(String datum) {
		this.pair_compl_end_th = Double.parseDouble(datum);
	}

	/**
	 * PRIMER_PICK_ANYWAY
	 * @param datum
	 */
	public void p3_set_pa_pick_anyway(String datum) {
		if(Integer.parseInt(datum) == 0 )
			this.pick_anyway = false;
		else
			this.pick_anyway = true;
	}

	/**
	 * PRIMER_GC_CLAMP
	 * @param datum
	 */
	public void p3_set_pa_gc_clamp(String datum) {
		this.gc_clamp = Integer.parseInt(datum);
	}

	/**
	 * PRIMER_MAX_END_GC
	 * @param datum
	 */
	public void p3_set_pa_max_end_gc(String datum) {
		this.max_end_gc = Integer.parseInt(datum);		
	}

	/**
	 * PRIMER_LIBERAL_BASE
	 * @param datum
	 */
	public void p3_set_pa_liberal_base(String datum) {

		if(Integer.parseInt(datum) != 0)
			this.liberal_base = true;
		else
			this.liberal_base = false;
	}

	/**
	 * PRIMER_FIRST_BASE_INDEX
	 * @param datum
	 */
	public void p3_set_pa_first_base_index(String datum) {
		this.first_base_index = Integer.parseInt(datum);
	}

	/**
	 * PRIMER_NUM_RETURN
	 * @param datum
	 */
	public void p3_set_pa_num_return(String datum) {
		this.num_return = Integer.parseInt(datum);
	}

	/**
	 * PRIMER_QUALITY_RANGE_MIN
	 * @param value
	 */
	public void set_min_left_three_prime_distance(int value) {
		// TODO Auto-generated method stub

	}
	public void set_min_right_three_prime_distance(int value) {
		// TODO Auto-generated method stub

	}
	public void p3_set_pa_min_left_three_prime_distance(String  datum) {
		// TODO Auto-generated method stub

	}
	public void p3_set_pa_min_right_three_prime_distance(String datum) {
		// TODO Auto-generated method stub

	}
	public void p3_set_pa_quality_range_min(String datum) {
		this.quality_range_min = Integer.parseInt(datum);
	}

	/**
	 * PRIMER_QUALITY_RANGE_MAX
	 * @param datum
	 */
	public void p3_set_pa_quality_range_max(String datum) {
		this.quality_range_max = Integer.parseInt(datum);
	}

	/**
	 * PRIMER_PRODUCT_MAX_TM
	 * @param datum
	 */
	public void p3_set_pa_product_max_tm(String datum) {
		this.product_max_tm = Integer.parseInt(datum);		
	}

	/**
	 * PRIMER_PRODUCT_MIN_TM
	 * @param datum
	 */
	public void p3_set_pa_product_min_tm(String datum) {
		this.product_min_tm = Integer.parseInt(datum);
	}

	/**
	 * PRIMER_PRODUCT_OPT_TM
	 * @param datum
	 */
	public void p3_set_pa_product_opt_tm(String datum) {
		this.product_opt_tm = Integer.parseInt(datum);
	}

	/**
	 * PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION
	 * @param datum
	 */
	public void p3_set_pa_min_5_prime_overlap_of_junction(String datum) {
		// TODO Auto-generated method stub

	}
	public void p3_set_pa_min_3_prime_overlap_of_junction(String datum) {
		this.min_3_prime_overlap_of_junction = Integer.parseInt(datum);
	}

	/**
	 * PRIMER_PICK_RIGHT_PRIMER
	 * @param datum
	 */
	public void p3_set_pa_pick_right_primer(String datum) {
		if(Integer.parseInt(datum) == 0)
			this.pick_right_primer = false;
		else
			this.pick_right_primer = true;

	}

	/**
	 * PRIMER_PICK_INTERNAL_OLIGO
	 * @param datum
	 */
	public void p3_set_pa_pick_internal_oligo(String datum) {
		if(Integer.parseInt(datum) == 0)
			this.pick_internal_oligo = false;
		else
			this.pick_internal_oligo = true;

	}

	/**
	 * PRIMER_PICK_LEFT_PRIMER
	 * @param datum
	 */
	public void p3_set_pa_pick_left_primer(String datum) {
		if(Integer.parseInt(datum) == 0)
			this.pick_left_primer = false;
		else
			this.pick_left_primer = true;		
	}


	//	public void p3_set_gs_primer_self_end(String datum) {
	//		this.p_args.set_max_self_end(datum);// = Double.parseDouble(datum);		
	//	}



	/**
	 * PRIMER_INTERNAL_MAX_SELF_END - o_args.set_max_self_end
	 * @param datum
	 */
	public void p3_set_gs_primer_internal_oligo_self_end(String datum) {
		this.o_args.set_max_self_end(datum);
	}

	/**
	 * PRIMER_PAIR_MAX_LIBRARY_MISPRIMING
	 * @param datum
	 */
	public void p3_set_pa_pair_repeat_compl(String datum) {
		this.pair_repeat_compl = Double.parseDouble(datum);		
	}

	/**
	 * PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING
	 * @param datum
	 */

	public void p3_set_pa_pair_max_template_mispriming(String datum) {
		this.pair_max_template_mispriming = Double.parseDouble(datum);		
	}

	/**
	 * PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH
	 * @param datum
	 */

	public void p3_set_pa_pair_max_template_mispriming_th(String datum) {
		this.pair_max_template_mispriming_th = Double.parseDouble(datum);
	}

	/**
	 * PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS
	 * @param datum
	 */

	public void p3_set_pa_lib_ambiguity_codes_consensus(String datum) {
		if(Integer.parseInt(datum) != 0)
			this.lib_ambiguity_codes_consensus = true;
		else 
			this.lib_ambiguity_codes_consensus = false;
	}

	/**
	 * PRIMER_INSIDE_PENALTY
	 * @param datum
	 */
	public void p3_set_pa_inside_penalty(String datum) {
		this.inside_penalty = Double.parseDouble(datum);
	}

	/**
	 * PRIMER_OUTSIDE_PENALTY
	 * @param datum
	 */
	public void p3_set_pa_outside_penalty(String datum) {
		this.outside_penalty = Double.parseDouble(datum);
	}

	/**
	 * PRIMER_MAX_END_STABILITY
	 * @param datum
	 */
	public void p3_set_pa_max_end_stability(String datum) {
		this.max_end_stability = Double.parseDouble(datum);
	}



	/**
	 * PRIMER_LOWERCASE_MASKING
	 */
	public void p3_set_pa_lowercase_masking(String datum) {
		if(Integer.parseInt(datum)==0)
			this.lowercase_masking = false;
		else
			this.lowercase_masking = true;
	}
	/**
	 * PRIMER_LOWERCASE_MASKING
	 */
	public void p3_set_pa_lowercase_masking(boolean value) {
		this.lowercase_masking = value;
	}
	/**
	 * 
	 * @param datum
	 */
	public boolean p3_set_pa_mask_template(String datum) {
		if(Integer.parseInt(datum)==0)
			this.mask_template = false;
		else
			this.mask_template = true;

		return this.mask_template;
	}


	/**
	 * PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT
	 */
	public void p3_set_pa_thermodynamic_oligo_alignment(String datum) {
		this.thermodynamic_oligo_alignment = Integer.parseInt(datum);
	}

	/**
	 * PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT
	 * @param datum
	 */

	public void p3_set_pa_thermodynamic_template_alignment(String datum) {
		this.thermodynamic_template_alignment = Integer.parseInt(datum);		
	}

	/**
	 * PRIMER_MASK_KMERLIST_PREFIX
	 */
	public void set_masking_parameters_KmerLis_Prefix(String datum) {
		this.mp.set_list_prefix(datum);
		this.masking_parameters_changed = true;
	}



	public boolean _pr_need_template_mispriming() {
		return this.p_args.max_template_mispriming >= 0
				|| this.p_args.weights.template_mispriming > 0.0
				|| this._pr_need_pair_template_mispriming();
	}

	public boolean _pr_need_pair_template_mispriming() {
		return this.pair_max_template_mispriming >= 0
				|| this.pr_pair_weights.template_mispriming > 0.0;
	}

	public boolean _pr_need_template_mispriming_thermod() {
		return
				this.p_args.max_template_mispriming_th >= 0
				|| this.p_args.weights.template_mispriming_th > 0.0
				|| this._pr_need_pair_template_mispriming_thermod();
	}

	public boolean _pr_need_pair_template_mispriming_thermod()
	{   
		return this.pair_max_template_mispriming_th >= 0
				|| this.pr_pair_weights.template_mispriming_th > 0.0;
	}

	void p3_set_gs_primer_task(String task_tmp)
	{
		if (task_tmp.equalsIgnoreCase( "pick_pcr_primers")) {
			this.primer_task = task.generic;
			this.pick_left_primer = true;
			this.pick_right_primer = true;
			this.pick_internal_oligo = false;
		} else if (task_tmp.equalsIgnoreCase( "pick_pcr_primers_and_hyb_probe")) {
			this.primer_task = task.generic;
			this.pick_left_primer = true;
			this.pick_right_primer = true;
			this.pick_internal_oligo = true;
		} else if (task_tmp.equalsIgnoreCase( "pick_left_only")) {
			this.primer_task = task.generic;
			this.pick_left_primer = true;
			this.pick_right_primer = false;
			this.pick_internal_oligo = false;
		} else if (task_tmp.equalsIgnoreCase( "pick_right_only")) {
			this.primer_task = task.generic;
			this.pick_left_primer = false;
			this.pick_right_primer = true;
			this.pick_internal_oligo = false;
		} else if (task_tmp.equalsIgnoreCase( "pick_hyb_probe_only")) {
			this.primer_task = task.generic;
			this.pick_left_primer = false;
			this.pick_right_primer = false;
			this.pick_internal_oligo = true;
		} else if (task_tmp.equalsIgnoreCase( "generic")) {
			this.primer_task = task.generic;
		} else if (task_tmp.equalsIgnoreCase( "pick_detection_primers")) {
			this.primer_task = task.generic; /* Deliberate duplication for
					    backward compatibility. */
		} else if (task_tmp.equalsIgnoreCase( "pick_cloning_primers")) {
			this.primer_task = task.pick_cloning_primers;
		} else if (task_tmp.equalsIgnoreCase( "pick_discriminative_primers")) {
			this.primer_task = task.pick_discriminative_primers;
		} else if (task_tmp.equalsIgnoreCase( "pick_sequencing_primers")) {
			this.primer_task = task.pick_sequencing_primers;
		} else if (task_tmp.equalsIgnoreCase( "pick_primer_list")) {
			this.primer_task = task.pick_primer_list;
		} else if (task_tmp.equalsIgnoreCase("check_primers")) {
			this.primer_task = task.check_primers;
		}
	}



	public void p3_print_args() {
		p3_global_settings p = this;
		int i;
		System.out.format("=============\n");
		System.out.format("BEGIN GLOBAL ARGS\n") ;
		System.out.format("  primer_task %s\n", p.primer_task);
		System.out.format("  pick_left_primer %s\n", p.pick_left_primer);
		System.out.format("  pick_right_primer %s\n", p.pick_right_primer);
		System.out.format("  pick_internal_oligo %s\n", p.pick_internal_oligo);
		System.out.format("  file_flag %s\n", p.file_flag) ;
		System.out.format("  first_base_index %s\n", p.first_base_index);
		System.out.format("  liberal_base %s\n", p.liberal_base );
		System.out.format("  num_return %s\n", p.num_return) ;
		System.out.format("  pick_anyway %s\n", p.pick_anyway);
		System.out.format("  lib_ambiguity_codes_consensus %s\n", p.lib_ambiguity_codes_consensus) ;
		System.out.format("  quality_range_min %s\n", p.quality_range_min) ;
		System.out.format("  quality_range_max %s\n", p.quality_range_max) ;

		System.out.format("  tm_santalucia %s\n", p.tm_santalucia) ;
		System.out.format("  salt_corrections %s\n", p.salt_corrections) ;
		System.out.format("  max_end_stability %f\n", p.max_end_stability) ;
		System.out.format("  gc_clamp %s\n", p.gc_clamp) ;
		System.out.format("  max_end_gc %s\n", p.max_end_gc);
		System.out.format("  lowercase_masking %s\n", p.lowercase_masking) ;
		System.out.format("  thermodynamic_oligo_alignment %s\n", p.thermodynamic_oligo_alignment);
		System.out.format("  thermodynamic_template_alignment %s\n", p.thermodynamic_template_alignment);
		System.out.format("  outside_penalty %f\n", p.outside_penalty) ;
		System.out.format("  inside_penalty %f\n", p.inside_penalty) ;
		System.out.format("  number of product size ranges: %d\n", p.num_intervals);
		System.out.format("  product size ranges:\n");
		for (i = 0; i < p.num_intervals; i++) {
			System.out.format("  %d - %d \n", p.pr_min[i], p.pr_max[i]);
		}
		System.out.format("  product_opt_size %s\n", p.product_opt_size) ;
		System.out.format("  product_max_tm %f\n", p.product_max_tm) ;
		System.out.format("  product_min_tm %f\n", p.product_min_tm) ;
		System.out.format("  product_opt_tm %f\n", p.product_opt_tm) ;
		System.out.format("  pair_max_template_mispriming %f\n", p.pair_max_template_mispriming) ;
		System.out.format("  pair_max_template_mispriming_th %f\n", p.pair_max_template_mispriming_th) ;
		System.out.format("  pair_repeat_compl %f\n", p.pair_repeat_compl) ;
		System.out.format("  pair_compl_any %f\n", p.pair_compl_any) ;
		System.out.format("  pair_compl_end %f\n", p.pair_compl_end) ;
		System.out.format("  pair_compl_any_th %f\n", p.pair_compl_any_th) ;
		System.out.format("  pair_compl_end_th %f\n", p.pair_compl_end_th) ;

		System.out.format("  min_left_three_prime_distance %s\n", p.min_left_three_prime_distance) ;
		System.out.format("  min_right_three_prime_distance %s\n", p.min_right_three_prime_distance) ;
		System.out.format("  min_5_prime_overlap_of_junction %s\n", p.min_5_prime_overlap_of_junction);
		System.out.format("  min_3_prime_overlap_of_junction %s\n", p.min_3_prime_overlap_of_junction);
		System.out.format("  mask_template %s\n", p.mask_template);
		System.out.format("  failure_rate %f\n", p.mp.failure_rate);
		System.out.format("  nucl_masked_in_5p_direction %s\n", p.mp.nucl_masked_in_5p_direction);
		System.out.format("  nucl_masked_in_3p_direction %s\n", p.mp.nucl_masked_in_3p_direction);
		System.out.format("  list_prefix %s\n", p.mp.list_prefix);
		System.out.format("  dump %s\n", p.dump);

		p.pr_pair_weights.p3_print_args();
		System.out.format("\n\n");
		System.out.format("=============\n");
		System.out.format("BEGIN primer_args\n");
		System.out.format("begin oligo_weights\n");
		System.out.format("temp_gt %f\n", p.p_args.weights.temp_gt ) ;
		System.out.format("temp_gt %f\n", p.p_args.weights.temp_gt) ;
		System.out.format("temp_lt %f\n", p.p_args.weights.temp_lt) ;
		System.out.format("gc_content_gt %f\n", p.p_args.weights.gc_content_gt) ;
		System.out.format("gc_content_lt %f\n", p.p_args.weights.gc_content_lt) ;
		System.out.format("compl_any %f\n", p.p_args.weights.compl_any) ;
		System.out.format("compl_end %f\n", p.p_args.weights.compl_end) ;
		System.out.format("compl_any_th %f\n", p.p_args.weights.compl_any_th) ;
		System.out.format("compl_end_th %f\n", p.p_args.weights.compl_end_th) ;
		System.out.format("hairpin %f\n", p.p_args.weights.hairpin_th) ;
		System.out.format("repeat_sim %f\n", p.p_args.weights.repeat_sim) ;
		System.out.format("length_lt %f\n", p.p_args.weights.length_lt) ;
		System.out.format("length_gt %f\n", p.p_args.weights.length_gt) ;
		System.out.format("seq_quality %f\n", p.p_args.weights.seq_quality) ;
		System.out.format("end_quality %f\n", p.p_args.weights.end_quality) ;
		System.out.format("pos_penalty %f\n", p.p_args.weights.pos_penalty) ;
		System.out.format("end_stability %f\n", p.p_args.weights.end_stability) ;
		System.out.format("num_ns %f\n", p.p_args.weights.num_ns) ;
		System.out.format("template_mispriming %f\n", p.p_args.weights.template_mispriming) ;
		System.out.format("template_mispriming_th %f\n", p.p_args.weights.template_mispriming_th) ;
		System.out.format("failure_rate %f\n", p.p_args.weights.failure_rate) ;
		System.out.format("end oligo_weights\n") ;

		System.out.format("opt_tm %f\n", p.p_args.opt_tm) ;
		System.out.format("min_tm %f\n", p.p_args.min_tm) ;
		System.out.format("max_tm %f\n", p.p_args.max_tm) ;
		System.out.format("opt_gc_content %f\n", p.p_args.opt_gc_content) ;
		System.out.format("max_gc %f\n", p.p_args.max_gc) ;
		System.out.format("min_gc %f\n", p.p_args.min_gc) ;
		System.out.format("divalent_conc %f\n", p.p_args.divalent_conc) ;
		System.out.format("dntp_conc %f\n", p.p_args.dntp_conc) ;
		System.out.format("dna_conc %f\n", p.p_args.dna_conc) ;
		System.out.format("num_ns_accepted %s\n", p.p_args.num_ns_accepted) ;
		System.out.format("opt_size %s\n", p.p_args.opt_size) ;
		System.out.format("min_size %s\n", p.p_args.min_size) ;
		System.out.format("max_size %s\n", p.p_args.max_size) ;
		System.out.format("max_poly_x %s\n", p.p_args.max_poly_x) ;
		System.out.format("min_end_quality %s\n", p.p_args.min_end_quality) ;
		System.out.format("min_quality %s\n", p.p_args.min_quality) ;
		System.out.format("max_self_any %f\n", p.p_args.max_self_any) ;
		System.out.format("max_self_end %f\n", p.p_args.max_self_end) ;
		System.out.format("max_self_any_th %f\n", p.p_args.max_self_any_th) ;
		System.out.format("max_self_end_th %f\n", p.p_args.max_self_end_th) ;
		System.out.format("max_hairpin %f\n", p.p_args.max_hairpin_th) ;
		System.out.format("max_repeat_compl %f\n", p.p_args.max_repeat_compl) ;
		System.out.format("max_template_mispriming %f\n", p.p_args.max_template_mispriming) ;
		System.out.format("max_template_mispriming_th %f\n", p.p_args.max_template_mispriming_th) ;
		System.out.format("end primer args\n") ;

		System.out.format("begin internal oligo args (p.o_args.)\n") ;

		System.out.format("  begin internal oligo_weights (p.o_args.weights.)\n") ;
		System.out.format("    temp_gt %f\n", p.o_args.weights.temp_gt) ;
		System.out.format("    temp_lt %f\n", p.o_args.weights.temp_lt) ;
		System.out.format("    gc_content_gt %f\n", p.o_args.weights.gc_content_gt) ;
		System.out.format("    gc_content_lt %f\n", p.o_args.weights.gc_content_lt) ;
		System.out.format("    compl_any %f\n", p.o_args.weights.compl_any) ;
		System.out.format("    compl_end %f\n", p.o_args.weights.compl_end) ;
		System.out.format("    compl_any_th %f\n", p.o_args.weights.compl_any_th) ;
		System.out.format("    compl_end_th %f\n", p.o_args.weights.compl_end_th) ;
		System.out.format("    hairpin %f\n", p.o_args.weights.hairpin_th) ;
		System.out.format("    repeat_sim %f\n", p.o_args.weights.repeat_sim) ;
		System.out.format("    length_lt %f\n", p.o_args.weights.length_lt) ;
		System.out.format("    length_gt %f\n", p.o_args.weights.length_gt) ;
		System.out.format("    seq_quality %f\n", p.o_args.weights.seq_quality) ;
		System.out.format("    end_quality %f\n", p.o_args.weights.end_quality) ;
		System.out.format("    pos_penalty %f\n", p.o_args.weights.pos_penalty) ;
		System.out.format("    end_stability %f\n", p.o_args.weights.end_stability) ;
		System.out.format("    num_ns %f\n", p.o_args.weights.num_ns) ;
		System.out.format("  end internal oligo_weights\n") ;

		System.out.format("  opt_tm %f\n", p.o_args.opt_tm) ;
		System.out.format("  min_tm %f\n", p.o_args.min_tm) ;
		System.out.format("  max_tm %f\n", p.o_args.max_tm) ;
		System.out.format("  opt_gc_content %f\n", p.o_args.opt_gc_content) ;
		System.out.format("  max_gc %f\n", p.o_args.max_gc) ;
		System.out.format("  min_gc %f\n", p.o_args.min_gc) ;
		System.out.format("  divalent_conc %f\n", p.o_args.divalent_conc) ;
		System.out.format("  dntp_conc %f\n", p.o_args.dntp_conc) ;
		System.out.format("  dna_conc %f\n", p.o_args.dna_conc) ;
		System.out.format("  num_ns_accepted %s\n", p.o_args.num_ns_accepted) ;
		System.out.format("  opt_size %s\n", p.o_args.opt_size) ;
		System.out.format("  min_size %s\n", p.o_args.min_size) ;
		System.out.format("  max_size %s\n", p.o_args.max_size) ;
		System.out.format("  max_poly_x %s\n", p.o_args.max_poly_x) ;
		System.out.format("  min_end_quality %s\n", p.o_args.min_end_quality) ;
		System.out.format("  min_quality %s\n", p.o_args.min_quality) ;
		System.out.format("  max_self_any %f\n", p.o_args.max_self_any) ;
		System.out.format("  max_self_end %f\n", p.o_args.max_self_end) ;
		System.out.format("  max_repeat_compl %f\n", p.o_args.max_repeat_compl) ;
		System.out.format("  end internal oligo args\n");
		System.out.format("\n");
		System.out.format("END GLOBAL ARGS\n");
		System.out.format("=============\n");
		System.out.format("\n");
	}

	public boolean PR_POSITION_PENALTY_IS_NULL() {
		return (LibPrimer3.PR_DEFAULT_INSIDE_PENALTY == this.inside_penalty 
				&& LibPrimer3.PR_DEFAULT_OUTSIDE_PENALTY == this.outside_penalty);
	}
}