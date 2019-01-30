package org.primer3.libprimer3;

import java.util.List;

import org.primer3.boulder;
import org.primer3.p3_seq_lib.seq_lib;
import org.primer3.search.P3BasicPairFinder;
import org.primer3.search.P3OptimzedFinder;
import org.primer3.search.Primer3Finder;
import org.primer3.sequence.Sequence;

/**
 * The return value for for primer3. 
 * After use, free memory with destroy_p3retval().
 */
public class P3RetVal {

	/* Arrays of oligo (primer) records. */
	// left primers 
	public OligoArray fwd;

	// internal/ probes
	public OligoArray intl;

	// right primers
	public OligoArray rev;

	/* Array of best primer pairs */
	public PairArrayT best_pairs ;

	/* Enum to store type of output */
	public P3OutputType output_type;

	
	// seq args associated with this result
	public SeqArgs sa;
	// global setting associated with this result 
	// This is goping to change after this result is collected
	// TODO :: side effect here refactor this to keep it's own copy of the args ??
	public P3GlobalSettings pa;
	
	
	/* Place for error messages */
	// Originally was pr_append_str -. changed to StringBuilder
	public StringBuilder glob_err = new StringBuilder();

	public StringBuilder per_sequence_err = new StringBuilder();
	public StringBuilder warnings = new StringBuilder();

	/* 
	 * An optional _output_, meaninful if a
	 * start_codon_pos is "not null".  The position of
	 * the intial base of the leftmost stop codon that
	 * is to the right of sa.start_codon_pos.
	 */
	int stop_codon_pos;

	int upstream_stop_codon;  /* TO DO needs docs */

	// this changed to constructor
	//  static public p3retval create_p3retval(){
	//	  
	//  }

	protected P3RetVal ()
	{
		this.fwd = new OligoArray(OligoType.OT_LEFT);
		this.intl = new OligoArray( OligoType.OT_INTL);
		this.rev = new OligoArray(OligoType.OT_RIGHT);

		this.fwd.type  = OligoType.OT_LEFT;
		this.intl.type = OligoType.OT_INTL;
		this.rev.type  = OligoType.OT_RIGHT;
		
		best_pairs = new PairArrayT() ;
	}

	public P3RetVal(P3GlobalSettings pa, SeqArgs sa) {
		this();
		/* Set the general output type */
		if (pa.isPickLeftPrimer() && pa.isPickRightPrimer()) {
			this.output_type = P3OutputType.primer_pairs;
		} else {
			this.output_type = P3OutputType.primer_list;
		}
		if (	pa.getPrimerTask() == P3Task.PICK_PRIMER_LIST ||
				pa.getPrimerTask() == P3Task.PICK_SEQUENCING_PRIMERS) {
			this.output_type = P3OutputType.primer_list;
		}
		this.pa = pa;
		this.sa = sa;
	}

	public OligoArray p3_get_rv_rev()
	{
		return rev;
	}

	public OligoArray p3_get_rv_fwd()
	{
		return fwd;
	}
	public OligoArray p3_get_rv_intl()
	{
		return intl;
	}

	public PairArrayT p3_get_rv_best_pairs()
	{
		return best_pairs;
	}


	public String get_glob_err() {
		return glob_err.toString();
	}
	public String get_per_sequence_err() {
		return per_sequence_err.toString();
	}
	public String get_warnings() {
		return warnings.toString();
	}


	public void print_boulder(int io_version, boolean   explain_flag) {

		P3RetVal retval = this;
		/* The pointers to warning tag */
		String warning;

		/* A place to put a string containing all error messages */
		StringBuilder combined_retval_err = null;

		/* A small spacer; WARNING this is a fixed size
	     buffer, but plenty bigger than
	     log(2^64, 10), the longest character
	     string that is needed for a 64 bit integer. */
		String suffix ;

		/* Pointers for the primer set just printing */
		PrimerRecord fwd, rev, intl;

		/* Variables only used for Primer Lists */
		int num_fwd, num_rev, num_int, num_pair, num_print;
		int print_fwd = 0;
		int print_rev = 0;
		int print_int = 0;

		/* Switches for printing this primer */
		int go_fwd = 0;
		int go_rev = 0;
		int go_int = 0;

		/* The number of loop cycles */
		int loop_max;

		/* That links to the included region */
		int i, incl_s = sa.getIncludedRegionStart();

		/* This deals with the renaming of the internal oligo */
		String new_oligo_name = "INTERNAL";
		String int_oligo = new_oligo_name;

		/* Check: are all pointers linked to something*/
		//	  PR_ASSERT(NULL != pa);
		//	  PR_ASSERT(NULL != sa);

		/* Check if there are warnings and print them */
		warning = retval.p3_get_rv_and_gs_warnings(pa);
		if (!warning.isEmpty()) { 
			System.out.format("PRIMER_WARNING=%s\n", warning);
		}

		/* Check if a settings file was read an print its id 
		 * Not needed anymore, we do this when we read the setting file. */
		/*
	  if (pa.settings_file_id != NULL) { 
	    System.out.format("P3_FILE_ID=%s\n", pa.settings_file_id);
	    free(warning);
	  }
		 */

		combined_retval_err = new StringBuilder();
		//	  if (NULL == combined_retval_err) exit(-2); /* Out of memory */

		combined_retval_err.append( retval.glob_err.toString());
		combined_retval_err.append( retval.per_sequence_err.toString());



		/* Check if there are errors, print and return */
		if (!combined_retval_err.toString().isEmpty()) {
			boulder.print_boulder_error(combined_retval_err.toString());
			combined_retval_err = null;
			return;
		}
		combined_retval_err = null;
		/* Get how many primers are in the array */
		num_fwd = retval.fwd.num_elem;
		num_rev = retval.rev.num_elem;
		num_int = retval.intl.num_elem;
		
		num_pair = retval.best_pairs.num_pairs;
		//num_pair = pa.getNumReturn();
		/* Prints out statistics about the primers */
		if (explain_flag) retval.print_all_explain(pa, sa, io_version);

		/* Print out the stop codon if a reading frame was specified */
		if (!sa.PR_START_CODON_POS_IS_NULL())
			System.out.format("PRIMER_STOP_CODON_POSITION=%d\n", retval.stop_codon_pos);

		/* How often has the loop to be done? */
		if (retval.output_type == P3OutputType.primer_list) {
			/* For Primer Lists: Figure out how many primers are in
			 * the array that can be printed. If more than needed,
			 * set it to the number requested. */

			/* Get how may primers should be printed */
			num_print = pa.getNumReturn();
			/* Set how many primers will be printed */
			print_fwd = (num_print < num_fwd) ? num_print : num_fwd;
			print_rev = (num_print < num_rev) ? num_print : num_rev;
			print_int = (num_print < num_int) ? num_print : num_int;
			/* Get which list has to print most primers */
			loop_max = 0;
			if (loop_max < print_fwd) {
				loop_max = print_fwd;
			}
			if (loop_max < print_rev) {
				loop_max = print_rev;
			}
			if (loop_max < print_int) {
				loop_max = print_int;
			}
			/* Now the vars are there how often we have to go
			 * through the loop and how many of each primer can
			 * be printed. */
			num_pair = 0;
		} else {
			loop_max = num_pair;
			/* Set how many primers will be printed */
			print_fwd = num_pair;
			print_rev = num_pair;
			if (num_int != 0) {
				print_int = num_pair;
			}
		}

		if (io_version == 4) {
			System.out.format("PRIMER_LEFT_NUM_RETURNED=%d\n", print_fwd);
			System.out.format("PRIMER_RIGHT_NUM_RETURNED=%d\n",  print_rev);
			System.out.format("PRIMER_%s_NUM_RETURNED=%d\n", int_oligo, print_int);
			System.out.format("PRIMER_PAIR_NUM_RETURNED=%d\n", num_pair);
		}

		/* --------------------------------------- */
		/* Start of the loop printing all pairs or primers or oligos */
		for(i=0; i<loop_max; i++) {
			/* What needs to be printed */
			/* The conditions for primer lists */

			if (retval.output_type == P3OutputType.primer_list) {
				
				// bug access empty list
				fwd  = null;
				rev  = null;
				intl = null;
				
				
				/* Attach the selected primers to the pointers */
				//fwd = retval.fwd.oligo.get(i);
				//rev = retval.rev.oligo.get(i);
				//intl = retval.intl.oligo.get(i);
				/* Do fwd oligos have to be printed? */
				if ((pa.isPickLeftPrimer()) && (i < print_fwd)) {
					go_fwd = 1;
					fwd = retval.fwd.oligo.get(i);
				} else {
					go_fwd = 0;
				}
				/* Do rev oligos have to be printed? */
				if ((pa.isPickRightPrimer()) && (i < print_rev)) {
					go_rev = 1;
					rev = retval.rev.oligo.get(i);
				} else {
					go_rev = 0;
				}
				/* Do int oligos have to be printed? */
				if ((pa.isPickInternalOligo()) && (i < print_int)) {
					go_int = 1;
					intl = retval.intl.oligo.get(i);
				} else {
					go_int = 0;
				}
			}  else {
				/* We will print primer pairs or pairs plus internal oligos */
				/* Get pointers to the primer_rec's that we will print */
				fwd  = retval.best_pairs.pairs.get(i).left;
				rev  = retval.best_pairs.pairs.get(i).right;
				
				// potential null exp bug here 
				// intl = retval.best_pairs.pairs.get(i).intl;
				intl = null;
				/* Pairs must have fwd and rev primers */
				go_fwd = 1;
				go_rev = 1;
				/* Do hyb oligos have to be printed? */
				if (pa.isPickInternalOligo() ) {
					go_int = 1;
					intl = retval.best_pairs.pairs.get(i).intl;
				} else {
					go_int = 0;
				}
			}

			/* Get the number for pimer counting in suffix[0] */
			suffix = "_"+ i;

			/* Print out the Pair Penalties */
			if (retval.output_type == P3OutputType.primer_pairs) {
				System.out.format("PRIMER_PAIR%s_PENALTY=%f\n", suffix,
						retval.best_pairs.pairs.get(i).pair_quality);
			}

			/* Print single primer penalty */
			if (go_fwd == 1)
				System.out.format("PRIMER_LEFT%s_PENALTY=%f\n", suffix, fwd.quality);
			if (go_rev == 1)
				System.out.format("PRIMER_RIGHT%s_PENALTY=%f\n", suffix, rev.quality);
			if (go_int == 1)
				System.out.format("PRIMER_%s%s_PENALTY=%f\n", int_oligo, suffix, intl.quality);

			/* Print the oligo_problems */
			if (io_version == 4) {
				if (go_fwd == 1 && fwd.p3_ol_has_any_problem())
					System.out.format("PRIMER_LEFT%s_PROBLEMS=%s\n", suffix, fwd.p3_get_ol_problem_string());
				if (go_rev == 1 && rev.p3_ol_has_any_problem())
					System.out.format("PRIMER_RIGHT%s_PROBLEMS=%s\n", suffix, rev.p3_get_ol_problem_string());
				if (go_int == 1 && intl.p3_ol_has_any_problem())
					System.out.format("PRIMER_%s%s_PROBLEMS=%s\n", int_oligo, suffix, intl.p3_get_ol_problem_string());
			}

			/* Print primer sequences. */
			if (go_fwd == 1)
				System.out.format("PRIMER_LEFT%s_SEQUENCE=%s\n", suffix,
						//						pr_oligo_sequence(sa, fwd));
						LibPrimer3.	string(fwd.pr_oligo_sequence(sa)));
			if (go_rev == 1)
				System.out.format("PRIMER_RIGHT%s_SEQUENCE=%s\n", suffix,
						LibPrimer3.	string(rev.pr_oligo_rev_c_sequence(sa )));
			if(go_int == 1)
				System.out.format("PRIMER_%s%s_SEQUENCE=%s\n", int_oligo, suffix,
					LibPrimer3.	string(intl.pr_oligo_sequence(sa)));

			/* Print primer start and length */
			if (go_fwd == 1)
				System.out.format("PRIMER_LEFT%s=%d,%d\n", suffix,
						fwd.start + incl_s + pa.getFirstBaseIndex(),
						fwd.length);
			if (go_rev == 1)
				System.out.format("PRIMER_RIGHT%s=%d,%d\n", suffix,
						rev.start + incl_s + pa.getFirstBaseIndex(),
						rev.length);
			if (go_int == 1)
				System.out.format("PRIMER_%s%s=%d,%d\n", int_oligo, suffix,
						intl.start + incl_s + pa.getFirstBaseIndex(),
						intl.length);

			/* Print primer Tm */
			if (go_fwd == 1)
				System.out.format("PRIMER_LEFT%s_TM=%.3f\n", suffix, fwd.temp);
			if (go_rev == 1)
				System.out.format("PRIMER_RIGHT%s_TM=%.3f\n", suffix, rev.temp);
			if (go_int == 1)
				System.out.format("PRIMER_%s%s_TM=%.3f\n", int_oligo, suffix, intl.temp);

			/* Print primer GC content */
			if (go_fwd == 1)
				System.out.format("PRIMER_LEFT%s_GC_PERCENT=%.3f\n", suffix, fwd.gc_content);
			if (go_rev == 1)
				System.out.format("PRIMER_RIGHT%s_GC_PERCENT=%.3f\n", suffix, rev.gc_content);
			if (go_int == 1)
				System.out.format("PRIMER_%s%s_GC_PERCENT=%.3f\n", int_oligo, suffix,
						intl.gc_content);

			/* Print primer self_any */
			if (go_fwd == 1 && pa.isThermodynamicOligoAlignment()==false)
				System.out.format("PRIMER_LEFT%s_SELF_ANY=%.2f\n", suffix,
						fwd.self_any);
			if (go_rev == 1 && pa.isThermodynamicOligoAlignment()==false)
				System.out.format("PRIMER_RIGHT%s_SELF_ANY=%.2f\n", suffix,
						rev.self_any);
			if (go_int == 1 && pa.isThermodynamicOligoAlignment()==false)
				System.out.format("PRIMER_%s%s_SELF_ANY=%.2f\n", int_oligo, suffix,
						intl.self_any);
			if (go_int == 1 && pa.isThermodynamicOligoAlignment()==true)
				System.out.format("PRIMER_%s%s_SELF_ANY_TH=%.2f\n", int_oligo, suffix,
						intl.self_any);
			/* Print primer self_any thermodynamical approach */
			if (go_fwd == 1 && pa.isThermodynamicOligoAlignment()==true)
				System.out.format("PRIMER_LEFT%s_SELF_ANY_TH=%.2f\n", suffix,
						fwd.self_any);
			if (go_rev == 1 && pa.isThermodynamicOligoAlignment()==true)
				System.out.format("PRIMER_RIGHT%s_SELF_ANY_TH=%.2f\n", suffix,
						rev.self_any);
			/* Print primer self_end*/
			if (go_fwd == 1 && pa.isThermodynamicOligoAlignment()==false)
				System.out.format("PRIMER_LEFT%s_SELF_END=%.2f\n", suffix,
						fwd.self_end);
			if (go_rev == 1 && pa.isThermodynamicOligoAlignment()==false)
				System.out.format("PRIMER_RIGHT%s_SELF_END=%.2f\n", suffix,
						rev.self_end);
			if (go_int == 1 && pa.isThermodynamicOligoAlignment()==false)
				System.out.format("PRIMER_%s%s_SELF_END=%.2f\n", int_oligo, suffix,
						intl.self_end);
			if (go_int == 1 && pa.isThermodynamicOligoAlignment()==true)
				System.out.format("PRIMER_%s%s_SELF_END_TH=%.2f\n", int_oligo, suffix,
						intl.self_end);
			/* Print primer self_end thermodynamical approach */
			if (go_fwd == 1 && pa.isThermodynamicOligoAlignment()==true)
				System.out.format("PRIMER_LEFT%s_SELF_END_TH=%.2f\n", suffix,
						fwd.self_end);
			if (go_rev == 1 && pa.isThermodynamicOligoAlignment()==true)
				System.out.format("PRIMER_RIGHT%s_SELF_END_TH=%.2f\n", suffix,
						rev.self_end);
			/* Print primer hairpin */
			if (go_fwd == 1 && pa.isThermodynamicOligoAlignment()==true)
				System.out.format("PRIMER_LEFT%s_HAIRPIN_TH=%.2f\n", suffix,
						fwd.hairpin_th);
			if (go_rev == 1 && pa.isThermodynamicOligoAlignment()==true)
				System.out.format("PRIMER_RIGHT%s_HAIRPIN_TH=%.2f\n", suffix,
						rev.hairpin_th);
			if (go_int == 1 && pa.isThermodynamicOligoAlignment()==true)
				System.out.format("PRIMER_%s%s_HAIRPIN_TH=%.2f\n", int_oligo, suffix,
						intl.hairpin_th);
			/*Print out primer mispriming scores */
			if (pa.primersArgs.repeat_lib != null) {
				if (go_fwd == 1)
					System.out.format("PRIMER_LEFT%s_LIBRARY_MISPRIMING=%.2f, %s\n", suffix,
							fwd.repeat_sim.score.get(fwd.repeat_sim.max),
							fwd.repeat_sim.name);
				if (go_rev == 1)
					System.out.format("PRIMER_RIGHT%s_LIBRARY_MISPRIMING=%.2f, %s\n", suffix,
							rev.repeat_sim.score.get(rev.repeat_sim.max),
							rev.repeat_sim.name);
				if (retval.output_type == P3OutputType.primer_pairs)
					System.out.format("PRIMER_PAIR%s_LIBRARY_MISPRIMING=%.2f, %s\n", suffix,
							retval.best_pairs.pairs.get(i).repeat_sim,
							retval.best_pairs.pairs.get(i).rep_name);
			}

			/* Print out internal oligo mispriming scores */
			if (go_int == 1 && pa.oligosArgs.repeat_lib != null)
				System.out.format("PRIMER_%s%s_LIBRARY_MISHYB=%.2f, %s\n", int_oligo, suffix,
						intl.repeat_sim.score.get(intl.repeat_sim.max),
						intl.repeat_sim.name);

			/* If a sequence quality was provided, print it*/
			if (null != sa.getSequenceQuality()){
				if (go_fwd == 1)
					System.out.format("PRIMER_LEFT%s_MIN_SEQ_QUALITY=%d\n", suffix,
							fwd.seq_quality);
				if (go_rev == 1) 
					System.out.format("PRIMER_RIGHT%s_MIN_SEQ_QUALITY=%d\n", suffix,
							rev.seq_quality);
				if (go_int == 1 && (retval.output_type == P3OutputType.primer_list)) 
					System.out.format("PRIMER_%s%s_MIN_SEQ_QUALITY=%d\n", int_oligo, suffix,
							intl.seq_quality);
				/* Has to be here and in primer pairs for backward compatibility */
			}

			/* Print position penalty, this is for backward compatibility */
			if (!pa.isDefaultPositionPenalties()
					|| !sa.PR_START_CODON_POS_IS_NULL()){
				System.out.format("PRIMER_LEFT%s_POSITION_PENALTY=%f\n", suffix,
						fwd.position_penalty);
				System.out.format("PRIMER_RIGHT%s_POSITION_PENALTY=%f\n", suffix,
						rev.position_penalty);
			}

			/* Print primer end stability */
			if (go_fwd == 1)
				System.out.format("PRIMER_LEFT%s_END_STABILITY=%.4f\n",
						suffix, fwd.end_stability);
			if (go_rev == 1)
				System.out.format("PRIMER_RIGHT%s_END_STABILITY=%.4f\n",
						suffix, rev.end_stability);

			/* Print primer template mispriming */
			if ( (!pa.isThermodynamicTemplateAlignment() ) && (go_fwd == 1) && 
					(fwd.oligo_max_template_mispriming() != LibPrimer3.ALIGN_SCORE_UNDEF))
				System.out.format("PRIMER_LEFT%s_TEMPLATE_MISPRIMING=%.4f\n", suffix,
						fwd.oligo_max_template_mispriming());
			if ( (!pa.isThermodynamicTemplateAlignment() ) && (go_rev == 1) && 
					(rev.oligo_max_template_mispriming() != LibPrimer3.ALIGN_SCORE_UNDEF))
				System.out.format("PRIMER_RIGHT%s_TEMPLATE_MISPRIMING=%.4f\n", suffix,
						rev.oligo_max_template_mispriming());

			/* Print primer template mispriming, thermodynamical approach*/
			if ( (pa.isThermodynamicTemplateAlignment()) && (go_fwd == 1) &&
					(fwd.oligo_max_template_mispriming_thermod() != LibPrimer3.ALIGN_SCORE_UNDEF)) {
				System.out.format("PRIMER_LEFT%s_TEMPLATE_MISPRIMING_TH=%.4f\n", suffix,
						fwd.oligo_max_template_mispriming_thermod());
			}

			if ( (pa.isThermodynamicTemplateAlignment() ) && (go_rev == 1) &&
					(rev.oligo_max_template_mispriming_thermod() != LibPrimer3.ALIGN_SCORE_UNDEF)) {
				System.out.format("PRIMER_RIGHT%s_TEMPLATE_MISPRIMING_TH=%.4f\n", suffix,
						rev.oligo_max_template_mispriming_thermod());
				//	#if 0
				//	       System.out.format("DEBUG_PRIMER_RIGHT%s_TEMPLATE_MISPRIMING_TOP_TH=%.4f\n", suffix,
				//		      rev.template_mispriming);
				//	       System.out.format("DEBUG_PRIMER_RIGHT%s_TEMPLATE_MISPRIMING_R_TH=%.4f\n", suffix,
				//		      rev.template_mispriming_r);
				//	#endif
			}
			/************************************************************************************/
			/* Print the pair parameters*/
			if (retval.output_type == P3OutputType.primer_pairs) {
				if (go_int == 1 && null != sa.getSequenceQuality()) /* FIX ME - Uptate the tests */
					System.out.format("PRIMER_%s%s_MIN_SEQ_QUALITY=%d\n", int_oligo,
							suffix, intl.seq_quality);
				/* Print pair comp_any */
				if(pa.isThermodynamicOligoAlignment()==false)
					System.out.format("PRIMER_PAIR%s_COMPL_ANY=%.2f\n", suffix,
							retval.best_pairs.pairs.get(i).compl_any);
				if(pa.isThermodynamicOligoAlignment()==true)
					System.out.format("PRIMER_PAIR%s_COMPL_ANY_TH=%.2f\n", suffix,
							retval.best_pairs.pairs.get(i).compl_any);
				/* Print pair comp_end */
				if(pa.isThermodynamicOligoAlignment()==false)
					System.out.format("PRIMER_PAIR%s_COMPL_END=%.2f\n", suffix,
							retval.best_pairs.pairs.get(i).compl_end);
				if(pa.isThermodynamicOligoAlignment()==true)
					System.out.format("PRIMER_PAIR%s_COMPL_END_TH=%.2f\n", suffix,
							retval.best_pairs.pairs.get(i).compl_end);
				/* Print product size */
				System.out.format("PRIMER_PAIR%s_PRODUCT_SIZE=%d\n", suffix,
						retval.best_pairs.pairs.get(i).product_size);
				/* Print the product Tm if a Tm range is defined */
				if (pa.getProductMaxTM() != LibPrimer3.PR_DEFAULT_PRODUCT_MAX_TM ||
						pa.getProductMinTM() != LibPrimer3.PR_DEFAULT_PRODUCT_MIN_TM) {
					System.out.format("PRIMER_PAIR%s_PRODUCT_TM=%.4f\n", suffix,
							retval.best_pairs.pairs.get(i).product_tm);

					System.out.format("PRIMER_PAIR%s_PRODUCT_TM_OLIGO_TM_DIFF=%.4f\n", suffix,
							retval.best_pairs.pairs.get(i).product_tm_oligo_tm_diff);

					System.out.format("PRIMER_PAIR%s_T_OPT_A=%.4f\n", suffix,
							retval.best_pairs.pairs.get(i).t_opt_a);
				}

				/* Print the primer pair template mispriming */
				if ((!pa.isThermodynamicTemplateAlignment() ) && (retval.best_pairs.pairs.get(i).template_mispriming != LibPrimer3.ALIGN_SCORE_UNDEF))
					System.out.format("PRIMER_PAIR%s_TEMPLATE_MISPRIMING=%.2f\n", suffix,
							retval.best_pairs.pairs.get(i).template_mispriming);
				/* Print the primer pair template mispriming. Thermodynamic approach.  */
				if ((pa.isThermodynamicTemplateAlignment() ) && (retval.best_pairs.pairs.get(i).template_mispriming != LibPrimer3.ALIGN_SCORE_UNDEF))
					System.out.format("PRIMER_PAIR%s_TEMPLATE_MISPRIMING_TH=%.2f\n", suffix,
							retval.best_pairs.pairs.get(i).template_mispriming);

			} /* End of print parameters of primer pairs */

		} /* End of the big loop printing all data */

		/* End the print with newline and flush all buffers */
		System.out.format("=\n");

	}

	private String p3_get_rv_and_gs_warnings(P3GlobalSettings pa) {

		String warning = "";

		if(pa.primersArgs.repeat_lib != null)
			warning += pa.primersArgs.repeat_lib.seq_lib_warning_data();

		if(pa.oligosArgs.repeat_lib != null && !pa.oligosArgs.repeat_lib.seq_lib_warning_data().isEmpty()) {
			warning += pa.oligosArgs.repeat_lib.seq_lib_warning_data() + " (for internal oligo)";
		}

		if(this.warnings.length() > 0)
			warning += this.warnings.toString();

		return warning;
	}

	private void print_all_explain(P3GlobalSettings pa, SeqArgs sa,
			int io_version) {
		if (pa.isPickLeftPrimer()
				&& !(pa.isPickAnyway() && sa.getLeftInput() != null))
			System.out.format("PRIMER_LEFT_EXPLAIN=%s\n",
					p3_get_rv_fwd().p3_get_oligo_array_explain_string());

		if (pa.isPickRightPrimer() 
				&& !(pa.isPickAnyway() && sa.getRightInput() != null))
			System.out.format("PRIMER_RIGHT_EXPLAIN=%s\n",
					p3_get_rv_rev().p3_get_oligo_array_explain_string());

		if ( pa.isPickInternalOligo() 
				&& !(pa.isPickAnyway() && sa.getInternalInput() != null)) 
			System.out.format("PRIMER_INTERNAL_EXPLAIN=%s\n",
					p3_get_rv_intl().p3_get_oligo_array_explain_string());

		if (pa.isPickRightPrimer()  
				&& pa.isPickLeftPrimer() ) {
			System.out.format("PRIMER_PAIR_EXPLAIN=%s\n", 
					p3_get_rv_best_pairs().p3_get_pair_array_explain_string());
		}
	}



	


	public void add_must_use_warnings(String text, OligoArray oarray) {


		OligoStats stats = oarray.expl;

		String sep = "/";
		StringBuilder s =  new StringBuilder();



		if (stats.size_min != 0) 
			s.append( sep + "Too short");
		if (stats.size_max != 0) 
			s.append(sep+ "Too long");
		if (stats.ns != 0) s.append(sep+ "Too many Ns");
		if (stats.target != 0) s.append(sep+ "Overlaps Target");
		if (stats.excluded != 0) s.append(sep+ "Overlaps Excluded Region");
		if (stats.gc != 0) s.append(sep+ "Unacceptable GC content");
		if (stats.gc_clamp != 0) s.append(sep+ "No GC clamp");
		if (stats.temp_min != 0) s.append(sep+ "Tm too low");
		if (stats.temp_max != 0) s.append(sep+ "Tm too high");
		if (stats.compl_any != 0) s.append(sep+ "High self complementarity");
		if (stats.compl_end != 0)
			s.append(sep+ "High end self complementarity");
		if (stats.hairpin_th != 0) s.append(sep+ "High hairpin stability (thermod. approach)");
		if (stats.repeat_score != 0)
			s.append(sep+ "High similarity to mispriming or mishyb library");
		if (stats.poly_x != 0) s.append(sep+ "Long poly-X");
		if (stats.seq_quality != 0) s.append(sep+ "Low sequence quality");
		if (stats.stability != 0) s.append(sep+ "High 3' stability");
		if (stats.no_orf != 0) s.append(sep+ "Would not amplify any ORF");
		if (stats.not_in_any_left_ok_region != 0) s.append(sep+ "Not in any ok left region");
		if (stats.not_in_any_right_ok_region != 0) s.append(sep+ "Not in any ok right region");

		/* edited by T. Koressaar for lowercase masking: */
		if (stats.gmasked != 0)
			s.append(sep+ "Masked with lowercase letter");

		if (stats.must_match_fail != 0)
			s.append(sep+ "Failed must_match requirements");

		if (s.length() > 0) {
			warnings.append(text + " is unacceptable: " + s.toString());
		}
	}

	/**
	 *  The position of the intial base of the rightmost stop codon that
	 * is to the left of sa.start_codon_pos; valid only if
	 * sa.start_codon_pos is "not null".  We will not want to include
	 * a stop codon to the right of the the start codon in the
	 * amplicon. 
	 */
	public void set_retval_both_stop_codons() {
		
		  this.upstream_stop_codon = Sequence.find_stop_codon(sa.getTrimmedSequence(),
		                                                sa.getStartCodonPos(), -1);
		  this.upstream_stop_codon += sa.getIncludedRegionStart();
		  this.stop_codon_pos = Sequence.find_stop_codon(sa.getTrimmedSequence(),
		                                             sa.getStartCodonPos(),  1);
		  this.stop_codon_pos += sa.getIncludedRegionStart();		
	}
	public Primer3Finder p3Finder;
	public void choose_pairs(DPAlArgHolder dpal_arg_to_use, THAlArgHolder thal_arg_to_use,
			THAlArgHolder thal_oligo_arg_to_use) throws Exception {
		
		p3Finder = new P3BasicPairFinder(this,dpal_arg_to_use, thal_arg_to_use, thal_oligo_arg_to_use);
//		Primer3Finder p3Finder = new P3OptimzedFinder(this,dpal_arg_to_use, thal_arg_to_use, thal_oligo_arg_to_use);
		p3Finder.getNextResult();
//		p3Finder.getNextResult();
	}

	
	public void inisSearch(DPAlArgHolder dpal_arg_to_use, THAlArgHolder thal_arg_to_use,
			THAlArgHolder thal_oligo_arg_to_use)
	{
		p3Finder = new P3OptimzedFinder(this,dpal_arg_to_use, thal_arg_to_use, thal_oligo_arg_to_use);
	}
	
	public PairArrayT best_pairs_All ;
	public List<PrimerPair> getNextRusult() throws Exception
	{
//		// cache current set
//		best_pairs.cacheCurrent();
//		p3Finder.getNextResult();
//		List<PrimerPair> newPairs = best_pairs.pairs;
//		best_pairs.mergeBests();
//		return newPairs;
		return p3Finder.getNextResult();
	}


}