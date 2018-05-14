package com.primer3;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

import com.primer3.libprimer3.p3_global_settings;
import com.primer3.libprimer3.p3retval;
import com.primer3.libprimer3.seq_args;
import com.primer3.libprimer3.task;
import com.primer3.p3_seq_lib.seq_lib;

public class boulder {
	/* 
	 * Read data from file_input until a "=" line occurs.  Assign
	 * parameter values for primer picking to pa and sarg. Perform initial
	 * data checking. Return 0 if no records or no _more_ records were
	 * found, and 1 otherwise.  If nonfatal_err.data is not NULL or
	 * fatal_err.data is not NULL then the data is erroneous and should
	 * not be processed. Echo the input lines to stdout.
	 */
	static public boolean read_boulder_record(Scanner scan,
	                        boolean strict_tags,
	                        int io_version,
	                        boolean echo_output,
	                        p3_file_type file_type,
	                        p3_global_settings pa, 
	                        seq_args sa,
	                        StringBuilder glob_err,
	                        StringBuilder non_fatal_err,
	                        StringBuilder warnings,
	                        read_boulder_record_results  res){
		StringBuilder parse_err = non_fatal_err;
		boolean data_found = false;
		String task_tmp = null ;
		
		boolean min_3_prime_distance_specific = false;
		boolean min_3_prime_distance_global = false;
		String line = "";
		String repeat_file_path = null;
		String int_repeat_file_path = null;
		while(scan.hasNext())
		{
			line = scan.nextLine().trim();
			if(line.equals("="))
				break;
			
			if (file_type == p3_file_type.settings && ! line.startsWith("PRIMER_") && ! line.startsWith("P3_FILE_ID") )
				continue;
			data_found = true;
			if (echo_output) 
				System.out.format("%s\n", line);
			
			String[] lineTokens = line.split("=");
			if(lineTokens.length <= 1 )
			{
				glob_err.append("Input line with no '=': " + line);
				continue;
			}
			else
			{
				String datum = lineTokens[1];
				String key = lineTokens[0];
				try
				{

					
					if (key.equals("SEQUENCE_TEMPLATE")) {   /* NEW WAY */
						if (sa.sequence != null) {
							parse_err.append("Duplicate tag: ");
							parse_err.append("SEQUENCE_TEMPLATE"); 
						} else {
							sa.p3_set_sa_sequence(datum); 
						}
					}
					else if (key.equals("SEQUENCE_QUALITY")) {
	//					 = parse_seq_quality(datum, sa);
						if (!sa.set_n_quality(datum)) {
							parse_err.append("Error in sequence quality data");
						}
					}
					else if (key.equals("SEQUENCE_ID"))
					{
						sa.p3_set_sa_sequence_name(datum);
					} 
					else if (key.equals("SEQUENCE_PRIMER")) {
						sa.p3_set_sa_left_input(datum);
					} 
					else if (key.equals("SEQUENCE_PRIMER_REVCOMP")) {
						sa.p3_set_sa_right_input(datum);	
					}
					else if (key.equals("SEQUENCE_INTERNAL_OLIGO")) { 
						sa.p3_set_sa_internal_input(datum);
					}
					else if (key.equals("SEQUENCE_PRIMER_PAIR_OK_REGION_LIST")) { 
						sa.p3_set_sa_ok_regions(datum);
					}
					else if (key.equals("SEQUENCE_TARGET")) { 
						sa.p3_set_sa_tar2(datum);
					}
					else if (key.equals("SEQUENCE_EXCLUDED_REGION")) { 
						sa.p3_set_sa_excl2(datum);
					}
					else if (key.equals("SEQUENCE_INTERNAL_EXCLUDED_REGION")) {
						sa.p3_set_sa_excl_internal2(datum);
					}
					else if (key.equals("SEQUENCE_OVERLAP_JUNCTION_LIST")) {
						if(!sa.p3_set_sa_primer_overlap_junctions(datum))
						{
							parse_err.append("Error in SEQUENCE_PRIMER_OVERLAP_JUNCTION_LIST");
						}
	//					if (parse_intron_list(datum, sa.p3_set_sa_primer_overlap_junctions, 
	//						      &sa.p3_set_sa_primer_overlap_junctions_count) == 0) {
	//			          pr_append_new_chunk(parse_err,
	//						      "Error in SEQUENCE_PRIMER_OVERLAP_JUNCTION_LIST");
	//			        }
					} 
					else if (key.equals("SEQUENCE_INCLUDED_REGION")) {
				        
						if(!sa.p3_set_sa_incl_sl(datum))
						{
							tag_syntax_error("SEQUENCE_INCLUDED_REGION", datum,			
			                             parse_err);
						}
					} 
					else if (key.equals("SEQUENCE_START_CODON_POSITION")) { 
						sa.p3_set_sa_start_codon_pos(datum);
					}
					else if (key.equals("SEQUENCE_FORCE_LEFT_START")) { 
						sa.p3_set_sa_force_left_start(datum);
					}
					else if (key.equals("SEQUENCE_FORCE_LEFT_END")) { 
						sa.p3_set_sa_force_left_end(datum);
					}
					else if (key.equals("SEQUENCE_FORCE_RIGHT_START")) { 
						sa.p3_set_sa_force_right_start(datum);
					}
					else if (key.equals("SEQUENCE_FORCE_RIGHT_END")) { 
						sa.p3_set_sa_force_right_end(datum);
						continue;
					}
				
					/* Process "Global" Arguments (those that persist between boulder
				       * records). */
					parse_err = glob_err;  /* These errors are considered fatal. */
					
					
					if (key.equals("PRIMER_PRODUCT_SIZE_RANGE")) {
						pa.p3_set_pa_product_size(datum);
					}
					else if (key.equals("PRIMER_OPT_SIZE")) { 
						pa.p_args.setOptSize(datum);
					}
					else if (key.equals("PRIMER_MIN_SIZE")) { 
						pa.p_args.setMinSize(datum);
					}
					else if (key.equals("PRIMER_MAX_SIZE")) { 
						pa.p_args.setMaxSize(datum);
					}
					else if (key.equals("PRIMER_MAX_POLY_X")) { 
						pa.p_args.setMaxPolyX(datum);
					}
					else if (key.equals("PRIMER_OPT_TM")) { 
						pa.p_args.setOptTm(datum);
					}
					else if (key.equals("PRIMER_OPT_GC_PERCENT")) { 
						pa.p_args.setOptGCContent(datum);
					}
					else if (key.equals("PRIMER_MIN_TM")) {
						pa.p_args.setMinTm(datum);
					}
					else if (key.equals("PRIMER_MAX_TM")) { 
						pa.p_args.setMaxTm(datum);
					}
					else if (key.equals("PRIMER_PAIR_MAX_DIFF_TM")) { 
						pa.p3_set_pa_max_diff_tm(datum);
					}
					else if (key.equals("PRIMER_TM_FORMULA")) { 
						pa.p3_set_pa_tm_method_type(datum);
					}
					else if (key.equals("PRIMER_SALT_CORRECTIONS")) { 
						pa.p3_set_pa_salt_corrections(datum);
					}
					else if (key.equals("PRIMER_MIN_GC")) { 
						//pa.p_args.set_min_gc
						pa.p_args.setMinGC(datum);
					}
					else if (key.equals("PRIMER_MAX_GC")) { 
						//pa.p_args.set_max_gc
						pa.p_args.setMaxGC(datum);
					}
					else if (key.equals("PRIMER_SALT_MONOVALENT")) { 
						//pa.p_args.set_salt_conc
						pa.p_args.setSaltConcentration(datum);
					}
					else if (key.equals("PRIMER_SALT_DIVALENT")) { 
						// pa.p_args.set_divalent_conc
						pa.p_args.setDivalentConcentration(datum);
					}
					else if (key.equals("PRIMER_DNTP_CONC")) { 
						//pa.p_args.set_dntp_conc
						pa.p_args.setDntpConcentration(datum);
					}
					else if (key.equals("PRIMER_DNA_CONC")) { 
						//pa.p_args.set_dna_conc
						pa.p_args.setDnaConcentration(datum);
					}
					else if (key.equals("PRIMER_MAX_NS_ACCEPTED")) { 
						// pa.p_args.set_num_ns_accepted
						pa.p_args.setMaxNumOfNsAccepted(datum);
					}
					else if (key.equals("PRIMER_PRODUCT_OPT_SIZE")) { 
						//pa.p3_set_pa_product_opt_size
						pa.p3_set_pa_product_opt_size(datum);
					}
					else if (key.equals("PRIMER_MAX_SELF_ANY")) { 
						// pa.p_args.set_max_self_any
						pa.p_args.setMaxSelfAny(datum);
					}
					else if (key.equals("PRIMER_MAX_SELF_END")) { 
						// p3_set_gs_primer_self_end
						pa.p_args.setMaxSelfEnd(datum);
					}
					else if (key.equals("PRIMER_MAX_SELF_ANY_TH")) { 
						//pa.p_args.set_max_self_any_th
						pa.p_args.setMaxSelfAnyTH(datum);
					}
					else if (key.equals("PRIMER_MAX_SELF_END_TH")) { 
						//pa.p_args.set_max_self_end_th
						pa.p_args.setMaxSelfEndTH(datum);
					}
					else if (key.equals("PRIMER_MAX_HAIRPIN_TH")) { 
						//pa.p_args.set_max_hairpin_th
						pa.p_args.setMaxHairPinTH(datum);
					}
					else if (key.equals("PRIMER_PAIR_MAX_COMPL_ANY")) { 
						//pa.p3_set_pa_pair_compl_any
						pa.p3_set_pa_pair_compl_any(datum);
					}
					else if (key.equals("PRIMER_PAIR_MAX_COMPL_END")) { 
						//pa.p3_set_pa_pair_compl_end
						pa.p3_set_pa_pair_compl_end(datum);
					}
					else if (key.equals("PRIMER_PAIR_MAX_COMPL_ANY_TH")) { 
						//pa.p3_set_pa_pair_compl_any_th
						pa.p3_set_pa_pair_compl_any_th(datum);
					}
					else if (key.equals("PRIMER_PAIR_MAX_COMPL_END_TH")) { 
						//pa.p3_set_pa_pair_compl_end_th
						pa.p3_set_pa_pair_compl_end_th(datum);
					}
					else if (key.equals("PRIMER_PICK_ANYWAY")) { 
						//pa.p3_set_pa_pick_anyway
						pa.p3_set_pa_pick_anyway(datum);
					}
					else if (key.equals("PRIMER_GC_CLAMP")) { 
						//pa.p3_set_pa_gc_clamp
						pa.p3_set_pa_gc_clamp(datum);
					}
					else if (key.equals("PRIMER_MAX_END_GC")) { 
						//pa.p3_set_pa_max_end_gc
						pa.p3_set_pa_max_end_gc(datum);
					}
					else if (key.equals("P3_FILE_FLAG")) { 
						res.file_flag = Integer.parseInt(datum);
					}
					else if (key.equals("PRIMER_EXPLAIN_FLAG")) { 
						res.explain_flag = Integer.parseInt(datum);
					}
					else if (key.equals("PRIMER_LIBERAL_BASE")) { 
						//pa.p3_set_pa_liberal_base
						pa.p3_set_pa_liberal_base(datum);
					}
					else if (key.equals("PRIMER_FIRST_BASE_INDEX")) { 
						//pa.p3_set_pa_first_base_index
						pa.p3_set_pa_first_base_index(datum);
					}
					else if (key.equals("PRIMER_NUM_RETURN")) { 
						// pa.p3_set_pa_num_return
						pa.p3_set_pa_num_return(datum);
					}
					else if (key.equals("PRIMER_MIN_QUALITY")) { 
						// pa.p_args.set_min_quality
						pa.p_args.setMinQuality(datum);
					}
					else if (key.equals("PRIMER_MIN_END_QUALITY")) {
						
	//					COMPARE_INT("PRIMER_MIN_END_QUALITY", pa.p_args.set_min_end_quality);
						pa.p_args.setMinEndQuality(datum);
					}
					else if (key.equals("PRIMER_MIN_THREE_PRIME_DISTANCE")) { 
						int min_three_prime_distance = Integer.parseInt(datum);
						if (min_3_prime_distance_specific == true) {
							glob_err.append("Both PRIMER_MIN_THREE_PRIME_DISTANCE and PRIMER_{LEFT/RIGHT}_MIN_THREE_PRIME_DISTANCE specified");
						} else {
							min_3_prime_distance_global = true;
							/* set up individual flags */
							pa.set_min_left_three_prime_distance(min_three_prime_distance);
							pa.set_min_right_three_prime_distance(min_three_prime_distance);
						}
					}
					else if (key.equals("PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE")) { 
						if (min_3_prime_distance_global) {
							  glob_err.append("Both PRIMER_MIN_THREE_PRIME_DISTANCE and PRIMER_{LEFT/RIGHT}_MIN_THREE_PRIME_DISTANCE specified");
						} else {
							pa.p3_set_pa_min_left_three_prime_distance(datum);
							min_3_prime_distance_specific = true;
						}
					}
					else if (key.equals("PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE")) { 
						if (min_3_prime_distance_global) {
							  glob_err.append("Both PRIMER_MIN_THREE_PRIME_DISTANCE and PRIMER_{LEFT/RIGHT}_MIN_THREE_PRIME_DISTANCE specified");
						} else {
							pa.p3_set_pa_min_right_three_prime_distance(datum);
							min_3_prime_distance_specific = true;
						}
					}
					else if (key.equals("PRIMER_QUALITY_RANGE_MIN")) { 
						//pa.p3_set_pa_quality_range_min
						pa.p3_set_pa_quality_range_min(datum);
					}
					else if (key.equals("PRIMER_QUALITY_RANGE_MAX")) { 
						//pa.p3_set_pa_quality_range_max
						pa.p3_set_pa_quality_range_max(datum);
					}
					else if (key.equals("PRIMER_PRODUCT_MAX_TM")) { 
						//pa.p3_set_pa_product_max_tm
						pa.p3_set_pa_product_max_tm(datum);
					}
					else if (key.equals("PRIMER_PRODUCT_MIN_TM")) { 
						//pa.p3_set_pa_product_min_tm
						pa.p3_set_pa_product_min_tm(datum);
					}
					else if (key.equals("PRIMER_PRODUCT_OPT_TM")) { 
						//pa.p3_set_pa_product_opt_tm
						pa.p3_set_pa_product_opt_tm(datum);
					}
					else if (key.equals("PRIMER_SEQUENCING_LEAD")) { 
						//pa.p3_set_pa_sequencing.lead
						pa.sequencing.set_lead(datum);
					}
					else if (key.equals("PRIMER_SEQUENCING_SPACING")) { 
						//pa.p3_set_pa_sequencing.spacing
						pa.sequencing.set_spacing(datum);
					}
					else if (key.equals("PRIMER_SEQUENCING_INTERVAL")) { 
						//pa.p3_set_pa_sequencing.interval
						pa.sequencing.set_interval(datum);
					}
					else if (key.equals("PRIMER_SEQUENCING_ACCURACY")) { 
						//pa.p3_set_pa_sequencing.accuracy
						pa.sequencing.set_accuracy(datum);
					}
					else if (key.equals("PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION")) { 
						//pa.p3_set_pa_min_5_prime_overlap_of_junction
						pa.p3_set_pa_min_5_prime_overlap_of_junction(datum);
					}
					else if (key.equals("PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION")) { 
						//pa.p3_set_pa_min_3_prime_overlap_of_junction
						pa.p3_set_pa_min_3_prime_overlap_of_junction(datum);
					}
					else if (key.equals("PRIMER_TASK")) { 
						task_tmp = datum ;
//						pa.p3_set_pa_(datum);
					}
					
					else if (key.equals("PRIMER_PICK_RIGHT_PRIMER")) { 
						pa.p3_set_pa_pick_right_primer(datum);
					}
					else if (key.equals("PRIMER_PICK_INTERNAL_OLIGO")) { 
						pa.p3_set_pa_pick_internal_oligo(datum);
					}
					else if (key.equals("PRIMER_PICK_LEFT_PRIMER")) { 
						pa.p3_set_pa_pick_left_primer(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_OPT_SIZE")) {
						pa.o_args.setOptSize(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_MAX_SIZE")) { 
						pa.o_args.setMaxSize(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_MIN_SIZE")) { 
						pa.o_args.setMinSize(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_MAX_POLY_X")) { 
						pa.o_args.setMaxPolyX(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_OPT_TM")) { 
						pa.o_args.setOptTm(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_OPT_GC_PERCENT")) {
						pa.o_args.setOptGCContent(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_MAX_TM")) { 
						pa.o_args.setMaxTm(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_MIN_TM")) { 
						pa.o_args.setMinTm(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_MIN_GC")) { 
						pa.o_args.setMinGC(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_MAX_GC")) { 
						pa.o_args.setMaxGC(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_SALT_MONOVALENT")) {
						pa.o_args.setSaltConcentration(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_SALT_DIVALENT")) {
						pa.o_args.setDivalentConcentration(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_DNTP_CONC")) {
						pa.o_args.setDntpConcentration(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_DNA_CONC")) { 
						pa.o_args.setDnaConcentration(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_MAX_NS_ACCEPTED")) { 
						pa.o_args.setMaxNumOfNsAccepted(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_MIN_QUALITY")) { 
						pa.o_args.setMinQuality(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_MAX_SELF_ANY")) {
						pa.o_args.setMaxSelfAny(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_MAX_SELF_END")) { 
						pa.p3_set_gs_primer_internal_oligo_self_end(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_MAX_SELF_ANY_TH")) {
						pa.o_args.setMaxSelfAnyTH(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_MAX_SELF_END_TH")) {
						pa.o_args.setMaxSelfEndTH(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_MAX_HAIRPIN_TH")) {
						pa.o_args.setMaxHairPinTH(datum);
					}
					else if (key.equals("PRIMER_MAX_LIBRARY_MISPRIMING")) {
						pa.p_args.setMaxRepeatCompl(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_MAX_LIBRARY_MISHYB")) {
						pa.o_args.setMaxRepeatCompl(datum);
					}
					else if (key.equals("PRIMER_PAIR_MAX_LIBRARY_MISPRIMING")) {
						pa.p3_set_pa_pair_repeat_compl(datum);
					}
				      /* Mispriming / mishybing in the template. */
					else if (key.equals("PRIMER_MAX_TEMPLATE_MISPRIMING")) {
						pa.p_args.setMaxTemplateMispriming(datum);
					}
					else if (key.equals("PRIMER_MAX_TEMPLATE_MISPRIMING_TH")) {
						pa.p_args.setMaxTemplateMisprimingTH(datum);
					}
					else if (key.equals("PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING")) {
						pa.p3_set_pa_pair_max_template_mispriming(datum);
					}
					else if (key.equals("PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH")) {
						pa.p3_set_pa_pair_max_template_mispriming_th(datum);
					}
				       /* Control interpretation of ambiguity codes in mispriming
				          and mishyb libraries. */
					else if (key.equals("PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS")) {
						pa.p3_set_pa_lib_ambiguity_codes_consensus(datum);
					}
					else if (key.equals("PRIMER_INSIDE_PENALTY")) {
						pa.p3_set_pa_inside_penalty(datum);
					}
					else if (key.equals("PRIMER_OUTSIDE_PENALTY")) {
						pa.p3_set_pa_outside_penalty(datum);
					}
					else if (key.equals("PRIMER_MISPRIMING_LIBRARY")) {
						repeat_file_path = datum;
					}
					else if (key.equals("PRIMER_INTERNAL_MISHYB_LIBRARY")) {
						int_repeat_file_path = datum;
						//				        TODO :: 
//						if (int_repeat_file_path != NULL) {
//				          pr_append_new_chunk(glob_err,
//				                              "Duplicate PRIMER_INTERNAL_MISHYB_LIBRARY tag");
//				          free(int_repeat_file_path);
//				          int_repeat_file_path = NULL;
//				        } else {
//				          int_repeat_file_path = (char*) _rb_safe_malloc(strlen(datum) + 1);
//				          strcpy(int_repeat_file_path, datum);
//				        }
					}
					else if (key.equals("P3_COMMENT")) {

					}
					else if (key.equals("PRIMER_MAX_END_STABILITY")) {
						pa.p3_set_pa_max_end_stability(datum);
					}
					else if (key.equals("PRIMER_LOWERCASE_MASKING")) {
						pa.p3_set_pa_lowercase_masking(datum);
					}
				      /* added by T. Koressaar */
					else if (key.equals("PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT")) { 
						pa.p3_set_pa_thermodynamic_oligo_alignment(datum);
					}
					else if (key.equals("PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT")) {
						pa.p3_set_pa_thermodynamic_template_alignment(datum);
					}
					else if (key.equals("PRIMER_THERMODYNAMIC_PARAMETERS_PATH")) {
//				        if (thermodynamic_params_path == NULL) {
//				          thermodynamic_params_path = (char*) _rb_safe_malloc(datum_len + 1);
//				          strcpy(thermodynamic_params_path, datum);
//				          thermodynamic_path_changed = 1;
//				        }
//				        /* check if path changes */
//					else if (strcmp(thermodynamic_params_path, datum)) {
//					  free(thermodynamic_params_path);
//					  thermodynamic_params_path = (char*) _rb_safe_malloc(datum_len + 1); 
//					  strcpy(thermodynamic_params_path, datum);
//					  thermodynamic_path_changed = 1;
//				        }
//				        continue;

					}
					else if (key.equals("PRIMER_MASK_KMERLIST_PATH")) {
//				           if (kmer_lists_path == NULL) {
//				               kmer_lists_path = (char*) _rb_safe_malloc(datum_len + 1);
//				               strcpy(kmer_lists_path, datum);
//				           }
					}
					/* added by A. Untergasser */
					else if (key.equals("PRIMER_MUST_MATCH_FIVE_PRIME")) { 
						pa.p_args.set_must_match_five_prime(datum);
					}
					else if (key.equals("PRIMER_MUST_MATCH_THREE_PRIME")) { 
						pa.p_args.set_must_match_three_prime(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME")) { 
						pa.o_args.set_must_match_five_prime(datum);
					}
					else if (key.equals("PRIMER_INTERNAL_MUST_MATCH_THREE_PRIME")) { 
						pa.o_args.set_must_match_three_prime(datum);
					}
				      /* weights for objective functions  */
				      /* CHANGE TEMP/temp . TM/tm */
					else if (key.equals("PRIMER_WT_TM_GT")) { 
						pa.p_args.weights.set_temp_gt(datum);
					}
					else if (key.equals("PRIMER_WT_TM_LT")) { 
						pa.p_args.weights.set_temp_lt(datum);
					}
				    else if (key.equals("PRIMER_WT_GC_PERCENT_GT")) { 
				    	pa.p_args.weights.set_gc_content_gt(datum);
				    }
				    else if (key.equals("PRIMER_WT_GC_PERCENT_LT")) { 
				    	pa.p_args.weights.set_gc_content_lt(datum);
				    }
				    else if (key.equals("PRIMER_WT_SIZE_LT")) { 
				    	pa.p_args.weights.set_length_lt(datum);
				    }
				    else if (key.equals("PRIMER_WT_SIZE_GT")) { 
				    	pa.p_args.weights.set_length_gt(datum);
				    }
				    else if (key.equals("PRIMER_WT_SELF_ANY")) { 
				    	pa.p_args.weights.set_compl_any(datum);
				    }
				    else if (key.equals("PRIMER_WT_SELF_END")) { 
				    	pa.p_args.weights.set_compl_end(datum);
				    }
				    else if (key.equals("PRIMER_WT_SELF_ANY_TH")) { 
				    	pa.p_args.weights.set_compl_any_th(datum);
				    }
				    else if (key.equals("PRIMER_WT_SELF_END_TH")) { 
				    	pa.p_args.weights.set_compl_end_th(datum);
				    }
				    else if (key.equals("PRIMER_WT_HAIRPIN_TH")) { 
				    	pa.p_args.weights.set_hairpin_th(datum);
				    }
				    else if (key.equals("PRIMER_WT_NUM_NS")) { 
				    	pa.p_args.weights.set_num_ns(datum);
				    }
				    else if (key.equals("PRIMER_WT_LIBRARY_MISPRIMING")) { 
				    	pa.p_args.weights.set_repeat_sim(datum);
				    }
				    else if (key.equals("PRIMER_WT_SEQ_QUAL")) { 
				    	pa.p_args.weights.set_seq_quality(datum);
				    }
				    else if (key.equals("PRIMER_WT_END_QUAL")) { 
				    	pa.p_args.weights.set_end_quality(datum);
				    }
				    else if (key.equals("PRIMER_WT_POS_PENALTY")) { 
				    	pa.p_args.weights.set_pos_penalty(datum);
				    }
				    else if (key.equals("PRIMER_WT_END_STABILITY")) { 
				    	pa.p_args.weights.set_end_stability(datum);
				    }
				    else if (key.equals("PRIMER_WT_TEMPLATE_MISPRIMING")) { 
				    	pa.p_args.weights.set_template_mispriming(datum);
				    }
				    else if (key.equals("PRIMER_WT_TEMPLATE_MISPRIMING_TH")) { 
				    	pa.p_args.weights.set_template_mispriming_th(datum);
				    }
				    else if (key.equals("PRIMER_WT_MASK_FAILURE_RATE")) { 
				    	pa.p_args.weights.set_failure_rate(datum);                        
				    }
				    else if (key.equals("PRIMER_INTERNAL_WT_TM_GT")) { 
				    	pa.o_args.weights.set_temp_gt(datum);
				    }
				    else if (key.equals("PRIMER_INTERNAL_WT_TM_LT")) { 
				    	pa.o_args.weights.set_temp_lt(datum);
				    }
				    else if (key.equals("PRIMER_INTERNAL_WT_GC_PERCENT_GT")) { 
				    	pa.o_args.weights.set_gc_content_gt(datum);
				    }
				    else if (key.equals("PRIMER_INTERNAL_WT_GC_PERCENT_LT")) { 
				    	pa.o_args.weights.set_gc_content_lt(datum);
				    }
				    else if (key.equals("PRIMER_INTERNAL_WT_SIZE_LT")) { 
				    	pa.o_args.weights.set_length_lt(datum);
				    }
				    else if (key.equals("PRIMER_INTERNAL_WT_SIZE_GT")) { 
				    	pa.o_args.weights.set_length_gt(datum);
				    }
				    else if (key.equals("PRIMER_INTERNAL_WT_SELF_ANY")) { 
				    	pa.o_args.weights.set_compl_any(datum);
				    }
				    else if (key.equals("PRIMER_INTERNAL_WT_SELF_END")) { 
				    	pa.o_args.weights.set_compl_end(datum);
				    }
				    else if (key.equals("PRIMER_INTERNAL_WT_SELF_ANY_TH")) { 
				    	pa.o_args.weights.set_compl_any_th(datum);
				    }
				    else if (key.equals("PRIMER_INTERNAL_WT_SELF_END_TH")) { 
				    	pa.o_args.weights.set_compl_end_th(datum);
				    }
				    else if (key.equals("PRIMER_INTERNAL_WT_HAIRPIN_TH")) { 
				    	pa.o_args.weights.set_hairpin_th(datum);
				    }
				    else if (key.equals("PRIMER_INTERNAL_WT_NUM_NS")) { 
				    	pa.o_args.weights.set_num_ns(datum);
				    }
				    else if (key.equals("PRIMER_INTERNAL_WT_LIBRARY_MISHYB")) { 
				    	pa.o_args.weights.set_repeat_sim(datum);
				    }
				    else if (key.equals("PRIMER_INTERNAL_WT_SEQ_QUAL")) { 
				    	pa.o_args.weights.set_seq_quality(datum);
				    }
				    else if (key.equals("PRIMER_INTERNAL_WT_END_QUAL")) { 
				    	pa.o_args.weights.set_end_quality(datum);
				    }
				    else if (key.equals("PRIMER_WT_TEMPLATE_MISPRIMING_TH")) { 
				    	pa.o_args.weights.set_template_mispriming_th(datum);
				    }
				    else if (key.equals("PRIMER_PAIR_WT_PR_PENALTY")) {
				    	pa.pr_pair_weights.set_primer_quality(datum);
				    }
				    else if (key.equals("PRIMER_PAIR_WT_IO_PENALTY")) { 
				    	pa.pr_pair_weights.set_io_quality(datum);
				    }
				    else if (key.equals("PRIMER_PAIR_WT_DIFF_TM")) {
				    	pa.pr_pair_weights.set_diff_tm(datum);
				    }
				    else if (key.equals("PRIMER_PAIR_WT_COMPL_ANY")) { 
				    	pa.pr_pair_weights.set_compl_any(datum);
				    }
				    else if (key.equals("PRIMER_PAIR_WT_COMPL_END")) { 
				    	pa.pr_pair_weights.set_compl_end(datum);
				    }
				    else if (key.equals("PRIMER_PAIR_WT_COMPL_ANY_TH")) { 
				    	pa.pr_pair_weights.set_compl_any_th(datum);
				    }
				    else if (key.equals("PRIMER_PAIR_WT_COMPL_END_TH")) { 
				    	pa.pr_pair_weights.set_compl_end_th(datum);
				    }
				    else if (key.equals("PRIMER_PAIR_WT_PRODUCT_TM_LT")) { 
				    	pa.pr_pair_weights.set_product_tm_lt(datum);
				    }
				    else if (key.equals("PRIMER_PAIR_WT_PRODUCT_TM_GT")) { 
				    	pa.pr_pair_weights.set_product_tm_gt(datum);
				    }
				    else if (key.equals("PRIMER_PAIR_WT_PRODUCT_SIZE_GT")) { 
				    	pa.pr_pair_weights.set_product_size_gt(datum);
				    }
				    else if (key.equals("PRIMER_PAIR_WT_PRODUCT_SIZE_LT")) {  
				    	pa.pr_pair_weights.set_product_size_lt(datum);
				    }
				    else if (key.equals("PRIMER_PAIR_WT_LIBRARY_MISPRIMING")) { 
				    	pa.pr_pair_weights.set_repeat_sim(datum);
				    }
				    else if (key.equals("PRIMER_PAIR_WT_TEMPLATE_MISPRIMING")) { 
				    	pa.pr_pair_weights.set_template_mispriming(datum);
				    }
				    else if (key.equals("PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH")) { 
				    	pa.pr_pair_weights.set_template_mispriming_th(datum);
				    }
				    else if (key.equals("PRIMER_MASK_TEMPLATE")) {

				    	boolean retValue = pa.p3_set_pa_mask_template(datum);
				    	pa.p3_set_pa_lowercase_masking(retValue);
				    	pa.masking_parameters_changed = retValue;
//	                    parse_int("PRIMER_MASK_TEMPLATE", datum, &pa.mask_template, parse_err);  
//	                    pa.lowercase_masking = pa.mask_template;
//	                    pa.masking_parameters_changed = pa.mask_template;
//	                    continue;
					}
				    else if (key.equals("PRIMER_MASK_FAILURE_RATE")) {
				    	pa.mp.set_failure_rate(datum);
				    }
				    else if (key.equals("PRIMER_MASK_5P_DIRECTION")) {
				    	pa.mp.set_nucl_masked_in_5p_direction(datum);
				    }
				    else if (key.equals("PRIMER_MASK_3P_DIRECTION")) {
				    	pa.mp.set_nucl_masked_in_3p_direction(datum);
				    }
				    else if (key.equals("PRIMER_MASK_KMERLIST_PREFIX")) {
				    	pa.set_masking_parameters_KmerLis_Prefix(datum);
				    }
				    else
				    {
				    	if(strict_tags) {
							glob_err.append("Unrecognized tag: " + line);
							System.err.println("Unrecognized tag: " + line);
							
						}
				    }
					
					
				} // try
				catch(Exception ex)
				{
					
				}
				// catch  here
				
			} // main else
			
			
		} // while
		
		
		
	
		
		/* Figure out the right settings for the tasks*/
		if (task_tmp != null ) {
			if (task_tmp.equals("pick_pcr_primers")) {
				pa.primer_task = task.generic;
				pa.pick_left_primer = true;
				pa.pick_right_primer = true;
				pa.pick_internal_oligo = false;
			} 
			else if (task_tmp.equals("pick_pcr_primers_and_hyb_probe")) {
				pa.primer_task = task.generic; 
				pa.pick_left_primer = true;
				pa.pick_right_primer = true;
				pa.pick_internal_oligo = true;
		    } 
			else if (task_tmp.equals( "pick_left_only")) {
				pa.primer_task = task.generic;
				pa.pick_left_primer = true;
				pa.pick_right_primer = false;
				pa.pick_internal_oligo = false;
		    } 
			else if (task_tmp.equals( "pick_right_only")) {
				pa.primer_task = task.generic;
				pa.pick_left_primer = false;
				pa.pick_right_primer = true;
				pa.pick_internal_oligo = false;
		    } 
			else if (task_tmp.equals( "pick_hyb_probe_only")) {
				pa.primer_task = task.generic;
				pa.pick_left_primer = false;
				pa.pick_right_primer = false;
				pa.pick_internal_oligo = true;
		    } 
			else if (task_tmp.equals( "generic")) {
				pa.primer_task = task.generic;
		    } 
			else if (task_tmp.equals( "pick_detection_primers")) {
				pa.primer_task = task.generic; /* Deliberate duplication for
						    backward compatibility. */
		    } 
			else if (task_tmp.equals( "pick_cloning_primers")) {
				pa.primer_task = task.pick_cloning_primers;
		    } 
			else if (task_tmp.equals( "pick_discriminative_primers")) {
				pa.primer_task = task.pick_discriminative_primers;
		    } 
			else if (task_tmp.equals( "pick_sequencing_primers")) {
		    	pa.primer_task = task.pick_sequencing_primers;
		    } 
			else if (task_tmp.equals( "pick_primer_list")) {
				pa.primer_task = task.pick_primer_list;
		    } else if (task_tmp.equals( "check_primers")) {
		    	pa.primer_task = task.check_primers;
		    	/* check_primers sets the picking flags itself */
		    	pa.pick_left_primer = false;
		    	pa.pick_right_primer = false;
		    	pa.pick_internal_oligo = false;
		    	if (sa.left_input != null){
		    		pa.pick_left_primer = true;
		    	}
		    	if (sa.right_input != null){
		    		pa.pick_right_primer = true;
		    	}
		    	if (sa.internal_input != null){
		    		pa.pick_internal_oligo = true;
		    	}
		    } 
		    else 
		    	glob_err.append("Unrecognized PRIMER_TASK");
		}
		
		
		/* Read in the repeat libraries */
		if (null != repeat_file_path) {
			pa.p_args.repeat_lib =  null;
			if (!repeat_file_path.isEmpty()) {
				try {
					pa.p_args.repeat_lib = seq_lib.read_and_create_seq_lib(repeat_file_path, 
				                              "mispriming library");
					if(!pa.p_args.repeat_lib.getErrors().isEmpty()) {
						glob_err.append(pa.p_args.repeat_lib.getErrors());
					}
				} catch (Exception e) {
					glob_err.append("Can not read mispriming library " + repeat_file_path);
				}
			}
			repeat_file_path = null;
		}
		
		  /* Read in the repeat libraries for internal oligo */
		if (null != int_repeat_file_path) {
			pa.o_args.repeat_lib =  null;
			if (!int_repeat_file_path.isEmpty()) {
				try {
					pa.o_args.repeat_lib = seq_lib.read_and_create_seq_lib(int_repeat_file_path, 
				                              "internal oligo mishyb library");
					if(!pa.o_args.repeat_lib.getErrors().isEmpty()) {
						glob_err.append(pa.o_args.repeat_lib.getErrors());
					}
				} catch (Exception e) {
					glob_err.append("Can not read internal oligo mishyb library" + int_repeat_file_path);
				}
			}
			int_repeat_file_path = null;
		}
		
		if(!scan.hasNext())
		{
			if(data_found)
			{
				glob_err.append("Final record not terminated by '='");
				return true;
			}
			return false;
		}
		
		return true;
		
	}

	

	// FIXME :: no need for this
	private static void tag_syntax_error(String string, String datum,
			StringBuilder parse_err) {
		// TODO Auto-generated method stub
		
	}

	/* Return null on error. */
	/* pr_append_str is an append-only string ADT. */
	static public boolean read_p3_file(
					 String file_name,
	                 p3_file_type expected_file_type,
	                 boolean echo_output,
	                 boolean strict_tags,
	                 p3_global_settings pa, 
	                 seq_args sarg,
	                 StringBuilder fatal_err,
	                 StringBuilder nonfatal_err,
	                 StringBuilder warnings,
	                 read_boulder_record_results res){
		
		String line1 = null;
		String line2 = null;
		String line3 = null;
		Scanner scan = null;
		int io_version = 4;
		p3_file_type file_type = p3_file_type.all_parameters;
		try {
			scan = new Scanner(new File(file_name));
			if(scan.hasNext())
				line1 = scan.nextLine();
			if(line1 == null)
			{
				fatal_err.append("Settings file is empty: " + file_name);
				return false;
			}
			
			
			
			
//			if ((strcmp(line1,"Primer3 File - http://primer3.org") != 0) &&
//				      (strcmp(line1,"Primer3 File - http://primer3.sourceforge.net") != 0)) {
//				      pr_append2(fatal_err,
//				                 "First line must be \"Primer3 File - http://primer3.org\" in ",
//						 file_name);
//				      return ret_par;
//				  }
			
			
			
			
			
			
			
			if(scan.hasNext())
				line2 = scan.nextLine();
			if(line2 == null)
			{
				fatal_err.append("Incorrect file format (too few lines) in " + file_name);
				return false;
			}
			
			
			if (line2.equals("P3_FILE_TYPE=all_parameters")) {
			    file_type = p3_file_type.all_parameters;
			  } else if (line2.equals("P3_FILE_TYPE=sequence")) {
			    file_type = p3_file_type.sequence;
			  } else if (line2.equals("P3_FILE_TYPE=settings")) {
			    file_type = p3_file_type.settings;
			  } else {
				  fatal_err.append( "Unknown file type in at line 2 (" + line2 + ") in " + file_name);
				  return false;
			  }
			  if (echo_output) {
			    System.out.format("P3_SETTINGS_FILE_USED=%s\n", file_name);
			    System.out.format("%s\n", line2);
			  }
			
			
			
			
			
			if(scan.hasNext())
				line3 = scan.nextLine();
			if(line3 == null)
			{
				fatal_err.append("Incorrect file format (too few lines) in " + file_name);
				return false;
			}
			
			if (!line3.isEmpty()) {
				  fatal_err.append("Line 3 must be empty in "+file_name);
				  return false;
			}
			
			
			
			
			/* read the file */
			
			 /* Check if the file type matches the expected type */
			if (file_type != expected_file_type){
			   nonfatal_err.append( "Unexpected P3 file type parsed");
			}
			
			boolean ret_par = read_boulder_record(scan, strict_tags, io_version, 
							echo_output, expected_file_type,
							pa, sarg, fatal_err, 
							nonfatal_err, warnings, 
							res);
			
			if (echo_output) 
				System.out.println("P3_SETTINGS_FILE_END=");

			
			
			return ret_par;
		} catch (FileNotFoundException e) {
			fatal_err.append("Cannot open "+ file_name);
		}
		finally{
			if(scan !=  null )
				scan.close();
		}
		
		return false;
		
	}

	
	


	public static void format_warning(String sequence_name,
			String warnings) {
		if(sequence_name != null && !sequence_name.isEmpty())
			System.out.format( "WARNINGS FOR %s\n\n", sequence_name);
		if(warnings != null && !warnings.isEmpty())
			System.out.format( "INPUT PROBLEM: %s\n\n", warnings);
	}

	

	public static void format_error(String sequence_name, String err) {
		if(sequence_name != null && !sequence_name.isEmpty())
			System.out.format( "PRIMER PICKING RESULTS FOR %s\n\n", sequence_name);
		if(err != null && !err.isEmpty())
			System.out.format( "INPUT PROBLEM: %s\n\n", err);		
	}

	public static void print_boulder_error(String err) {
		  System.out.format("PRIMER_ERROR=%s\n=\n", err);
	}
	public static void print_boulder_warning(String err) {
		  System.out.format("PRIMER_WARNING=%s\n=\n", err);
	}



	public static void print_format_output(int io_version,
			p3_global_settings global_pa, seq_args sarg, p3retval retval,
			String pr_release, int explain_flag) {
		// TODO Auto-generated method stub
		
	}
	
	
	
}
