package org.primer3.libprimer3;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.primer3.dpal.AlignmentException;
import org.primer3.dpal.DPAlignmentArgs;
import org.primer3.dpal.DPAlignmentResults;
import org.primer3.dpal.DPAlignment;
import org.primer3.masker.input_sequence;
import org.primer3.masker.masker;
import org.primer3.masker.masking_direction;
import org.primer3.masker.output_sequence;
import org.primer3.oligotm.OligoTMCalculator;
import org.primer3.p3_seq_lib.seq_lib;
import org.primer3.sequence.Sequence;
import org.primer3.thal.ThermodynamicAlignmentException;
import org.primer3.thal.ThermodynamicAlignment;
import org.primer3.thal.ThermodynamicAlignmentArguments;
import org.primer3.thal.ThermodynamicAlignmentResult;
import org.primer3.thal.ThAl;


/**
 * ReImplemntation to libprimer3 
 * Part of primer3 software 
 * @author ahafez
 *
 */
public class LibPrimer3 {

	
	
	// Debug info 
	static int thal_trace = 0;
	static boolean choose_pair_or_triple_trace_me = false;

	

	/* ALIGN_SCORE_UNDEF is used only libprimer3 and clients, not in dpal */
	// #define ALIGN_SCORE_UNDEF            -DBL_MAX
	static final public double   ALIGN_SCORE_UNDEF = -Double.MAX_VALUE;

	/* These next 5 are exposed for format_output.c -- probabaly should be reviewed. */
	//	#define PR_INFINITE_POSITION_PENALTY -1.0
	static final public double PR_INFINITE_POSITION_PENALTY = -1.0;
	//	#define PR_DEFAULT_INSIDE_PENALTY     PR_INFINITE_POSITION_PENALTY
	static final public double PR_DEFAULT_INSIDE_PENALTY = PR_INFINITE_POSITION_PENALTY;
	//	#define PR_DEFAULT_OUTSIDE_PENALTY    0.0
	static final public double PR_DEFAULT_OUTSIDE_PENALTY = 0.0;
	//	#define PR_DEFAULT_PRODUCT_MAX_TM     1000000.0
	static final public double PR_DEFAULT_PRODUCT_MAX_TM = 1000000.0;
	//	#define PR_DEFAULT_PRODUCT_MIN_TM     -1000000.0
	static final public double PR_DEFAULT_PRODUCT_MIN_TM =  -1000000.0;

	static final public int PR_NULL_FORCE_POSITION  = -1000000;

	/*  Exposed in the read_boulder input routine.... */
	static final public int PR_NULL_START_CODON_POS     =  -1000000;
	static final public int PR_DEFAULT_START_CODON_POS  =  -2000000 ;
	//   static final public int PR_START_CODON_POS_IS_NULL(SA) ((SA).start_codon_pos <= PR_NULL_START_CODON_POS)
	//	
	//	#define _PR_DEFAULT_POSITION_PENALTIES(PA) \
	//	    (PR_DEFAULT_INSIDE_PENALTY == pa.inside_penalty \
	//	     && PR_DEFAULT_OUTSIDE_PENALTY == pa.outside_penalty)
	//	
	static final public double PR_ALIGN_SCORE_PRECISION  = 100.0;
	//	
	//	#define MACRO_STRING(X) #X

	/* pr_progam_name must be set in main(). */
	//	#define PR_ASSERT(COND)                                  \
	//	if (!(COND)) {                                           \
	//	    fprintf(stderr, "%s:%s:%d, assertion (%s) failed\n", \
	//	           pr_program_name, __FILE__, __LINE__,          \
	//	           MACRO_STRING(COND));                          \
	//	    abort();                                             \
	//	}





	/* Maxima needed for interface data structures. */
	// #define PR_MAX_INTERVAL_ARRAY 200
	static public int PR_MAX_INTERVAL_ARRAY = 200;
	/* 
	 * Maximum number of input intervals
	 * supported; used for targets, excluded
	 * regions, product-size intervals, etc.
	 */

	// this might be not used
	static int[][] interval_array_t = new int[PR_MAX_INTERVAL_ARRAY][2];; // typedef int interval_array_t[PR_MAX_INTERVAL_ARRAY][2];



	static DPAlArgHolder dpal_arg_to_use =  null;
	static THAlArgHolder thal_arg_to_use =  null;
	static THAlArgHolder thal_oligo_arg_to_use = null;



	/*
	 * OPTIMIZE_OK_REGIONS 1 allows _optimize_ok_regions_list() to use the
	 * max/min product size info and the max/min oligo length to reduce
	 * the sizes of the ok regions (while still generating the same primer
	 * pairs in the same order).  Set OPTIMIZE_OK_REGIONS to 0 confirm
	 * that the results do not change.  (The output tags
	 * PRIMER_{LEFT,RIGHT,PAIR}_EXPLAIN _will_ likely change.)
	 */
	static public boolean OPTIMIZE_OK_REGIONS = true;

	//	#ifndef MAX_PRIMER_LENGTH
	static public int MAX_PRIMER_LENGTH = 36;
	//	#endif
	//	#if (MAX_PRIMER_LENGTH > DPAL_MAX_ALIGN)
	//	#error "MAX_PRIMER_LENGTH must be <= DPAL_MAX_ALIGN"
	//	#endif
	//	#if (MAX_PRIMER_LENGTH > THAL_MAX_ALIGN)
	//	# error "MAX_PRIMER_LENGTH must be <= THAL_MAX_ALIGN"
	//	#endif
	public static int MAX_NN_TM_LENGTH = 36; /* The maxium length for which to use the
	                               nearest neighbor model when calculating
	                               oligo Tms. */

	//	#define MACRO_CAT_2(A,B) A##B
	//	#define MACRO_VALUE_AS_STRING(A) MACRO_STRING(A)

	//	#define PR_POSITION_PENALTY_IS_NULL(PA) \
	//	(PR_DEFAULT_INSIDE_PENALTY == (PA).inside_penalty \
	//	 && PR_DEFAULT_OUTSIDE_PENALTY == (PA).outside_penalty)

	static public int INITIAL_LIST_LEN  =   2000 ; /* Initial size of oligo lists. */
	static public int INITIAL_NUM_RETURN   = 5;    /* Initial space to allocate for pairs to
	                                     return. */



	//	#define OK_OR_MUST_USE(H) (!p3_ol_has_any_problem(H) || (H).must_use)

	static public int  PR_UNDEFINED_INT_OPT    =      Integer.MIN_VALUE ;
	static public double  PR_UNDEFINED_DBL_OPT    =      Double.MIN_VALUE ;

	/* Undefined value for alignment score (meaning do not check) used for
	   maximum template mispriming or mishyb. */
	static public double  PR_UNDEFINED_ALIGN_OPT   =     -100.0;

	//	#define TRIMMED_SEQ_LEN(X) ((X).incl_l)


	static public int use_end_for_th_template_mispriming = 1;

	// TODO :: check this value
	static public double DEFAULT_OPT_GC_PERCENT = 55;// PR_UNDEFINED_INT_OPT;





	/** 
	 * Choose individual primers or oligos, or primer pairs, or primer
	 * pairs with internal oligos. On ENOMEM return NULL and set errno. 
	 * Otherwise return retval (updated).  Errors are returned in 
	 * in retval.
	 */
	public static P3RetVal choose_primers(P3GlobalSettings pa,  SeqArgs sa){

		/* Create retval and set were to find the results */
		P3RetVal retval = new P3RetVal();

		// 		TODO :: check input parameters
		//		PR_ASSERT(NULL != pa);
		//		PR_ASSERT(NULL != sa);
		if (pa.isDump()) {
			System.out.println("Start of choose_primers:\n");
			p3_print_args(pa, sa) ;
		}

		/* Set the general output type */
		if (pa.isPickLeftPrimer() && pa.isPickRightPrimer()) {
			retval.output_type = P3OutputType.primer_pairs;
		} else {
			retval.output_type = P3OutputType.primer_list;
		}
		if (	pa.getPrimerTask() == P3Task.PICK_PRIMER_LIST ||
				pa.getPrimerTask() == P3Task.PICK_SEQUENCING_PRIMERS) {
			retval.output_type = P3OutputType.primer_list;
		}





		/* 
		 * TODO :: Error checking will be converted to Exception handling
		 * For catching ENOMEM.  WARNING: We can only use longjmp to escape
		 * from errors that have been called through choose_primers().
		 * Therefore, if we subsequently update other static functions in
		 * this file to have external linkage then we need to check whether
		 * they call (or use functions that in turn call) longjmp to handle
		 * ENOMEM.
		 */
		//		 if (setjmp(_jmp_buf) != 0) {
		//		    /* Check if this was a thermodynamic alignment length error. */
		//		    if (thermodynamic_alignment_length_error == 1) {
		//		      thermodynamic_alignment_length_error = 0;
		//		      /* Set per sequence error */
		//		      pr_append_new_chunk(&retval.per_sequence_err, 
		//					  thermodynamic_alignment_length_error_msg);
		//		      free(thermodynamic_alignment_length_error_msg);
		//		      thermodynamic_alignment_length_error_msg = NULL;
		//		      /* Other necessary cleanup. */
		//		      free_pair_memory(retval.rev.num_elem);
		//		      return retval;
		//		    }
		//		    /* This was a memory error. */
		//		    destroy_p3retval(retval);
		//		    return NULL;  /* If we get here, that means errno should be ENOMEM. */
		//		  }
		// List of possible Exceptions :
		// thermodynamic_alignment_length_error
		//
		try
		{
			/* Change some parameters to fit the task */
			sa._adjust_seq_args(pa,  retval);
			if (pa.isDump()) {
				System.out.println("After _adjust_seq_args\n");
				p3_print_args(pa, sa) ;
			}


			if(!retval.get_per_sequence_err().isEmpty()) 
				return retval;

			/* Check if the input in sa and pa makes sense */
			if (_pr_data_control(pa, sa,retval.glob_err,retval.per_sequence_err,retval.warnings)) {
				return retval;
			}



			retval.set_retval_both_stop_codons(sa);
			/* Set the parameters for alignment functions
		     if dpal_arg_to_use, a static variable that has 'file'
		     scope, has not yet been initialized. */
			if (dpal_arg_to_use == null)
				dpal_arg_to_use = DPAlArgHolder.create_dpal_arg_holder();// create_dpal_arg_holder();

			// TODO :: Refactor and clean
			if(thal_arg_to_use == null) {
				thal_arg_to_use = THAlArgHolder.create_thal_arg_holder(pa.primersArgs);// create_thal_arg_holder(&pa.p_args);
			} else {
				// destroy_thal_arg_holder(thal_arg_to_use);
				thal_arg_to_use = THAlArgHolder.create_thal_arg_holder(pa.primersArgs);// create_thal_arg_holder(&pa.p_args);
			}
			if(thal_oligo_arg_to_use == null) {
				thal_oligo_arg_to_use = THAlArgHolder.create_thal_arg_holder(pa.oligosArgs);
			} else {
				// destroy_thal_arg_holder(thal_oligo_arg_to_use);
				thal_oligo_arg_to_use = THAlArgHolder.create_thal_arg_holder(pa.oligosArgs);
			} 
			if (pa.getPrimerTask() == P3Task.PICK_PRIMER_LIST) {
				make_complete_primer_lists(retval, pa, sa,
						dpal_arg_to_use,thal_arg_to_use,thal_oligo_arg_to_use);
			} else if (pa.getPrimerTask() == P3Task.PICK_SEQUENCING_PRIMERS) {
				pick_sequencing_primer_list(retval, pa, sa, dpal_arg_to_use,thal_arg_to_use);
			} else if (pa.getPrimerTask() == P3Task.CHECK_PRIMERS) {
				add_primers_to_check(retval, pa, sa, dpal_arg_to_use, thal_arg_to_use, thal_oligo_arg_to_use);
			} else { 
				/* The general way to pick primers */
				/* Populate the forward and reverse primer lists */
				if (make_detection_primer_lists(retval, pa, sa, dpal_arg_to_use,thal_arg_to_use) != 0) {
					/* There was an error */
					return retval;
				}
				/* Populate the internal oligo lists */
				if ( pa.isPickInternalOligo()) {
					if (make_internal_oligo_list(retval, pa, sa, dpal_arg_to_use,thal_oligo_arg_to_use) != 0) {
						/* There was an error*/
						return retval;
					}
				}
			}


			if (pa.isPickRightPrimer() && ( pa.getPrimerTask() != P3Task.PICK_SEQUENCING_PRIMERS))
				retval.rev.sort_primer_array();
			if (pa.isPickLeftPrimer() && ( pa.getPrimerTask() != P3Task.PICK_SEQUENCING_PRIMERS))
				retval.fwd.sort_primer_array();

			/** 
			 * If we are returning a list of internal oligos, sort them by their 'goodness'. We do not care if these are sorted if we end up in
			 * choose_pair_or_triple(), since this calls
			 * choose_internal_oligo(), which selects the best internal oligo
			 * for a given primer pair. 
			 */
			if (retval.output_type == P3OutputType.primer_list && pa.isPickInternalOligo())
				retval.intl.sort_primer_array();

			/* Select primer pairs if needed */
			if (retval.output_type == P3OutputType.primer_pairs) {
				choose_pair_or_triple(retval, pa, sa, dpal_arg_to_use, thal_arg_to_use,
						thal_oligo_arg_to_use);
			}

			if (pa.isDump()) {
				System.out.println("End of choose_primers:\n");
				p3_print_args(pa, sa) ;
			}

		}
		catch (Exception ex)
		{
			ex.printStackTrace();
			return null;
		}



		return retval;
	}


	/** ============================================================ */
	/** BEGIN choose_pair_or_triple
	 * 
	 * This function uses retval.fwd and retval.rev and
	 * updates the oligos in these array.
	 * 
	 * This function posibly uses retval.intl and updates its
	 * elements via choose_internal_oligo().
	 * 
	 * This function examines primer pairs or triples to find
	 * pa.num_return pairs or triples to return.
	 * 
	 *  Results are returned in best_pairs and in
	 *   retval.best_pairs.expl
	 * ============================================================ 
	 * @throws Exception */
	private static void choose_pair_or_triple(
			P3RetVal retval,
			P3GlobalSettings pa,
			SeqArgs sa,
			DPAlArgHolder dpal_arg_to_use,
			THAlArgHolder thal_arg_to_use,
			THAlArgHolder thal_oligo_arg_to_use
			) throws Exception {
		PairArrayT best_pairs = retval.best_pairs;


//		int i,j;   /* Loop index. */
		int n_int; /* Index of the internal oligo */
		/*int *max_j_seen; */  /* The maximum value of j (loop index for forward primers)
		                            that has been examined for every reverse primer
		                            index (i) -- global variable now */
		boolean update_stats = true;  /* Flag to indicate whether pair_stats
		                            should be updated. */
		PrimerPair h;             /* The current pair which is being evaluated. */
		PrimerPair the_best_pair = new PrimerPair(); /* The best pair is being "remembered". */
		PairStats pair_expl = retval.best_pairs.expl; /* For statistics */

		int product_size_range_index = 0;
		int the_best_i, the_best_j;

		/* Hash maps used to store pairs that were computed */

		/* std::hash_map<int, primer_pair*> **pairs; */
		/* pairs is an array of pointers to hash maps.  It will be indexed
		     by the indices of the reverse primers in retval.rev. -- global var now */
		HashMap<Integer, PrimerPair> hmap = null, best_hmap = null;
		/* hmap and best_hmap will be pointers to hash maps also pointed to
		by elements of pairs. */

		//		  std::hash_map<int, primer_pair*>::iterator it;
		PrimerPair pp = null, best_pp = null;
		boolean pair_found = false;

		HashMap<Integer, PrimerPair>[] pairs = new HashMap[retval.rev.num_elem];
		int[] max_j_seen = new int[retval.rev.num_elem];
		for (int i = 0; i < max_j_seen.length; i++) max_j_seen[i] = -1;

		while(true) {
			the_best_i = -1;
			the_best_j = -1;
			/* To start put penalty to the maximum */
			the_best_pair = new PrimerPair();
			the_best_pair.pair_quality = Double.MAX_VALUE;

			for (int i = 0; i < retval.rev.num_elem; i++) {
				
				// keep retval.rev.oligo.get(i) in right
				PrimerRecord right = retval.rev.oligo.get(i);
				hmap = pairs[i];  
				/* Pairs[i] is NULL if there has never been an assignment to
				 pairs[i] because pairs was allocated by calloc, which
				 sets the allocated memory to 0. */

				/* Only use a primer that *might be* legal or that the caller
		         has provided and specified as "must use".  Primers are *NOT*
		         FULLY ASSESSED until the call to characterize_pair(), in
		         order to avoid expensive computations (mostly alignments)
		         unless necessary. */
				if (!right.OK_OR_MUST_USE()) {
					/* Can free the memory used by the hmap associated to this reverse primer */
					if (hmap != null) {
//						 	for (it=hmap.begin(); it!=hmap.end(); it++) {
						/* it.second is the second element (i.e. the 'value', as opposed to the 'key'). */
//						 		pp = it.second;
//						 		delete pp;
//						 	}
						if (hmap == best_hmap) 
							best_hmap = null;
//						 	delete hmap;
						hmap.clear();
						hmap = null;
						pairs[i] = null;
					}
					continue;
				}

				/* If the pair cannot be better than the one already 
				 * selected, then we can skip remaining reverse primers */
				if (pa.getPrPairWeights().primer_quality *
						(right.quality + retval.fwd.oligo.get(0).quality)
						> the_best_pair.pair_quality) {
					break;
				}

				if (right.overlaps) {
					/* The stats will not keep track of the pair correctly
			            after the first pass, because an oligo might
			            have been legal on one pass but become illegal on
			            a subsequent pass. */
					if (update_stats) {
						if (choose_pair_or_triple_trace_me)
							System.err.format( "i=%d, j=%d, overlaps_oligo_in_better_pair++\n", i, 0); // this was j
						pair_expl.overlaps_oligo_in_better_pair++;
					}
					/* Can free the memory used by the hmap associated to this reverse primer */
					if (hmap != null ) {
						//						for (it=hmap.begin(); it!=hmap.end(); it++) {
						//							/* it.second is the second element (i.e. the 'value', as opposed to the 'key'). */
						//							pp = it.second;
						//							delete pp;
						//						}
						if (hmap == best_hmap) 
							best_hmap = null;
						//						delete hmap;
						hmap = null;
						pairs[i] = null;
					}
					continue;
				}


				/* Loop over forward primers */
				for (int j=0; j<retval.fwd.num_elem; j++) {
					PrimerRecord left = retval.fwd.oligo.get(j);

					/* We check the reverse oligo again, because we may
			           have determined that it is "not ok", even though
			           (as a far as we knew), it was ok above. */
					if ( ! right.OK_OR_MUST_USE()) {
						/* Can free the memory used by the hmap associated to this reverse primer */
						if (hmap != null) {
							//							for (it=hmap.begin(); it!=hmap.end(); it++) {
							//								/* it.second is the second element (i.e. the 'value', as opposed to the 'key'). */
							//								pp = it.second;
							//								delete pp;
							//							}
							if (hmap == best_hmap) 
								best_hmap = null;
							//							delete hmap;
							hmap = null;
							pairs[i] = null;
						}
						break;
					}

					/* Only use a primer that is legal, or that the caller
			           has provided and specified as "must use". */
					if (!left.OK_OR_MUST_USE()) continue;

					/* If the pair cannot be better than the one already 
					 * selected, then we can skip remaining forward primers 
					 * for this reverse primer */
					if (pa.getPrPairWeights().primer_quality *
							(left.quality + right.quality) 
							> the_best_pair.pair_quality) {
						break;
					}

					/* Need to have this here because if we break just above, then,
			           at a later iteration, we may need to examine the oligo
			           pair with reverse oligo at i and forward oligo at j. */
					update_stats = false;
					if (j > max_j_seen[i]) {
						if (choose_pair_or_triple_trace_me)
							System.err.format("updates ON: i=%d, j=%d, max_j_seen[%d]=%d\n",  i, j, i, max_j_seen[i]);
						max_j_seen[i] = j;
						if (choose_pair_or_triple_trace_me)
							System.err.format("max_j_seen[%d] -. %d\n", i, max_j_seen[i]);
						if (choose_pair_or_triple_trace_me) System.err.format( "updates on\n");
						update_stats = true;
					}

					if (left.overlaps) {
						/* The stats will not keep track of the pair correctly
			             after the first pass, because an oligo might
			             have been legal on one pass but become illegal on
			             a subsequent pass. */
						if (update_stats) {
							if (choose_pair_or_triple_trace_me)
								System.err.format("i=%d, j=%d, overlaps_oligo_in_better_pair++\n", i, j);
							pair_expl.overlaps_oligo_in_better_pair++;
						}
						continue;
					}

					/* Some simple checks first, before searching the hashmap */
					boolean must_use = false;
					if ((pa.getPrimerTask() == P3Task.CHECK_PRIMERS) || 
							((left.must_use != false) &&
									(right.must_use != false))) {
						must_use = true;
					}

					/* Determine if overlap with an overlap point is required, and
				   if so, whether one of the primers in the pairs overlaps
				   that point. */
					if ((sa.getPrimerOverlapJunctionsList().size() > 0)
							&& !(right.overlaps_overlap_position
									|| left.overlaps_overlap_position)
							) {
						if (update_stats) { 
							pair_expl.considered++;
							pair_expl.does_not_overlap_a_required_point++; 
						}
						if (!must_use) continue;
					}

					/* Check product size now */
					double product_size  = right.start - left.start+1;

					
					//  pa.pr_min[product_size_range_index]  pa.pr_max[product_size_range_index]
					if (product_size <  pa.getProductSizeRange(product_size_range_index).getLeft()  ||
							product_size > pa.getProductSizeRange(product_size_range_index).getRight()) {
						if (update_stats) {
							/* This line NEW */ 
							if (!must_use)
								pair_expl.considered++;
							pair_expl.product++; 
						}
						if (!must_use) continue;
					}

					/* Check if pair was already computed */
					pair_found = false;
					if (hmap != null) {
						if (hmap.containsKey(j)) {
							pair_found = true;

							/* it.second is the second element (i.e. the 'value', as opposed to the 'key'). */
							pp = hmap.get(j); 
							if (pp != null) { 
								/* The pair was computed, it isn't illegal and it wasn't selected yet */
								if (update_stats) {
									pair_expl.considered++;
									if (choose_pair_or_triple_trace_me)
										System.err.format("ok++\n");
									pair_expl.ok++;
								}
								/* Check if this is a better pair */
								if (compare_primer_pair(pp, the_best_pair) < 0) {
									the_best_pair = pp;
									the_best_i = i;
									the_best_j = j;
									best_hmap = hmap;
									best_pp = pp;
								}

								/* There cannot be a better pair */
								if (the_best_pair.pair_quality == 0) {
									break;
								} 
							} /* else - pp is NULL - it's illegal or already selected */
						}
					} else {
						/* Create this hashmap */
						hmap = new HashMap<Integer, PrimerPair>();
						pairs[i] = hmap;
					}

					if (!pair_found) {
						/* Characterize the pair. h is initialized by this call. */
						h = new PrimerPair();
						int tmp =  h.characterize_pair(retval, pa, sa, j, i,
								product_size_range_index, dpal_arg_to_use,
								thal_arg_to_use,
								update_stats);
						if (tmp == PrimerPair.PAIR_OK) {

							/* Choose internal oligo if needed */
							if (pa.isPickRightPrimer() && pa.isPickLeftPrimer()
									&& pa.isPickInternalOligo()) {
								n_int = choose_internal_oligo(retval, h.left, h.right, sa, pa, dpal_arg_to_use, thal_oligo_arg_to_use);
								if ( n_int == -1) {

									/* We were UNable to choose an internal oligo. */
									if (update_stats) { 
										pair_expl.internal++;
									}

									/* Mark the pair as not good - the entry in the hash map will be a NULL */
									hmap.put( j, null);
									continue;
								} else {
									/* We DID choose an internal oligo, and we
			                   set h.intl to point to it. */
									h.intl = retval.intl.oligo.get(n_int);
								}
							}

							if (update_stats) { 
								if (choose_pair_or_triple_trace_me)
									System.err.format("ok++\n");
								pair_expl.ok++;
							}

							/* Calculate the pair penalty */
							h.pair_quality = h.obj_fn(pa);
							//							PR_ASSERT(h.pair_quality >= 0.0);

							/* Save the pair */
//							pp = new primer_pair();

							pp = h;
							hmap.put(j,pp);

							/* The current pair (h) is the new best pair if it is better than the best pair so far. */
							if (compare_primer_pair(h, the_best_pair) < 0) {
								the_best_pair = h;
								the_best_i = i;
								the_best_j = j;
								best_hmap = hmap;
								best_pp = pp;
							}

							/* There cannot be a better pair */
							if (the_best_pair.pair_quality == 0) {
								break;
							}
						} else if (tmp == PrimerPair.PAIR_FAILED) {
							/* Illegal pair */
							hmap.put(j,null);
						} 
					}
				}  /* for (j=0; j<retval.fwd.num_elem; j++) -- inner loop */
				/* Check if there cannot be a better pair than best found */
			      
				if (the_best_pair.pair_quality == 0) {
			        break;
				}

			} /* for (i = 0; i < retval.rev.num_elem; i++) --- outer loop */
			

			if (the_best_pair.pair_quality == Double.MAX_VALUE ) {
				/* No pair was found. Try another product-size-range,
			       if one exists. */

				product_size_range_index++;
				/* Re-set the high-water marks for the indices
			         for reverse and forward primers:  */
				for (int i = 0; i < retval.rev.num_elem; i++) max_j_seen[i] = -1;

				if (!(product_size_range_index < pa.getProductSizeRangesNumber())) {
					/* We ran out of product-size-ranges. Exit the while loop. */
					break;

					/* Our bookkeeping was incorrect unless the assertion below is
				   true. If num_intervals > 1 or min_three_*_prime_distance >
				   -1 the assertion below might not be true. */
					//			        PR_ASSERT(!((pa.num_intervals == 1) &&  ((pa.min_left_three_prime_distance == -1) || (pa.min_right_three_prime_distance == -1))) || (best_pairs.num_pairs == pair_expl.ok));
				}

			} else {
				/* Store the best primer for output */

				if (choose_pair_or_triple_trace_me)
					System.err.format("ADD pair i=%d, j=%d\n", the_best_i, the_best_j);

				best_pairs.add_pair(the_best_pair);

				/* Mark the pair as already selected */
				//			      delete best_pp;
				best_hmap.put(the_best_j, null);

				/* Update the overlaps flags */
				for (int i = 0; i < retval.rev.num_elem; i++) {
					PrimerRecord right = retval.rev.oligo.get(i);
					if (right_oligo_in_pair_overlaps_used_oligo(right,
							the_best_pair,
							pa.getMinRight3PrimeDistance())) {
						right.overlaps = true;
					}
				}
				for (int j = 0; j < retval.fwd.num_elem; j++) {
					PrimerRecord left = retval.fwd.oligo.get(j);
					if (left_oligo_in_pair_overlaps_used_oligo(left,
							the_best_pair,
							pa.getMinLeft3PrimeDistance())) {
						left.overlaps = true;
					}
				}

				/* If we have enough then stop the while loop */
				if (pa.getNumReturn() == best_pairs.num_pairs) {
					break;
				}
			}

		}
	}

	


	private static boolean left_oligo_in_pair_overlaps_used_oligo(
			PrimerRecord left, PrimerPair best_pair,
			int min_dist) {
		int best_pos, pair_pos;

		if (min_dist == -1)
			return false;

		best_pos =  best_pair.left.start + best_pair.left.length - 1;

		pair_pos = left.start + left.length -  1;

		if ((Math.abs(best_pos - pair_pos) < min_dist)
				&& (min_dist != 0)) { return true; }

		if ((best_pair.left.length == left.length)
				&& (best_pair.left.start == left.start)
				&& (min_dist == 0)) {
			return true;
		}

		return false;
	}


	private static boolean right_oligo_in_pair_overlaps_used_oligo(
			PrimerRecord right, PrimerPair best_pair,
			int min_dist) {
		int best_pos, pair_pos;

		if (min_dist == -1)
			return false;

		best_pos = best_pair.right.start - best_pair.right.length + 1;

		pair_pos = right.start - right.length + 1;

		if ((Math.abs(best_pos - pair_pos) < min_dist)
				&& (min_dist != 0)) { return true; }

		if ((best_pair.right.length == right.length)
				&& (best_pair.right.start == right.start)
				&& (min_dist == 0)) {
			return true;
		}

		return false;
	}


	/**
	 * Choose best internal oligo for given pair of left and right primers.
	 * return -1 if it did not found one
	 */
	private static int choose_internal_oligo(P3RetVal retval,
			PrimerRecord left,
			PrimerRecord right,
			SeqArgs sa,
			P3GlobalSettings pa,
			DPAlArgHolder dpal_arg_to_use,
			THAlArgHolder thal_oligo_arg_to_use) throws AlignmentException, ThermodynamicAlignmentException {
		//		int i;
		int k;
		double min;
		char[] oligo_seq, revc_oligo_seq;
		PrimerRecord h;
		min = 1000000.0;
		//		i = -1;
		int nm = -1;
		for (k=0; k < retval.intl.num_elem; k++) {
			h = retval.intl.oligo.get(k);  /* h is the record for the oligo currently
			                                    under consideration */

			if ((h.start > (left.start + (left.length-1)))
					&& ((h.start + (h.length-1))
							< (right.start-right.length+1))
							&& (h.quality < min)
							&& (h.OK_OR_MUST_USE())) {

				if (h.self_any == ALIGN_SCORE_UNDEF && !pa.isThermodynamicOligoAlignment()) {

					oligo_seq =  Sequence._pr_substr(sa.getTrimmedSequence(), h.start, h.length);
					revc_oligo_seq = Sequence.p3_reverse_complement(oligo_seq);

					h.oligo_compl( pa.oligosArgs, retval.intl.expl,
							dpal_arg_to_use, oligo_seq, revc_oligo_seq);
					if (!h.OK_OR_MUST_USE()) continue;
				}

				if (h.self_any == ALIGN_SCORE_UNDEF && pa.isThermodynamicOligoAlignment()) {
					oligo_seq=  Sequence._pr_substr(sa.getTrimmedSequence(), h.start, h.length);
					revc_oligo_seq = Sequence.p3_reverse_complement(oligo_seq);

					h.oligo_compl_thermod( pa.oligosArgs, retval.intl.expl,
							thal_oligo_arg_to_use, oligo_seq, oligo_seq);
					if (!h.OK_OR_MUST_USE()) continue;
				}
				if(h.hairpin_th == ALIGN_SCORE_UNDEF && pa.isThermodynamicOligoAlignment()) {
					oligo_seq = Sequence._pr_substr(sa.getTrimmedSequence(), h.start, h.length);
					h.oligo_hairpin( pa.oligosArgs,
							retval.intl.expl, thal_oligo_arg_to_use,
							oligo_seq);
					if (!h.OK_OR_MUST_USE()) continue;
				}

				if (h.repeat_sim.score == null) {
					h.oligo_repeat_library_mispriming( pa, sa, OligoType.OT_INTL, retval.intl.expl,
							dpal_arg_to_use, retval.glob_err);
					if (!h.OK_OR_MUST_USE()) continue;
				}

				min = h.quality;
				nm=k;

			} /* if ((h.start.... */

		}  /* for (k=0;..... */

		//		nm = i;
		return nm;



	}

	static final double epsilon = 1e-6;

	/**
	 *  Compare function for sorting primer records.
	 */
	private static int compare_primer_pair(PrimerPair a1 , PrimerPair a2) {
		int y1, y2;

		if ((a1.pair_quality + epsilon) < a2.pair_quality) return -1;
		if (a1.pair_quality > (a2.pair_quality + epsilon)) return 1;

		/*
		 * The following statements ensure that we get a stable order that
		 * is the same on all systems.
		 */

		y1 = a1.left.start;
		y2 = a2.left.start;
		if (y1 > y2) return -1; /* prefer left primers to the right. */
		if (y1 < y2) return 1;

		y1 = a1.right.start;
		y2 = a2.right.start;
		if (y1 < y2) return -1; /* prefer right primers to the left. */
		if (y1 > y2) return 1;

		y1 = a1.left.length;
		y2 = a2.left.length;
		if (y1 < y2) return -1; /* prefer shorter primers. */
		if (y1 > y2) return 1;

		y1 = a1.right.length;
		y2 = a2.right.length;
		if (y1 < y2) return -1; /* prefer shorter primers. */
		if (y1 > y2) return 1;

		return 0;
	}


	



	private static int make_internal_oligo_list(P3RetVal retval,
			P3GlobalSettings pa, SeqArgs sa,
			DPAlArgHolder dpal_arg_to_use,
			THAlArgHolder thal_arg_to_use) throws PrimerRecordException, AlignmentException, ThermodynamicAlignmentException {


		int ret;
		int left = 0;
		retval.intl.extreme = 0;
		  /* Use the primer provided */
		  if ((sa.getInternalInput() != null ) || (pa.getPrimerTask() == P3Task.CHECK_PRIMERS )){
		         ret = add_one_primer(sa.getInternalInput(),  retval.intl,
		                               pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
		  }
		  else {
		    /* Pick all good in the given range */

		    /* Use the settings to select a proper range */
		   int length = sa.getTrimmedSequence().length - pa.getOligosArgs().getMinSize();
		   int start = pa.getOligosArgs().getMinSize() - 1;
		   left = 0;

		    ret = pick_primer_range(start, length, retval.intl,
		                            pa, sa, dpal_arg_to_use, thal_arg_to_use,
		                            retval);
		  }
		  return ret;
	}



	/**
	 * Make lists of acceptable left and right primers.  After return, the
	 * lists are stored in retval.fwd.oligo and retval.rev.oligo and the
	 * coresponding list sizes are stored in retval.fwd.num_elem and
	 * retval.rev.num_elem.  Return 1 if one of lists is empty or if
	 * leftmost left primer and rightmost right primer do not provide
	 * sufficient product size.
	 * @throws PrimerRecordException 
	 * @throws ThermodynamicAlignmentException 
	 * @throws AlignmentException 
	 */
	private static int make_detection_primer_lists(
			P3RetVal retval,
			P3GlobalSettings pa,
			SeqArgs sa,
			DPAlArgHolder dpal_arg_to_use,
			THAlArgHolder thal_arg_to_use) throws AlignmentException, ThermodynamicAlignmentException, PrimerRecordException 
			{
		int left, right;
		int length, start;
		int i,n, pr_min;
		int tar_l, tar_r, f_b, r_b;
		PairStats pair_expl = retval.best_pairs.expl; /* To store the statistics for pairs */

		/* Var to save the very left and very right primer */
		left = right = 0;

		/* Set pr_min to the very smallest
		     allowable product size. */
		pr_min = Integer.MAX_VALUE;
		for (i=0; i < pa.getProductSizeRangesNumber(); i++)
			if(pa.getProductSizeRange(i).getLeft() < pr_min)
				pr_min = pa.getProductSizeRange(i).getLeft();

		n = sa.getTrimmedSequence().length;
		//		  PR_ASSERT(INT_MAX > (n=strlen(sa.trimmed_seq)));
		tar_r = 0; /* Target start position */
		tar_l = n; /* Target length */

		/* Iterate over target array */
		for (i=0; i < sa.getTargetRegions().getCount(); i++) {

			/* Select the rightmost target start */
			if (sa.getTargetRegions().getInterval(i)[0] > tar_r)
				tar_r = sa.getTargetRegions().getInterval(i)[0];

			/* Select the rightmost target end */
			if (sa.getTargetRegions().getInterval(i)[0] + sa.getTargetRegions().getInterval(i)[1] - 1 < tar_l)
				tar_l = sa.getTargetRegions().getInterval(i)[0] + sa.getTargetRegions().getInterval(i)[1] - 1;
		}

		if (pa.isDefaultPositionPenalties()) {
			if (0 == tar_r) tar_r = n;
			if (tar_l == n) tar_l = 0;
		} else {
			tar_r = n;
			tar_l = 0;
		}

		/* We use some global information to restrict the region
		     of the input sequence in which we generate candidate
		     oligos. */
		if (retval.output_type == P3OutputType.primer_list && pa.isPickLeftPrimer())
			f_b = n - 1;
		else if (tar_r - 1 < n - pr_min + pa.primersArgs.getMaxSize() - 1
				&& !(pa.isPickAnyway() && sa.getLeftInput() != null ))
			f_b=tar_r - 1;
		else
			f_b = n - pr_min + pa.primersArgs.getMaxSize()-1;

		if (pa.isPickLeftPrimer()) {
			/* We will need a left primer. */
			left=n; right=0;
			length = f_b - pa.primersArgs.getMinSize() + 1;
			start = pa.primersArgs.getMinSize() - 1;

			// just in case

			retval.fwd.extreme =  left;


			/* Use the primer provided */
			if (sa.getLeftInput() != null) {

				add_one_primer(sa.getLeftInput(), retval.fwd,
						pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
				//				left = retval.fwd.extreme;

			}
			/* Pick primers at one position */
			else if(sa.getForceLeftStart() > -1 ||
					sa.getForceLeftEnd() > -1) {
				pick_primers_by_position(sa.getForceLeftStart(), sa.getForceLeftEnd(), retval.fwd, pa, sa,
						dpal_arg_to_use, thal_arg_to_use, retval);
				//				left = retval.fwd.extreme;
			}
			/* Or pick all good in the given range */
			else {
				pick_primer_range(start, length,  retval.fwd,
						pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
				//				left = retval.fwd.extreme;
			}
			left = retval.fwd.extreme;

		}  /* if (pa.pick_left_primer) */
		if (retval.output_type == P3OutputType.primer_list && pa.isPickRightPrimer())
			r_b = 0;
		else if (tar_l+1>pr_min - pa.primersArgs.getMaxSize()
				&& !(pa.isPickAnyway() && sa.getRightInput() != null))
			r_b = tar_l+1;
		else
			r_b = pr_min - pa.primersArgs.getMaxSize();

		if ( pa.isPickRightPrimer() ) {
			// init extreme
			retval.fwd.extreme = right;
			/* We will need a right primer */
			length = n-pa.primersArgs.getMinSize() - r_b + 1;
			start = r_b;

			/* Use the primer provided */
			if (sa.getRightInput() != null) {
				add_one_primer(sa.getRightInput(), retval.rev,
						pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
				/*  pick_right_primers(start, length, &right, &retval.rev,
		          pa, sa, dpal_arg_to_use, retval);*/

			}
			/* Pick primers at one position */
			else if(sa.getForceRightStart() > -1 ||
					sa.getForceRightEnd() > -1) {
				pick_primers_by_position(sa.getForceRightStart(), sa.getForceRightEnd(),
						retval.rev, pa, sa,
						dpal_arg_to_use, thal_arg_to_use, retval);
			}
			/* Or pick all good in the given range */
			else {
				pick_primer_range(start, length,  retval.rev,
						pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
			}

			right = retval.rev.extreme;
		}

		/*
		 * Return 1 if either the left primer list or the right primer
		 * list is empty or if leftmost left primer and
		 * rightmost right primer do not provide sufficient product size.
		 */

		if ((pa.isPickLeftPrimer() && 0 == retval.fwd.num_elem)
				|| ((pa.isPickRightPrimer())  && 0 == retval.rev.num_elem)) {
			return 1;
		} else if (!((sa.getRightInput() != null)
				&& (sa.getLeftInput() != null))
				&& pa.isPickLeftPrimer()
				&& pa.isPickRightPrimer()
				&& (right - left) < (pr_min - 1)) {
			pair_expl.product    = 1;
			pair_expl.considered = 1;
			return 1;
		} else 
			return 0;


			}



	private static int pick_primers_by_position(int start,
			int end, 
			OligoArray oligo,
			P3GlobalSettings pa,
			SeqArgs sa,
			DPAlArgHolder dpal_arg_to_use,
			THAlArgHolder thal_arg_to_use,
			P3RetVal retval) throws AlignmentException, ThermodynamicAlignmentException, PrimerRecordException {
		int found_primer, length, j, ret, new_start;
		found_primer = 1;
		ret = 1;

		if(start > -1 && end > -1) {
			if (oligo.type != OligoType.OT_RIGHT) {
				length = end - start + 1;
			} else {
				length = start - end + 1;
			}

			found_primer = add_one_primer_by_position(start, length, oligo,
					pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
			return found_primer;
		} else if (start > -1) {
			/* Loop over possible primer lengths, from min to max */
			ret = 0;
			for (j = pa.primersArgs.getMinSize(); j <= pa.primersArgs.getMaxSize(); j++) {
				ret += add_one_primer_by_position(start, j, oligo,
						pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
				if (ret == 0) {
					found_primer = 0;
				}
			}
			return found_primer;
		} else if (end > -1) {
			/* Loop over possible primer lengths, from min to max */
			ret = 0;
			for (j = pa.primersArgs.getMinSize(); j <= pa.primersArgs.getMaxSize(); j++) {
				if (oligo.type != OligoType.OT_RIGHT) {
					new_start = end - j + 1;
				} else {
					new_start = end + j - 1;
				}
				ret += add_one_primer_by_position(new_start, j, oligo,
						pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
				if (ret == 0) {
					found_primer = 0;
				}
			}
			return found_primer;
		} else {
			/* Actually this should never happen */
			retval.warnings.append("Calculation error in forced primer position calculation");
			return 1;
		}		
	}



	private static int add_one_primer_by_position(int start, int length,
			OligoArray oligo, P3GlobalSettings pa, SeqArgs sa,
			DPAlArgHolder dpal_arg_to_use, THAlArgHolder thal_arg_to_use,
			P3RetVal retval) throws AlignmentException, ThermodynamicAlignmentException, PrimerRecordException {
		/* Variables for the loop */
		int i;
		int n, found_primer;

		/* Array to store one primer sequences in */
		char[] oligo_seq;

		/* Struct to store the primer parameters in */
		PrimerRecord h = new PrimerRecord();
		//		  memset(&h, 0, sizeof(primer_rec));

		/* Retun 1 for no primer found */
		found_primer = 1;

		n = sa.getTrimmedSequence().length;
		//		  PR_ASSERT(INT_MAX > (n=strlen(sa.trimmed_seq)));

		/* Just to be sure */
		if (start < 0) {
			return 1;
		}
		if (start >= n) {
			return 1;
		}
		if (oligo.type != OligoType.OT_RIGHT) {
			if ((start + length) > n) {
				return 1;
			}
		} else {
			if ((start - length + 1) < 0) {
				return 1;
			}
		}

		//		  oligo_seq[0] = '\0';

		/* Set the length of the primer */
		h.length = length;

		/* Figure out positions for forward primers */
		if (oligo.type != OligoType.OT_RIGHT) {
			/* Set the start of the primer */
			h.start = start;

			/* Put the real primer sequence in s */
			oligo_seq = Sequence._pr_substr(sa.getTrimmedSequence(), h.start, length);
		}
		/* Figure out positions for reverse primers */
		else {
			i = start - length + 1;
			/* Set the start of the primer */
			h.start = start;

			/* Put the real primer sequence in s */
			oligo_seq = Sequence._pr_substr(sa.getTrimmedSequence(), i, length);
		}

		/* Force primer3 to use this oligo */
		h.must_use = (true && pa.isPickAnyway());

		h.overlaps = false;

		/* Add it to the considered statistics */
		oligo.expl.considered++;

		/* Calculate all the primer parameters */
		h.calc_and_check_oligo_features(pa,  oligo.type, dpal_arg_to_use, thal_arg_to_use,
				sa, oligo.expl, retval, oligo_seq);

		/* If primer has to be used or is OK */
		if (h.OK_OR_MUST_USE()) {
			/* Calculate the penalty */
			h.quality = h.p_obj_fn(pa, oligo.type);
			/* Save the primer in the array */
			oligo.add_oligo_to_oligo_array( h);
			found_primer = 0;
			/* Update the most extreme primer variable */
			//		      if (( h.start < *extreme) &&
			//		          (oligo.type != OT_RIGHT))
			//		        *extreme =  h.start;
			/* Update the most extreme primer variable */
			//		      if (( h.start > *extreme) &&
			//		          (oligo.type == OT_RIGHT))
			//		        *extreme =  h.start;
			/* Update the number of primers */
		} else {
			/* Free memory used by this primer. */
			h.free_primer_repeat_sim_score();
		}
		/* Update array with how many primers are good */
		/* Update statistics with how many primers are good */
		oligo.expl.ok = oligo.num_elem;

		/* return 0 for success */
		return found_primer;
	}



	/**
	 * add_one_primer finds one primer in the trimmed sequence and stores
	 * it in *oligo The main difference to the general fuction is that it
	 * calculates its length and it will add a primer of any length to the
	 * list 
	 * @throws ThermodynamicAlignmentException 
	 * @throws AlignmentException 
	 * @throws PrimerRecordException 
	 */
	private static int add_one_primer(char[] primer, 
			//			int extreme,
			OligoArray oligo,
			P3GlobalSettings pa,
			SeqArgs sa,
			DPAlArgHolder dpal_arg_to_use,
			THAlArgHolder thal_arg_to_use,
			P3RetVal retval) throws AlignmentException, ThermodynamicAlignmentException, PrimerRecordException {

		/* Variables for the loop */
		int i, j;
		int n;

		/* Array to store one primer sequences in */
		char[] oligo_seq , test_oligo;

		
		//		  memset(&h, 0, sizeof(primer_rec));

		/* Copy *primer into test_oligo */
		//		  test_oligo[0] = '\0';
		if (oligo.type != OligoType.OT_RIGHT) {
			test_oligo  = Sequence.subSeq(primer,0, primer.length);
		} else {
			test_oligo = Sequence.p3_reverse_complement(primer );
		}
		n = sa.getTrimmedSequence().length;
		//		  PR_ASSERT(INT_MAX > (n)));

		/* This time we already know the size of the primer */
		j = (primer.length);

		/* Loop over the whole sequence */
		for(i = sa.getTrimmedSequence().length; i >= 0; i--) {
			//		    oligo_seq[0] = '\0';
			/* Struct to store the primer parameters in */
			PrimerRecord h = new PrimerRecord();
			/* Set the length of the primer */
			h.length = j;

			/* Figure out positions for forward primers */
			if (oligo.type != OligoType.OT_RIGHT) {
				/* Break if the primer is bigger than the sequence left*/
				if(i-j <= -1) continue;

				/* Set the start of the primer */
				h.start = i - j;//  +1;

				/* Put the real primer sequence in s */
//				if(h.start+ j > sa.trimmed_seq.length)
//					j = j - ((h.start+ j)- sa.trimmed_seq.length);
				oligo_seq =  Sequence._pr_substr(sa.getTrimmedSequence(), h.start, j );
			}
			/* Figure out positions for reverse primers */
			else {
				/* Break if the primer is bigger than the sequence left*/
				if(i+j>n) continue;

				/* Set the start of the primer */
				h.start=i+j-1;

				/* Put the real primer sequence in s */
				oligo_seq =  Sequence._pr_substr(sa.getTrimmedSequence(),  i, j);
			}

			/* Compare the primer with the sequence */
			if (strcmp_nocase(test_oligo, oligo_seq) != 0)
				continue;

			/* Force primer3 to use this oligo */
			h.must_use = (true && pa.isPickAnyway());

			h.overlaps = false;

			/* Add it to the considered statistics */
			oligo.expl.considered++;

			/* Calculate all the primer parameters */
			h.calc_and_check_oligo_features(pa,  oligo.type, dpal_arg_to_use, thal_arg_to_use,
					sa, oligo.expl, retval, oligo_seq);

			/* If primer has to be used or is OK */
			if (h.OK_OR_MUST_USE()) {
				/* Calculate the penalty */
				h.quality = h.p_obj_fn(pa, oligo.type);
				/* Save the primer in the array */
				oligo.add_oligo_to_oligo_array( h);
				/* Update the most extreme primer variable */
				//				if ((h.start < extreme) && (oligo.type != oligo_type.OT_RIGHT))
				//					extreme = h.start;
				/* Update the most extreme primer variable */
				//				if ((h.start > extreme) && (oligo.type == oligo_type.OT_RIGHT))
				//					extreme = h.start;
			} else {
				/* Free memory used by this primer. */
				h.free_primer_repeat_sim_score();
			}
		} /* i: Loop over the sequence */
		/* Update array with how many primers are good */
		/* Update statistics with how many primers are good */
		oligo.expl.ok = oligo.num_elem;

		if (oligo.num_elem == 0) return 1;
		else {
			if (oligo.num_elem > 1) {
				retval.warnings.append("More than one position in template for input oligo " + Arrays.toString(primer));
				//		      pr_append(&retval.warnings, );
			}
			return 0; /* Success */
		}

	}


	/**
	 * Add caller-specified primers to retval.fwd, retval.rev, and/or
	 * retval.intl.  The primers do not have to be "legal". It is possible
	 * to add more than one primer to retval.fwd, etc, if the input
	 * primer is found in multiple locations in the template. Similar to
	 * make_complete_primer_lists.
	 * @throws PrimerRecordException 
	 * @throws ThermodynamicAlignmentException 
	 * @throws AlignmentException 
	 */
	private static int add_primers_to_check(P3RetVal retval, P3GlobalSettings pa,
			SeqArgs sa, DPAlArgHolder dpal_arg_to_use,
			THAlArgHolder thal_arg_to_use,
			THAlArgHolder thal_oligo_arg_to_use) throws AlignmentException, ThermodynamicAlignmentException, PrimerRecordException {
	      

		if (sa.getLeftInput() != null) {
			add_one_primer(sa.getLeftInput(),  retval.fwd,
					pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
		}

		if (sa.getRightInput()!= null) {
			add_one_primer(sa.getRightInput(),  retval.rev,
					pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
		}

		if (sa.getInternalInput()!= null) {
			add_one_primer(sa.getInternalInput(), retval.intl,
					pa, sa, dpal_arg_to_use, thal_oligo_arg_to_use, retval);
		}

		return 0;
	}


	/**
	 * Make lists of acceptable left and right primers.  After return, the
	 * lists are stored in retval.fwd.oligo and retval.rev.oligo and the
	 * coresponding list sizes are stored in retval.fwd.num_elem and
	 * retval.rev.num_elem.  Return 1 if one of lists is empty or if
	 * leftmost left primer and rightmost right primer do not provide
	 * sufficient product size.
	 * @throws PrimerRecordException 
	 * @throws ThermodynamicAlignmentException 
	 * @throws AlignmentException 
	 */
	private static void pick_sequencing_primer_list(P3RetVal retval,
			P3GlobalSettings pa, SeqArgs sa,
			DPAlArgHolder dpal_arg_to_use,
			THAlArgHolder thal_arg_to_use) throws AlignmentException, ThermodynamicAlignmentException, PrimerRecordException {
		int length, start;
		int n, rest_accuracy;
		int primer_nr; /* number of primers we need to pick */
		int tar_n; /* counter for the targets */
		int step_nr; /* counter to step through the targets */
		int sequenced_len; /* bp sequenced in good quality */
		int extra_seq; /* bp sequenced additionally on both sides */

		/* best location for the 3' end of the fwd primer: */
		int pr_position_f;

		/* best location for the 3' end of the rev primer: */
		int pr_position_r;

		/* Get the length of the sequence */
		n= sa.getTrimmedSequence().length;
		//		  PR_ASSERT(INT_MAX > (N));

		/* For each target needed loop*/
		for (tar_n=0; tar_n < sa.getTargetRegions().getCount(); tar_n++) {

			/* Calculate the amount of primers needed */
			primer_nr = 1;
			if ((pa.isPickLeftPrimer()) && (pa.isPickRightPrimer())){
				sequenced_len = pa.getSequencingParameters().getInterval();
				while(sequenced_len < sa.getTargetRegions().getInterval(tar_n)[1]) {
					primer_nr++;
					sequenced_len = pa.getSequencingParameters().getSpacing() * (primer_nr - 1)
							+ pa.getSequencingParameters().getInterval();
				}
			} else {
				sequenced_len = pa.getSequencingParameters().getSpacing();
				while(sequenced_len < sa.getTargetRegions().getInterval(tar_n)[1]) {
					primer_nr++;
					sequenced_len = pa.getSequencingParameters().getSpacing() * primer_nr;
				}
			}
			/* Calculate the overlap on the sides */
			extra_seq = (sequenced_len - sa.getTargetRegions().getInterval(tar_n)[1]) / 2;

			/* Pick primers for each position */
			for ( step_nr = 0 ; step_nr < primer_nr ; step_nr++ ) {
				pr_position_f = sa.getTargetRegions().getInterval(tar_n)[0] - extra_seq
						+ ( pa.getSequencingParameters().getSpacing() * step_nr )
						- pa.getSequencingParameters().getLead();
				if ((pa.isPickLeftPrimer()) && (pa.isPickRightPrimer())) {
					pr_position_r = sa.getTargetRegions().getInterval(tar_n)[0] - extra_seq
							+ ( pa.getSequencingParameters().getSpacing() * step_nr )
							+ pa.getSequencingParameters().getInterval()
							+ pa.getSequencingParameters().getLead();
				} else {
					pr_position_r = sa.getTargetRegions().getInterval(tar_n)[0] - extra_seq
							+ ( pa.getSequencingParameters().getSpacing() * (step_nr+1))
							+ pa.getSequencingParameters().getLead();        
				}
				/* Check if calculated positions make sense */
				/* position_f cannot be outside included region */
				if (pr_position_f < (pa.primersArgs.getMinSize() -1)) {
					pr_position_f = pa.primersArgs.getMinSize() - 1;
				}
				if (pr_position_f > (n - pa.primersArgs.getMinSize() - 1)) {
					pr_position_f = n - pa.primersArgs.getMinSize() - 1;
					/* Actually this should never happen */
					retval.warnings.append("Calculation error in forward sequencing position calculation");
				}
				/* position_r cannot be outside included region */
				if (pr_position_r < (pa.primersArgs.getMinSize() - 1)) {
					pr_position_r = pa.primersArgs.getMinSize() - 1;
					/* Actually this should never happen */
					retval.warnings.append("Calculation error in reverse sequencing position calculation");
				}
				if (pr_position_r > (n - pa.primersArgs.getMinSize() - 1)) {
					pr_position_r = n - pa.primersArgs.getMinSize() - 1;
				}
				/* Now all pr_positions are within the sequence */
				if (pa.isPickLeftPrimer()) {
					/* Set the start and length for the regions */
					start = pr_position_f - pa.getSequencingParameters().getAccuracy();
					if (start < 0) {
						rest_accuracy = pr_position_f + 1;
						start = 0;
					} else {
						rest_accuracy = pa.getSequencingParameters().getAccuracy();
					}
					length = rest_accuracy + pa.getSequencingParameters().getAccuracy() ;
					if ((start + length) > n) {
						length = n - start;
					}
					/* Pick all good in the given range */
					pick_only_best_primer(start, length, retval.fwd,
							pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
				}
				if (pa.isPickRightPrimer()) {
					start = pr_position_r - pa.getSequencingParameters().getAccuracy();
					if (start < 0) {
						rest_accuracy = pr_position_r + 1;
						start = 0;
					} else {
						rest_accuracy = pa.getSequencingParameters().getAccuracy();
					}
					length = rest_accuracy + pa.getSequencingParameters().getAccuracy() ;
					if ((start + length) > n) {
						length = n - start;
					}
					/* Pick all good in the given range */
					pick_only_best_primer(start, length, retval.rev,
							pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
				}
			}

		} /* End of Target Loop */

		/* Print an error if not all primers will be printed */
		if (retval.fwd.num_elem > pa.getNumReturn() 
				|| retval.rev.num_elem > pa.getNumReturn()) {
			retval.warnings.append("Increase PRIMER_NUM_RETURN to obtain all sequencing primers");
		}

		//		  return 0;
	}

	/**
	 * pick_primer_range picks all primers which have their 3' end in the range
	 * from start to start+length and store only the best oligo  
	 * @throws ThermodynamicAlignmentException 
	 * @throws AlignmentException 
	 * @throws PrimerRecordException 
	 */
	private static int pick_only_best_primer(int start, int length,
			OligoArray oligo, P3GlobalSettings pa, SeqArgs sa,
			DPAlArgHolder dpal_arg_to_use, THAlArgHolder thal_arg_to_use,
			P3RetVal retval) throws AlignmentException, ThermodynamicAlignmentException, PrimerRecordException {
		/* Variables for the loop */
		int i, j, primer_size_small, primer_size_large;
		int n, found_primer;
		//		  char p_number = number[0];
		//		  int temp_value;

		/* Array to store one primer sequences in */
		char[] oligo_seq ;//[MAX_PRIMER_LENGTH+1];

		/* Struct to store the primer parameters in */
		PrimerRecord h;
		PrimerRecord best = null;
		//		  memset(&h, 0, sizeof(primer_rec));
		//		  memset(&best, 0, sizeof(primer_rec));
		//		  best.quality = 1000.00;
		double best_quality = 1000.00;
		found_primer = 0;

		/* Set n to the length of included region */
		n=sa.getTrimmedSequence().length;
		//		  PR_ASSERT(INT_MAX > (n));

		/* Conditions for primer length */
		if (oligo.type == OligoType.OT_INTL) {
			primer_size_small=pa.oligosArgs.getMinSize();
			primer_size_large=pa.oligosArgs.getMaxSize();
		}
		else {
			primer_size_small=pa.primersArgs.getMinSize();
			primer_size_large=pa.primersArgs.getMaxSize();
		}

		/* Loop over locations in the sequence */
		for(i = start + length - 1; i >= start; i--) {
			//		    oligo_seq[0] = '\0';

			/* Loop over possible primer lengths, from min to max */
			for (j = primer_size_small; j <= primer_size_large; j++) {
				h = new PrimerRecord();  
				/* Set the length of the primer */
				h.length = j;

				/* Set repeat_sim to NULL as indicator that the repeat_sim
		         struct is not initialized. */
				h.repeat_sim.score = null;

				/* Figure out positions for left primers and internal oligos */
				if (oligo.type != OligoType.OT_RIGHT) {
					/* Break if the primer is bigger than the sequence left */
					if(i-j < -1) continue;

					/* Set the start of the primer */
					h.start = i - j + 1;

					/* Put the real primer sequence in oligo_seq */
					oligo_seq =  Sequence._pr_substr(sa.getTrimmedSequence(), h.start, j);
				} else {
					/* Figure out positions for reverse primers */
					/* Break if the primer is bigger than the sequence left*/
					if(i+j > n) continue;

					/* Set the start of the primer */
					h.start = i+j-1;

					/* Put the real primer sequence in s */
					oligo_seq =  Sequence._pr_substr(sa.getTrimmedSequence(),  i, j);

				}

				/* Do not force primer3 to use this oligo */
				h.must_use = false;

				h.overlaps = false;

				/* Add it to the considered statistics */
				oligo.expl.considered++;
				/* Calculate all the primer parameters */
				h.calc_and_check_oligo_features(pa, oligo.type, dpal_arg_to_use, thal_arg_to_use,
						sa, oligo.expl, retval, oligo_seq);

				/* If primer has to be used or is OK */
				if (h.OK_OR_MUST_USE()) {
					/* Calculate the penalty */
					h.quality = h.p_obj_fn(pa, oligo.type);
					/* Save the primer in the array */
					if (h.quality < best_quality) {
						/* Free memory used by previous best primer. */
						//		        	best.free_primer_repeat_sim_score();
						best_quality =  h.quality;
						best = h;
						found_primer = 1;
					} else {
						/* Free memory used by this primer. */
						h.free_primer_repeat_sim_score();
					}
				}
				else {
					/* Free memory used by this primer. */
					h.free_primer_repeat_sim_score();
					if (h.any_5_prime_ol_extension_has_problem()) {
						/* Break from the inner for loop, because there is no
			     legal longer oligo with the same 3' sequence. */
						break;
					}
				}
			} /* j: Loop over possible primer length from min to max */
		} /* i: Loop over the sequence */

		if (found_primer == 1) {
			/* Add the best to the array */
			oligo.add_oligo_to_oligo_array( best);
			/* Update statistics with how many primers are good */
			oligo.expl.ok = oligo.expl.ok + 1;
		} else {
			if (oligo.type == OligoType.OT_RIGHT) {
				retval.warnings.append( "No right primer found in range ");
			} else {
				retval.warnings.append("No left primer found in range ");
			}

			retval.warnings.append((start + pa.getFirstBaseIndex()) + " - " + (start + length + pa.getFirstBaseIndex()) ); 
		}
		
		if (oligo.num_elem == 0) 
			return 1;
		else 
			return 0;		
	}


	/**
	 * Make lists of acceptable left and right primers.  After return, the
	 * lists are stored in retval.fwd.oligo and retval.rev.oligo and the
	 * coresponding list sizes are stored in retval.fwd.num_elem and
	 * retval.rev.num_elem.  Return 1 if one of lists is empty or if
	 * leftmost left primer and rightmost right primer do not provide
	 * sufficient product size.
	 * @throws PrimerRecordException 
	 * @throws ThermodynamicAlignmentException 
	 * @throws AlignmentException 
	 */
	private static int make_complete_primer_lists(P3RetVal retval,
			P3GlobalSettings pa, SeqArgs sa,
			DPAlArgHolder dpal_arg_to_use2, THAlArgHolder thal_arg_to_use2,
			THAlArgHolder thal_oligo_arg_to_use2) throws PrimerRecordException, AlignmentException, ThermodynamicAlignmentException {


		// extreme moved to oligo_array

		//		int extreme;
		int length, start;
		int n;
		n=sa.getTrimmedSequence().length;
		/* Get the length of the sequence */
		//		PR_ASSERT(INT_MAX > (n=strlen(sa.trimmed_seq)));
		if (pa.isPickLeftPrimer()) {
			/* We will need a left primer. */
			retval.fwd.extreme = 0;
			length = n - pa.primersArgs.getMinSize();
			start = pa.primersArgs.getMinSize() - 1;

			/* Pick all good in the given range */
			pick_primer_range(start, length,retval.fwd,
					pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);

		}  /* if (pa.pick_left_primer) */
		if ( pa.isPickRightPrimer() ) {
			/* We will need a right primer */
			retval.rev.extreme = n;
			length = n - pa.primersArgs.getMinSize() + 1;
			start = 0;

			/* Pick all good in the given range */
			pick_primer_range(start, length,  retval.rev,
					pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
		}

		if ( pa.isPickInternalOligo() ) {
			/* We will need a internal oligo */
			length = n - pa.oligosArgs.getMinSize();
			start = pa.oligosArgs.getMinSize() - 1;
			retval.intl.extreme = 0;

			/* Pick all good in the given range */
			pick_primer_range(start, length,  retval.intl,
					pa, sa, dpal_arg_to_use, thal_oligo_arg_to_use, retval);
		}

		return 0;		
	}


	/**
	 * pick_primer_range picks all legal primers in the range from start
	 * to start+length and stores them in *oligo  
	 * @throws PrimerRecordException 
	 * @throws ThermodynamicAlignmentException 
	 * @throws AlignmentException 
	 */
	private static int pick_primer_range(int start, int length,
			OligoArray oligo,
			P3GlobalSettings pa,
			SeqArgs sa,
			DPAlArgHolder dpal_arg_to_use,
			THAlArgHolder thal_arg_to_use,
			P3RetVal retval) throws PrimerRecordException, AlignmentException, ThermodynamicAlignmentException {
		/* Variables for the loop */
		int i, j;
		int primer_size_small, primer_size_large;
		int pr_min, n;

		/* Array to store one primer sequences in */
		char[] oligo_seq = null ; //  new char[MAX_PRIMER_LENGTH+1];

		/* Struct to store the primer parameters in */
		PrimerRecord h = new PrimerRecord();
		//		memset(&h, 0, sizeof(primer_rec));

		/* Set pr_min to the very smallest
			allowable product size. */
		pr_min = Integer.MAX_VALUE;
		for (i=0; i < pa.getProductSizeRangesNumber(); i++)
			if(pa.getProductSizeRange(i).getLeft() < pr_min)
				pr_min = pa.getProductSizeRange(i).getLeft();

		/* Set n to the length of included region */
		n = sa.getTrimmedSequence().length;
		//		PR_ASSERT(INT_MAX > (n=strlen(sa.trimmed_seq)));

		if (oligo.type == OligoType.OT_INTL) {
			primer_size_small=pa.oligosArgs.getMinSize(); 
			primer_size_large=pa.oligosArgs.getMaxSize();
		}
		else {
			primer_size_small=pa.primersArgs.getMinSize();
			primer_size_large=pa.primersArgs.getMaxSize();
		}

		/* Loop over locations in the sequence */
		for(i = start + length; i >= start; i--) {
			//			oligo_seq[0] = '\0';

			/* Loop over possible primer lengths, from min to max */
			for (j = primer_size_small; j <= primer_size_large; j++) {
				h = new PrimerRecord();
				/* Set the length of the primer */
				h.length = j;

				
				/* Figure out positions for left primers and internal oligos */
				if (oligo.type != OligoType.OT_RIGHT) {
					/* Check if the product is of sufficient size -- could be optimized herer */
					if (i-j > n-pr_min-1 && retval.output_type == P3OutputType.primer_pairs
							&& oligo.type == OligoType.OT_LEFT) 
						continue;

					/* Break if the primer is bigger than the sequence left */
					if (i-j < -1) 
						break;

					/* Set the start of the primer */
					h.start = i - j + 1;

					/* Put the real primer sequence in oligo_seq */
					oligo_seq = Sequence._pr_substr(sa.getTrimmedSequence(), h.start, j);
				} else {
					/* Figure out positions for reverse primers */
					/* Check if the product is of sufficient size */
					if (i+j < pr_min && retval.output_type == P3OutputType.primer_pairs) continue;

					/* Break if the primer is bigger than the sequence left*/
					if (i+j > n) 
						break;

					/* Set the start of the primer */
					h.start = i+j-1;

					/* Put the real primer sequence in s */
					oligo_seq = Sequence._pr_substr(sa.getTrimmedSequence(),  i, j);

				}
				
				
				/* Do not force primer3 to use this oligo */
				h.must_use = false;

				h.overlaps = false;

				/* Add it to the considered statistics */
				oligo.expl.considered++;
				/* Calculate all the primer parameters */
				h.calc_and_check_oligo_features(pa, oligo.type, dpal_arg_to_use, thal_arg_to_use,
						sa, oligo.expl, retval, oligo_seq);
				/* If primer has to be used or is OK */
				if (h.OK_OR_MUST_USE()) {
					/* Calculate the penalty */
					h.quality = h.p_obj_fn(pa,oligo.type);
					/* Save the primer in the array */
					oligo.add_oligo_to_oligo_array( h);
					/* Update the most extreme primer variable */
					if (( h.start < oligo.extreme) && (oligo.type != OligoType.OT_RIGHT))
						oligo.extreme = h.start;
					/* Update the most extreme primer variable */
					if ((h.start > oligo.extreme) && (oligo.type == OligoType.OT_RIGHT))
						oligo.extreme = h.start;
				} else {
					/* Free memory used by this primer. */
					h.free_primer_repeat_sim_score();
					if (h.any_5_prime_ol_extension_has_problem()) {
						/* Break from the inner for loop, because there is no legal longer oligo with the same 3' sequence. */
						break;
					}
				}
			} /* j: Loop over possible primer length from min to max */
		} /* i: Loop over the sequence */

		/* Update statistics with how many primers are good */
		oligo.expl.ok = oligo.num_elem;

		if (oligo.num_elem == 0) return 1;
		else return 0;		
	}







	/**
	 * 
	 * Return 1 on error, 0 on success.  Set sa.trimmed_seq and possibly modify
	 * sa.tar.  Upcase and check all bases in sa.trimmed_seq.
	 * TO DO -- this would probably be cleaner if it only
	 * checked, rather than updated, sa.
	 * Check if the input in sa and pa makes sense 
	 * @param pa
	 * @param sa
	 * @param glob_err
	 * @param nonfatal_err
	 * @param warning
	 * @return false 
	 */
	private static boolean _pr_data_control(P3GlobalSettings pa, SeqArgs sa,
			StringBuilder glob_err,
			StringBuilder nonfatal_err,
			StringBuilder warning) {
		//		  char s1[MAX_PRIMER_LENGTH+1];
		int i, pr_min, seq_len;
		char offending_char = '\0';

		seq_len = sa.getSequence().length;

		/* If sequence quality is provided, is it as long as the sequence? */
		if (sa.getSequenceQuality() != null && sa.getSequenceQuality().length != seq_len)
			nonfatal_err.append(
					"Error in sequence quality data");

		if ((pa.getMinLeft3PrimeDistance() < -1) ||
				(pa.getMinRight3PrimeDistance() < -1))
			nonfatal_err.append( "Minimum 3' distance must be >= -1 (min_*_three_prime_distance)");

		if ((pa.primersArgs.getMinQuality() != 0 || pa.oligosArgs.getMinQuality() != 0)
				&& (sa.getSequenceQuality() == null || sa.getSequenceQuality().length == 0 ) )
			nonfatal_err.append( "Sequence quality data missing");

		if (pa.getFirstBaseIndex() < PR_NULL_FORCE_POSITION) {
			glob_err.append( "Value too small at tag PRIMER_FIRST_BASE_INDEX");
			return true;
		}

		if (pa.primersArgs.getMaxTemplateMispriming() > Short.MAX_VALUE && !pa.isThermodynamicTemplateAlignment() ) {
			glob_err.append( "Value too large at tag PRIMER_MAX_TEMPLATE_MISPRIMING");
			return true;
		}

		if (pa.getPairMaxTemplateMispriming() > Short.MAX_VALUE && !pa.isThermodynamicTemplateAlignment()) {
			glob_err.append( "Value too large at tag PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING");
			return true;
		}

		if (pa.primersArgs.getMaxRepeatCompl() > Short.MAX_VALUE && !pa.isThermodynamicOligoAlignment() ) {
			glob_err.append( "Value too large at tag PRIMER_MAX_LIBRARY_MISPRIMING");
			return true;
		}

		if (pa.oligosArgs.getMaxRepeatCompl() > Short.MAX_VALUE && !pa.isThermodynamicOligoAlignment() ) {
			glob_err.append( "Value too large at tag PRIMER_INTERNAL_MAX_LIBRARY_MISHYB");
			return true;
		}

		if (pa.getPairRepeatCompl() > Short.MAX_VALUE && !pa.isThermodynamicOligoAlignment() ) {
			glob_err.append( "Value too large at tag PRIMER_PAIR_MAX_LIBRARY_MISPRIMING");
			return true;
		}

		if (pa.oligosArgs.getMaxTemplateMispriming() >= 0 && !pa.isThermodynamicTemplateAlignment())
			glob_err.append("PRIMER_INTERNAL_MAX_TEMPLATE_MISHYB is not supported");
		if (pa.oligosArgs.getMaxTemplateMisprimingTH() >= 0 && pa.isThermodynamicTemplateAlignment())
			glob_err.append("PRIMER_INTERNAL_MAX_TEMPLATE_MISHYB_TH is not supported");
		if (pa.primersArgs.getMinSize() < 1)
			glob_err.append( "PRIMER_MIN_SIZE must be >= 1");

		if (pa.primersArgs.getMaxSize() > MAX_PRIMER_LENGTH) {
			glob_err.append( "PRIMER_MAX_SIZE exceeds built-in maximum of " + MAX_PRIMER_LENGTH);
			return true;
		}

		if (pa.primersArgs.getOptSize() > pa.primersArgs.getMaxSize()) {
			glob_err.append( "PRIMER_{OPT,DEFAULT}_SIZE > PRIMER_MAX_SIZE");
			return true;
		}

		if (pa.primersArgs.getOptSize() < pa.primersArgs.getMinSize()) {
			glob_err.append("PRIMER_{OPT,DEFAULT}_SIZE < PRIMER_MIN_SIZE");
			return true;
		}

		if (pa.oligosArgs.getMaxSize() > MAX_PRIMER_LENGTH) {
			glob_err.append("PRIMER_INTERNAL_MAX_SIZE exceeds built-in maximum");
			return true;
		}

		if (pa.oligosArgs.getOptSize() > pa.oligosArgs.getMaxSize()) {
			glob_err.append("PRIMER_INTERNAL_{OPT,DEFAULT}_SIZE > MAX_SIZE");
			return true;
		}

		if (pa.oligosArgs.getOptSize() < pa.oligosArgs.getMinSize()) {
			glob_err.append("PRIMER_INTERNAL_{OPT,DEFAULT}_SIZE < MIN_SIZE");
			return true;
		}

		/* A GC clamp can not be bigger then the primer */
		if (pa.getGcClamp() > pa.primersArgs.getMinSize()) {
			glob_err.append("PRIMER_GC_CLAMP > PRIMER_MIN_SIZE");
			return true;
		}

		/* Must be >= 0 and <= 5 */
		if ((pa.getMaxEndGC() < 0) 
				|| (pa.getMaxEndGC() > 5)) {
			glob_err.append("PRIMER_MAX_END_GC must be between 0 to 5");
			return true;
		}

		/* Product size must be provided */
		if (0 == pa.getProductSizeRangesNumber()) {
			glob_err.append("Empty value for PRIMER_PRODUCT_SIZE_RANGE");
			return true;
		}
		for (i = 0; i < pa.getProductSizeRangesNumber() ; i++) {
			if (pa.getProductSizeRange(i).getLeft() > pa.getProductSizeRange(i).getRight() || pa.getProductSizeRange(i).getLeft() < 0) {
				glob_err.append( "Illegal element in PRIMER_PRODUCT_SIZE_RANGE");
				return true;
			}
		}

		pr_min = Integer.MAX_VALUE;
		/* Check if the primer is bigger then the product */
		for(i=0;i<pa.getProductSizeRangesNumber();i++)
			if(pa.getProductSizeRange(i).getLeft() < pr_min) pr_min=pa.getProductSizeRange(i).getLeft();

		if (pa.primersArgs.getMaxSize() > pr_min) {
			glob_err.append("PRIMER_MAX_SIZE > min PRIMER_PRODUCT_SIZE_RANGE");
			return true;
		}

		if ((pa.isPickInternalOligo())
				&& pa.oligosArgs.getMaxSize() > pr_min) {
			glob_err.append( "PRIMER_INTERNAL_MAX_SIZE > min PRIMER_PRODUCT_SIZE_RANGE");
			return true;
		}

		/* There must be at least one primer returned */
		if (pa.getNumReturn() < 1) {
			glob_err.append("PRIMER_NUM_RETURN < 1");
			return true;
		}

		if ((pa.primersArgs.must_match_five_prime != null) && ((pa.primersArgs.must_match_five_prime.length) != 5)) {
			glob_err.append( "PRIMER_MUST_MATCH_FIVE_PRIME must have 5 characters");
			return true;
		}
		if ((pa.primersArgs.must_match_three_prime != null) && ((pa.primersArgs.must_match_three_prime.length) != 5)) {
			glob_err.append(
					"PRIMER_MUST_MATCH_THREE_PRIME must have 5 characters");
			return true;
		}
		if ((pa.oligosArgs.must_match_five_prime != null) && ((pa.oligosArgs.must_match_five_prime.length) != 5)) {
			glob_err.append(
					"PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME must have 5 characters");
			return true;
		}
		if ((pa.oligosArgs.must_match_three_prime != null) && ((pa.oligosArgs.must_match_three_prime.length) != 5)) {
			glob_err.append(
					"PRIMER_INTERNAL_MUST_MATCH_THREE_PRIME must have 5 characters");
			return true;
		}

		if (sa.getIncludedRegionLength() >= Integer.MAX_VALUE) {
			nonfatal_err.append( "Value for SEQUENCE_INCLUDED_REGION too large");
			return true;
		}

		if (sa.getIncludedRegionStart() < 0 || sa.getIncludedRegionLength() < 0
				|| sa.getIncludedRegionStart() + sa.getIncludedRegionLength() > seq_len) {
			nonfatal_err.append( "Illegal value for SEQUENCE_INCLUDED_REGION");
			return true;
		}

		/* The product must fit in the included region */
		if (sa.getIncludedRegionLength() < pr_min && pa.isPickLeftPrimer() 
				&& pa.isPickRightPrimer() ) {
			if (pa.getPrimerTask() == P3Task.CHECK_PRIMERS) {
				warning.append("SEQUENCE_INCLUDED_REGION length < min PRIMER_PRODUCT_SIZE_RANGE");
			} else if (pa.getPrimerTask() != P3Task.PICK_PRIMER_LIST) {
				nonfatal_err.append("SEQUENCE_INCLUDED_REGION length < min PRIMER_PRODUCT_SIZE_RANGE");
			}

			if (pa.getPrimerTask() == P3Task.GENERIC) {
				return true;
			}
		}

		if (pa.getMaxEndStability() < 0) {
			nonfatal_err.append(
					"PRIMER_MAX_END_STABILITY must be non-negative");
			return true;
		}

		/* Is the startodon ATG and in the incl. region */
		if (!sa.PR_START_CODON_POS_IS_NULL()) {
			if (!pa.isDefaultPositionPenalties()) {
				nonfatal_err.append( 	"Cannot accept both SEQUENCE_START_CODON_POSITION and non-default " + 
						"arguments for PRIMER_INSIDE_PENALTY or PRIMER_OUTSIDE_PENALTY");
			}
			if (sa.getStartCodonPos()  > (sa.getIncludedRegionStart() + sa.getIncludedRegionLength() - 3)) {
				nonfatal_err.append("Start codon position not contained in SEQUENCE_INCLUDED_REGION");
				return true;
			} else {
				if (sa.getStartCodonPos() >= 0
						&& ((sa.getSequence()[sa.getStartCodonPos()] != 'A'
						&& sa.getSequence()[sa.getStartCodonPos()] != 'a')
						|| (sa.getSequence()[sa.getStartCodonPos() + 1] != 'T'
						&& sa.getSequence()[sa.getStartCodonPos() + 1] != 't')
						|| (sa.getSequence()[sa.getStartCodonPos() + 2] != 'G'
						&& sa.getSequence()[sa.getStartCodonPos() + 2] != 'g'))) {
					nonfatal_err.append("No start codon at SEQUENCE_START_CODON_POSITION");
					return true;
				}
			}
		}

		if (null != sa.getSequenceQuality()) {
			if(pa.primersArgs.getMinQuality() != 0 && pa.primersArgs.getMinQuality() < pa.getQualityRangeMin()) {
				glob_err.append(
						"PRIMER_MIN_QUALITY < PRIMER_QUALITY_RANGE_MIN");
				return true;
			}
			if (pa.primersArgs.getMinQuality() != 0 && pa.primersArgs.getMinQuality() > pa.getQualityRangeMax()) {
				glob_err.append(
						"PRIMER_MIN_QUALITY > PRIMER_QUALITY_RANGE_MAX");
				return true;
			}
			if (pa.oligosArgs.getMinQuality() != 0 && pa.oligosArgs.getMinQuality() <pa.getQualityRangeMin()) {
				glob_err.append(
						"PRIMER_INTERNAL_MIN_QUALITY < PRIMER_QUALITY_RANGE_MIN");
				return true;
			}
			if (pa.oligosArgs.getMinQuality() != 0 && pa.oligosArgs.getMinQuality() > pa.getQualityRangeMax()) {
				glob_err.append(
						"PRIMER_INTERNAL_MIN_QUALITY > PRIMER_QUALITY_RANGE_MAX");
				return true;
			}
			for(i=0; i < sa.getSequenceQuality().length; i++) {
				if(sa.getSequenceQuality()[i] < pa.getQualityRangeMin() ||
						sa.getSequenceQuality()[i] > pa.getQualityRangeMax()) {
					nonfatal_err.append(
							"Sequence quality score out of range");
					return true;
				}
			}
		}
		else if (pa.primersArgs.weights.seq_quality != 0 || pa.oligosArgs.weights.seq_quality != 0) {
			nonfatal_err.append(
					"Sequence quality is part of objective function but sequence quality is not defined");
			return true;
		}
		offending_char = dna_to_upper(sa.getTrimmedSequence(), 0);
		if (offending_char != '\0' ) {
			if (pa.isLiberalBase()) {
				warning.append( "Unrecognized base in input sequence");
			}
			else {
				nonfatal_err.append(  "Unrecognized base in input sequence");
				return true;
			}
		}

		if (pa.primersArgs.getOptTm() < pa.primersArgs.getMinTm()
				|| pa.primersArgs.getOptTm() > pa.primersArgs.getMaxTm()) {
			glob_err.append("Optimum primer Tm lower than minimum or higher than maximum");
			return true;
		}

		if (pa.oligosArgs.getOptTm() < pa.oligosArgs.getMinTm()
				|| pa.oligosArgs.getOptTm() > pa.oligosArgs.getMaxTm()) {
			glob_err.append("Optimum internal oligo Tm lower than minimum or higher than maximum");
			return true;
		}

		if (pa.primersArgs.getMinGC() > pa.primersArgs.getMaxGC()
				|| pa.primersArgs.getMinGC() > 100
				|| pa.primersArgs.getMaxGC() < 0){
			glob_err.append("Illegal value for PRIMER_MAX_GC and PRIMER_MIN_GC");
			return true;
		}

		if (pa.oligosArgs.getMinGC() > pa.oligosArgs.getMaxGC()
				|| pa.oligosArgs.getMinGC() > 100
				|| pa.oligosArgs.getMaxGC() < 0) {
			glob_err.append(
					"Illegal value for PRIMER_INTERNAL_OLIGO_GC");
			return true;
		}
		if (pa.primersArgs.getMaxNumOfNsAccepted() < 0) {
			glob_err.append(
					"Illegal value for PRIMER_MAX_NS_ACCEPTED");
			return true;
		}
		if (pa.oligosArgs.getMaxNumOfNsAccepted() < 0){
			glob_err.append(
					"Illegal value for PRIMER_INTERNAL_MAX_NS_ACCEPTED");
			return true;
		}
		if (pa.primersArgs.getMaxSelfAny() < 0 || pa.primersArgs.getMaxSelfAny() > Short.MAX_VALUE
				|| pa.primersArgs.getMaxSelfEnd() < 0 || pa.primersArgs.getMaxSelfEnd() > Short.MAX_VALUE
				|| pa.getPairComplAny() < 0 || pa.getPairComplAny() > Short.MAX_VALUE 
				|| pa.getPairComplEnd() < 0 || pa.getPairComplEnd() > Short.MAX_VALUE) {
			glob_err.append( "Illegal value for primer complementarity restrictions");
			return true;
		}

		if (pa.primersArgs.getMaxSelfAnyTH() < 0
				|| pa.primersArgs.getMaxSelfEndTH() < 0 || pa.primersArgs.getMaxHairPinTH() < 0
				|| pa.getPairComplAnyTH() < 0 || pa.getPairComplEndTH() < 0) {
			glob_err.append(
					"Illegal value for primer complementarity restrictions (thermod. approach)");
			return true;
		}

		if( pa.oligosArgs.getMaxSelfAny() < 0 || pa.oligosArgs.getMaxSelfAny() > Short.MAX_VALUE
				|| pa.oligosArgs.getMaxSelfEnd() < 0 || pa.oligosArgs.getMaxSelfEnd() > Short.MAX_VALUE) {
			glob_err.append(
					"Illegal value for internal oligo complementarity restrictions");
			return true;
		}

		if( pa.oligosArgs.getMaxSelfAnyTH() < 0
				|| pa.oligosArgs.getMaxSelfEndTH() < 0 || pa.oligosArgs.getMaxHairPinTH() < 0) {
			glob_err.append(
					"Illegal value for internal oligo complementarity restrictions");
			return true;
		}

		if (pa.primersArgs.getSaltConcentration() <= 0 || pa.primersArgs.getDnaConcentration()<=0){
			glob_err.append(
					"Illegal value for primer salt or dna concentration");
			return true;
		}

		if((pa.primersArgs.getDntpConcentration() < 0 && pa.primersArgs.getDivalentConcentration() !=0 )
				|| pa.primersArgs.getDivalentConcentration()<0){ /* added by T.Koressaar */
			glob_err.append( "Illegal value for primer divalent salt or dNTP concentration");
			return true;
		}

		if(pa.oligosArgs.getSaltConcentration()<=0||pa.oligosArgs.getDnaConcentration()<=0){
			glob_err.append(
					"Illegal value for internal oligo salt or dna concentration");
			return true;
		}

		if((pa.oligosArgs.getDntpConcentration()<0 && pa.oligosArgs.getDivalentConcentration()!=0)
				|| pa.oligosArgs.getDivalentConcentration() < 0) { /* added by T.Koressaar */
			glob_err.append(
					"Illegal value for internal oligo divalent salt or dNTP concentration");
			return true;
		}

		if (!pa.isDefaultPositionPenalties() && sa.getTargetRegions().getCount() > 1) {
			nonfatal_err.append(
					"Non-default inside penalty or outside penalty ");
			nonfatal_err.append("is valid only when number of targets <= 1");
		}
		if (!pa.isDefaultPositionPenalties() && 0 == sa.getTargetRegions().getCount()) {
			warning.append(
					"Non-default inside penalty or outside penalty ");
			warning.append("has no effect when number of targets is 0");
		}
		if (pa.isPickInternalOligo() != true && sa.getInternalInput() != null) {
			nonfatal_err.append("Not specified to pick internal oligos");
			nonfatal_err.append(" but a specific internal oligo is provided");
		}
		if (sa.getInternalInput() != null) {
			if ((sa.getInternalInput().length) > MAX_PRIMER_LENGTH) {
				glob_err.append("Specified internal oligo exceeds built-in maximum of " + MAX_PRIMER_LENGTH);
				return true;
			}
			if ((sa.getInternalInput().length) >  pa.oligosArgs.getMaxSize())
				warning.append("Specified internal oligo > PRIMER_INTERNAL_MAX_SIZE");

			if ((sa.getInternalInput().length) <  pa.oligosArgs.getMinSize())
				warning.append("Specified internal oligo < PRIMER_INTERNAL_MIN_SIZE");

			if (!strstr_nocase(sa.getSequence(), sa.getInternalInput()))
				nonfatal_err.append( "Specified internal oligo not in sequence");
			else if (!strstr_nocase(sa.getTrimmedSequence(), sa.getInternalInput()))
				nonfatal_err.append( "Specified internal oligo not in Included Region");
		}
		if (sa.getLeftInput() != null) {
			if ((sa.getLeftInput().length) > MAX_PRIMER_LENGTH) {
				glob_err.append(
						"Specified left primer exceeds built-in maximum of " + MAX_PRIMER_LENGTH);
				return true;
			}
			if ((sa.getLeftInput().length) >  pa.primersArgs.getMaxSize())
				warning.append("Specified left primer > PRIMER_MAX_SIZE");
			if ((sa.getLeftInput().length) < pa.primersArgs.getMinSize())
				warning.append( "Specified left primer < PRIMER_MIN_SIZE");
			if (!strstr_nocase(sa.getSequence(), sa.getLeftInput()))
				nonfatal_err.append(
						"Specified left primer not in sequence");
			else if (!strstr_nocase(sa.getTrimmedSequence(), sa.getLeftInput()))
				nonfatal_err.append(
						"Specified left primer not in Included Region");
		}
		if (sa.getRightInput() != null) {
			if ((sa.getRightInput().length) > MAX_PRIMER_LENGTH) {
				glob_err.append(
						"Specified right primer exceeds built-in maximum of " + MAX_PRIMER_LENGTH);
				return true;
			}
			if ((sa.getRightInput().length) < pa.primersArgs.getMinSize())
				warning.append( "Specified right primer < PRIMER_MIN_SIZE");
			if ((sa.getRightInput().length) >  pa.primersArgs.getMaxSize()) {
				warning.append( "Specified right primer > PRIMER_MAX_SIZE");
			} else { /* We do not want to overflow s1. */
				char[] s1 =  Sequence.p3_reverse_complement(sa.getRightInput());
				if (!strstr_nocase(sa.getSequence(), s1))
					nonfatal_err.append(
							"Specified right primer not in sequence");
				else if (!strstr_nocase(sa.getTrimmedSequence(), s1))
					nonfatal_err.append( "Specified right primer not in Included Region");
			}
		}

		if ((pa.getPrPairWeights().product_tm_lt != 0 ||
				pa.getPrPairWeights().product_tm_gt != 0)
				&& pa.getProductOptTM() == PR_UNDEFINED_DBL_OPT) {
			glob_err.append("Product temperature is part of objective function while optimum temperature is not defined");
			return true;
		}

		if((pa.getPrPairWeights().product_size_lt != 0 ||
				pa.getPrPairWeights().product_size_gt != 0)
				&& pa.getProductOptSize() == PR_UNDEFINED_INT_OPT){
			glob_err.append( "Product size is part of objective function while optimum size is not defined");
			return true;
		}

		if ((pa.primersArgs.weights.gc_content_lt != 0 ||
				pa.primersArgs.weights.gc_content_gt != 0)
				&& pa.primersArgs.getOptGC() == DEFAULT_OPT_GC_PERCENT) {
			glob_err.append( "Primer GC content is part of objective function while optimum gc_content is not defined");
			return true;
		}

		if ((pa.oligosArgs.weights.gc_content_lt != 0 ||
				pa.oligosArgs.weights.gc_content_gt != 0 )
				&& pa.oligosArgs.getOptGC() == DEFAULT_OPT_GC_PERCENT) {
			glob_err.append("Hyb probe GC content is part of objective function while optimum gc_content is not defined");
			return true;
		}

		if ((pa.isPickInternalOligo() != true ) &&
				(pa.getPrPairWeights().io_quality > 0)) {
			glob_err.append( "Internal oligo quality is part of objective function while internal oligo choice is not required");
			return true;
		}

		if (pa.primersArgs.weights.repeat_sim > 0 && (seq_lib.seq_lib_num_seq(pa.primersArgs.repeat_lib) == 0)) {
			glob_err.append( "Mispriming score is part of objective function, but mispriming library is not defined");
			return true;
		}

		if (pa.oligosArgs.weights.repeat_sim > 0 && ( seq_lib.seq_lib_num_seq(pa.oligosArgs.repeat_lib) == 0)) {
			glob_err.append( "Internal oligo mispriming score is part of objective function while mishyb library is not defined");
			return true;
		}

		if (pa.getPrPairWeights().repeat_sim > 0 && ((seq_lib.seq_lib_num_seq(pa.primersArgs.repeat_lib) == 0))) {
			glob_err.append("Mispriming score is part of objective function, but mispriming library is not defined");
			return true;
		}

		if(pa.getPrPairWeights().io_quality > 0
				&& pa.isPickInternalOligo() == false ) {
			glob_err.append( "Internal oligo quality is part of objective function while internal oligo choice is not required");
			return true;
		}

		if (pa.getSequencingParameters().getLead() < 0) {
			glob_err.append(
					"Illegal value for PRIMER_SEQUENCING_LEAD");
			return true;
		}

		if (pa.getSequencingParameters().getInterval() < 0) {
			glob_err.append(
					"Illegal value for PRIMER_SEQUENCING_INTERVAL");
			return true;
		}

		if (pa.getSequencingParameters().getAccuracy() < 0) {
			glob_err.append(
					"Illegal value for PRIMER_SEQUENCING_ACCURACY");
			return true;
		}

		if (pa.getSequencingParameters().getSpacing() < 0) {
			glob_err.append(
					"Illegal value for PRIMER_SEQUENCING_SPACING");
			return true;
		}

		if(pa.getSequencingParameters().getInterval() > pa.getSequencingParameters().getSpacing()) {
			glob_err.append(
					"PRIMER_SEQUENCING_INTERVAL > PRIMER_SEQUENCING_SPACING");
			return true;
		}

		if(pa.getSequencingParameters().getAccuracy() > pa.getSequencingParameters().getSpacing()) {
			glob_err.append(
					"PRIMER_SEQUENCING_ACCURACY > PRIMER_SEQUENCING_SPACING");
			return true;
		}

		if(pa.getSequencingParameters().getLead() > pa.getSequencingParameters().getSpacing()) {
			glob_err.append(
					"PRIMER_SEQUENCING_LEAD > PRIMER_SEQUENCING_SPACING");
			return true;
		}

		if((sa.getForceLeftStart() > -1) && (sa.getForceLeftEnd() > -1)
				&& (sa.getForceLeftStart() > sa.getForceLeftEnd())) {
			glob_err.append(
					"SEQUENCE_FORCE_LEFT_START > SEQUENCE_FORCE_LEFT_END");
			return true;
		}

		if((sa.getForceRightEnd() > -1) && (sa.getForceRightStart() > -1)
				&& (sa.getForceRightEnd() > sa.getForceRightStart())) {
			glob_err.append(
					"SEQUENCE_FORCE_RIGHT_END > SEQUENCE_FORCE_RIGHT_START");
			return true;
		}

		if (pa.getMin5PrimeOverlapOfJunction() < 1) {
			glob_err.append("Illegal value for PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION");
			return true;
		}

		if (pa.getMin3PrimeOverlapOfJunction() < 1) {
			glob_err.append("Illegal value for PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION");
			return true;
		}

		if ((sa.getPrimerOverlapJunctionsList().size() > 0) &&
				(pa.getMin5PrimeOverlapOfJunction() > (pa.primersArgs.getMaxSize() / 2))) {
			glob_err.append("PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION > PRIMER_MAX_SIZE / 2");
			return true;
		}

		if ((sa.getPrimerOverlapJunctionsList().size() > 0) && 
				(pa.getMin3PrimeOverlapOfJunction() > (pa.primersArgs.getMaxSize() / 2))) {
			glob_err.append("PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION > PRIMER_MAX_SIZE / 2");
			return true;
		}

		if (pa.primersArgs.getDivalentConcentration() > 0.0 && pa.primersArgs.getDntpConcentration() <= 0.0) {
			warning.append( "PRIMER_SALT_DIVALENT > 0.0 but PRIMER_DNTP_CONC <= 0.0; use reasonable value for PRIMER_DNTP_CONC");
		}

		if ((pa.primersArgs.must_match_five_prime != null) &&
				(test_must_match_parameters(pa.primersArgs.must_match_five_prime))) {
			glob_err.append("Illegal values for PRIMER_MUST_MATCH_FIVE_PRIME");
			return true;
		}

		if ((pa.primersArgs.must_match_three_prime != null) &&
				(test_must_match_parameters(pa.primersArgs.must_match_three_prime))) {
			glob_err.append( "Illegal values for PRIMER_MUST_MATCH_THREE_PRIME");
			return true;
		}

		if ((pa.oligosArgs.must_match_five_prime != null) &&
				(test_must_match_parameters(pa.oligosArgs.must_match_five_prime))) {
			glob_err.append( "Illegal values for PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME");
			return true;
		}

		if ((pa.oligosArgs.must_match_three_prime != null) &&
				(test_must_match_parameters(pa.oligosArgs.must_match_three_prime))) {
			glob_err.append("Illegal values for PRIMER_INTERNAL_MUST_MATCH_THREE_PRIME");
			return true;
		}


		return (0 == nonfatal_err.length() && 0 == glob_err.length() ) ? false : true;

	}






	// TODO FIXME Find a better way to do this 
	private static boolean strstr_nocase(char[] s1, char[] s2) {

		if(s1 == null || s2 == null) return false;


		// s2 should be a substring of s
		String s1Str = new String(s1).toLowerCase();
		String s2Str = new String(s2).toLowerCase();
		if(s1Str.contains(s2Str))
		{
			return true;
		}
		else
			return false;
	}

	/**
	 *  return false if both are equal 
	 */
	static int strcmp_nocase(char[] s1 , char[] s2){
		if(s1 == null || s2 == null) return 1;
		if(s1.length != s2.length) return 1;

		for(int i =0;i< s1.length;i++)
		{
			if(Character.toLowerCase(s1[i]) != Character.toLowerCase(s2[i])  )
			{
				return 1;
			}
		}

		return 0;

	}

	private static boolean test_must_match_parameters(
			char[] seq) {


		if(seq.length != 5 )
			return true;

		for(int i = 0 ; i < seq.length ; i++ )
		{
			seq[i] =  Character.toUpperCase(seq[i]);
			char x = seq[i];
			//			Character.isAlphabetic(X);
			/* Check it is NACTGRYWSMKBHDV */
			if (	(x == 'N') ||
					(x == 'A') || (x == 'C') ||
					(x == 'T') || (x == 'G') ||
					(x == 'R') || (x == 'Y') ||
					(x == 'W') || (x == 'S') ||
					(x == 'M') || (x == 'K') ||
					(x == 'B') || (x == 'H') ||
					(x == 'D') || (x == 'V'))
			{
				// do nothing
			}
			else
			{
				return true;
			}
		}

		return false;
	}



	static char dna_to_upper(char[] seq, int ambiguity_code_ok) {

		char unrecognized_base = '\0';
		for(int i = 0; i < seq.length ; i++)
		{
			switch (seq[i]) {
			case 'a': case 'A': seq[i]='A'; break;
			case 'c': case 'C': seq[i]='C'; break;
			case 'g': case 'G': seq[i]='G'; break;
			case 't': case 'T': seq[i]='T'; break;
			case 'n': case 'N': seq[i]='N'; break;
			default:
				if (ambiguity_code_ok == 1) {
					switch (seq[i]) {
					case 'r': case 'R': seq[i] = 'R'; break;
					case 'y': case 'Y': seq[i] = 'Y'; break;
					case 'm': case 'M': seq[i] = 'M'; break;
					case 'w': case 'W': seq[i] = 'W'; break;
					case 's': case 'S': seq[i] = 'S'; break;
					case 'k': case 'K': seq[i] = 'K'; break;
					case 'd': case 'D': seq[i] = 'D'; break;
					case 'h': case 'H': seq[i] = 'H'; break;
					case 'v': case 'V': seq[i] = 'V'; break;
					case 'b': case 'B': seq[i] = 'B'; break;
					}
				} else {
					if (unrecognized_base == '\0' ) unrecognized_base = seq[i];
					seq[i] = 'N';
				}
				break;
			}

		}
		return unrecognized_base;
	}



	private static void p3_print_args(P3GlobalSettings pa, SeqArgs sa) {

		pa.p3_print_args();
		sa.p3_print_args();
	}



	public static String libprimer3_release() {
		// TODO Auto-generated method stub
		return null;
	}



	public static void p3_set_program_name(String pr_program_name) {
		// TODO Auto-generated method stub

	}



	public static void p3_set_gs_primer_file_flag(P3GlobalSettings global_pa,
			int file_flag) {
		// TODO Auto-generated method stub

	}



	public static double align_thermod(char[] s1, char[] s2,
			ThermodynamicAlignmentArguments a) throws ThermodynamicAlignmentException {
		
		ThermodynamicAlignmentResult r = null;
		ThermodynamicAlignment thal = new ThermodynamicAlignment(s1, s2, a);
		try {
			r = thal.thAlign();

			if (thal_trace != 0) {
				System.out.format(  "thal, thal_args, type=%s maxLoop=%d mv=%f dv=%f dntp=%f dna_conc=%f, temp=%f, temponly=%b dimer=%d\n",
						a.getAlignmentType(), a.getMaxLoop(), a.getMonovalentConc(), a.getDivalentConc(), a.getDntpConc(), a.getDnaConc(), 
						a.getTemperature(), a.isTempOnly(), a.getCalcDimer());
				System.out.format(  "thal: s1=%s s2=%s temp=%f msg=%s end1=%d end2=%d\n", 
					string(s1), string(s2), r.temp, r.msg, r.align_end_1, r.align_end_2);
			}

			// TODO :: add this to calc_thal
			//			PR_ASSERT(r.temp <= DBL_MAX);

			if(r.temp == thal.THAL_ERROR_SCORE )
				throw new ThermodynamicAlignmentException(r.msg);


			return ((r.temp < 0.0) ? 0.0 : (double) (r.temp));

		} catch (ThermodynamicAlignmentException e) {
			String msg = e.getMessage();
			if(r != null )
				msg +=  " : " + r.msg;
			e.printStackTrace();
			throw new ThermodynamicAlignmentException(msg);
		}


	}



	static String string(char[] s2) {
		if(s2 == null)
			return "";
		return new String(s2);
	}


	public static double align(char[] s1, char[] s2, DPAlignmentArgs a) throws AlignmentException {


		DPAlignmentResults r =  null;

		if(a.flag == DPAlignment.DPAL_LOCAL || a.flag == DPAlignment.DPAL_LOCAL_END) {
			if (s2.length < 3) {
				/* For extremely short alignments we simply
		         max out the score, because the dpal subroutines
		         for these cannot handle this case.
		         TO DO: this can probably be corrected in dpal. */
				return s2.length;
			}
		}
		//		  dpal((const unsigned char *) s1, (const unsigned char *) s2, a, &r);
		r = DPAlignment.dpAlign(s1, s2, a);
		//		  PR_ASSERT(r.score <= Short.MAX_VALUE);
		if (r.score == DPAlignment.DPAL_ERROR_SCORE) {
			throw new AlignmentException(r.msg)	;
		}
		return ((r.score < 0.0) ? 0.0 : r.score / LibPrimer3.PR_ALIGN_SCORE_PRECISION);
	}




	public static void primer_mispriming_to_template(PrimerRecord h,
			P3GlobalSettings pa,
			SeqArgs sa,
			OligoType l,
			OligoStats ostats,
			int first,
			int last,
			char[] s,
			char[] s_r,
			DPAlignmentArgs align_args) throws AlignmentException {

		char[] oseq;
		char[] target, target_r;
		int tmp, seqlen;
		boolean debug = false;
		int first_untrimmed, last_untrimmed;
		/* Indexes of first and last bases of the oligo in sa.seq,
		                     that is, WITHIN THE TOTAL SEQUENCE INPUT. */

		/* first, last are indexes of first and last bases of the oligo in
		     sa.trimmed_seq, that is, WITHIN THE INCLUDED REGION. */

		//		char   tmp_char;
		double tmp_score;

		seqlen = sa.getUpcasedSeq().length;
		first_untrimmed = sa.getIncludedRegionStart() + first;
		last_untrimmed = sa.getIncludedRegionStart() + last;


		if (l == OligoType.OT_LEFT) {
			oseq = s;
			target = sa.getUpcasedSeq();
			target_r = sa.getUpcasedSeqRev();
		} else {  /* l == OT_RIGHT */
			if (debug)
				System.err.format("first_untrimmed = %d, last_untrimmed = %d\n",
						first_untrimmed, last_untrimmed);
			oseq = s_r;
			target = sa.getUpcasedSeqRev();
			target_r = sa.getUpcasedSeq();
			/* We need to adjust first_untrimmed and last_untrimmed so that
			       they are correct in the reverse-complemented
			       sequence.
			 */
			tmp = (seqlen - last_untrimmed) - 1;
			last_untrimmed  = (seqlen - first_untrimmed) - 1;
			first_untrimmed = tmp;
		}

		/* 1. Align to the template 5' of the oligo. */
		//		tmp_char = target[first_untrimmed];
		//		target[first_untrimmed] = '\0';

		tmp_score = align(oseq, Sequence.subSeqRange(target, 0, first_untrimmed-1)  , align_args);

		if (debug) {
			if (l == OligoType.OT_LEFT) 
				System.err.format("\n************ OLIGO = LEFT\n");
			else  
				System.err.format("\n************ OLIGO = RIGHT\n");
			System.err.format("first_untrimmed = %d, last_untrimmed = %d,first = %d, last = %d\n",
					first_untrimmed, last_untrimmed, first, last);

			System.err.format("5' of oligo: Score %f aligning %s against %s\n\n", tmp_score,
					oseq, target);
		}

		//		target[first_untrimmed] = tmp_char;

		/* 2. Align to the template 3' of the oligo. */
		//		h.template_mispriming = align(oseq, target[0] + last_untrimmed + 1, align_args);
		h.template_mispriming = align(oseq, Sequence.subSeqRange(target,last_untrimmed + 1 , target.length-1), align_args);

		if (debug)
			System.err.format("3' of oligo Score %f aligning %s against %s\n\n",
					h.template_mispriming, oseq, Sequence.subSeqRange(target,last_untrimmed + 1 , target.length-1));

		/* 3. Take the max of 1. and 2. */
		if (tmp_score > h.template_mispriming)
			h.template_mispriming = tmp_score;

		/* 4. Align to the reverse strand of the template. */
		h.template_mispriming_r
		= align(oseq, target_r, align_args);

		if (debug)
			System.err.format("other strand Score %f aligning %s against %s\n\n",
					h.template_mispriming_r, oseq, target_r);

		if (pa.primersArgs.getMaxTemplateMispriming() >= 0) {
			if (h.oligo_max_template_mispriming()
					> pa.primersArgs.getMaxTemplateMispriming()) {
				h.op_set_high_similarity_to_multiple_template_sites();
				if (l == OligoType.OT_LEFT  || l == OligoType.OT_RIGHT  ) {
					ostats.template_mispriming++;
					ostats.ok--;
				} else {
					//		    	  PR_ASSERT(0); 
					/* Should not get here. */
				}
			} else {
				/* This oligo is ok, mark it so we do not do this check again. */
				h.template_mispriming_ok = 1;
			}
		}

	}



	public static void primer_mispriming_to_template_thermod(
			PrimerRecord h,
			P3GlobalSettings pa,
			SeqArgs sa,
			OligoType l,
			OligoStats ostats,
			int first,
			int last,
			char[] s,
			char[] s_r,
			ThermodynamicAlignmentArguments align_args) throws ThermodynamicAlignmentException {

		char[] oseq;
		char[] target, target_r;
		int tmp, seqlen;
		boolean debug = false;
		int first_untrimmed, last_untrimmed;
		/* Indexes of first and last bases of the oligo in sa.seq,
		 *                      that is, WITHIN THE TOTAL SEQUENCE INPUT. */

		/* first, last are indexes of first and last bases of the oligo in
		 *      sa.trimmed_seq, that is, WITHIN THE INCLUDED REGION. */

		//		   char   tmp_char;
		double  tmp_score;

		seqlen = sa.getUpcasedSeq().length;
		first_untrimmed = sa.getIncludedRegionStart() + first;
		last_untrimmed = sa.getIncludedRegionStart() + last;

		if (l == OligoType.OT_RIGHT) {
			oseq = s_r;
			target = sa.getUpcasedSeq();
			target_r = sa.getUpcasedSeqRev();
		} else { /* l == OT_RIGHT */
			if (debug)
				System.err.format( "first_untrimmed = %d, last_untrimmed = %d\n",
						first_untrimmed, last_untrimmed);
			oseq = s;
			target = sa.getUpcasedSeqRev();
			target_r = sa.getUpcasedSeq();
			/* We need to adjust first_untrimmed and last_untrimmed so that
			 * they are correct in the reverse-complemented
			 * sequence.
			 *     */
			tmp = (seqlen - last_untrimmed) - 1;
			last_untrimmed  = (seqlen - first_untrimmed) - 1;
			first_untrimmed = tmp;
		}

		/* 1. Align to the template 5' of the oligo. */
		/* tmp_char = target_r[first_untrimmed]; Corrected 2012-05-30 */
		//		   tmp_char = target[first_untrimmed];
		//		   target[first_untrimmed] = '\0';

		tmp_score = align_thermod(oseq, Sequence.subSeqRange(target,0, first_untrimmed-1), align_args);

		if (debug) {
			if (l == OligoType.OT_LEFT) 
				System.err.format( "\n************ OLIGO = LEFT\n");
			else 
				System.err.format( "\n************ OLIGO = RIGHT\n");
			System.err.format( "first_untrimmed = %d, last_untrimmed = %d, first = %d, last = %d\n",
					first_untrimmed, last_untrimmed, first, last);

			System.err.format( "5' of oligo: Score %f aligning %s against %s\n\n", tmp_score,
					oseq, target);
		}
		//		   target[first_untrimmed] = tmp_char;

		/* 2. Align to the template 3' of the oligo. &target[0] + last_untrimmed + 1 */ 
		h.template_mispriming
		= align_thermod(oseq, Sequence.subSeqRange(target,last_untrimmed + 1 , target.length-1) , align_args);

		if (debug)
			System.err.format("3' of oligo Score %f aligning %s against %s\n\n",
					h.template_mispriming, oseq, Sequence.subSeqRange(target,last_untrimmed + 1 , target.length-1));

		/* 3. Take the max of 1. and 2. */
		if (tmp_score > h.template_mispriming)
			h.template_mispriming = tmp_score;

		/* 4. Align to the reverse strand of the template. */
		h.template_mispriming_r
		= align_thermod(oseq, target_r, align_args);

		if (debug)
			System.err.format( "In primer_mispriming_to_template_thermod\n other strand Score %f aligning %s against %s\n\n",
					h.template_mispriming_r, oseq, target_r);

		if (pa.primersArgs.getMaxTemplateMisprimingTH() >= 0) {
			if (h.oligo_max_template_mispriming_thermod()
					> pa.primersArgs.getMaxTemplateMisprimingTH()) {
				h.op_set_high_similarity_to_multiple_template_sites();
				if (l == OligoType.OT_LEFT|| l == OligoType.OT_RIGHT  ) {
					ostats.template_mispriming++;
					ostats.ok--;
				}
				else {
					//		    	   PR_ASSERT(0); 
					/* Should not get here. */
				}
			} else {
				/* This oligo is ok, mark it so we do not do this check again. */
				h.template_mispriming_ok = 1;
			}
		}

	}



	public static void p3_print_oligo_lists(P3RetVal retval, SeqArgs sarg,
			P3GlobalSettings global_pa, StringBuilder per_sequence_err,
			String sequence_name) {
		// TODO Auto-generated method stub

	}

}
