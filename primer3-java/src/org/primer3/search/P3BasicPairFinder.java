package org.primer3.search;

import java.util.HashMap;

import org.primer3.libprimer3.DPAlArgHolder;
import org.primer3.libprimer3.DebugInfo;
import org.primer3.libprimer3.LibPrimer3;
import org.primer3.libprimer3.P3RetVal;
import org.primer3.libprimer3.P3Task;
import org.primer3.libprimer3.PairArrayT;
import org.primer3.libprimer3.PairStats;
import org.primer3.libprimer3.PrimerPair;
import org.primer3.libprimer3.PrimerRecord;
import org.primer3.libprimer3.THAlArgHolder;

public class P3BasicPairFinder extends Primer3Finder {

	public P3BasicPairFinder(P3RetVal retval, DPAlArgHolder dpal_arg_to_use, THAlArgHolder thal_arg_to_use,
			THAlArgHolder thal_oligo_arg_to_use) {
		super(retval, dpal_arg_to_use, thal_arg_to_use, thal_oligo_arg_to_use);
		// TODO Auto-generated constructor stub
	}

	@Override
	public void getLocalNextResult() throws Exception {
		// TODO Auto-generated method stub
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
				if (retval.pa.getPrPairWeights().primer_quality *
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
						if (DebugInfo.choose_pair_or_triple_trace_me)
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
					if (retval.pa.getPrPairWeights().primer_quality *
							(left.quality + right.quality) 
							> the_best_pair.pair_quality) {
						break;
					}

					/* Need to have this here because if we break just above, then,
			           at a later iteration, we may need to examine the oligo
			           pair with reverse oligo at i and forward oligo at j. */
					update_stats = false;
					if (j > max_j_seen[i]) {
						if (DebugInfo.choose_pair_or_triple_trace_me)
							System.err.format("updates ON: i=%d, j=%d, max_j_seen[%d]=%d\n",  i, j, i, max_j_seen[i]);
						max_j_seen[i] = j;
						if (DebugInfo.choose_pair_or_triple_trace_me)
							System.err.format("max_j_seen[%d] -. %d\n", i, max_j_seen[i]);
						if (DebugInfo.choose_pair_or_triple_trace_me) System.err.format( "updates on\n");
						update_stats = true;
					}

					if (left.overlaps) {
						/* The stats will not keep track of the pair correctly
			             after the first pass, because an oligo might
			             have been legal on one pass but become illegal on
			             a subsequent pass. */
						if (update_stats) {
							if (DebugInfo.choose_pair_or_triple_trace_me)
								System.err.format("i=%d, j=%d, overlaps_oligo_in_better_pair++\n", i, j);
							pair_expl.overlaps_oligo_in_better_pair++;
						}
						continue;
					}

					/* Some simple checks first, before searching the hashmap */
					boolean must_use = false;
					if ((retval.pa.getPrimerTask() == P3Task.CHECK_PRIMERS) || 
							((left.must_use != false) &&
									(right.must_use != false))) {
						must_use = true;
					}

					/* Determine if overlap with an overlap point is required, and
				   if so, whether one of the primers in the pairs overlaps
				   that point. */
					if ((retval.sa.getPrimerOverlapJunctionsList().size() > 0)
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
					if (product_size <  retval.pa.getProductSizeRange(product_size_range_index).getLeft()  ||
							product_size > retval.pa.getProductSizeRange(product_size_range_index).getRight()) {
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
									if (DebugInfo.choose_pair_or_triple_trace_me)
										System.err.format("ok++\n");
									pair_expl.ok++;
								}
								/* Check if this is a better pair */
								if (LibPrimer3.compare_primer_pair(pp, the_best_pair) < 0) {
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
						int tmp =  h.characterize_pair(retval,
//								retval.pa, retval.sa,
								left, right,
								dpal_arg_to_use,
								thal_arg_to_use,
								update_stats);
						if (tmp == PrimerPair.PAIR_OK) {

							/* Choose internal oligo if needed */
							if (retval.pa.isPickRightPrimer() && retval.pa.isPickLeftPrimer()
									&& retval.pa.isPickInternalOligo()) {
								n_int = LibPrimer3.choose_internal_oligo(retval, h.left, h.right, retval.sa, retval.pa, dpal_arg_to_use, thal_oligo_arg_to_use);
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
								if (DebugInfo.choose_pair_or_triple_trace_me)
									System.err.format("ok++\n");
								pair_expl.ok++;
							}

							/* Calculate the pair penalty */
							h.pair_quality = h.obj_fn(retval.pa);
							//							PR_ASSERT(h.pair_quality >= 0.0);

							/* Save the pair */
//							pp = new primer_pair();

							pp = h;
							hmap.put(j,pp);

							/* The current pair (h) is the new best pair if it is better than the best pair so far. */
							if (LibPrimer3.compare_primer_pair(h, the_best_pair) < 0) {
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

				if (!(product_size_range_index < retval.pa.getProductSizeRangesNumber())) {
					/* We ran out of product-size-ranges. Exit the while loop. */
					break;

					/* Our bookkeeping was incorrect unless the assertion below is
				   true. If num_intervals > 1 or min_three_*_prime_distance >
				   -1 the assertion below might not be true. */
					//			        PR_ASSERT(!((pa.num_intervals == 1) &&  ((pa.min_left_three_prime_distance == -1) || (pa.min_right_three_prime_distance == -1))) || (best_pairs.num_pairs == pair_expl.ok));
				}

			} else {
				/* Store the best primer for output */

				if (DebugInfo.choose_pair_or_triple_trace_me)
					System.err.format("ADD pair i=%d, j=%d\n", the_best_i, the_best_j);

				best_pairs.add_pair(the_best_pair);

				/* Mark the pair as already selected */
				//			      delete best_pp;
				best_hmap.put(the_best_j, null);

				/* Update the overlaps flags */
				for (int i = 0; i < retval.rev.num_elem; i++) {
					PrimerRecord right = retval.rev.oligo.get(i);
					if (LibPrimer3.right_oligo_in_pair_overlaps_used_oligo(right,
							the_best_pair,
							retval.pa.getMinRight3PrimeDistance())) {
						right.overlaps = true;
					}
				}
				for (int j = 0; j < retval.fwd.num_elem; j++) {
					PrimerRecord left = retval.fwd.oligo.get(j);
					if (LibPrimer3.left_oligo_in_pair_overlaps_used_oligo(left,
							the_best_pair,
							retval.pa.getMinLeft3PrimeDistance())) {
						left.overlaps = true;
					}
				}

				/* If we have enough then stop the while loop */
				if (retval.pa.getNumReturn() <= best_pairs.num_pairs) {
					break;
				}
			}

		}
	}

}
