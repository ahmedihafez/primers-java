package org.primer3.search;

import java.util.Comparator;
import java.util.TreeSet;

import org.primer3.libprimer3.DPAlArgHolder;
import org.primer3.libprimer3.DebugInfo;
import org.primer3.libprimer3.LibPrimer3;
import org.primer3.libprimer3.P3GlobalSettings;
import org.primer3.libprimer3.P3RetVal;
import org.primer3.libprimer3.P3Task;
import org.primer3.libprimer3.PairArrayT;
import org.primer3.libprimer3.PairStats;
import org.primer3.libprimer3.SeqArgs;
import org.primer3.libprimer3.THAlArgHolder;
import org.primer3.primer.PrimerPair;
import org.primer3.primer.PrimerRecord;

// Stateful search of primers pairs 
public class P3OptimzedFinder extends Primer3Finder {

	public P3OptimzedFinder(P3RetVal retval, DPAlArgHolder dpal_arg_to_use, THAlArgHolder thal_arg_to_use,
			THAlArgHolder thal_oligo_arg_to_use) {
		super(retval, dpal_arg_to_use, thal_arg_to_use, thal_oligo_arg_to_use);
		
		
		max_j_seen = new int[retval.rev.num_elem];
		
		resetSearch();
		
	}
	
	
	public void resetSearch()
	{
		for (int i = 0; i < max_j_seen.length; i++) max_j_seen[i] = -1;
		
		// in reset do not clear the cache calcPairs
	}
	class ijPairs {
		int i ;
		int j ;
		double estimatedQuilty = Double.MAX_VALUE;
		public ijPairs(int i , int j , double q) {
			this.i = i;
			this.j = j;
			this.estimatedQuilty =q;
		}
		
//		public double quality = Double.MAX_VALUE;
		public PrimerPair pairRecord =  null;
	}
	
	// Objects used in the Stateful search 
	int[] max_j_seen;
	TreeSet<ijPairs> calcPairs = new TreeSet<ijPairs>(new Comparator<ijPairs>() {

		@Override
		public int compare(ijPairs o1, ijPairs o2) {
			int cRes2 =  LibPrimer3.compare_primer_pair(o1.pairRecord,o2.pairRecord);
			return cRes2;
		}
	});
	/** ============================================================ */
	/** BEGIN choose_pair_or_triple_optimized
	 * 
	 * change :
	 * Product checking of if it satisfy any length
	 * use A TreeSet to hash calculated pairs instead of the old array of hashmaps used 
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
	@Override
	public void getLocalNextResult()  throws Exception {
		
		
		// who will clear the result ??
		PairArrayT best_pairs = retval.best_pairs;


		int n_int; /* Index of the internal oligo */
		boolean update_stats = true;  /* Flag to indicate whether pair_stats
		                            should be updated. */
		PrimerPair h;             /* The current pair which is being evaluated. */
		PrimerPair the_best_pair = new PrimerPair(); /* The best pair is being "remembered". */
		PairStats pair_expl = retval.best_pairs.expl; /* For statistics */

		int the_best_i, the_best_j;



//		HashMap<Integer, PrimerPair>[] pairs = new HashMap[retval.rev.num_elem];
		
		

		
		// #############
		
		
		
		
		//#############
		
		if(retval.fwd.num_elem == 0 || retval.rev.num_elem == 0)
			return;
		
		while(true) {
			/* If we have enough then stop the while loop */
			if (retval.pa.getNumReturn() <= best_pairs.num_pairs) {
				break;
			}
			the_best_i = -1;
			the_best_j = -1;
			/* To start put penalty to the maximum */
			if(calcPairs.size() > 0 )
			{
				ijPairs pair = calcPairs.first();
				the_best_pair = pair.pairRecord;
				the_best_i = pair.i;
				the_best_j = pair.j;
			}
			else {
				the_best_pair = new PrimerPair();
				the_best_pair.pair_quality = Double.MAX_VALUE;
			}

			for (int i = 0; i < retval.rev.num_elem; i++) {
				// keep retval.rev.oligo.get(i) in right
				PrimerRecord right = retval.rev.oligo.get(i);
				/* Pairs[i] is NULL if there has never been an assignment to
				 pairs[i] because pairs was allocated by calloc, which
				 sets the allocated memory to 0. */

				/* Only use a primer that *might be* legal or that the caller
		         has provided and specified as "must use".  Primers are *NOT*
		         FULLY ASSESSED until the call to characterize_pair(), in
		         order to avoid expensive computations (mostly alignments)
		         unless necessary. */
				if (!right.OK_OR_MUST_USE()) {
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
					continue;
				}


				/* Loop over forward primers */
				for (int j= max_j_seen[i]+1; j<retval.fwd.num_elem; j++) {
					PrimerRecord left = retval.fwd.oligo.get(j);

					/* We check the reverse oligo again, because we may
			           have determined that it is "not ok", even though
			           (as a far as we knew), it was ok above. */
					if ( ! right.OK_OR_MUST_USE()) {
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
					if (j > max_j_seen[i]) {  // this is dead code
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
					int product_size  = right.start - left.start+1;

					
					//  pa.pr_min[product_size_range_index]  pa.pr_max[product_size_range_index]
					if (!retval.pa.chackProductSizeRanges( product_size)) {
						if (update_stats) {
							/* This line NEW */ 
							if (!must_use)
								pair_expl.considered++;
							pair_expl.product++; 
						}
						if (!must_use) continue;
					}


					{
						/* Characterize the pair. h is initialized by this call. */
						h = new PrimerPair();
						int pairStatus =  h.characterize_pair(retval, 
//								retval.pa, retval.sa, 
								left, right,
								dpal_arg_to_use,
								thal_arg_to_use,
								update_stats);
						if (pairStatus == PrimerPair.PAIR_OK) {

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
//									hmap.put( j, null);
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



							ijPairs newPairInfo = new ijPairs(i, j, retval.pa.getPrPairWeights().primer_quality * 
									(left.quality + right.quality));
//							newPairInfo.quality = h.pair_quality;
							newPairInfo.pairRecord = h;
							boolean resAdded = calcPairs.add(newPairInfo);

							
							/* The current pair (h) is the new best pair if it is better than the best pair so far. */
							if (LibPrimer3.compare_primer_pair(h, the_best_pair) < 0) {
								the_best_pair = h;
								the_best_i = i;
								the_best_j = j;

							}

							/* There cannot be a better pair */
							if (the_best_pair.pair_quality == 0) {
								break;
							}
						} 
 
					}
				}  /* for (j=0; j<retval.fwd.num_elem; j++) -- inner loop */
				/* Check if there cannot be a better pair than best found */
			      
				if (the_best_pair.pair_quality == 0) {
			        break;
				}

			} /* for (i = 0; i < retval.rev.num_elem; i++) --- outer loop */
			

			if (the_best_pair.pair_quality == Double.MAX_VALUE ) {
				// search exhausted  no more to explore
				break;
			} 
			else 
			{
				/* Store the best primer for output */

				if (DebugInfo.choose_pair_or_triple_trace_me)
					System.err.format("ADD pair i=%d, j=%d\n", the_best_i, the_best_j);
				
				// remove the first one from the set
				calcPairs.remove(calcPairs.first());
				
				best_pairs.add_pair(the_best_pair);

				/* Update the overlaps flags */   // do not update all -- just as needed ??
				int minRight3PrimeDistance = retval.pa.getMinRight3PrimeDistance();
				
				if(minRight3PrimeDistance > 0 )
				{
					for (int i = 0; i < retval.rev.num_elem; i++) {
						PrimerRecord right = retval.rev.oligo.get(i);
						if ( LibPrimer3.right_oligo_in_pair_overlaps_used_oligo(right,
								the_best_pair,
								minRight3PrimeDistance)) {
							right.overlaps = true;
						}
					}
				}
				else if (minRight3PrimeDistance == 0)
				{
					// mark current only so we do not select it again
					the_best_pair.right.overlaps = true;
				}
				int minLeft3PrimeDistance = retval.pa.getMinLeft3PrimeDistance();
				if (minLeft3PrimeDistance > 0 )
				{
					for (int j = 0; j < retval.fwd.num_elem; j++) {
						PrimerRecord left = retval.fwd.oligo.get(j);
						if (LibPrimer3.left_oligo_in_pair_overlaps_used_oligo(left,
								the_best_pair,
								minLeft3PrimeDistance)) {
							left.overlaps = true;
						}
					}
				}
				else if (minLeft3PrimeDistance == 0)
				{
					the_best_pair.left.overlaps = true;
				}
				
			}

		}
	}





	
	
	
}
