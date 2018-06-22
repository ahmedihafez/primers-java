package org.primer3.libprimer3;

import org.primer3.dpal.AlignmentException;
import org.primer3.oligotm.OligoTMCalculator;
import org.primer3.p3_seq_lib.seq_lib;
import org.primer3.sequence.Sequence;
import org.primer3.thal.ThermodynamicAlignmentException;
import org.primer3.thal.ThermodynamicAlignmentArguments;

public class PrimerPair {



	static public int PAIR_OK = 1;
	static public int PAIR_FAILED = 0;


	/**
	 *  Penalty value of the primer pair 
	 */
	double pair_quality;  

	/** Absolute value of the difference between melting
	 * 	temperatures for left and right primers. 
	 */
	double diff_tm;       

	/* Estimated melting temperature of the product. */
	public double product_tm;    
	
	/* 
	 * Difference in Tm between the primer with lowest Tm
	 * the product Tm.
	 */
	double product_tm_oligo_tm_diff;
	
	double t_opt_a;

	/**
	 * Local complementarity score between left and right
	 * primers (* 100).
	 */
	double compl_any;     

	/** 
	 * 3'-anchored global complementatory score between *
	 * left and right primers (* 100).
	 */
	double compl_end;     


	
	/**
	 *  Maximum total mispriming score of both primers
	 *  to ectopic sites in the template, on "same"
	 *  strand (* 100). 
	 */
	double template_mispriming;
	
	/** Maximum total similarity of both primers to the
	 *  sequence from given file in fasta format.
	 */
	double repeat_sim;    
	public PrimerRecord left;     /* Left primer. */
	public PrimerRecord right;    /* Right primer. */
	public PrimerRecord intl;     /* Internal oligo. */

	boolean   must_use;

	int    product_size;  /* Product size. */
	int    target;        /* 
	 * 1 if there is a target between the right and left
	 * primers.
	 */
	// char *
	String rep_name;


	/**
	 * Return the value of the objective function value for given primer pair.
	 * We must know that the pair is acceptable before calling obj_fn.
	 * @throws Exception 
	 */
	public double obj_fn(P3GlobalSettings pa) throws Exception {
		double sum;
		double lower_tm;

		sum = 0.0;
		lower_tm = this.right.temp;
		if(this.left.temp < this.right.temp) lower_tm = this.left.temp;

		if(pa.getPrPairWeights().primer_quality > 0)   /*  HERE 1 */
			sum += pa.getPrPairWeights().primer_quality * (this.left.quality + this.right.quality);

		if(pa.getPrPairWeights().io_quality > 0 && pa.isPickRightPrimer()
				&& pa.isPickLeftPrimer() && pa.isPickInternalOligo())
			sum += pa.getPrPairWeights().io_quality * this.intl.quality;

		if(pa.getPrPairWeights().diff_tm > 0)
			sum += pa.getPrPairWeights().diff_tm * this.diff_tm;

		if (!pa.isThermodynamicOligoAlignment()) {
			if (pa.getPrPairWeights().compl_any > 0)
				sum += pa.getPrPairWeights().compl_any * this.compl_any;
			if (pa.getPrPairWeights().compl_end > 0)
				sum += pa.getPrPairWeights().compl_end * this.compl_end;
		} else if (pa.isThermodynamicOligoAlignment() ) {

			if (pa.getPrPairWeights().compl_any_th > 0 && ((lower_tm - pa.getPrPairWeights().temp_cutoff) <= this.compl_any))
				sum += pa.getPrPairWeights().compl_any_th * (this.compl_any - (lower_tm - pa.getPrPairWeights().temp_cutoff - 1.0));
			if (pa.getPrPairWeights().compl_any_th  > 0 && ((lower_tm - pa.getPrPairWeights().temp_cutoff) > this.compl_any))
				sum += pa.getPrPairWeights().compl_any_th * (1/(lower_tm - pa.getPrPairWeights().temp_cutoff + 1.0 - this.compl_any));

			if (pa.getPrPairWeights().compl_end_th > 0 && ((lower_tm - pa.getPrPairWeights().temp_cutoff) <= this.compl_end))
				sum += pa.getPrPairWeights().compl_end_th * (this.compl_end - (lower_tm - pa.getPrPairWeights().temp_cutoff - 1.0));
			if (pa.getPrPairWeights().compl_end_th > 0 && ((lower_tm - pa.getPrPairWeights().temp_cutoff) > this.compl_end))
				sum += pa.getPrPairWeights().compl_end_th * (1/(lower_tm - pa.getPrPairWeights().temp_cutoff + 1.0 - this.compl_end));
		}

		if(pa.getPrPairWeights().product_tm_lt > 0 && this.product_tm < pa.getProductOptTM())
			sum += pa.getPrPairWeights().product_tm_lt *
			(pa.getProductOptTM() - this.product_tm);

		if(pa.getPrPairWeights().product_tm_gt  > 0 && this.product_tm > pa.getProductOptTM())
			sum += pa.getPrPairWeights().product_tm_gt *
			(this.product_tm - pa.getProductOptTM());

		if(pa.getPrPairWeights().product_size_lt > 0 &&
				this.product_size < pa.getProductOptSize())
			sum += pa.getPrPairWeights().product_size_lt *
			(pa.getProductOptSize() - this.product_size);

		if(pa.getPrPairWeights().product_size_gt > 0 &&
				this.product_size > pa.getProductOptSize())
			sum += pa.getPrPairWeights().product_size_gt *
			(this.product_size - pa.getProductOptSize());

		if(pa.getPrPairWeights().repeat_sim > 0)
			sum += pa.getPrPairWeights().repeat_sim * this.repeat_sim;

		if (pa.getPrPairWeights().template_mispriming > 0 && !pa.isThermodynamicTemplateAlignment()) {
			//			PR_ASSERT(pa.pr_pair_weights.template_mispriming >= 0.0);
			//			TODO :: PR_ASSERT(this.template_mispriming >= 0.0);
			sum += pa.getPrPairWeights().template_mispriming * this.template_mispriming;
		}
		if (pa.getPrPairWeights().template_mispriming_th > 0 && pa.isThermodynamicTemplateAlignment()) {
			//			PR_ASSERT(pa.pr_pair_weights.template_mispriming_th >= 0.0);
			//			TODO :: PR_ASSERT(this.template_mispriming >= 0.0);
			if((lower_tm - pa.getPrPairWeights().temp_cutoff) <= this.template_mispriming)
				sum += pa.getPrPairWeights().template_mispriming_th * 
				(this.template_mispriming - (lower_tm - pa.getPrPairWeights().temp_cutoff - 1.0));
			if((lower_tm - pa.getPrPairWeights().temp_cutoff) > this.template_mispriming)
				sum += pa.getPrPairWeights().template_mispriming_th *
				(1/(lower_tm - pa.getPrPairWeights().temp_cutoff + 1.0 - this.template_mispriming));
		}
		if( !(sum >= 0.0) ) 
			throw new Exception("calculation error in obj_fn@primer_pair");
		this.pair_quality = sum;
		return sum;

	}

	/**
	 * Return 1 if 'pair' spans any target (in sa->tar); otherwise return 0; An
	 * oligo pair spans a target, t, if the last base of the left primer is to
	 * left of the last base of t and the first base of the right primer is to
	 * the right of the first base of t.  Of course the primers must
	 * still be in a legal position with respect to each other.
	 */
	public boolean pair_spans_target(SeqArgs sa) {
		int i;
		int last_of_left = this.left.start + this.left.length - 1;
		int first_of_right = this.right.start - this.right.length + 1;
		int target_first, target_last;
		for (i = 0; i < sa.targetRegions.getCount(); i++) {
			target_first = sa.targetRegions.getInterval(i)[0];
			target_last = target_first + sa.targetRegions.getInterval(i)[1] - 1;
			if (last_of_left <= target_last
					&& first_of_right >= target_first
					&& last_of_left < first_of_right)
				return true;
		}
		return false;
	}

	/**
	 * Defines parameter values for given primer pair. Returns PAIR_OK if the pair is
	 * acceptable; PAIR_FAILED otherwise.  This function sets the
	 * various elements of *ppair, and also calculates some key EXPENSIVE
	 * parameters for individual primers if they have not been set yet.
	 * @throws ThermodynamicAlignmentException 
	 * @throws AlignmentException 
	 */
	public  int characterize_pair(P3RetVal retval,
			P3GlobalSettings pa,
			SeqArgs sa,
			int m,
			int n,
			int int_num,
			DPAlArgHolder dpal_arg_to_use,
			THAlArgHolder thal_arg_to_use,
			boolean update_stats) throws ThermodynamicAlignmentException, AlignmentException {
		
		PrimerPair ppair = this;

		char[] s1, s2,s1_rev, s2_rev;
		double compl_end;
		PairStats pair_expl = retval.best_pairs.expl;
		boolean must_use = false;
		double min_oligo_tm;
		int i;
		ThermodynamicAlignmentArguments thal_args_for_template_mispriming  = LibPrimer3.use_end_for_th_template_mispriming != 0  ? 
				thal_arg_to_use.end1 : thal_arg_to_use.any;
		ppair.left = retval.fwd.oligo.get( m);
		ppair.right = retval.rev.oligo.get(n);
		ppair.product_size = retval.rev.oligo.get(n).start - retval.fwd.oligo.get(m).start+1;

		ppair.target = 0;
		ppair.compl_any = ppair.compl_end = 0;
		if (update_stats) { 
			pair_expl.considered++;
		}

		if (pa.getPrimerTask() == P3Task.CHECK_PRIMERS) {
			must_use = true;
		}

		/* We must use the pair if the caller specifed
		     both the left and the right primer. */
		if (ppair.left.must_use != false && ppair.right.must_use != false) {

			/* But if caller specified primers are in
		       reversed order we cannot analyze them
		       as a pair. */
			if (ppair.product_size < 1) {
				pair_expl.reversed++;
				return PAIR_FAILED;
			}

			must_use = true;
		}

		ppair.must_use = must_use;

		if (sa.targetRegions.getCount() > 0) {
			if (ppair.pair_spans_target(sa)) {
				ppair.target = 1;
			} else {
				ppair.target = -1;
				if (update_stats) { pair_expl.target++; }
				if (!must_use) return PAIR_FAILED;
			}
		}

		/* ============================================================= */

		/* Check that the pair is included in one of the specified ok regions */
		if ((sa.ok_regions.getCount() > 0) && (!sa.ok_regions.any_pair)) {
			
			
			boolean included = sa.ok_regions.checkIncludedInAny(ppair);
			// TODO :: Check if it work
//			boolean included = false;
//			int l_start = ppair.left.start, l_end = ppair.left.start + ppair.left.length - 1;
//			int r_start = ppair.right.start - ppair.right.length + 1, r_end = ppair.right.start;
//			for (i=0; i<sa.ok_regions.count; i++) {
//				if (sa.ok_regions.left_pairs[i][0] == -1) {
//					/* any left primer, check right primer */
//					if ((r_start >= sa.ok_regions.right_pairs[i][0]) &&
//							(r_end <= sa.ok_regions.right_pairs[i][0] + sa.ok_regions.right_pairs[i][1] - 1)) {
//						included = true;
//						break;
//					}
//				} else if (sa.ok_regions.right_pairs[i][0] == -1) {
//					/* any right primer, check the left primer */
//					if ((l_start >= sa.ok_regions.left_pairs[i][0]) && 
//							(l_end <= sa.ok_regions.left_pairs[i][0] + sa.ok_regions.left_pairs[i][1] - 1)) {
//						included = true;
//						break;
//					}
//				} else {
//					/* check both primers */
//					if ((l_start >= sa.ok_regions.left_pairs[i][0]) && 
//							(l_end <= sa.ok_regions.left_pairs[i][0] + sa.ok_regions.left_pairs[i][1] - 1) &&
//							(r_start >= sa.ok_regions.right_pairs[i][0]) &&
//							(r_end <= sa.ok_regions.right_pairs[i][0] + sa.ok_regions.right_pairs[i][1] - 1)) {
//						included = true;
//						break;
//					}  
//				}
//			}
			
			
			if (!included) {
				if (update_stats) { pair_expl.not_in_any_ok_region++; }
				if (!must_use) return PAIR_FAILED;
			}
		}
		/* ============================================================= */
		/* Compute product Tm and related parameters; check constraints. */

		if (ppair.right.start - ppair.left.start + 1 <= 0) {
			System.err.print("temporary");
		}
		//		  PR_ASSERT(ppair.right.start - ppair.left.start + 1 > 0)
		ppair.product_tm =  OligoTMCalculator.longSeqTM(sa.trimmed_seq, ppair.left.start,
				ppair.right.start - ppair.left.start + 1,
				/* TO DO -- skewed, it would be better to not use p_args elements here */
				pa.primersArgs.getSaltConcentration(),
				pa.primersArgs.getDivalentConcentration(),
				pa.primersArgs.getDntpConcentration());

		//		  PR_ASSERT(ppair.product_tm != oligo_tm.OLIGOTM_ERROR);

		min_oligo_tm
		= ppair.left.temp > ppair.right.temp ? ppair.right.temp : ppair.left.temp;
		ppair.product_tm_oligo_tm_diff = ppair.product_tm - min_oligo_tm;
		ppair.t_opt_a  = 0.3 * min_oligo_tm + 0.7 * ppair.product_tm - 14.9;

		if (pa.getProductMinTM() != LibPrimer3.PR_DEFAULT_PRODUCT_MIN_TM
				&& ppair.product_tm < pa.getProductMinTM()) {
			if (update_stats) { pair_expl.low_tm++; }
			if (!must_use) return PAIR_FAILED;
		}

		if (pa.getProductMaxTM() != LibPrimer3.PR_DEFAULT_PRODUCT_MAX_TM
				&& ppair.product_tm > pa.getProductMaxTM()) {
			if (update_stats) { pair_expl.high_tm++; }
			if (!must_use) return PAIR_FAILED;
		}

		ppair.diff_tm =  Math.abs(retval.fwd.oligo.get(m).temp - retval.rev.oligo.get(n).temp);
		if (ppair.diff_tm > pa.getMaxDiffTm()) {
			if (update_stats) { pair_expl.temp_diff++; }
			if (!must_use) return PAIR_FAILED;
		}

		/* End of product-temperature related computations. */
		/* ============================================================= */

		/* ============================================================= */
		/* BEGIN "EXPENSIVE" computations on _individual_ primers
		     in the pair.  These have been postponed until
		     this point in the interests of efficiency.
		 */

		/* ============================================================= */
		/* Estimate secondary structure and primer-dimer of _individual_
		     primers if not already calculated. */

		/* s1 is the forward oligo. */
		s1 = Sequence._pr_substr(sa.trimmed_seq, retval.fwd.oligo.get(m).start, retval.fwd.oligo.get(m).length );

		/* s2 is the reverse oligo. */
		s2 = Sequence._pr_substr(sa.trimmed_seq, retval.rev.oligo.get(n).start - retval.rev.oligo.get(n).length + 1, retval.rev.oligo.get(n).length);

		s1_rev = Sequence.p3_reverse_complement(s1);
		s2_rev = Sequence.p3_reverse_complement(s2);


		if (retval.fwd.oligo.get(m).self_any == LibPrimer3.ALIGN_SCORE_UNDEF 
				&& !pa.isThermodynamicOligoAlignment()) {
			/* We have not yet computed the 'self_any' paramter,
			       which is an estimate of self primer-dimer and secondary
			       structure propensity. */
			retval.fwd.oligo.get(m).oligo_compl( pa.primersArgs, retval.fwd.expl, dpal_arg_to_use, s1, s1_rev);

			if (!retval.fwd.oligo.get(m).OK_OR_MUST_USE()) {
				pair_expl.considered--;
				if (!must_use) return PAIR_FAILED;
			}
		}
		/* Thermodynamic approach, fwd-primer */
		if (retval.fwd.oligo.get(m).self_any == LibPrimer3.ALIGN_SCORE_UNDEF 
				&& pa.isThermodynamicOligoAlignment()) {
			retval.fwd.oligo.get(m).oligo_compl_thermod(pa.primersArgs,
					retval.fwd.expl, thal_arg_to_use, s1, s1); /* ! s1, s1_rev */

			if (!retval.fwd.oligo.get(m).OK_OR_MUST_USE()) {
				pair_expl.considered--;
				if (!must_use) return PAIR_FAILED;
			}
		}   
		if (retval.fwd.oligo.get(m).hairpin_th == LibPrimer3.ALIGN_SCORE_UNDEF 
				&& pa.isThermodynamicOligoAlignment()) {
			retval.fwd.oligo.get(m).oligo_hairpin( pa.primersArgs,
					retval.fwd.expl, thal_arg_to_use, s1);
			if (!retval.fwd.oligo.get(m).OK_OR_MUST_USE()) {
				pair_expl.considered--;
				if (!must_use) return PAIR_FAILED;
			}
		}

		/* End of thermodynamic approach */
		
		
		
		
		  if (retval.rev.oligo.get(n).self_any == LibPrimer3.ALIGN_SCORE_UNDEF && !pa.isThermodynamicOligoAlignment()) {
			  retval.rev.oligo.get(n).oligo_compl( pa.primersArgs, retval.rev.expl, dpal_arg_to_use, s2_rev, s2);

			    if (!retval.rev.oligo.get(n).OK_OR_MUST_USE()) {
			      pair_expl.considered--;
			      if (!must_use) return PAIR_FAILED;
			    }
			  }
			   /* Thermodynamic approach */
			   if (retval.rev.oligo.get(n).self_any == LibPrimer3.ALIGN_SCORE_UNDEF && pa.isThermodynamicOligoAlignment()) {
				   retval.rev.oligo.get(n).oligo_compl_thermod( pa.primersArgs,
			                          retval.rev.expl, thal_arg_to_use, s2_rev, s2_rev); /* s2_rev, s2 */
			      
			      if (!retval.rev.oligo.get(n).OK_OR_MUST_USE()) {
			         pair_expl.considered--;
			         if (!must_use) return PAIR_FAILED;
			      }  
			   }
			   if (retval.rev.oligo.get(n).hairpin_th == LibPrimer3.ALIGN_SCORE_UNDEF && pa.isThermodynamicOligoAlignment()) {
				   retval.rev.oligo.get(n).oligo_hairpin( pa.primersArgs,
			                    retval.rev.expl, thal_arg_to_use, s2_rev);
			      if (!retval.rev.oligo.get(n).OK_OR_MUST_USE()) {
			         pair_expl.considered--;
			         if (!must_use) return PAIR_FAILED;
			      }
			   }
			   
			   /* End of thermodynamic approach */
			  /* End of secondary structure and primer-dimer of _individual_ 
			     primers. */
			  /* ============================================================= */


			  /* ============================================================= */
			  /* Mispriming of _individual_ primers to template and 
			     to repeat libraries. */

			  if (retval.fwd.oligo.get(m).repeat_sim.score == null) {
			    /* We have not yet checked the oligo against the repeat library. */
				  retval.fwd.oligo.get(m).oligo_repeat_library_mispriming( pa, sa, OligoType.OT_LEFT,
			                                    retval.fwd.expl,dpal_arg_to_use, retval.glob_err);
			    if (retval.fwd.oligo.get(m).OK_OR_MUST_USE()) {
			    	retval.fwd.oligo.get(m).oligo_template_mispriming( pa, sa, OligoType.OT_LEFT,
							retval.fwd.expl,
							dpal_arg_to_use.local_end,
							thal_args_for_template_mispriming);
			    }
			    if (!retval.fwd.oligo.get(m).OK_OR_MUST_USE()) {
			      pair_expl.considered--;
			      if (!must_use) return PAIR_FAILED;
			    }
			  }
			   
			  if (retval.rev.oligo.get(n).repeat_sim.score == null) {
				  retval.rev.oligo.get(n).oligo_repeat_library_mispriming(pa, sa, OligoType.OT_RIGHT,
			                                    retval.rev.expl, dpal_arg_to_use,retval.glob_err);
			    if (retval.rev.oligo.get(n).OK_OR_MUST_USE()) {
			    	retval.rev.oligo.get(n).oligo_template_mispriming( pa, sa, OligoType.OT_RIGHT, retval.rev.expl, 
							dpal_arg_to_use.local_end,
							thal_args_for_template_mispriming);
			    }
			    if (!retval.rev.oligo.get(n).OK_OR_MUST_USE()) {
			      pair_expl.considered--;
			      if (!must_use) return PAIR_FAILED;
			    }
			  }
			   
			  /* End of mispriming of _individual_ primers to template and
			     to repeat libraries. */
			  /* ============================================================= */


			  /* ============================================================= */
			  /*
			   * Similarity between s1 and s2 is equivalent to complementarity between
			   * s2's complement and s1.  (Both s1 and s2 are taken from the same strand.)
			   */
			   if(!pa.isThermodynamicOligoAlignment()) {
			      ppair.compl_any = LibPrimer3.align(s1, s2, dpal_arg_to_use.local);
			      if (ppair.compl_any > pa.getPairComplAny()) {
			         if (update_stats) { pair_expl.compl_any++; }
			         if (!must_use) return PAIR_FAILED;
			      }
			      
			      ppair.compl_end = LibPrimer3.align(s1, s2, dpal_arg_to_use.end);
			      if (ppair.compl_end > pa.getPairComplEnd()) {
			         if (update_stats) { pair_expl.compl_end++; }
			         if (!must_use) return PAIR_FAILED;
			      }
			   } else {
			      /* thermodynamical approach */
			      ppair.compl_any = LibPrimer3.align_thermod(s1, s2_rev, thal_arg_to_use.any);
			      if (ppair.compl_any > pa.getPairComplAnyTH()) {
			         if (update_stats) {
			            pair_expl.compl_any++; 
			         }
			         if (!must_use) return PAIR_FAILED;
			      }
			      ppair.compl_end = LibPrimer3.align_thermod(s1, s2_rev, thal_arg_to_use.end1);
			      compl_end        = LibPrimer3.align_thermod(s1, s2_rev, thal_arg_to_use.end2); /* Triinu Please check */
			      if (ppair.compl_end < compl_end) {
			         ppair.compl_end = compl_end;
			      }
			      if (ppair.compl_end > pa.getPairComplEndTH()) {
			         if (update_stats) {
			            pair_expl.compl_end++; 
			         }
			         if (!must_use) return PAIR_FAILED;
			      }
			   }
			   
			  /*
			   * It is conceivable (though unlikely) that
			   * align(s2_rev, s1_rev, end_args) > align(s1,s2,end_args).
			   */
			  if (!pa.isThermodynamicOligoAlignment() && (compl_end = LibPrimer3.align(s2_rev, s1_rev, dpal_arg_to_use.end))
			      > ppair.compl_end) {
			    if (compl_end > pa.primersArgs.getMaxSelfEnd()) {
			      if (update_stats) { pair_expl.compl_end++; }
			      if (!must_use) return PAIR_FAILED;
			    }
			    ppair.compl_end = compl_end;
			  }

			  if ((ppair.repeat_sim = ppair.pair_repeat_sim( pa))
			      > pa.getPairRepeatCompl()) {
			    if (update_stats) { pair_expl.repeat_sim++; }
			    if (!must_use) return PAIR_FAILED;
			  }
			   /* thermodynamic approach */
			   if (pa.isThermodynamicOligoAlignment() ) {
				   // I think this is better than prev inline check
				   // maybe the second call will return a greater value
				   double compl_end1 = LibPrimer3.align_thermod(s2, s1_rev, thal_arg_to_use.end1) ; // ) > ppair.compl_end || 
				   double compl_end2 = LibPrimer3.align_thermod(s2, s1_rev, thal_arg_to_use.end2) ; // ) > ppair.compl_end)
				   if(compl_end1 > ppair.compl_end ||  compl_end2 > ppair.compl_end)
				   { 
					   compl_end = compl_end1;
					   if(compl_end1 < compl_end2 )
					   {
						   compl_end = compl_end2;
					   }
					  
				      if (compl_end > pa.primersArgs.getMaxSelfEndTH()) {
				         if (update_stats) {
				            pair_expl.compl_end++; 
				         }
				         if (!must_use) return PAIR_FAILED;
				      }
				      ppair.compl_end = compl_end;
				   }
			   }
			   
			  /* ============================================================= */


			  /* ============================================================= */
			  /* Calculate _pair_ mispriming, if necessary. */
			 
			   if (!pa.isThermodynamicTemplateAlignment() ) {
			     if (!pa.needPairTemplateMispriming())
			        ppair.template_mispriming = LibPrimer3.ALIGN_SCORE_UNDEF;
			     else {
//			        PR_ASSERT(ppair.left.template_mispriming != LibPrimer3.ALIGN_SCORE_UNDEF);
//			        PR_ASSERT(ppair.left.template_mispriming_r != LibPrimer3.ALIGN_SCORE_UNDEF);
//			        PR_ASSERT(ppair.right.template_mispriming != LibPrimer3.ALIGN_SCORE_UNDEF);
//			        PR_ASSERT(ppair.right.template_mispriming_r != LibPrimer3.ALIGN_SCORE_UNDEF);
			        ppair.template_mispriming =
			           ppair.left.template_mispriming + ppair.right.template_mispriming_r;
			        if ((ppair.left.template_mispriming_r + ppair.right.template_mispriming)
			             > ppair.template_mispriming)
			           ppair.template_mispriming
			           = ppair.left.template_mispriming_r + ppair.right.template_mispriming;

			        if (pa.getPairMaxTemplateMispriming() >= 0.0
			             && ppair.template_mispriming > pa.getPairMaxTemplateMispriming()) {
			            if (update_stats) { pair_expl.template_mispriming++; }
			            if (!must_use) return PAIR_FAILED;
			        }
			     }
			   } else { /* thermodynamic approach */
			     if (!pa.needPairTemplateMisprimingTH())
			       ppair.template_mispriming = LibPrimer3.ALIGN_SCORE_UNDEF;
			     else {
//			       PR_ASSERT(ppair.left.template_mispriming != LibPrimer3.ALIGN_SCORE_UNDEF);
//			       PR_ASSERT(ppair.left.template_mispriming_r != LibPrimer3.ALIGN_SCORE_UNDEF);
//			       PR_ASSERT(ppair.right.template_mispriming != LibPrimer3.ALIGN_SCORE_UNDEF);
//			       PR_ASSERT(ppair.right.template_mispriming_r != LibPrimer3.ALIGN_SCORE_UNDEF);
			       ppair.template_mispriming =
			           ppair.left.template_mispriming + ppair.right.template_mispriming_r;
			       if ((ppair.left.template_mispriming_r + ppair.right.template_mispriming)
			             > ppair.template_mispriming)
			           ppair.template_mispriming
			           = ppair.left.template_mispriming_r + ppair.right.template_mispriming;

			       if (pa.getPairMaxTemplateMisprimingTH() > 0  && ppair.template_mispriming > pa.getPairMaxTemplateMisprimingTH()) {
			            if (update_stats) {
			               pair_expl.template_mispriming++;
			            }
			            if (!must_use) return PAIR_FAILED;
			       }
			     }
			   }
			      
			   /* End of calculating _pair_ mispriming if necessary. */
			   /* ============================================================= */

			  return PAIR_OK;
	}


	private double pair_repeat_sim(P3GlobalSettings pa) {
		int i, n, max, w;
		  PrimerRecord fw, rev;

		  fw = this.left;
		  rev = this.right;

		  max = 0;
		  n = seq_lib.seq_lib_num_seq(pa.primersArgs.repeat_lib);
		  if(n == 0) return 0;
		  this.rep_name =  pa.primersArgs.repeat_lib.names.get(0) ;
		  for (i = 0; i < n; i++) {
			  w = (int) (fw.repeat_sim.score[i] + rev.repeat_sim.score[i]);
			  if (w > max) {
				  max = w;
				  this.rep_name =  pa.primersArgs.repeat_lib.names.get(i) ;
		    }
		  }
		  return max;
	}


}