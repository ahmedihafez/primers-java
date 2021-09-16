/*
    This file is part of Primer3 porting to java (https://github.com/primer3-org/primer3)


	Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008
	Whitehead Institute for Biomedical Research, Steve Rozen
	(http://purl.com/STEVEROZEN/), Andreas Untergasser and Helen Skaletsky
	All rights reserved to Primer3 authors.

    Primer3 and the libprimer3 library are free software;
    you can redistribute them and/or modify them under the terms
    of the GNU General Public License as published by the Free
    Software Foundation; either version 2 of the License, or (at
    your option) any later version.

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this software (file gpl-2.0.txt in the source
    distribution); if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
package org.primer3.primer;

import org.primer3.boulder;
import org.primer3.dpal.AlignmentException;
import org.primer3.libprimer3.DPAlArgHolder;
import org.primer3.libprimer3.LibPrimer3;
import org.primer3.libprimer3.OligoType;
import org.primer3.libprimer3.P3GlobalSettings;
import org.primer3.libprimer3.P3OutputType;
import org.primer3.libprimer3.P3RetVal;
import org.primer3.libprimer3.P3Task;
import org.primer3.libprimer3.PairStats;
import org.primer3.libprimer3.SeqArgs;
import org.primer3.libprimer3.THAlArgHolder;
import org.primer3.oligotm.OligoTMCalculator;
import org.primer3.p3_seq_lib.seq_lib;
import org.primer3.sequence.Sequence;
import org.primer3.thal.ThermodynamicAlignmentException;
import org.primer3.thal.ThermodynamicAlignmentArguments;

public class PrimerPair {



	final static public int PAIR_OK = 1;
	final static public int PAIR_FAILED = 0;
	final static public int PAIR_UNCHARACTRIZED = -1;
	
	
	public int pairStatus = PAIR_UNCHARACTRIZED;
	
	
	
	/**
	 *  Penalty value of the primer pair 
	 */
	public double pair_quality;  

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
	public double product_tm_oligo_tm_diff;
	
	public double t_opt_a;

	/**
	 * Local complementarity score between left and right
	 * primers (* 100).
	 */
	public double compl_any;     

	/** 
	 * 3'-anchored global complementatory score between *
	 * left and right primers (* 100).
	 */
	public double compl_end;     


	
	/**
	 *  Maximum total mispriming score of both primers
	 *  to ectopic sites in the template, on "same"
	 *  strand (* 100). 
	 */
	public double template_mispriming;
	
	/** Maximum total similarity of both primers to the
	 *  sequence from given file in fasta format.
	 */
	public double repeat_sim;    
	public PrimerRecord left;     /* Left primer. */
	public PrimerRecord right;    /* Right primer. */
	public PrimerRecord intl;     /* Internal oligo. */

	public boolean   must_use;

	public int    product_size;  /* Product size. */
	public int    target;        /* 
	 * 1 if there is a target between the right and left
	 * primers.
	 */
	// char *
	public String rep_name;


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
		for (i = 0; i < sa.getTargetRegions().getCount(); i++) {
			target_first = sa.getTargetRegions().getInterval(i)[0];
			target_last = target_first + sa.getTargetRegions().getInterval(i)[1] - 1;
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
//			P3GlobalSettings pa,
//			SeqArgs sa,
			PrimerRecord left, // int m,
			PrimerRecord right, // int n,
			//int int_num,
			DPAlArgHolder dpal_arg_to_use,
			THAlArgHolder thal_arg_to_use,
			boolean update_stats) throws ThermodynamicAlignmentException, AlignmentException {
		
		// TODO :: check if left or right is not good then return fail without continue
		
//		PrimerPair ppair = this;

		P3GlobalSettings pa = retval.pa;
		SeqArgs sa= retval.sa;
		
		char[] s1, s2,s1_rev, s2_rev;
		double compl_end;
		PairStats pair_expl = retval.best_pairs.expl;
		boolean must_use = false;
		double min_oligo_tm;
		int i;
		ThermodynamicAlignmentArguments thal_args_for_template_mispriming  = LibPrimer3.use_end_for_th_template_mispriming != 0  ? 
				thal_arg_to_use.end1 : thal_arg_to_use.any;
		this.left = left;//retval.fwd.oligo.get(m);
		this.right = right;//retval.rev.oligo.get(n);
		this.product_size = right.start - left.start+1;

		this.target = 0;
		this.compl_any = this.compl_end = 0;
		if (update_stats) { 
			pair_expl.considered++;
		}

		if (pa.getPrimerTask() == P3Task.CHECK_PRIMERS) {
			must_use = true;
		}

		/* We must use the pair if the caller specifed
		     both the left and the right primer. */
		if (this.left.must_use != false && this.right.must_use != false) {

			/* But if caller specified primers are in
		       reversed order we cannot analyze them
		       as a pair. */
			if (this.product_size < 1) {
				pair_expl.reversed++;
				this.pairStatus = PAIR_FAILED;
				return PAIR_FAILED;
			}

			must_use = true;
		}

		this.must_use = must_use;

		if (sa.getTargetRegions().getCount() > 0) {
			if (this.pair_spans_target(sa)) {
				this.target = 1;
			} else {
				this.target = -1;
				if (update_stats) { pair_expl.target++; }
				if (!must_use) {
					this.pairStatus = PAIR_FAILED;
					return PAIR_FAILED;
				}
			}
		}

		/* ============================================================= */

		/* Check that the pair is included in one of the specified ok regions */
		if ((sa.getOkRegions().getCount() > 0) && (!sa.getOkRegions().any_pair)) {
			
			
			boolean included = sa.getOkRegions().checkIncludedInAny(this);
			// TODO :: Check if it work
//			boolean included = false;
//			int l_start = this.left.start, l_end = this.left.start + this.left.length - 1;
//			int r_start = this.right.start - this.right.length + 1, r_end = this.right.start;
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
				if (!must_use) {
					this.pairStatus = PAIR_FAILED;
					return PAIR_FAILED;
				}
			}
		}
		/* ============================================================= */
		/* Compute product Tm and related parameters; check constraints. */

		if (this.right.start - this.left.start + 1 <= 0) {
			System.err.print("temporary");
		}
		//		  PR_ASSERT(this.right.start - this.left.start + 1 > 0)
		this.product_tm =  OligoTMCalculator.longSeqTM(sa.getTrimmedSequence(), this.left.start,
				this.right.start - this.left.start + 1,
				/* TO DO -- skewed, it would be better to not use p_args elements here */
				pa.primersArgs.getSaltConcentration(),
				pa.primersArgs.getDivalentConcentration(),
				pa.primersArgs.getDntpConcentration());

		//		  PR_ASSERT(this.product_tm != oligo_tm.OLIGOTM_ERROR);

		min_oligo_tm
		= this.left.temp > this.right.temp ? this.right.temp : this.left.temp;
		this.product_tm_oligo_tm_diff = this.product_tm - min_oligo_tm;
		this.t_opt_a  = 0.3 * min_oligo_tm + 0.7 * this.product_tm - 14.9;

		if (pa.getProductMinTM() != LibPrimer3.PR_DEFAULT_PRODUCT_MIN_TM
				&& this.product_tm < pa.getProductMinTM()) {
			if (update_stats) { pair_expl.low_tm++; }
			if (!must_use) {
				this.pairStatus = PAIR_FAILED;
				return PAIR_FAILED;
			}
		}

		if (pa.getProductMaxTM() != LibPrimer3.PR_DEFAULT_PRODUCT_MAX_TM
				&& this.product_tm > pa.getProductMaxTM()) {
			if (update_stats) { pair_expl.high_tm++; }
			if (!must_use) {
				this.pairStatus = PAIR_FAILED;
				return PAIR_FAILED;
			}
		}

		this.diff_tm =  Math.abs(left.temp - right.temp);
		if (this.diff_tm > pa.getMaxDiffTm()) {
			if (update_stats) { pair_expl.temp_diff++; }
			if (!must_use) {
				this.pairStatus = PAIR_FAILED;
				return PAIR_FAILED;
			}
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
//		s1 = Sequence._pr_substr(sa.getTrimmedSequence(), left.start, left.length );
		s1 = this.left.getOligoSeq();
		/* s2 is the reverse oligo. */
//		s2 = Sequence._pr_substr(sa.getTrimmedSequence(), right.start - right.length + 1, right.length);
		s2 = this.right.getOligoSeq();
		
//		s1_rev = Sequence.p3_reverse_complement(s1);
//		s2_rev = Sequence.p3_reverse_complement(s2);

		s1_rev = this.left.getOligoRevSeq();
		s2_rev = this.right.getOligoRevSeq();
				

		if (left.self_any == LibPrimer3.ALIGN_SCORE_UNDEF 
				&& !pa.isThermodynamicOligoAlignment()) {
			/* We have not yet computed the 'self_any' paramter,
			       which is an estimate of self primer-dimer and secondary
			       structure propensity. */
			left.oligo_compl( pa.primersArgs, retval.fwd.expl, dpal_arg_to_use, s1, s1_rev);

			if (!left.OK_OR_MUST_USE()) {
				pair_expl.considered--;
				if (!must_use) {
					this.pairStatus = PAIR_FAILED;
					return PAIR_FAILED;
				}
			}
		}
		/* Thermodynamic approach, fwd-primer */
		if (left.self_any == LibPrimer3.ALIGN_SCORE_UNDEF 
				&& pa.isThermodynamicOligoAlignment()) {
			left.oligo_compl_thermod(pa.primersArgs,
					retval.fwd.expl, thal_arg_to_use, s1, s1); /* ! s1, s1_rev */

			if (!left.OK_OR_MUST_USE()) {
				pair_expl.considered--;
				if (!must_use) {
					this.pairStatus = PAIR_FAILED;
					return PAIR_FAILED;
				}
			}
		}   
		if (left.hairpin_th == LibPrimer3.ALIGN_SCORE_UNDEF 
				&& pa.isThermodynamicOligoAlignment()) {
			left.oligo_hairpin( pa.primersArgs,
					retval.fwd.expl, thal_arg_to_use, s1);
			if (!left.OK_OR_MUST_USE()) {
				pair_expl.considered--;
				if (!must_use) {
					this.pairStatus = PAIR_FAILED;
					return PAIR_FAILED;
				}
			}
		}

		/* End of thermodynamic approach */
		
		
		
		
		  if (right.self_any == LibPrimer3.ALIGN_SCORE_UNDEF && !pa.isThermodynamicOligoAlignment()) {
			  right.oligo_compl( pa.primersArgs, retval.rev.expl, dpal_arg_to_use, s2_rev, s2);

			    if (!right.OK_OR_MUST_USE()) {
			      pair_expl.considered--;
			      if (!must_use) {
						this.pairStatus = PAIR_FAILED;
						return PAIR_FAILED;
					}
			    }
			  }
			   /* Thermodynamic approach */
			   if (right.self_any == LibPrimer3.ALIGN_SCORE_UNDEF && pa.isThermodynamicOligoAlignment()) {
				   right.oligo_compl_thermod( pa.primersArgs,
			                          retval.rev.expl, thal_arg_to_use, s2_rev, s2_rev); /* s2_rev, s2 */
			      
			      if (!right.OK_OR_MUST_USE()) {
			         pair_expl.considered--;
			         if (!must_use) {
							this.pairStatus = PAIR_FAILED;
							return PAIR_FAILED;
						}
			      }  
			   }
			   if (right.hairpin_th == LibPrimer3.ALIGN_SCORE_UNDEF && pa.isThermodynamicOligoAlignment()) {
				   right.oligo_hairpin( pa.primersArgs,
			                    retval.rev.expl, thal_arg_to_use, s2_rev);
			      if (!right.OK_OR_MUST_USE()) {
			         pair_expl.considered--;
			         if (!must_use) {
							this.pairStatus = PAIR_FAILED;
							return PAIR_FAILED;
						}
			      }
			   }
			   
			   /* End of thermodynamic approach */
			  /* End of secondary structure and primer-dimer of _individual_ 
			     primers. */
			  /* ============================================================= */


			  /* ============================================================= */
			  /* Mispriming of _individual_ primers to template and 
			     to repeat libraries. */

			  if (left.repeat_sim.score.size() == 0) {
			    /* We have not yet checked the oligo against the repeat library. */
				  left.oligo_repeat_library_mispriming( pa, sa, OligoType.OT_LEFT,
			                                    retval.fwd.expl,dpal_arg_to_use, retval.glob_err);
			    if (left.OK_OR_MUST_USE()) {
			    	left.oligo_template_mispriming( pa, sa, OligoType.OT_LEFT,
							retval.fwd.expl,
							dpal_arg_to_use.local_end,
							thal_args_for_template_mispriming);
			    }
			    if (!left.OK_OR_MUST_USE()) {
			      pair_expl.considered--;
			      if (!must_use) {
						this.pairStatus = PAIR_FAILED;
						return PAIR_FAILED;
					}
			    }
			  }
			   
			  if (right.repeat_sim.score.size() == 0 ) {
				  right.oligo_repeat_library_mispriming(pa, sa, OligoType.OT_RIGHT,
			                                    retval.rev.expl, dpal_arg_to_use,retval.glob_err);
			    if (right.OK_OR_MUST_USE()) {
			    	right.oligo_template_mispriming( pa, sa, OligoType.OT_RIGHT, retval.rev.expl, 
							dpal_arg_to_use.local_end,
							thal_args_for_template_mispriming);
			    }
			    if (!right.OK_OR_MUST_USE()) {
			      pair_expl.considered--;
			      if (!must_use) {
						this.pairStatus = PAIR_FAILED;
						return PAIR_FAILED;
					}
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
			      this.compl_any = LibPrimer3.align(s1, s2, dpal_arg_to_use.local);
			      if (this.compl_any > pa.getPairComplAny()) {
			         if (update_stats) { pair_expl.compl_any++; }
			         if (!must_use) {
							this.pairStatus = PAIR_FAILED;
							return PAIR_FAILED;
						}
			      }
			      
			      this.compl_end = LibPrimer3.align(s1, s2, dpal_arg_to_use.end);
			      if (this.compl_end > pa.getPairComplEnd()) {
			         if (update_stats) { pair_expl.compl_end++; }
			         if (!must_use) {
							this.pairStatus = PAIR_FAILED;
							return PAIR_FAILED;
						}
			      }
			   } else {
			      /* thermodynamical approach */
			      this.compl_any = LibPrimer3.align_thermod(s1, s2_rev, thal_arg_to_use.any);
			      if (this.compl_any > pa.getPairComplAnyTH()) {
			         if (update_stats) {
			            pair_expl.compl_any++; 
			         }
			         if (!must_use) {
							this.pairStatus = PAIR_FAILED;
							return PAIR_FAILED;
						}
			      }
			      this.compl_end = LibPrimer3.align_thermod(s1, s2_rev, thal_arg_to_use.end1);
			      compl_end        = LibPrimer3.align_thermod(s1, s2_rev, thal_arg_to_use.end2); /* Triinu Please check */
			      if (this.compl_end < compl_end) {
			         this.compl_end = compl_end;
			      }
			      if (this.compl_end > pa.getPairComplEndTH()) {
			         if (update_stats) {
			            pair_expl.compl_end++; 
			         }
			         if (!must_use) {
							this.pairStatus = PAIR_FAILED;
							return PAIR_FAILED;
						}
			      }
			   }
			   
			  /*
			   * It is conceivable (though unlikely) that
			   * align(s2_rev, s1_rev, end_args) > align(s1,s2,end_args).
			   */
			  if (!pa.isThermodynamicOligoAlignment() && (compl_end = LibPrimer3.align(s2_rev, s1_rev, dpal_arg_to_use.end))
			      > this.compl_end) {
			    if (compl_end > pa.primersArgs.getMaxSelfEnd()) {
			      if (update_stats) { pair_expl.compl_end++; }
			      if (!must_use) {
						this.pairStatus = PAIR_FAILED;
						return PAIR_FAILED;
					}
			    }
			    this.compl_end = compl_end;
			  }

			  if ((this.repeat_sim = this.pair_repeat_sim( pa))
			      > pa.getPairRepeatCompl()) {
			    if (update_stats) { pair_expl.repeat_sim++; }
			    if (!must_use) {
					this.pairStatus = PAIR_FAILED;
					return PAIR_FAILED;
				}
			  }
			   /* thermodynamic approach */
			   if (pa.isThermodynamicOligoAlignment() ) {
				   // I think this is better than prev inline check
				   // maybe the second call will return a greater value
				   double compl_end1 = LibPrimer3.align_thermod(s2, s1_rev, thal_arg_to_use.end1) ; // ) > this.compl_end || 
				   double compl_end2 = LibPrimer3.align_thermod(s2, s1_rev, thal_arg_to_use.end2) ; // ) > this.compl_end)
				   if(compl_end1 > this.compl_end ||  compl_end2 > this.compl_end)
				   { 
					  compl_end =  Math.max(compl_end1, compl_end2);
				      if (compl_end > pa.primersArgs.getMaxSelfEndTH()) {
				         if (update_stats) {
				            pair_expl.compl_end++; 
				         }
				         if (!must_use) {
								this.pairStatus = PAIR_FAILED;
								return PAIR_FAILED;
							}
				      }
				      this.compl_end = compl_end;
				   }
			   }
			   
			  /* ============================================================= */


			  /* ============================================================= */
			  /* Calculate _pair_ mispriming, if necessary. */
			 
			   if (!pa.isThermodynamicTemplateAlignment() ) {
			     if (!pa.needPairTemplateMispriming())
			        this.template_mispriming = LibPrimer3.ALIGN_SCORE_UNDEF;
			     else {
//			        PR_ASSERT(this.left.template_mispriming != LibPrimer3.ALIGN_SCORE_UNDEF);
//			        PR_ASSERT(this.left.template_mispriming_r != LibPrimer3.ALIGN_SCORE_UNDEF);
//			        PR_ASSERT(this.right.template_mispriming != LibPrimer3.ALIGN_SCORE_UNDEF);
//			        PR_ASSERT(this.right.template_mispriming_r != LibPrimer3.ALIGN_SCORE_UNDEF);
			        this.template_mispriming =
			           this.left.template_mispriming + this.right.template_mispriming_r;
			        if ((this.left.template_mispriming_r + this.right.template_mispriming)
			             > this.template_mispriming)
			           this.template_mispriming
			           = this.left.template_mispriming_r + this.right.template_mispriming;

			        if (pa.getPairMaxTemplateMispriming() >= 0.0
			             && this.template_mispriming > pa.getPairMaxTemplateMispriming()) {
			            if (update_stats) { pair_expl.template_mispriming++; }
			            if (!must_use) return this.pairStatus = PAIR_FAILED;
			        }
			     }
			   } else { /* thermodynamic approach */
			     if (!pa.needPairTemplateMisprimingTH())
			       this.template_mispriming = LibPrimer3.ALIGN_SCORE_UNDEF;
			     else {
//			       PR_ASSERT(this.left.template_mispriming != LibPrimer3.ALIGN_SCORE_UNDEF);
//			       PR_ASSERT(this.left.template_mispriming_r != LibPrimer3.ALIGN_SCORE_UNDEF);
//			       PR_ASSERT(this.right.template_mispriming != LibPrimer3.ALIGN_SCORE_UNDEF);
//			       PR_ASSERT(this.right.template_mispriming_r != LibPrimer3.ALIGN_SCORE_UNDEF);
			       this.template_mispriming =
			           this.left.template_mispriming + this.right.template_mispriming_r;
			       if ((this.left.template_mispriming_r + this.right.template_mispriming)
			             > this.template_mispriming)
			           this.template_mispriming
			           = this.left.template_mispriming_r + this.right.template_mispriming;

			       if (pa.getPairMaxTemplateMisprimingTH() > 0  && this.template_mispriming > pa.getPairMaxTemplateMisprimingTH()) {
			            if (update_stats) {
			               pair_expl.template_mispriming++;
			            }
			            if (!must_use) return this.pairStatus = PAIR_FAILED;
			       }
			     }
			   }
			      
			   /* End of calculating _pair_ mispriming if necessary. */
			   /* ============================================================= */

			  return this.pairStatus = PAIR_OK;
	}


	private double pair_repeat_sim(P3GlobalSettings pa) {
		if (pa.primersArgs.repeat_lib == null)  return 0 ;
		int i, n;
		double max, w;
		PrimerRecord fw, rev;

		fw = this.left;
		rev = this.right;
		// TO
		if (fw.repeat_sim.score.size() != rev.repeat_sim.score.size() )
			System.err.println("Stop");
		
		max = 0;
		n = pa.primersArgs.repeat_lib.seq_lib_num_seq();
		if(n == 0) return 0;
		this.rep_name =  pa.primersArgs.repeat_lib.getName(0) ;
		for (i = 0; i < n; i++) {
			w = (double) (fw.repeat_sim.score.get(i) + rev.repeat_sim.score.get(i));
			if (w > max) {
				max = w;
				this.rep_name =  pa.primersArgs.repeat_lib.getName(i) ;
			}
		}
		return max;
	}

	
	public void print_boulder(P3RetVal retval , String i, int io_version) {

		/* A place to put a string containing all error messages */
//		StringBuilder combined_retval_err = null;

		/* A small spacer; WARNING this is a fixed size
	     buffer, but plenty bigger than
	     log(2^64, 10), the longest character
	     string that is needed for a 64 bit integer. */
		String suffix ;

		/* Pointers for the primer set just printing */
		PrimerRecord fwd, rev, intl;



		/* Switches for printing this primer */
		int go_fwd = 0;
		int go_rev = 0;
		int go_int = 0;

		/* The number of loop cycles */
//		int loop_max;

		/* That links to the included region */
		int  incl_s = retval.sa.getIncludedRegionStart();

		/* This deals with the renaming of the internal oligo */
		String new_oligo_name = "INTERNAL";
		String int_oligo = new_oligo_name;


		{
			/* We will print primer pairs or pairs plus internal oligos */
			/* Get pointers to the primer_rec's that we will print */
			fwd  = this.left;
			rev  = this.right;

			// potential null exp bug here 
			// intl = retval.best_pairs.pairs.get(i).intl;
			intl = null;
			/* Pairs must have fwd and rev primers */
			go_fwd = 1;
			go_rev = 1;
			/* Do hyb oligos have to be printed? */
			if (this.intl != null ) {
				go_int = 1;
				intl = this.intl;
			} else {
				go_int = 0;
			}


			/* Get the number for pimer counting in suffix[0] */
			suffix = "_"+ i;

			/* Print out the Pair Penalties */
//			if (retval.output_type == P3OutputType.primer_pairs) {
				System.out.format("PRIMER_PAIR%s_PENALTY=%f\n", suffix,
						this.pair_quality);
//			}

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
						LibPrimer3.	string(fwd.getOligoSeq()));
			if (go_rev == 1)
				System.out.format("PRIMER_RIGHT%s_SEQUENCE=%s\n", suffix,
						LibPrimer3.	string(rev.getOligoRevSeq()));
			if(go_int == 1)
				System.out.format("PRIMER_%s%s_SEQUENCE=%s\n", int_oligo, suffix,
					LibPrimer3.	string(intl.getOligoSeq()));

			/* Print primer start and length */
			if (go_fwd == 1)
				System.out.format("PRIMER_LEFT%s=%d,%d\n", suffix,
						fwd.start + incl_s + retval.pa.getFirstBaseIndex(),
						fwd.length);
			if (go_rev == 1)
				System.out.format("PRIMER_RIGHT%s=%d,%d\n", suffix,
						rev.start + incl_s + retval.pa.getFirstBaseIndex(),
						rev.length);
			if (go_int == 1)
				System.out.format("PRIMER_%s%s=%d,%d\n", int_oligo, suffix,
						intl.start + incl_s + retval.pa.getFirstBaseIndex(),
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
			if (go_fwd == 1 && retval.pa.isThermodynamicOligoAlignment()==false)
				System.out.format("PRIMER_LEFT%s_SELF_ANY=%.2f\n", suffix,
						fwd.self_any);
			if (go_rev == 1 && retval.pa.isThermodynamicOligoAlignment()==false)
				System.out.format("PRIMER_RIGHT%s_SELF_ANY=%.2f\n", suffix,
						rev.self_any);
			if (go_int == 1 && retval.pa.isThermodynamicOligoAlignment()==false)
				System.out.format("PRIMER_%s%s_SELF_ANY=%.2f\n", int_oligo, suffix,
						intl.self_any);
			if (go_int == 1 && retval.pa.isThermodynamicOligoAlignment()==true)
				System.out.format("PRIMER_%s%s_SELF_ANY_TH=%.2f\n", int_oligo, suffix,
						intl.self_any);
			/* Print primer self_any thermodynamical approach */
			if (go_fwd == 1 && retval.pa.isThermodynamicOligoAlignment()==true)
				System.out.format("PRIMER_LEFT%s_SELF_ANY_TH=%.2f\n", suffix,
						fwd.self_any);
			if (go_rev == 1 && retval.pa.isThermodynamicOligoAlignment()==true)
				System.out.format("PRIMER_RIGHT%s_SELF_ANY_TH=%.2f\n", suffix,
						rev.self_any);
			/* Print primer self_end*/
			if (go_fwd == 1 && retval.pa.isThermodynamicOligoAlignment()==false)
				System.out.format("PRIMER_LEFT%s_SELF_END=%.2f\n", suffix,
						fwd.self_end);
			if (go_rev == 1 && retval.pa.isThermodynamicOligoAlignment()==false)
				System.out.format("PRIMER_RIGHT%s_SELF_END=%.2f\n", suffix,
						rev.self_end);
			if (go_int == 1 && retval.pa.isThermodynamicOligoAlignment()==false)
				System.out.format("PRIMER_%s%s_SELF_END=%.2f\n", int_oligo, suffix,
						intl.self_end);
			if (go_int == 1 && retval.pa.isThermodynamicOligoAlignment()==true)
				System.out.format("PRIMER_%s%s_SELF_END_TH=%.2f\n", int_oligo, suffix,
						intl.self_end);
			/* Print primer self_end thermodynamical approach */
			if (go_fwd == 1 && retval.pa.isThermodynamicOligoAlignment()==true)
				System.out.format("PRIMER_LEFT%s_SELF_END_TH=%.2f\n", suffix,
						fwd.self_end);
			if (go_rev == 1 && retval.pa.isThermodynamicOligoAlignment()==true)
				System.out.format("PRIMER_RIGHT%s_SELF_END_TH=%.2f\n", suffix,
						rev.self_end);
			/* Print primer hairpin */
			if (go_fwd == 1 && retval.pa.isThermodynamicOligoAlignment()==true)
				System.out.format("PRIMER_LEFT%s_HAIRPIN_TH=%.2f\n", suffix,
						fwd.hairpin_th);
			if (go_rev == 1 && retval.pa.isThermodynamicOligoAlignment()==true)
				System.out.format("PRIMER_RIGHT%s_HAIRPIN_TH=%.2f\n", suffix,
						rev.hairpin_th);
			if (go_int == 1 && retval.pa.isThermodynamicOligoAlignment()==true)
				System.out.format("PRIMER_%s%s_HAIRPIN_TH=%.2f\n", int_oligo, suffix,
						intl.hairpin_th);
			/*Print out primer mispriming scores */
			if (retval.pa.primersArgs.repeat_lib != null) {
				if (go_fwd == 1)
					System.out.format("PRIMER_LEFT%s_LIBRARY_MISPRIMING=%.2f, %s\n", suffix,
							fwd.repeat_sim.score.get(fwd.repeat_sim.max),
							fwd.repeat_sim.name);
				if (go_rev == 1)
					System.out.format("PRIMER_RIGHT%s_LIBRARY_MISPRIMING=%.2f, %s\n", suffix,
							rev.repeat_sim.score.get(rev.repeat_sim.max),
							rev.repeat_sim.name);
//				if (retval.output_type == P3OutputType.primer_pairs)
					System.out.format("PRIMER_PAIR%s_LIBRARY_MISPRIMING=%.2f, %s\n", suffix,
							this.repeat_sim,
							this.rep_name);
			}

			/* Print out internal oligo mispriming scores */
			if (go_int == 1 && retval.pa.oligosArgs.repeat_lib != null)
				System.out.format("PRIMER_%s%s_LIBRARY_MISHYB=%.2f, %s\n", int_oligo, suffix,
						intl.repeat_sim.score.get(intl.repeat_sim.max),
						intl.repeat_sim.name);

			/* If a sequence quality was provided, print it*/
			if (null != retval.sa.getSequenceQuality()){
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
			if (!retval.pa.isDefaultPositionPenalties()
					|| !retval.sa.PR_START_CODON_POS_IS_NULL()){
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
			if ( (!retval.pa.isThermodynamicTemplateAlignment() ) && (go_fwd == 1) && 
					(fwd.oligo_max_template_mispriming() != LibPrimer3.ALIGN_SCORE_UNDEF))
				System.out.format("PRIMER_LEFT%s_TEMPLATE_MISPRIMING=%.4f\n", suffix,
						fwd.oligo_max_template_mispriming());
			if ( (!retval.pa.isThermodynamicTemplateAlignment() ) && (go_rev == 1) && 
					(rev.oligo_max_template_mispriming() != LibPrimer3.ALIGN_SCORE_UNDEF))
				System.out.format("PRIMER_RIGHT%s_TEMPLATE_MISPRIMING=%.4f\n", suffix,
						rev.oligo_max_template_mispriming());

			/* Print primer template mispriming, thermodynamical approach*/
			if ( (retval.pa.isThermodynamicTemplateAlignment()) && (go_fwd == 1) &&
					(fwd.oligo_max_template_mispriming_thermod() != LibPrimer3.ALIGN_SCORE_UNDEF)) {
				System.out.format("PRIMER_LEFT%s_TEMPLATE_MISPRIMING_TH=%.4f\n", suffix,
						fwd.oligo_max_template_mispriming_thermod());
			}

			if ( (retval.pa.isThermodynamicTemplateAlignment() ) && (go_rev == 1) &&
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
				if (go_int == 1 && null != retval.sa.getSequenceQuality()) /* FIX ME - Uptate the tests */
					System.out.format("PRIMER_%s%s_MIN_SEQ_QUALITY=%d\n", int_oligo,
							suffix, intl.seq_quality);
				/* Print pair comp_any */
				if(retval.pa.isThermodynamicOligoAlignment()==false)
					System.out.format("PRIMER_PAIR%s_COMPL_ANY=%.2f\n", suffix,
							this.compl_any);
				if(retval.pa.isThermodynamicOligoAlignment()==true)
					System.out.format("PRIMER_PAIR%s_COMPL_ANY_TH=%.2f\n", suffix,
							this.compl_any);
				/* Print pair comp_end */
				if(retval.pa.isThermodynamicOligoAlignment()==false)
					System.out.format("PRIMER_PAIR%s_COMPL_END=%.2f\n", suffix,
							this.compl_end);
				if(retval.pa.isThermodynamicOligoAlignment()==true)
					System.out.format("PRIMER_PAIR%s_COMPL_END_TH=%.2f\n", suffix,
							this.compl_end);
				/* Print product size */
				System.out.format("PRIMER_PAIR%s_PRODUCT_SIZE=%d\n", suffix,
						this.product_size);
				/* Print the product Tm if a Tm range is defined */
				if (retval.pa.getProductMaxTM() != LibPrimer3.PR_DEFAULT_PRODUCT_MAX_TM ||
						retval.pa.getProductMinTM() != LibPrimer3.PR_DEFAULT_PRODUCT_MIN_TM) {
					System.out.format("PRIMER_PAIR%s_PRODUCT_TM=%.4f\n", suffix,
							this.product_tm);

					System.out.format("PRIMER_PAIR%s_PRODUCT_TM_OLIGO_TM_DIFF=%.4f\n", suffix,
							this.product_tm_oligo_tm_diff);

					System.out.format("PRIMER_PAIR%s_T_OPT_A=%.4f\n", suffix,
							this.t_opt_a);
				}

				/* Print the primer pair template mispriming */
				if ((!retval.pa.isThermodynamicTemplateAlignment() ) && (this.template_mispriming != LibPrimer3.ALIGN_SCORE_UNDEF))
					System.out.format("PRIMER_PAIR%s_TEMPLATE_MISPRIMING=%.2f\n", suffix,
							this.template_mispriming);
				/* Print the primer pair template mispriming. Thermodynamic approach.  */
				if ((retval.pa.isThermodynamicTemplateAlignment() ) && (this.template_mispriming != LibPrimer3.ALIGN_SCORE_UNDEF))
					System.out.format("PRIMER_PAIR%s_TEMPLATE_MISPRIMING_TH=%.2f\n", suffix,
							this.template_mispriming);

			} /* End of print parameters of primer pairs */

		} /* End of the big loop printing all data */

		/* End the print with newline and flush all buffers */

	}
	

}