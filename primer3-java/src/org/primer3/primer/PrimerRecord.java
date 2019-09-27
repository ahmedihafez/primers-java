package org.primer3.primer;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import org.primer3.dpal.AlignmentException;
import org.primer3.dpal.DPAlignment;
import org.primer3.dpal.DPAlignmentArgs;
import org.primer3.dpal.DPAlignmentResults;
import org.primer3.libprimer3.DPAlArgHolder;
import org.primer3.libprimer3.LibPrimer3;
import org.primer3.libprimer3.OligoStats;
import org.primer3.libprimer3.OligoType;
import org.primer3.libprimer3.P3GlobalSettings;
import org.primer3.libprimer3.P3OutputType;
import org.primer3.libprimer3.P3RetVal;
import org.primer3.libprimer3.P3Task;
import org.primer3.libprimer3.PrimerRecordException;
import org.primer3.libprimer3.PrimersOligosArguments;
import org.primer3.libprimer3.RepSim;
import org.primer3.libprimer3.SeqArgs;
import org.primer3.libprimer3.THAlArgHolder;
import org.primer3.masker.masker;
import org.primer3.masker.oligo_pair;
import org.primer3.oligotm.OligoTMCalculator;
import org.primer3.p3_seq_lib.seq_lib;
import org.primer3.sequence.Sequence;
import org.primer3.thal.ThermodynamicAlignmentException;
import org.primer3.thal.ThermodynamicAlignmentArguments;

public class PrimerRecord  {
	// not used right now TODO :: use it as an option
 	static public boolean CACHE_SEQ = true;
	
	
 	int status = PrimerPair.PAIR_UNCHARACTRIZED;
 	
	/* This is 5'.3' on the template sequence: */
	protected char[] oligoSeq = null;
	// reverse of above
	protected char[] oligoRevSeq =  null;
	
	// Store the source sequence and input setting here
	public SeqArgs sa;
	
	public String getTargetName() {
		if(sa != null)
		return this.sa.getSequenceName();
		return null;
	}
	
	public void setOligoSeq() {
		
		if (rec_type == OligoType.OT_RIGHT )
		{
			oligoSeq = Sequence._pr_substr(sa.getTrimmedSequence(), this.start - this.length + 1, this.length);

		}
		else
		{
			oligoSeq =  Sequence._pr_substr(sa.getTrimmedSequence(), this.start, this.length );
		}
		oligoRevSeq = Sequence.p3_reverse_complement(oligoSeq);
		
	}
	
	
	
	/**
	 * OligoSeq AS 5' to 3' on the template sequence
	 * from the  forward strand in template sequence 
	 * @return
	 */
	public char[] getOligoSeq() {
		return oligoSeq;
	}
	
	/**
	 * rev complement
	 * this is 5-3 in the reverse template strand
	 * @return
	 */
	public char[] getOligoRevSeq() {
		return oligoRevSeq;
	}	
	
	
	public PrimerRecord(OligoType rec_type) {
		this.rec_type = rec_type;
	}
	
	
	// new add here and should be maintained
	public OligoType rec_type;

	public RepSim repeat_sim = new RepSim();
	/*
	 * Information on the best repeat library (mispriming library) match for
	 * this oligo (primer), plus additional scores.
	 */

	public double temp; /*
						 * The oligo melting temperature calculated for the
						 * primer.
						 */

	public double gc_content;

	public double position_penalty;
	/*
	 * Penalty for distance from "ideal" position as specified by inside_penalty
	 * and outside_penalty.
	 */

	public double quality; /* Part of objective function due to this primer. */

	public double end_stability;
	/* Delta G of disription of 5 3' bases. */

	public int start; /*
					 * Position of the 5'-most base within the primer WITH
					 * RESPECT TO THE seq_args FIELD trimmed_seq.
					 */

	public int seq_quality; /* Minimum quality score of bases included. */
	public int seq_end_quality; /* Minimum quality core of the 5 3' bases. */

	public double self_any; /* Self complementarity as local alignment * 100. */

	public double self_end; /* Self complementarity at 3' end * 100. */

	public double hairpin_th; /*
							 * hairpin, thermodynamical approach and calculated
							 * as any
							 */

	public double template_mispriming;
	/*
	 * Max 3' complementarity to any ectopic site in template on the given
	 * template strand.
	 */
	public double template_mispriming_r;
	/*
	 * Max 3' complementarity to any ectopic site in the template on the reverse
	 * complement of the given template strand.
	 */
	public int length; /* Length of the oligo. */
	public int num_ns; /* Number of Ns in the oligo. */

	public boolean must_use; /* Non-0 if the oligo must be used even if it is illegal. */
	public boolean overlaps; /*
					 * Non-0 if the oligo overlaps some oligo used in one of the
					 * best pairs.
					 */

	public long oligoProblems;
	public boolean overlaps_overlap_position;

	public char template_mispriming_ok; /*
								 * Non-0 if the oligo was checked for this
								 * already and it is ok.
								 */

	public double failure_rate; /* Primer failure rate due to non-specific priming */

	public void free_primer_repeat_sim_score() {

		if (this.repeat_sim.score != null) {
			this.repeat_sim.score = null;
		}
	}

	/* ============================================================ */
	/* START functions for setting and getting oligo problems */
	/* ============================================================ */

	public boolean OK_OR_MUST_USE() {
		return !p3_ol_has_any_problem() || this.must_use;
	}

	public void initialize_op() {
		this.oligoProblems = 0L; /* 0UL is the unsigned long zero */
	}

	/*
	 * We use bitfields to store all primer data. The idea is that a primer is
	 * uninitialized if all bits are 0. It was evaluated if OP_PARTIALLY_WRITTEN
	 * is true. If OP_COMPLETELY_WRITTEN is also true the primer was evaluated
	 * till the end - meaning if all OP_... are false, the primer is OK, if some
	 * are true primer3 was forced to use it (must_use).
	 */

	static long OP_PARTIALLY_WRITTEN = (1L << 0);
	static long OP_COMPLETELY_WRITTEN = (1L << 1);
	static long BF_OVERLAPS_TARGET = (1L << 2);
	static long BF_OVERLAPS_EXCL_REGION = (1L << 3);
	static long BF_INFINITE_POSITION_PENALTY = (1L << 4);
	/* Space for more bitfields */

	static long OP_NOT_IN_ANY_OK_REGION = (1L << 7);
	static long OP_TOO_MANY_NS = (1L << 8); /* 3prime problem */
	static long OP_OVERLAPS_TARGET = (1L << 9); /* 3prime problem */
	static long OP_HIGH_GC_CONTENT = (1L << 10);
	static long OP_LOW_GC_CONTENT = (1L << 11);
	static long OP_HIGH_TM = (1L << 12);
	static long OP_LOW_TM = (1L << 13);
	static long OP_OVERLAPS_EXCL_REGION = (1L << 14);/* 3prime problem */
	static long OP_HIGH_SELF_ANY = (1L << 15);/* 3prime problem */
	static long OP_HIGH_SELF_END = (1L << 16);
	static long OP_NO_GC_CLAMP = (1L << 17);/* 3prime problem */
	static long OP_HIGH_END_STABILITY = (1L << 18);/* 3prime problem */
	static long OP_HIGH_POLY_X = (1L << 19);/* 3prime problem */
	static long OP_LOW_SEQUENCE_QUALITY = (1L << 20);/* 3prime problem */
	static long OP_LOW_END_SEQUENCE_QUALITY = (1L << 21);/* 3prime problem */
	static long OP_HIGH_SIM_TO_NON_TEMPLATE_SEQ = (1L << 22);/* 3prime problem */
	static long OP_HIGH_SIM_TO_MULTI_TEMPLATE_SITES = (1L << 23);
	static long OP_OVERLAPS_MASKED_SEQ = (1L << 24);
	static long OP_TOO_LONG = (1L << 25);/* 3prime problem */
	static long OP_TOO_SHORT = (1L << 26);
	static long OP_DOES_NOT_AMPLIFY_ORF = (1L << 27);
	static long OP_TOO_MANY_GC_AT_END = (1L << 28);/* 3prime problem */
	static long OP_HIGH_HAIRPIN = (1L << 29);/* 3prime problem */
	static long OP_MUST_MATCH_ERR = (1L << 30);

	/* all bits 1 except bits 0 to 6 */
	/* (~0UL) ^ 127UL = 1111 1111 1111 1111 1111 1111 1000 0000 */
	static long any_problem = (~0L) ^ 127L;
	/* 310297344UL = 0001 0010 0111 1110 1100 0011 0000 0000 */
	static long five_prime_problem = 310297344L;

	boolean p3_ol_is_uninitialized() {
		return (this.oligoProblems == 0L);
	}

	boolean p3_ol_is_ok() {
		return (this.oligoProblems & OP_COMPLETELY_WRITTEN) != 0;
	}

	/**
	 * Return 1 iff the argument 'oligo' has any problems -- i.e. violations of
	 * design constraints.
	 */
	public boolean p3_ol_has_any_problem() {
		return (this.oligoProblems & any_problem) != 0;
	}

	public boolean any_5_prime_ol_extension_has_problem() {
		return (this.oligoProblems & five_prime_problem) != 0;
	}

	/**
	 * Return a string details the the problems in 'oligo', i.e. the constraints
	 * that 'oligo' violates. WARNING: Returns a pointer to static storage,
	 * which is over-written on next call. p3_primer_rec_problems_to_string
	 */
	public String p3_get_ol_problem_string() {
		long prob = this.oligoProblems;

		StringBuilder output = new StringBuilder();
		// (prob & OP_PARTIALLY_WRITTEN) && !(prob & OP_COMPLETELY_WRITTEN))
		if ((prob & OP_PARTIALLY_WRITTEN) != 0
				&& (prob & OP_COMPLETELY_WRITTEN) == 0) {
			output.append(" Not completely checked;");
		}

		if (0 != (prob & OP_TOO_MANY_NS))
			output.append(" Too many Ns;");
		if (0 != (prob & OP_OVERLAPS_TARGET))
			output.append(" Overlaps target;");
		if (0 != (prob & OP_HIGH_GC_CONTENT))
			output.append(" GC content too high;");
		if (0 != (prob & OP_LOW_GC_CONTENT))
			output.append(" GC content too low;");
		if (0 != (prob & OP_HIGH_TM))
			output.append(" Temperature too high;");
		if (0 != (prob & OP_LOW_TM))
			output.append(" Temperature too low;");
		if (0 != (prob & OP_OVERLAPS_EXCL_REGION))
			output.append(" Overlaps an excluded region;");
		if (0 != (prob & OP_NOT_IN_ANY_OK_REGION))
			output.append(" Not in any ok region;");
		if (0 != (prob & OP_HIGH_SELF_ANY))
			output.append(" Similarity to self too high;");
		if (0 != (prob & OP_HIGH_SELF_END))
			output.append(" Similary to 3' end of self too high;");
		if (0 != (prob & OP_HIGH_HAIRPIN))
			output.append(" Hairpin stability too high;");
		if (0 != (prob & OP_NO_GC_CLAMP))
			output.append(" No 3' GC clamp;");
		if (0 != (prob & OP_TOO_MANY_GC_AT_END))
			output.append(" Too many GCs at 3' end;");
		if (0 != (prob & OP_HIGH_END_STABILITY))
			output.append(" 3' end too stable (delta-G too high);");
		if (0 != (prob & OP_HIGH_POLY_X))
			output.append(" Contains too-long poly nucleotide tract;");
		if (0 != (prob & OP_LOW_SEQUENCE_QUALITY))
			output.append(" Template sequence quality too low;");
		if (0 != (prob & OP_LOW_END_SEQUENCE_QUALITY))
			output.append(" Template sequence quality at 3' end too low;");
		if (0 != (prob & OP_HIGH_SIM_TO_NON_TEMPLATE_SEQ))
			output.append(" Similarity to non-template sequence too high;");
		if (0 != (prob & OP_HIGH_SIM_TO_MULTI_TEMPLATE_SITES))
			output.append(" Similarity to multiple sites in template;");
		if (0 != (prob & OP_OVERLAPS_MASKED_SEQ))
			output.append(" 3' base overlaps masked sequence;");
		if (0 != (prob & OP_TOO_LONG))
			output.append(" Too long;");
		if (0 != (prob & OP_TOO_SHORT))
			output.append(" Too short;");
		if (0 != (prob & OP_DOES_NOT_AMPLIFY_ORF))
			output.append(" Would not amplify an open reading frame;");
		if (0 != (prob & OP_MUST_MATCH_ERR))
			output.append(" Failed must_match requirements;");
		return output.toString();
	}

	public void bf_set_overlaps_target(int val) {
		if (val == 0) {
			this.oligoProblems |= BF_OVERLAPS_TARGET;
			this.oligoProblems ^= BF_OVERLAPS_TARGET;
		} else {
			this.oligoProblems |= BF_OVERLAPS_TARGET;
		}
	}

	public boolean bf_get_overlaps_target() {
		return (this.oligoProblems & BF_OVERLAPS_TARGET) != 0;
	}

	public void bf_set_overlaps_excl_region(int val) {
		if (val == 0) {
			this.oligoProblems |= BF_OVERLAPS_EXCL_REGION;
			this.oligoProblems ^= BF_OVERLAPS_EXCL_REGION;
		} else {
			this.oligoProblems |= BF_OVERLAPS_EXCL_REGION;
		}
	}

	public boolean bf_get_overlaps_excl_region() {
		return (this.oligoProblems & BF_OVERLAPS_EXCL_REGION) != 0;
	}

	public void bf_set_infinite_pos_penalty(int val) {
		if (val == 0) {
			this.oligoProblems |= BF_INFINITE_POSITION_PENALTY;
			this.oligoProblems ^= BF_INFINITE_POSITION_PENALTY;
		} else {
			this.oligoProblems |= BF_INFINITE_POSITION_PENALTY;
		}
	}

	public boolean bf_get_infinite_pos_penalty() {
		return (this.oligoProblems & BF_INFINITE_POSITION_PENALTY) != 0;
	}

	public void op_set_does_not_amplify_orf() {
		this.oligoProblems |= OP_DOES_NOT_AMPLIFY_ORF;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	public void op_set_completely_written() {
		this.oligoProblems |= OP_COMPLETELY_WRITTEN;
	}

	public void op_set_must_match_err() {
		this.oligoProblems |= OP_MUST_MATCH_ERR;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	public void op_set_too_many_ns() {
		this.oligoProblems |= OP_TOO_MANY_NS;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	public void op_set_overlaps_target() {
		this.oligoProblems |= OP_OVERLAPS_TARGET;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	public void op_set_high_gc_content() {
		this.oligoProblems |= OP_HIGH_GC_CONTENT;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	public void op_set_low_gc_content() {
		this.oligoProblems |= OP_LOW_GC_CONTENT;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	public void op_set_high_tm() {
		this.oligoProblems |= OP_HIGH_TM;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	public void op_set_low_tm() {
		this.oligoProblems |= OP_LOW_TM;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	public void op_set_overlaps_excluded_region() {
		this.oligoProblems |= OP_OVERLAPS_EXCL_REGION;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	public void op_set_not_in_any_ok_region() {
		this.oligoProblems |= OP_NOT_IN_ANY_OK_REGION;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	public void op_set_high_self_any() {
		this.oligoProblems |= OP_HIGH_SELF_ANY;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	public void op_set_high_self_end() {
		this.oligoProblems |= OP_HIGH_SELF_END;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	public void op_set_high_hairpin_th() {
		this.oligoProblems |= OP_HIGH_HAIRPIN;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	public void op_set_no_gc_glamp() {
		this.oligoProblems |= OP_NO_GC_CLAMP;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	public void op_set_too_many_gc_at_end() {
		this.oligoProblems |= OP_TOO_MANY_GC_AT_END;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	/* Must not be called on a hybridization probe / internal oligo */
	public void op_set_high_end_stability() {
		this.oligoProblems |= OP_HIGH_END_STABILITY;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	public void op_set_high_poly_x() {
		this.oligoProblems |= OP_HIGH_POLY_X;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	public void op_set_low_sequence_quality() {
		this.oligoProblems |= OP_LOW_SEQUENCE_QUALITY;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	public void op_set_low_end_sequence_quality() {
		this.oligoProblems |= OP_LOW_END_SEQUENCE_QUALITY;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	public void op_set_high_similarity_to_non_template_seq() {
		this.oligoProblems |= OP_HIGH_SIM_TO_NON_TEMPLATE_SEQ;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	public void op_set_high_similarity_to_multiple_template_sites() {
		this.oligoProblems |= OP_HIGH_SIM_TO_MULTI_TEMPLATE_SITES;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	public void op_set_overlaps_masked_sequence() {
		this.oligoProblems |= OP_OVERLAPS_MASKED_SEQ;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	public void op_set_too_long() {
		this.oligoProblems |= OP_TOO_LONG;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	public void op_set_too_short() {
		this.oligoProblems |= OP_TOO_SHORT;
		this.oligoProblems |= OP_PARTIALLY_WRITTEN;
	}

	/**
	 * Calculate the part of the objective function due to one primer.
	 */
	public double p_obj_fn(P3GlobalSettings pa, OligoType type)
			throws PrimerRecordException {
		double sum = 0;

		if (type == OligoType.OT_LEFT || type == OligoType.OT_RIGHT) {
			if (pa.primersArgs.weights.getTemperatureQT() > 0
					&& this.temp > pa.primersArgs.getOptTm())
				sum += pa.primersArgs.weights.getTemperatureQT()
						* (this.temp - pa.primersArgs.getOptTm());
			if (pa.primersArgs.weights.getTemperatureLT() > 0
					&& this.temp < pa.primersArgs.getOptTm())
				sum += pa.primersArgs.weights.getTemperatureLT()
						* (pa.primersArgs.getOptTm() - this.temp);

			if (pa.primersArgs.weights.getGcContentQT() > 0
					&& this.gc_content > pa.primersArgs.getOptGC())
				sum += pa.primersArgs.weights.getGcContentQT()
						* (this.gc_content - pa.primersArgs.getOptGC());
			if (pa.primersArgs.weights.getGcContentLT() > 0
					&& this.gc_content < pa.primersArgs.getOptGC())
				sum += pa.primersArgs.weights.getGcContentLT()
						* (pa.primersArgs.getOptGC() - this.gc_content);

			if (pa.primersArgs.weights.getLengthLT() > 0
					&& this.length < pa.primersArgs.getOptSize())
				sum += pa.primersArgs.weights.getLengthLT()
						* (pa.primersArgs.getOptSize() - this.length);
			if (pa.primersArgs.weights.getLengthQT() > 0
					&& this.length > pa.primersArgs.getOptSize())
				sum += pa.primersArgs.weights.getLengthQT()
						* (this.length - pa.primersArgs.getOptSize());

			/* edited by M. Lepamets */
			if (pa.primersArgs.weights.getFailureRate() > 0) {
				sum += pa.primersArgs.weights.getFailureRate() * this.failure_rate;
			}

			/* BEGIN: secondary structures */
			// pa.thermodynamic_oligo_alignment ==
			// P3GlobalSettings.DONOT_USE_THERMODYNAMICS_ALIGNMENT
			if (!pa.isThermodynamicOligoAlignment()) {

				if (pa.primersArgs.weights.getComplAny() > 0)
					sum += pa.primersArgs.weights.getComplAny() * this.self_any;
				if (pa.primersArgs.weights.getComplEnd() > 0)
					sum += pa.primersArgs.weights.getComplEnd() * this.self_end;

			} else if (pa.isThermodynamicOligoAlignment()) { // ==
																// P3GlobalSettings.USE_THERMODYNAMICS_ALIGNMENT

				if (pa.primersArgs.weights.getComplAnyTh() > 0) {
					if ((this.temp - pa.primersArgs.weights.getTemperatureCutoff()) <= this.self_any)
						sum += pa.primersArgs.weights.getComplAnyTh()
								* (this.self_any - (this.temp
										- pa.primersArgs.weights.getTemperatureCutoff() - 1.0)); /*
																					 * -
																					 * 1.0
																					 * is
																					 * added
																					 * for
																					 * the
																					 * case
																					 * where
																					 * ==
																					 */
					else
						sum += pa.primersArgs.weights.getComplAnyTh()
								* (1 / (this.temp
										- pa.primersArgs.weights.getTemperatureCutoff()
										+ 1.0 - this.self_any));
				}

				if (pa.primersArgs.weights.getComplEndTh() > 0) {
					if ((this.temp - pa.primersArgs.weights.getTemperatureCutoff()) <= this.self_end)
						sum += pa.primersArgs.weights.getComplEndTh()
								* (this.self_end - (this.temp
										- pa.primersArgs.weights.getTemperatureCutoff() - 1.0));
					else
						sum += pa.primersArgs.weights.getComplEndTh()
								* (1 / (this.temp
										- pa.primersArgs.weights.getTemperatureCutoff()
										+ 1.0 - this.self_end));
				}

				if (pa.primersArgs.weights.getHairpinTh() > 0) {
					if ((this.temp - pa.primersArgs.weights.getTemperatureCutoff()) <= this.hairpin_th)
						sum += pa.primersArgs.weights.getHairpinTh()
								* (this.hairpin_th - (this.temp
										- pa.primersArgs.weights.getTemperatureCutoff() - 1.0));
					else
						sum += pa.primersArgs.weights.getHairpinTh()
								* (1 / (this.temp
										- pa.primersArgs.weights.getTemperatureCutoff()
										+ 1.0 - this.hairpin_th));
				}

			}
			// else {
			// throw new
			// PrimerRecordException("  illegal value for pa.thermodynamic_oligo_alignment ");
			// /* Programming error,
			// illegal value for pa.thermodynamic_oligo_alignment */
			// }
			/* END: secondary structures */

			if (pa.primersArgs.weights.getNumNs() > 0)
				sum += pa.primersArgs.weights.getNumNs() * this.num_ns;
			if (pa.primersArgs.weights.getRepeatSimilarity() > 0)
				sum += pa.primersArgs.weights.getRepeatSimilarity()
						* this.repeat_sim.score.get(this.repeat_sim.max);
			if (!bf_get_overlaps_target()) {
				/*
				 * We might be evaluating p_obj_fn with this.target if the
				 * client supplied 'pick_anyway' and specified a primer or
				 * oligo.
				 */
				if ((bf_get_infinite_pos_penalty()))
					throw new PrimerRecordException(
							"  bf_get_infinite_pos_penalty is set ");
				;
				if (pa.primersArgs.weights.getPosPenalty() > 0)
					sum += pa.primersArgs.weights.getPosPenalty()
							* this.position_penalty;
			}
			if (pa.primersArgs.weights.getEndStability() > 0)
				sum += pa.primersArgs.weights.getEndStability()
						* this.end_stability;

			/* FIX ME QUALITY WT Change here */
			if (pa.primersArgs.weights.getSeqQuality() > 0)
				sum += pa.primersArgs.weights.getSeqQuality()
						* (pa.getQualityRangeMax() - this.seq_quality); /*
																		 * Look
																		 * for
																		 * end
																		 * seq
																		 * quality
																		 */

			if (pa.primersArgs.weights.getTemplateMispriming() > 0
					&& !pa.isThermodynamicTemplateAlignment()) {
				if (!(this.oligo_max_template_mispriming() != LibPrimer3.ALIGN_SCORE_UNDEF))
					throw new PrimerRecordException(
							" Error in : oligo_max_template_mispriming exceed MaxValue ");
				sum += pa.primersArgs.weights.getTemplateMispriming()
						* this.oligo_max_template_mispriming();
			}

			if (pa.primersArgs.weights.getTemplateMisprimingTh() > 0
					&& pa.isThermodynamicTemplateAlignment()) {

				if (this.oligo_max_template_mispriming_thermod() == LibPrimer3.ALIGN_SCORE_UNDEF)
					throw new PrimerRecordException(
							" Error in : oligo_max_template_mispriming_thermod exceed MaxValue ");

				if ((this.temp - pa.primersArgs.weights.getTemperatureCutoff()) <= this
						.oligo_max_template_mispriming_thermod())
					sum += pa.primersArgs.weights.getTemplateMisprimingTh()
							* (this.oligo_max_template_mispriming_thermod() - (this.temp
									- pa.primersArgs.weights.getTemperatureCutoff() - 1.0));
				if ((this.temp - pa.primersArgs.weights.getTemperatureCutoff()) > this
						.oligo_max_template_mispriming_thermod())
					sum += pa.primersArgs.weights.getTemplateMisprimingTh()
							* (1 / (this.temp
									- pa.primersArgs.weights.getTemperatureCutoff() + 1.0 - this
										.oligo_max_template_mispriming_thermod()));
			}
			return sum;
		} else if (type == OligoType.OT_INTL) {
			if (pa.oligosArgs.weights.getTemperatureQT() > 0
					&& this.temp > pa.oligosArgs.getOptTm())
				sum += pa.oligosArgs.weights.getTemperatureQT()
						* (this.temp - pa.oligosArgs.getOptTm());
			if (pa.oligosArgs.weights.getTemperatureLT() > 0
					&& this.temp < pa.oligosArgs.getOptTm())
				sum += pa.oligosArgs.weights.getTemperatureLT()
						* (pa.oligosArgs.getOptTm() - this.temp);

			if (pa.oligosArgs.weights.getGcContentQT() > 0
					&& this.gc_content > pa.oligosArgs.getOptGC())
				sum += pa.oligosArgs.weights.getGcContentQT()
						* (this.gc_content - pa.oligosArgs.getOptGC());
			if (pa.oligosArgs.weights.getGcContentLT() > 0
					&& this.gc_content < pa.oligosArgs.getOptGC())
				sum += pa.oligosArgs.weights.getGcContentLT()
						* (pa.oligosArgs.getOptGC() - this.gc_content);

			if (pa.oligosArgs.weights.getLengthLT() > 0
					&& this.length < pa.oligosArgs.getOptSize())
				sum += pa.oligosArgs.weights.getLengthLT()
						* (pa.oligosArgs.getOptSize() - this.length);
			if (pa.oligosArgs.weights.getLengthQT() > 0
					&& this.length > pa.oligosArgs.getOptSize())
				sum += pa.oligosArgs.weights.getLengthQT()
						* (this.length - pa.oligosArgs.getOptSize());
			// Refactored
			// pa.thermodynamic_oligo_alignment==P3GlobalSettings.DONOT_USE_THERMODYNAMICS_ALIGNMENT

			if (!pa.isThermodynamicOligoAlignment()) {
				if (pa.oligosArgs.weights.getComplAny() > 0) // &&
															// !pa.thermodynamic_oligo_alignment
					sum += pa.oligosArgs.weights.getComplAny() * this.self_any;
				if (pa.oligosArgs.weights.getComplEnd() > 0) // &&
															// !pa.thermodynamic_oligo_alignment
					sum += pa.oligosArgs.weights.getComplEnd() * this.self_end;
			}
			/* begin thermodynamical approach */
			// == P3GlobalSettings.USE_THERMODYNAMICS_ALIGNMENT
			if (pa.isThermodynamicOligoAlignment()) {

				if (pa.oligosArgs.weights.getComplAnyTh() > 0) {
					if ((this.temp - pa.oligosArgs.weights.getTemperatureCutoff()) <= this.self_any)
						sum += pa.oligosArgs.weights.getComplAnyTh()
								* (this.self_any - (this.temp
										- pa.oligosArgs.weights.getTemperatureCutoff() - 1.0)); /*
																					 * -
																					 * 1.0
																					 * is
																					 * added
																					 * for
																					 * the
																					 * case
																					 * where
																					 * ==
																					 */
					else
						sum += pa.oligosArgs.weights.getComplAnyTh()
								* (1 / (this.temp
										- pa.oligosArgs.weights.getTemperatureCutoff()
										+ 1.0 - this.self_any));
				}

				if (pa.oligosArgs.weights.getComplEndTh() > 0) {
					if ((this.temp - pa.oligosArgs.weights.getTemperatureCutoff()) <= this.self_end)
						sum += pa.oligosArgs.weights.getComplEndTh()
								* (this.self_end - (this.temp
										- pa.oligosArgs.weights.getTemperatureCutoff() - 1.0));
					else
						sum += pa.oligosArgs.weights.getComplEndTh()
								* (1 / (this.temp
										- pa.oligosArgs.weights.getTemperatureCutoff()
										+ 1.0 - this.self_end));
				}

				if (pa.oligosArgs.weights.getHairpinTh() > 0) {
					if ((this.temp - pa.oligosArgs.weights.getTemperatureCutoff()) <= this.hairpin_th)
						sum += pa.oligosArgs.weights.getHairpinTh()
								* (this.hairpin_th - (this.temp
										- pa.oligosArgs.weights.getTemperatureCutoff() - 1.0));
					else
						sum += pa.oligosArgs.weights.getHairpinTh()
								* (1 / (this.temp
										- pa.oligosArgs.weights.getTemperatureCutoff()
										+ 1.0 - this.hairpin_th));
				}
			}
			/* end thermodynamical approach */

			if (pa.oligosArgs.weights.getNumNs() > 0)
				sum += pa.oligosArgs.weights.getNumNs() * this.num_ns;
			if (pa.oligosArgs.weights.getRepeatSimilarity() > 0)
				sum += pa.oligosArgs.weights.getRepeatSimilarity()
						* this.repeat_sim.score.get(this.repeat_sim.max);

			/* FIXME :: QUALITY WT */
			if (pa.oligosArgs.weights.getSeqQuality() > 0)
				sum += pa.oligosArgs.weights.getSeqQuality()
						* (pa.getQualityRangeMax() - this.seq_quality);

			this.quality = sum;
			return sum;
		} else {
			throw new PrimerRecordException("Using unknown Oligo Type"); /*
																		 * Programmig
																		 * error
																		 * .
																		 */
		}

	}

	/*
	 * Return max of this.template_mispriming and this.template_mispriming_r
	 * (max template mispriming on either strand).
	 */
	public double oligo_max_template_mispriming() {
		return this.template_mispriming > this.template_mispriming_r ? this.template_mispriming
				: this.template_mispriming_r;
	}

	public  double oligo_max_template_mispriming_thermod() {
		return this.template_mispriming > this.template_mispriming_r ? this.template_mispriming
				: this.template_mispriming_r;
	}

	/* ============================================================ */
	/* END functions for setting and getting oligo problems */
	/* ============================================================ */

	/**
	 * Compute various characteristics of the oligo, and determine if it is
	 * acceptable.
	 * 
	 * @throws AlignmentException
	 * @throws ThermodynamicAlignmentException
	 */
	public int calc_and_check_oligo_features(P3GlobalSettings pa,
			OligoType otype, DPAlArgHolder dpal_arg_to_use,
			THAlArgHolder thal_arg_to_use, SeqArgs sa, OligoStats stats,
			P3RetVal retval
			//,
			/* This is 5'.3' on the template sequence: */
			// char[] input_oligo_seq
			) throws AlignmentException,
			ThermodynamicAlignmentException {

		
		if(status != PrimerPair.PAIR_UNCHARACTRIZED)
			System.err.println("Stop");
		
		this.sa = sa;
		
		
		
		final double OUTSIDE_START_WT = 30.0;
		final double INSIDE_START_WT = 20.0;
		final double INSIDE_STOP_WT = 100.0;
		final double OUTSIDE_STOP_WT = 0.5;

		int i, j, k, for_i, gc_count;
		int three_prime_pos; /* position of 3' base of oligo */
		// oligo_type l = otype;
		int poly_x, max_poly_x;
		boolean must_use = this.must_use;
		boolean three_conditions = (must_use || pa.getFileFlag() != 0 || retval.output_type == P3OutputType.primer_list);
		char[] seq = sa.getTrimmedSequence();
		ThermodynamicAlignmentArguments thal_args_for_template_mispriming = LibPrimer3.use_end_for_th_template_mispriming == 1 ? thal_arg_to_use.end1
				: thal_arg_to_use.any;



		PrimersOligosArguments po_args;
		oligo_pair op = new oligo_pair();

		/* Initialize slots in h */
		this.initialize_op();
		this.overlaps = false;

		/*
		 * Set repeat_sim to NULL as indicator that the repeat_sim struct is not
		 * initialized.
		 */
		this.repeat_sim.score.clear();

		this.gc_content = this.num_ns = 0;
		this.overlaps_overlap_position = false;
		this.template_mispriming = this.template_mispriming_r = LibPrimer3.ALIGN_SCORE_UNDEF;
		this.template_mispriming_ok = 0;

		// no need for this
		// PR_ASSERT(OT_LEFT == l || OT_RIGHT == l || OT_INTL == l);



		if (OligoType.OT_INTL == otype) {
			po_args = pa.oligosArgs;
		} else {
			po_args = pa.primersArgs;
		}

		/* Set j and k, and sanity check */
		if (OligoType.OT_LEFT == otype || OligoType.OT_INTL == otype) {
			j = this.start;
			three_prime_pos = k = j + this.length - 1;
		} else {
			three_prime_pos = j = this.start - this.length + 1;
			k = this.start;
		}

		if (k >= 0)
			; // PR_ASSERT
		if (k < sa.getIncludedRegionLength())
			; // PR_ASSERT TRIMMED_SEQ_LEN

		if ((otype == OligoType.OT_LEFT)
				// FIXME :: do not depend on this values PR_NULL_START_CODON_POS
				&& !(sa.getStartCodonPos() <= LibPrimer3.PR_NULL_START_CODON_POS)
				/*
				 * Make sure the primer would amplify at least part of the ORF.
				 */
				&& (0 != (this.start - sa.getStartCodonPos()) % 3
						|| this.start <= retval.upstream_stop_codon 
						|| (retval.stop_codon_pos != -1 && this.start >= retval.stop_codon_pos))) {
			stats.no_orf++;
			this.op_set_does_not_amplify_orf();
			if (!pa.isPickAnyway())
				return status = PrimerPair.PAIR_FAILED;
		}

		/* edited by T. Koressaar and M. Lepamets for lowercase masking */
		if (pa.isLowercaseMasking()) {
			char[] sequence_check = sa.getTrimmedOrigSequence();
			if (pa.isMaskTemplate()) {
				if (otype == OligoType.OT_LEFT)
					sequence_check = sa.getTrimmedMaskedSeq();
				else if (otype == OligoType.OT_RIGHT)
					sequence_check = sa.getTrimmedMaskedSeqRev();
			}
			// is_lowercase_masked(three_prime_pos,sequence_check,this, stats)
			if (is_lowercase_masked(sequence_check[three_prime_pos], stats)) {
				if (!must_use)
					return status = PrimerPair.PAIR_FAILED;
			}
		} else if (pa.isLowercaseMasking() == false
				&& pa.isMaskTemplate() == true) {
			retval.warnings
					.append("Use PRIMER_LOWERCASE_MASKING=1 when using PRIMER_MASK_TEMPLATE=1.");
			return status = PrimerPair.PAIR_FAILED;
		}
		
		
		
		this.setOligoSeq();
		//char[] s1_rev;
		char[] oligo_seq;
		char[] revc_oligo_seq;
//		s1_rev = Sequence.p3_reverse_complement(input_oligo_seq);

		if (OligoType.OT_RIGHT == otype) {
			oligo_seq = this.oligoRevSeq;
			revc_oligo_seq = this.oligoSeq;
		} else {
			oligo_seq = this.oligoSeq;
			revc_oligo_seq = this.oligoRevSeq;
		}
		
		
		
		/* edited by A. Untergasser for forcing sequence use */
		if ((po_args.must_match_five_prime != null)
				|| (po_args.must_match_three_prime != null)) {
			if (this.primer_must_match(pa, stats, oligo_seq,
					po_args.must_match_three_prime,
					po_args.must_match_five_prime)) {
				if (!must_use) {
					op_set_must_match_err();
					return status = PrimerPair.PAIR_FAILED;
				}
			}
		}
		/* end A. Untergasser's changes */

		gc_and_n_content(j, k - j + 1, sa.getTrimmedSequence());

		if (this.num_ns > po_args.getMaxNumOfNsAccepted()) {
			this.op_set_too_many_ns();
			stats.ns++;
			if (!must_use)
				return status = PrimerPair.PAIR_FAILED;
		}

		/*
		 * Upstream error checking has ensured that we use non-default position
		 * penalties only when there is 0 or 1 target.
		 */
		// PR_ASSERT(sa.tar2.count <= 1 || pa.PR_DEFAULT_POSITION_PENALTIES());

		if (pa.getPrimerTask() == P3Task.PICK_SEQUENCING_PRIMERS) {
			this.position_penalty = 0.0;
		} else if (otype != OligoType.OT_INTL
				&& pa.isDefaultPositionPenalties()
				&& sa.getTargetRegions().oligoOverlapsInterval(j, k - j + 1)) {
			this.position_penalty = 0.0;
			this.bf_set_infinite_pos_penalty(1);
			this.bf_set_overlaps_target(1);
		} else if (otype != OligoType.OT_INTL
				&& !pa.isDefaultPositionPenalties() && 1 == sa.getTargetRegions().getCount()) {
			this.compute_position_penalty(pa, sa, otype);
			if (this.bf_get_infinite_pos_penalty()) {
				this.bf_set_overlaps_target(1);
			}
		} else {
			this.position_penalty = 0.0;
			this.bf_set_infinite_pos_penalty(0);
		}

		if (!sa.PR_START_CODON_POS_IS_NULL()) {
			if (OligoType.OT_LEFT == otype) {
				if (sa.getStartCodonPos() > this.start)
					this.position_penalty = (sa.getStartCodonPos() - this.start)
							* OUTSIDE_START_WT;
				else
					this.position_penalty = (this.start - sa.getStartCodonPos())
							* INSIDE_START_WT;
			} else if (OligoType.OT_RIGHT == otype) {
				if (-1 == retval.stop_codon_pos) {
					this.position_penalty = ((sa.getIncludedRegionLength()) - this.start - 1)
							* INSIDE_STOP_WT;
				} else if (retval.stop_codon_pos < this.start) {
					this.position_penalty = (this.start - retval.stop_codon_pos)
							* OUTSIDE_STOP_WT;
				} else {
					this.position_penalty = (retval.stop_codon_pos - this.start)
							* INSIDE_STOP_WT;
				}
			}
		}

		/* TO DO Simplify logic here */
		if (otype != OligoType.OT_INTL
				&& sa.getExcludedRegions().oligoOverlapsInterval(j, k - j + 1))
			bf_set_overlaps_excl_region(1);

		if (otype == OligoType.OT_INTL
				&& sa.getExcludedInternalRegions().oligoOverlapsInterval(j, k - j + 1))
			bf_set_overlaps_excl_region(1);

		if (otype != OligoType.OT_INTL && bf_get_overlaps_target()) {
			op_set_overlaps_target();
			stats.target++;
			if (!must_use)
				return  status = PrimerPair.PAIR_FAILED;
		}

		if (bf_get_overlaps_excl_region()) {
			op_set_overlaps_excluded_region();
			stats.excluded++;
			if (!must_use)
				return status = PrimerPair.PAIR_FAILED;
		}

		/* Check if the oligo is included in any ok region */
		
		boolean included = false;
		included = sa.getOkRegions().chechIncludedInAny(otype,j,k);
		if (!included) {
			op_set_not_in_any_ok_region();
			if(otype == OligoType.OT_LEFT)
				stats.not_in_any_left_ok_region++;
			else
				stats.not_in_any_right_ok_region++;
			if (!must_use)
				return status = PrimerPair.PAIR_FAILED;
		}
//		if ((otype == OligoType.OT_LEFT) && (sa.ok_regions.count > 0)
//				&& (!sa.ok_regions.any_left)) {
//			boolean included = false;
//			for (i = 0; i < sa.ok_regions.count; i++) {
//				if ((j >= sa.ok_regions.left_pairs[i][0])
//						&& (k <= sa.ok_regions.left_pairs[i][0]
//								+ sa.ok_regions.left_pairs[i][1] - 1)) {
//					included = true;
//					break;
//				}
//			}
//			if (!included) {
//				op_set_not_in_any_ok_region();
//				stats.not_in_any_left_ok_region++;
//				if (!must_use)
//					return;
//			}
//		}
//		
//		if ((otype == OligoType.OT_RIGHT) && (sa.ok_regions.count > 0)
//				&& (!sa.ok_regions.any_right)) {
//			boolean included = false;
//			for (i = 0; i < sa.ok_regions.count; i++) {
//				if ((j >= sa.ok_regions.right_pairs[i][0])
//						&& (k <= sa.ok_regions.right_pairs[i][0]
//								+ sa.ok_regions.right_pairs[i][1] - 1)) {
//					included = true;
//					break;
//				}
//			}
//			if (!included) {
//				op_set_not_in_any_ok_region();
//				stats.not_in_any_right_ok_region++;
//				if (!must_use)
//					return;
//			}
//		}

		if (this.gc_content < po_args.getMinGC()) {
			op_set_low_gc_content();
			stats.gc++;
			if (!must_use)
				return status = PrimerPair.PAIR_FAILED;
		} else if (this.gc_content > po_args.getMaxGC()) {
			op_set_high_gc_content();
			stats.gc++;
			if (!must_use)
				return status = PrimerPair.PAIR_FAILED;
		}

		if (OligoType.OT_LEFT == otype || OligoType.OT_RIGHT == otype) {
			/*
			 * gc_clamp is applicable only to primers (as opposed to primers and
			 * hybridzations oligos.
			 */
			for (i = 0; i < pa.getGcClamp(); i++) {
				/*
				 * We want to look at the 3' end of the oligo being assessed, so
				 * we look in the 5' end of its reverse-complement
				 */
				if (revc_oligo_seq[i] != 'G' && revc_oligo_seq[i] != 'C') {
					op_set_no_gc_glamp();
					stats.gc_clamp++;
					if (!must_use)
						return status = PrimerPair.PAIR_FAILED;
					else
						break;
				}
			}
		}

		if (pa.getMaxEndGC() < 5
				&& (OligoType.OT_LEFT == otype || OligoType.OT_RIGHT == otype)) {
			/*
			 * The CGs are only counted in the END of primers (as opposed to
			 * primers and hybridzations oligos.
			 */
			gc_count = 0;
			for (i = 0; i < 5; i++) {
				/*
				 * We want to look at the 3' end of the oligo being assessed, so
				 * we look in the 5' end of its reverse-complement
				 */
				if (revc_oligo_seq[i] == 'G' || revc_oligo_seq[i] == 'C') {
					gc_count++;
				}
			}
			if (gc_count > pa.getMaxEndGC()) {
				op_set_too_many_gc_at_end();
				stats.gc_end_high++;
				if (!must_use)
					return  status = PrimerPair.PAIR_FAILED;
			}
		}

		/*
		 * Tentative interface for generic oligo testing function
		 */
		if (!sequence_quality_is_ok(pa, otype, sa, j, k, stats, po_args)
				&& !must_use)
			return status = PrimerPair.PAIR_FAILED;

		max_poly_x = po_args.getMaxPolyX();
		if (max_poly_x > 0) {
			poly_x = 1;
			for (i = j + 1; i <= k; i++) {
				if (seq[i] == seq[i - 1] || seq[i] == 'N') {
					poly_x++;
					if (poly_x > max_poly_x) {
						op_set_high_poly_x();
						stats.poly_x++;
						if (!must_use)
							return status = PrimerPair.PAIR_FAILED;
						else
							break;
					}
				} else
					poly_x = 1;
			}
		}

		this.temp = OligoTMCalculator.sequenceTM(oligo_seq, po_args.getDnaConcentration(),
				po_args.getSaltConcentration(),
				po_args.getDivalentConcentration(),
				po_args.getDntpConcentration(), LibPrimer3.MAX_NN_TM_LENGTH,
				pa.getMeltingTemperatureMethod(), pa.getSaltCorrectionMethod());

		if (this.temp < po_args.getMinTm()) {
			op_set_low_tm();
			stats.temp_min++;
			if (!must_use)
				return status = PrimerPair.PAIR_FAILED;
		}

		if (this.temp > po_args.getMaxTm()) {
			op_set_high_tm();
			stats.temp_max++;
			if (!must_use)
				return status = PrimerPair.PAIR_FAILED;
		}

		/*
		 * End stability is applicable only to primers (not to oligos)
		 */
		if (OligoType.OT_LEFT == otype || OligoType.OT_RIGHT == otype) {
			if ((this.end_stability = OligoTMCalculator.end_oligodg(oligo_seq, 5,
					pa.getMeltingTemperatureMethod())) > pa.getMaxEndStability()) {
				/* Must not be called on a hybridization probe / internal oligo: */
				op_set_high_end_stability();
				stats.stability++;
				if (!must_use)
					return status = PrimerPair.PAIR_FAILED;
			}
		}

		if ((must_use || pa.getFileFlag() != 0
				|| retval.output_type == P3OutputType.primer_list
				|| po_args.weights.getComplAny() > 0 || po_args.weights.getComplEnd() > 0)
				&& !pa.isThermodynamicOligoAlignment()) {

			oligo_compl(po_args, stats, dpal_arg_to_use, oligo_seq,
					revc_oligo_seq);

			if ((!(p3_ol_is_uninitialized())) && !must_use) {
				// PR_ASSERT
				if (!p3_ol_is_ok())
					;
				return status = PrimerPair.PAIR_FAILED;
			}
		} else {
			/* Thermodynamical approach: for primers only */
			if ((must_use || pa.getFileFlag() != 0
					|| retval.output_type == P3OutputType.primer_list
					|| po_args.weights.getComplAnyTh() > 0 || po_args.weights.getComplEndTh() > 0)
					&& pa.isThermodynamicOligoAlignment()) {
				oligo_compl_thermod(po_args, stats, thal_arg_to_use, oligo_seq,
						oligo_seq);
				/* input_oligo_seq, input_oligo_seq); */
				/* oligo_seq, revc_oligo_seq); */
				if ((!(p3_ol_is_uninitialized())) && !must_use) {
					// PR_ASSERT
					if (!p3_ol_is_ok())
						;
					return status = PrimerPair.PAIR_FAILED;
				}
			} else {
				this.self_any = this.self_end = LibPrimer3.ALIGN_SCORE_UNDEF;
			}
		}

		if ((three_conditions || po_args.weights.getHairpinTh() > 0
		// #if 0
		// || po_args.weights.compl_any_th /* Triinu, is this needed? */
		// || po_args.weights.compl_end_th /* Triinu, is this needed? */
		// #endif
				)
				&& pa.isThermodynamicOligoAlignment()) {
			this.oligo_hairpin(po_args, stats, thal_arg_to_use,
			/* input_oligo_seq); */
			oligo_seq);
			if ((!(p3_ol_is_uninitialized())) && !must_use) {
				// PR_ASSERT
				if (!p3_ol_is_ok())
					;
				return status = PrimerPair.PAIR_FAILED;
			}
		} else {
			/*
			 * This will get calculated later if necessary, in
			 * characterize_pair.
			 */
			this.hairpin_th = LibPrimer3.ALIGN_SCORE_UNDEF;
		}
		/* end of thermod. approach */

		// Ignored block in source file removed

		// #else

		if (three_conditions || po_args.weights.getRepeatSimilarity() > 0) {
			this.oligo_repeat_library_mispriming(pa, sa, otype, stats,
					dpal_arg_to_use, retval.glob_err);
		}

		if (!OK_OR_MUST_USE())
			return status = PrimerPair.PAIR_FAILED;

		if (three_conditions
				|| ( /* Do we need template mispriming for the penalty function? */
				(OligoType.OT_RIGHT == otype || OligoType.OT_LEFT == otype) && ((pa.primersArgs.weights.getTemplateMispriming() > 0 && !pa
						.isThermodynamicTemplateAlignment()) || (pa.primersArgs.weights.getTemplateMisprimingTh() > 0 && pa
						.isThermodynamicTemplateAlignment())))) {
			if (OK_OR_MUST_USE()) {
				this.oligo_template_mispriming(pa, sa, otype, stats,
						dpal_arg_to_use.local_end,
						thal_args_for_template_mispriming);
			}
		}
		// #endif

		if (this.length > po_args.getMaxSize()) {
			op_set_too_long();
			stats.size_max++;
			if (!must_use)
				return status = PrimerPair.PAIR_FAILED;
		}

		if (this.length < po_args.getMinSize()) {
			op_set_too_short();
			stats.size_min++;
			if (!must_use)
				return  status = PrimerPair.PAIR_FAILED;
		}

		for (for_i = 0; for_i < sa.getPrimerOverlapJunctionsList().size(); for_i++) {
			int value_for_i = sa.getPrimerOverlapJunctionsList().get(for_i);
			if (OligoType.OT_LEFT == otype
					&& ((this.start + pa.getMin5PrimeOverlapOfJunction() - 1) <= value_for_i)
					&& ((this.start + this.length - pa
							.getMin3PrimeOverlapOfJunction())) > value_for_i) {
				this.overlaps_overlap_position = true;
				/* no need to continue checking */
				break;
			}
			if (OligoType.OT_RIGHT == otype
					&& ((this.start - this.length + pa
							.getMin3PrimeOverlapOfJunction()) <= value_for_i)
					&& ((this.start - pa.getMin5PrimeOverlapOfJunction() + 1)) > value_for_i) {
				this.overlaps_overlap_position = true;
				/* no need to continue checking */
				break;
			}
		}
		/* Calculate failure rate */
		/* Added by M. Lepamets */
		this.failure_rate = 0.0;
		if (pa.isMaskTemplate()
				&& this.length >= pa.getMaskingParameters().window_size) {
			op.fwd = masker.string_to_word(oligo_seq, this.length,
					pa.getMaskingParameters().window_size);
			op.rev = op.fwd; /* not used in this calculation */
			op.calculate_scores(pa.getMaskingParameters(),
					pa.getMaskingParameters().window_size);
			this.failure_rate = op.score_fwd;
		}
		/* FIX ME FIXME Steve, is this really needed? */
		op_set_completely_written();
		return  status = PrimerPair.PAIR_OK;
	}

	public void oligo_template_mispriming(P3GlobalSettings pa, SeqArgs sa,
			OligoType l, OligoStats ostats, DPAlignmentArgs d_align_args,
			ThermodynamicAlignmentArguments t_align_args) throws AlignmentException,
			ThermodynamicAlignmentException {
		/* Check if we already did this and the oligo was ok. */

		PrimerRecord h = this;
		if (h.template_mispriming_ok != 0) {
			return;
		}
		int first, last; /*
						 * Indexes of first and last bases of the oligo in
						 * sa.trimmed_seq, that is, WITHIN THE INCLUDED REGION.
						 */
		int[] cor_res = getIndeces(l);
		first = cor_res[0];
		last = cor_res[1];

		char[] s =  this.oligoSeq;//oligo_compute_sequence_and_reverse(sa, l);
		char[] s_r = this.oligoRevSeq ;//Sequence.p3_reverse_complement(s);

		/* Calculate maximum similarity to ectopic sites in the template. */
		if (l == OligoType.OT_RIGHT || l == OligoType.OT_LEFT) {
			if (!pa.isThermodynamicTemplateAlignment()
					&& pa.needTemplateMispriming())
				LibPrimer3.primer_mispriming_to_template(h, pa, sa, l, ostats,
						first, last, s, s_r, d_align_args);
			if (pa.isThermodynamicTemplateAlignment()
					&& pa.needTemplateMisprimingTH())
				LibPrimer3.primer_mispriming_to_template_thermod(h, pa, sa, l,
						ostats, first, last, s, s_r, t_align_args);
		}
	}

	int[] getIndeces(OligoType l) {
		/*
		 * Indexes of first and last bases of the oligo in sa.trimmed_seq, that
		 * is, WITHIN THE INCLUDED REGION.
		 */
		int first = (l == OligoType.OT_LEFT || l == OligoType.OT_INTL) ? this.start
				: this.start - this.length + 1;
		int last = (l == OligoType.OT_LEFT || l == OligoType.OT_INTL) ? this.start
				+ this.length - 1
				: this.start;
		return new int[] { first, last };
	}

	
	
	
	public void oligo_repeat_library_mispriming(P3GlobalSettings pa,
			seq_lib lib,
			OligoType l, OligoStats ostats,
			DPAlArgHolder dpal_arg_to_use, StringBuilder glob_err , Set<String> exculdeSeq) throws AlignmentException
	{
		

		PrimerRecord h = this;
		double w;
		int min, max;
		double max_lib_compl;
		boolean  max_lib_compl_is_percent = false;
		/* First, check the oligo against the repeat library. */
		if (l == OligoType.OT_INTL) {
			max_lib_compl =  pa.oligosArgs.getMaxRepeatCompl();
			max_lib_compl_is_percent = pa.oligosArgs.isMaxRepeatComplIsPercent();
		} else {
			max_lib_compl =  pa.primersArgs.getMaxRepeatCompl();
			max_lib_compl_is_percent = pa.primersArgs.isMaxRepeatComplIsPercent();

		}

		char[] s =  this.oligoSeq ;// oligo_compute_sequence_and_reverse(sa, l);
		char[] s_r =  this.oligoRevSeq;//Sequence.p3_reverse_complement(s);
//		h.repeat_sim.score = new ArrayList<Double>();
		h.repeat_sim.max = h.repeat_sim.min = 0;
		max = min = 0;
		h.repeat_sim.name = lib.getName(0);
		max = h.repeat_sim.max;
		for (int i = 0; i < lib.seq_lib_num_seq(); i++) {
			w = 0;
			if(exculdeSeq == null || !exculdeSeq.contains(lib.getName(i)))
			{
				if (l == OligoType.OT_LEFT)
					w = lib.weight.get(i)
							* LibPrimer3
									.align(s,
											lib.getSeq(i),
											(pa.isLibAmbiguityCodesConsensus() ? dpal_arg_to_use.local_end_ambig
													: dpal_arg_to_use.local_end) ,pa.getMispriming3EndScore() );
	
				else if (l == OligoType.OT_INTL)
					w = lib.weight.get(i)
							* LibPrimer3
									.align(s,
											lib.getSeq(i),
											(pa.isLibAmbiguityCodesConsensus() ? dpal_arg_to_use.local_ambig
													: dpal_arg_to_use.local));
	
				else
					w = lib.weight.get(i)
							* LibPrimer3
									.align(s_r,
											lib.getSeqRevCompl(i),
											(pa.isLibAmbiguityCodesConsensus() ? dpal_arg_to_use.local_end_ambig
													: dpal_arg_to_use.local));
				}
			// if (w > SHRT_MAX || w < SHRT_MIN) {
			// /* This check is necessary for the next 9 lines */
			// pr_append_new_chunk( error,
			// "Out of range error occured calculating match to repeat library");
			// return;
			// }
			
			// TODO :: calc w as percent ,, this is an optional value
			h.repeat_sim.score.add(w);
			if (max_lib_compl_is_percent)
				w = 100 *( w / this.length);

			if (w > max) {
				max = (int) w;
				h.repeat_sim.max = h.repeat_sim.score.size()-1;
				h.repeat_sim.name = lib.getName(i);
			}
			if (w < min) {
				min = (int) w;
				h.repeat_sim.min = (short) i;
			}

			if (w > max_lib_compl) {
				op_set_high_similarity_to_non_template_seq();
				ostats.repeat_score++;
				ostats.ok--;
				if (!h.must_use)
					return;
			} /* if w > max_lib_compl */
		} /* for */
	}
	
	
	
	
	
	public void oligo_repeat_library_mispriming(
			P3GlobalSettings pa,
			// TODO :: clean sa args no longer needed
			SeqArgs sa, 
			// TODO l is a member of h
			OligoType l, 
			OligoStats ostats,
			DPAlArgHolder dpal_arg_to_use, StringBuilder glob_err)
			throws AlignmentException {
		
		seq_lib lib;
//		int i;

//		int min, max;
//		double max_lib_compl;
//		boolean  max_lib_compl_is_percent = false;
		/* First, check the oligo against the repeat library. */
		if (l == OligoType.OT_INTL) {
			lib = pa.oligosArgs.repeat_lib;
//			max_lib_compl =  pa.oligosArgs.getMaxRepeatCompl();
//			max_lib_compl_is_percent = pa.oligosArgs.maxRepeatComplIsPercent;
		} else {
			lib = pa.primersArgs.repeat_lib;
//			max_lib_compl =  pa.primersArgs.getMaxRepeatCompl();
//			max_lib_compl_is_percent = pa.primersArgs.maxRepeatComplIsPercent;

		}

//		char[] s =  this.oligoSeq ;// oligo_compute_sequence_and_reverse(sa, l);
//		char[] s_r =  this.oligoRevSeq;//Sequence.p3_reverse_complement(s);

		/*
		 * Calculate maximum similarity to sequences from user defined repeat
		 * library. Compare it with maximum allowed repeat similarity.
		 */

		if (lib != null) {
			/* Library exists and is non-empty. */
			oligo_repeat_library_mispriming(pa,
					lib,
					l, ostats,
					dpal_arg_to_use, glob_err, null);

		} /* if library exists and is non-empty */
		/* End of checking against the repeat library */

	}

	// oligo_compute_sequence_and_reverse(primer_rec *h, const seq_args
	// *sa,oligo_type l, int *first, int *last,char *s, char *s_r)

	/**
	 * TODO :: revise this to the orginal verion return only sequence
	 * 
	 * @param sa
	 * @param l
	 * @return
	 */
	@Deprecated
	private char[] oligo_compute_sequence_and_reverse(SeqArgs sa, OligoType l) {

		int first = (l == OligoType.OT_LEFT || l == OligoType.OT_INTL) ? this.start
				: this.start - this.length + 1;
		int last = (l == OligoType.OT_LEFT || l == OligoType.OT_INTL) ? this.start
				+ this.length - 1
				: this.start;
		return Sequence._pr_substr(sa.getTrimmedSequence(), first, this.length);
	}

	/**
	 * Calculate self complementarity (both h->self_any and h->self_end) which
	 * we use as an approximation for both secondary structure and self
	 * primer-dimer.
	 */
	public void oligo_hairpin(PrimersOligosArguments po_args, OligoStats ostats,
			THAlArgHolder thal_arg_to_use, char[] oligo_seq)
			throws ThermodynamicAlignmentException {
		this.hairpin_th = LibPrimer3.align_thermod(oligo_seq, oligo_seq,
				thal_arg_to_use.hairpin_th);
		if (this.hairpin_th > po_args.getMaxHairPinTH()) {
			op_set_high_hairpin_th();
			ostats.hairpin_th++;
			ostats.ok--;
		}
	}

	public void oligo_compl_thermod(PrimersOligosArguments po_args,
			OligoStats ostats, THAlArgHolder thal_arg_to_use,
			char[] oligo_seq, char[] revc_oligo_seq)
			throws ThermodynamicAlignmentException {
		this.self_any = LibPrimer3.align_thermod(oligo_seq, revc_oligo_seq,
				thal_arg_to_use.any);
		if (this.self_any > po_args.getMaxSelfAnyTH()) {
			this.op_set_high_self_any();
			ostats.compl_any++;
			ostats.ok--;
			if (!this.must_use)
				return;
		}
		this.self_end = LibPrimer3.align_thermod(oligo_seq, revc_oligo_seq,
				thal_arg_to_use.end1);
		if (this.self_end > po_args.getMaxSelfEndTH()) {
			this.op_set_high_self_end();
			ostats.compl_end++;
			ostats.ok--;
			if (!this.must_use)
				return;
		}
	}

	public void oligo_compl(PrimersOligosArguments po_args, OligoStats ostats,
			DPAlArgHolder dpal_arg_to_use, char[] oligo_seq,
			char[] revc_oligo_seq) throws AlignmentException {
		this.self_any = LibPrimer3.align(oligo_seq, revc_oligo_seq,
				dpal_arg_to_use.local);
		if (this.self_any > po_args.getMaxSelfAny()) {
			this.op_set_high_self_any();
			ostats.compl_any++;
			ostats.ok--;
			if (!this.must_use)
				return;
		}

		this.self_end = LibPrimer3.align(oligo_seq, revc_oligo_seq,
				dpal_arg_to_use.end);
		if (this.self_end > po_args.getMaxSelfEnd()) {
			this.op_set_high_self_end();
			ostats.compl_end++;
			ostats.ok--;
			if (!this.must_use)
				return;
		}
	}

	/**
	 * Calculate the minimum sequence quality and the minimum sequence quality
	 * of the 3' end of a primer or oligo. Set h.seq_quality and
	 * h.seq_end_quality with these values. Return 1 (== ok) if sa.quality is
	 * undefined. Otherwise, return 1 (ok) if h.seq_quality and
	 * h.seq_end_quality are within range, or else return 0.
	 */
	private boolean sequence_quality_is_ok(P3GlobalSettings pa, OligoType l,
			SeqArgs sa, int j, int k, OligoStats global_oligo_stats,
			PrimersOligosArguments po_args) {

		int i, min_q, min_q_end, m, q;
		boolean retval = true;

		if (null == sa.getSequenceQuality()) {
			this.seq_end_quality = this.seq_quality = pa.getQualityRangeMax();
			return true;
		}

		q = pa.getQualityRangeMax();

		min_q = po_args.getMinQuality();
		if (OligoType.OT_LEFT == l || OligoType.OT_RIGHT == l) {
			min_q_end = po_args.getMinEndQuality();
		} else {
			min_q_end = min_q;
		}

		if (OligoType.OT_LEFT == l || OligoType.OT_INTL == l) {

			for (i = k - 4; i <= k; i++) {
				if (i < j)
					continue;
				m = sa.getSequenceQuality()[i + sa.getIncludedRegionStart()];
				if (m < q)
					q = m;
			}
			min_q_end = q;

			for (i = j; i <= k - 5; i++) {
				m = sa.getSequenceQuality()[i + sa.getIncludedRegionStart()];
				if (m < q)
					q = m;
			}
			min_q = q;

		} else if (OligoType.OT_RIGHT == l) {
			for (i = j; i < j + 5; i++) {
				if (i > k)
					break;
				m = sa.getSequenceQuality()[i + sa.getIncludedRegionStart()];
				if (m < q)
					q = m;
			}
			min_q_end = q;

			for (i = j + 5; i <= k; i++) {
				m = sa.getSequenceQuality()[i + sa.getIncludedRegionStart()];
				if (m < q)
					q = m;
			}
			min_q = q;
		}
		// else {
		// PR_ASSERT(0); /* Programming error. */
		// }

		this.seq_quality = min_q;
		this.seq_end_quality = min_q_end;

		if (this.seq_quality < po_args.getMinQuality()) {
			op_set_low_sequence_quality();
			global_oligo_stats.seq_quality++;
			// retval = 0;
			return false;
		}

		if (OligoType.OT_LEFT == l || OligoType.OT_RIGHT == l) {
			if (this.seq_end_quality < po_args.getMinEndQuality()) {
				op_set_low_end_sequence_quality();
				global_oligo_stats.seq_quality++;
				retval = false;
			}
		}
		return retval;
	}

	private void compute_position_penalty(P3GlobalSettings pa, SeqArgs sa,
			OligoType o_type) {
		int three_prime_base;
		boolean inside_flag = false;
		int target_begin, target_end;

		// PR_ASSERT(oligo_type.OT_LEFT == o_type || oligo_type.OT_RIGHT ==
		// o_type);
		// PR_ASSERT(1 == sa.tar2.count);
		target_begin = sa.getTargetRegions().getInterval(0)[0];
		target_end = target_begin + sa.getTargetRegions().getInterval(0)[1] - 1;

		three_prime_base = OligoType.OT_LEFT == o_type ? this.start
				+ this.length - 1 : this.start - this.length + 1;
		this.bf_set_infinite_pos_penalty(1);
		this.position_penalty = 0.0;

		if (OligoType.OT_LEFT == o_type) {
			if (three_prime_base <= target_end) {
				this.bf_set_infinite_pos_penalty(0);
				if (three_prime_base < target_begin)
					this.position_penalty = target_begin - three_prime_base - 1;
				else {
					this.position_penalty = three_prime_base - target_begin + 1;
					inside_flag = true;
				}
			}
		} else { /* OT_RIGHT == o_type */
			if (three_prime_base >= target_begin) {
				this.bf_set_infinite_pos_penalty(0);
				if (three_prime_base > target_end) {
					this.position_penalty = three_prime_base - target_end - 1;
				} else {
					this.position_penalty = target_end - three_prime_base + 1;
					inside_flag = true;
				}
			}
		}
		if (!inside_flag)
			this.position_penalty *= pa.getOutsidePenalty();
		else {
			this.position_penalty *= pa.getInsidePenalty();
		}
	}

	private void gc_and_n_content(int start, int len, char[] sequence) {

		double[] res = Sequence.calcGCandN(sequence, start, len);
		this.gc_content = res[0];
		this.num_ns = (int) res[1];

	}

	private boolean primer_must_match(P3GlobalSettings pa,
			OligoStats global_oligo_stats,
			/* This is 5'.3' on the template sequence: */
			char[] input_oligo_seq, char[] match_three_prime,
			char[] match_five_prime) {

		int seq = 0;
		int test = 0;
		int length = this.length - 5;
		if (match_five_prime != null) {
			// seq = input_oligo_seq;
			// test = match_five_prime;
			for (int i = 0; i < 5; i++) {
				if (!compare_nucleotides(input_oligo_seq[seq],
						match_five_prime[test])) {
					global_oligo_stats.must_match_fail++;
					return true;
				}
				seq++;
				test++;
			}
		}
		if (match_three_prime != null) {
			seq = 0;
			test = 0;
			// seq = input_oligo_seq;
			// test = match_three_prime;
			seq = seq + length;
			for (int i = 0; i < 5; i++) {
				if (!compare_nucleotides(input_oligo_seq[seq],
						match_three_prime[test])) {
					global_oligo_stats.must_match_fail++;
					return true;
				}
				seq++;
				test++;
			}
		}
		return false;
	}

	/**
	 * For a [NACTG] is allowed, for b [NACTGRYWSMKBHDV].
	 */
	private boolean compare_nucleotides(char a, char b) {

		char x = Character.toUpperCase(a);
		char y = Character.toUpperCase(b);
		/* Convert to uppercase */

		if (x == y) {
			return true;
		}
		if ((x == 'N') || (y == 'N')) {
			return true;
		}
		if (x == 'A') {
			if ((y == 'R') || (y == 'W') || (y == 'M') || (y == 'H')
					|| (y == 'D') || (y == 'V')) {
				return true;
			}
		}
		if (x == 'G') {
			if ((y == 'R') || (y == 'S') || (y == 'K') || (y == 'B')
					|| (y == 'D') || (y == 'V')) {
				return true;
			}
		}
		if (x == 'C') {
			if ((y == 'Y') || (y == 'S') || (y == 'M') || (y == 'B')
					|| (y == 'H') || (y == 'V')) {
				return true;
			}
		}
		if (x == 'T') {
			if ((y == 'Y') || (y == 'W') || (y == 'K') || (y == 'B')
					|| (y == 'H') || (y == 'D')) {
				return true;
			}
		}

		return false;
	}

	/**
	 * Edited by T. Koressaar for lowercase masking. This function checks if the
	 * 3' end of the primer has been masked by lowercase letter. Function
	 * created/Added by Eric Reppo, July 9, 2002
	 */
	private boolean is_lowercase_masked(char p, OligoStats global_oligo_stats) {

		if ('a' == p || 'c' == p || 'g' == p || 't' == p) {
			op_set_overlaps_masked_sequence();
			global_oligo_stats.gmasked++;
			return true;
		}
		return false;
	}

	// TODO :: change this
	public char[] pr_oligo_sequence(SeqArgs sa) {
		// int seq_len;

		// TODO :: add check here
		// PR_ASSERT(NULL != sa);
		// PR_ASSERT(NULL != oligo);
		// seq_len = strlen(sa.sequence);

		// PR_ASSERT(oligo.start + sa.incl_s >= 0);
		// PR_ASSERT(oligo.start + sa.incl_s + oligo.length <= seq_len);

		char[] s = Sequence._pr_substr(sa.getSequence(), sa.getIncludedRegionStart() + this.start,
				this.length);
		return s;
		// return null;
	}
	// TODO :: change this
	public char[] pr_oligo_rev_c_sequence(SeqArgs sa) {
		// TODO :: add check
		int seq_len, start;
		// PR_ASSERT(NULL != sa);
		// PR_ASSERT(NULL != o);
		// seq_len = strlen(sa.sequence);
		start = sa.getIncludedRegionStart() + this.start - this.length + 1;
		// PR_ASSERT(start >= 0);
		// PR_ASSERT(start + o.length <= seq_len);
		char[] s = Sequence._pr_substr(sa.getSequence(), start, this.length);
		return Sequence.p3_reverse_complement(s);
	}


	// Multi targets specfic primers : other targets that this primer can bind to
	// exact
	public HashMap<String,Integer> targetSpecificIndex = new HashMap<String, Integer>();
	public HashMap<String,Double> targetSpecificScore = new HashMap<String, Double>();
	// targets specfic primers
	/**
	 * primer is specific regarding its target list. it could have multiple targert
	 */
	protected boolean isTargetSpecific = true;
	
	
	public void calcSpecific(seq_lib targets_lib,String targetName, String reversePrefixTargets, DPAlignmentArgs dpAlignmentArgs) throws AlignmentException {
		
		HashSet<String> targetsName  = new HashSet<String>();
		
		targetsName.add(targetName);
		calcSpecific(targets_lib, targetsName  , reversePrefixTargets, dpAlignmentArgs);
	}
	public void calcSpecific(seq_lib targets_lib, String reversePrefixTargets, DPAlignmentArgs dpAlignmentArgs) throws AlignmentException {
		
		HashSet<String> targetsName  = new HashSet<String>();
		
		targetsName.add(this.sa.getSequenceName());
		calcSpecific(targets_lib, targetsName  , reversePrefixTargets, dpAlignmentArgs);
	}
	// TODO :: targetName should be internal
	public void calcSpecific(seq_lib targets_lib,Set<String> targetsNameToExculde, String reversePrefixTargets, DPAlignmentArgs dpAlignmentArgs) throws AlignmentException {
		
		
		for(int i = 0 ; i < targets_lib.seq_lib_num_seq(); i++)
		{
			
			DPAlignmentResults r = null;
			
			String oTargetName = targets_lib.getName(i);
			// skip for targets already we know them
			if( targetsNameToExculde.contains(oTargetName) || oTargetName.contains(reversePrefixTargets))
				continue;
			
			if( rec_type == OligoType.OT_RIGHT) {
				r = DPAlignment.dpAlign(this.getOligoRevSeq(), 
						targets_lib.getSeqRevCompl(i),dpAlignmentArgs);
			}else
			{
				r = DPAlignment.dpAlign(this.getOligoSeq(), 
						targets_lib.getSeq(i),dpAlignmentArgs);
			}
			// if score == len then is the same
			// if score is about 90% then see the last 3 and 
			
			if (r.score  ==  (this.length*100 ))
			{
				// TODO :: see if the this is the correct coordinates 
				if(this.rec_type ==OligoType.OT_RIGHT ) {
					targetSpecificIndex.put(oTargetName, (targets_lib.getSeqRevCompl(i).length -r.align_end_2 ) - this.length);
				}
				else
				{
					targetSpecificIndex.put(oTargetName, r.align_end_2-this.length);
				}
				//isTargetSpecific = false;
			}
			else
			{
				if(r.align_end_1 != this.length-1)
				{
					System.err.println("r.align_end_1 != this.length-1");
				}
				// TODO :: read from pa
				double maxScore = 80;
				int end3ScoreCalc = 4;
				
				double finalScore = r.score < 0.0 ? 0.0 : r.score / LibPrimer3.PR_ALIGN_SCORE_PRECISION;
				
				String last3End_Target = "" ,last3End_Primer = "";
				if( rec_type == OligoType.OT_RIGHT) {
					last3End_Target =  String.copyValueOf(targets_lib.getSeqRevCompl(i), r.align_end_2-(end3ScoreCalc-1), end3ScoreCalc);
					last3End_Primer =  String.copyValueOf(this.getOligoRevSeq(), r.align_end_1-(end3ScoreCalc-1), end3ScoreCalc);
				}
				else if ( rec_type == OligoType.OT_LEFT) {	
					last3End_Target =  String.copyValueOf(targets_lib.getSeq(i), r.align_end_2-(end3ScoreCalc-1), end3ScoreCalc);
					last3End_Primer =  String.copyValueOf(this.getOligoSeq(), r.align_end_1-(end3ScoreCalc-1), end3ScoreCalc);
				}
				if(     last3End_Primer.length() ==  last3End_Target.length() ) {
					int matches = 0;
					for(int j = 0; j < end3ScoreCalc;j++ )
					{
						if(last3End_Primer.charAt(j) == last3End_Target.charAt(j)) {
							finalScore++;
							matches++;
						}
					}
					if(matches == end3ScoreCalc )
						finalScore++;
				}
				
				
				
				double wScore = 100 * (finalScore/this.length);
				if(wScore > maxScore ) {
					isTargetSpecific = false;
					targetSpecificScore.put(oTargetName, finalScore);
				}
				
				
				String targetSeqMatch,thisSeqMatch;
//				if( Math.abs( (r.score - (this.length*100))) <= 1000   ) {
//					if( rec_type == OligoType.OT_RIGHT) {
//	
//						last3End_Target =  String.copyValueOf(targets_lib.getSeqRevCompl(i), r.align_end_2-2, 3);
//						last3End_Primer =  String.copyValueOf(this.getOligoRevSeq(), r.align_end_1-2, 3);
//						if(     last3End_Primer.length() == 3 && last3End_Target.length()== 3 )
//						{
//							
//							if(last3End_Primer.charAt(2) == last3End_Target.charAt(2) || 
//									( last3End_Primer.charAt(1) == last3End_Target.charAt(1) || 
//										last3End_Primer.charAt(0) == last3End_Target.charAt(0))) {
//								isTargetSpecific = false;
//								// one is enough so return
//								return;
//							}
//						}
//					}
//					else if ( rec_type == OligoType.OT_LEFT) {
//	
//						last3End_Target =  String.copyValueOf(targets_lib.getSeq(i), r.align_end_2-2, 3);
//						last3End_Primer =  String.copyValueOf(this.getOligoSeq(), r.align_end_1-2, 3);
//						if(     last3End_Primer.length() == 3 && last3End_Target.length()== 3 )
//						{
//							
//							if(last3End_Primer.charAt(0) == last3End_Target.charAt(0) ||
//									( last3End_Primer.charAt(1) == last3End_Target.charAt(1) || 
//										last3End_Primer.charAt(2) == last3End_Target.charAt(2))) {
//								isTargetSpecific = false;
//								// one is enough so return
//								return;
//							}
//						}
//					}
//					
//				}
			}
				
		}
		
	}



	int findSubSeq(char[] seq, char[] oligoSeq2, int maxMisMatch) {
		int subLen = oligoSeq2.length;
		int numMisMatch = 0;
		int offset3End = 0;
		if(rec_type != OligoType.OT_RIGHT)
		{
			offset3End = subLen - 4;
		}
		// last one should differ
		
		for(int i = 0 ; i < seq.length - subLen;i++)
		{
			numMisMatch = 0;
			for(int j =0 ; j < subLen ; j++ )
			{
				if(seq[i+j] != oligoSeq2[j])
				{
					numMisMatch++;
					if(numMisMatch > maxMisMatch)
						break;
				}
				
			}
			if(numMisMatch == 0) 
			{
				
				return i;
			}
			// not yet it could 
			// compare last 3 or first three
//			if(seq[i] ==  )
			
		}
		return -1;
	}

	public boolean getIsTargetSpecific()
	{
		return isTargetSpecific;
	}
	
	@Override
	public String toString() {
		return  String.copyValueOf(oligoSeq) + " Q: " + this.quality +  " " + p3_get_ol_problem_string();
	}

	public static int compare(PrimerRecord a1, PrimerRecord a2) {
		if(a1.quality < a2.quality) return -1;
		if (a1.quality > a2.quality) return 1;

		/*
		 * We want primer_rec_comp to always return a non-0 result, because
		 * different implementations of qsort differ in how they treat "equal"
		 * elements, making it difficult to compare test output on different
		 * systems.
		 */
		if(a1.start > a2.start) return -1;
		if(a1.start < a2.start) return 1;
		if(a1.length < a2.length) return -1;
		return 1;
	}
	
}