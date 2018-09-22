package org.primer3.libprimer3;

import java.util.ArrayList;
import java.util.List;

import org.primer3.masker.input_sequence;
import org.primer3.masker.masker;
import org.primer3.masker.masking_direction;
import org.primer3.masker.output_sequence;
import org.primer3.sequence.Sequence;


/*
 * Arguments relating to a single particular source sequence (for which
 * we will pick primer(s), etc.
 */
public class SeqArgs {

	/* The net next 3 slots are presented as
	 * indexes within the sequence slot, but
	 * they are recalculated to be indexes
	 * within trimmed_seq (i.e. within the
	 * "included region").
	 */

	private IntervalList targetRegions = new IntervalList(); /* The targets.  tar2.pairs[i][0] is the start
	 * of the ith target, tar2.pairs[i][1] its length.  */

	private IntervalList excludedRegions = new IntervalList();/* The number of excluded regions. */

	private IntervalList excludedInternalRegions = new IntervalList(); 
	/* Number of excluded regions for internal
	           oligo; similar to excl2.*/

	private PairIntervalList okRegions = new PairIntervalList();

	private List<Integer> primerOverlapJunctionsList = new ArrayList<Integer>();// int[LibPrimer3.PR_MAX_INTERVAL_ARRAY]; 
	/* List of overlap junction positions. */

	private int primer_overlap_junctions_count;

	/**
	 *  The 0-based start of included region. 
	 */
	private int includedRegionStart;             

	/**
	 * TRIMMED_SEQ_LEN
	 * The length of the included region, which is
	 * also the length of the trimmed_seq field.
	 */
	private int includedRegionLength;             
	/* Index of first base of the start codon. */
	private int  startCodonPos;   


	private int[] sequenceQuality;             /* Vector of quality scores. */
//	int  n_quality;            /* Number of valid elements in 'quality' */
//	int  quality_storage_size; /* Amount of storage quality points to. */

	private char[] sequence;         /* The template sequence itself as input, 
	           not trimmed, not up-cased. */
	
	private String sequenceName;    /* An identifier for the sequence. */
	/* Another identifer for the sequence. */
	private char[] sequence_file;    
	/* The included region only, _UPCASED_. */
	private char[] trimmedSequence;      

	/* Element add by T. Koressaar support lowercase masking: */
	/**
	 *  Trimmed version of the original, mixed-case sequence. 
	 */
	private char[] trimmedOrigSequence; 

	/* Added by M. Lepamets */
	/** Masked version of the trimmed seq */
	private char[] trimmedMaskedSeq;
	/** Masked version of the other strand of the trimmed seq */
	private char[] trimmedMaskedSeqRev;                            

	
	/**
	 *  Upper case version of sequence (_not_ trimmed). 
     */
	private char[] upcasedSeq;      

	/**
	 *  Upper case version of sequence, 
	 *  other strand (_not_ trimmed). 
	 */
	private char[] upcasedSeqRev;    

	
	/**
	 *  A left primer to check or design around. 
	 */
	private char[] leftInput;       

	/**
	 *  A right primer to check or design around. 
	 */
	private char[] rightInput;      

	/**
	 *  An internal oligo to check or design around. 
	 */
	private char[] internalInput;   

	private int forceLeftStart;   /* The 0-based forced 5' start left primer. */
	private int forceLeftEnd;     /* The 0-based forced 3' end left primer. */
	private int forceRightStart;  /* The 0-based forced 5' start right primer. */
	private int forceRightEnd;    /* The 0-based forced 3' end right primer. */


	/** 
	 * default init from create_seq_arg()
	 */
	public SeqArgs()
	{
		this.startCodonPos = LibPrimer3.PR_DEFAULT_START_CODON_POS;
		this.includedRegionLength = -1; /* Indicates logical NULL. */

		this.forceLeftStart = LibPrimer3.PR_NULL_FORCE_POSITION; /* Indicates logical NULL. */
		this.forceLeftEnd = LibPrimer3.PR_NULL_FORCE_POSITION; /* Indicates logical NULL. */
		this.forceRightStart = LibPrimer3.PR_NULL_FORCE_POSITION; /* Indicates logical NULL. */
		this.forceRightEnd = LibPrimer3.PR_NULL_FORCE_POSITION; /* Indicates logical NULL. */
		this.primer_overlap_junctions_count = 0;

//		this.n_quality = 0;
		this.sequenceQuality = null;
	}
	
	
	static public SeqArgs create_seq_arg() 
	{
		SeqArgs r = new SeqArgs();

//		r.start_codon_pos = LibPrimer3.PR_DEFAULT_START_CODON_POS;
//		r.incl_l = -1; /* Indicates logical NULL. */
//
//		r.force_left_start = LibPrimer3.PR_NULL_FORCE_POSITION; /* Indicates logical NULL. */
//		r.force_left_end = LibPrimer3.PR_NULL_FORCE_POSITION; /* Indicates logical NULL. */
//		r.force_right_start = LibPrimer3.PR_NULL_FORCE_POSITION; /* Indicates logical NULL. */
//		r.force_right_end = LibPrimer3.PR_NULL_FORCE_POSITION; /* Indicates logical NULL. */
//		r.primer_overlap_junctions_count = 0;
//
//		r.n_quality = 0;
//		r.quality = null;

		return r;
	}


	public boolean PR_START_CODON_POS_IS_NULL() {
		return this.startCodonPos <= LibPrimer3.PR_NULL_START_CODON_POS;
	}


	/**
	 * SEQUENCE_ID
	 * @param datum
	 */
	public void setSequenceName(String datum) {
		this.sequenceName = datum;
	}

	/**
	 * @return the sequenceName
	 */
	public String getSequenceName() {
		return sequenceName;
	}


	/**
	 * SEQUENCE_TEMPLATE
	 * @param datum
	 */
	public void setSequence(String sequence) {
		this.sequence = sequence.toCharArray();		
	}

	/**
	 * SEQUENCE_QUALITY
	 * @param datum
	 * @return true if parsed with no errors
	 */
	public boolean set_n_quality(String datum) {

		String[] qsStr =datum.trim().split(" ");
		this.sequenceQuality  = new int[qsStr.length];
		for(int i = 0 ; i< qsStr.length;i++)
		{
			this.sequenceQuality[i] = Integer.parseInt(qsStr[i]);
		}
//		this.n_quality = quality.length;
		return true;
	}

	/**
	 * SEQUENCE_PRIMER
	 * @param datum
	 */
	public void p3_set_sa_left_input(String datum) {
		this.leftInput = datum.toCharArray();		
	}

	/**
	 * SEQUENCE_PRIMER_REVCOMP
	 * @param datum
	 */
	public void p3_set_sa_right_input(String datum) {
		this.rightInput = datum.toCharArray();		
	}

	/**
	 * SEQUENCE_INTERNAL_OLIGO
	 * @param datum
	 */
	public void p3_set_sa_internal_input(String datum) {
		this.internalInput = datum.toCharArray();			
	}

	/**
	 * SEQUENCE_PRIMER_PAIR_OK_REGION_LIST
	 * @param datum
	 */
	public void p3_set_sa_ok_regions(String datum) {
		this.okRegions = new PairIntervalList();
		String numSep = ",";
		String intervalSep = ";";

		String[] intervalStrs = datum.split(intervalSep);
		for(String intervalStr : intervalStrs){
			String[] intervalNums = intervalStr.split(numSep);
			
			
			
			int i1,i2,i3,i4;
			i1 = i2 = i3 = i4 = -1;
			if(intervalNums.length > 0 && !intervalNums[0].trim().isEmpty())
				i1 = Integer.parseInt(intervalNums[0].trim());
			if(intervalNums.length > 1 && !intervalNums[1].trim().isEmpty())
				i2 = Integer.parseInt(intervalNums[1].trim());
			if(intervalNums.length > 2 && !intervalNums[2].trim().isEmpty())
				i3 = Integer.parseInt(intervalNums[2].trim());
			if(intervalNums.length > 3 && !intervalNums[3].trim().isEmpty())
				i4 = Integer.parseInt(intervalNums[3].trim());
			this.okRegions.addIntervalPair(i1, i2, i3, i4);
		}


	}
	/**
	 * SEQUENCE_TARGET
	 * @param datum
	 */
	public void p3_set_sa_tar2(String datum) {
//		this.tar2 = IntervalArrayT2.append_interval(datum);	
		this.targetRegions.append_interval(datum);	

	}



	/**
	 * SEQUENCE_EXCLUDED_REGION
	 * @param datum
	 */
	public void p3_set_sa_excl2(String datum) {
//		this.excl2 = IntervalArrayT2.append_interval(datum);	
		this.excludedRegions.append_interval(datum);	

	}
	/**
	 * SEQUENCE_INTERNAL_EXCLUDED_REGION
	 * @param datum
	 */
	public void p3_set_sa_excl_internal2(String datum) {
//		this.excl_internal2 = IntervalArrayT2.append_interval(datum);	
		this.excludedInternalRegions.append_interval(datum);	


	}
	/**
	 * SEQUENCE_OVERLAP_JUNCTION_LIST
	 * @param datum
	 * @return 
	 */
	public boolean p3_set_sa_primer_overlap_junctions(String datum) {
		//		parse_intron_list(datum, this.primer_overlap_junctions, 
		//			      &this.primer_overlap_junctions_count)
		String intervalSep = " ";
		String[] intervalStrs = datum.split(intervalSep);
		this. primer_overlap_junctions_count = 0;
		for(String item : intervalStrs) 
		{
//			this.primerOverlapJunctionsList[this. primer_overlap_junctions_count] = Integer.parseInt(item);
			this.primerOverlapJunctionsList.add(Integer.parseInt(item));
			this. primer_overlap_junctions_count++;
		}
		return true;
	}
	/**
	 * SEQUENCE_INCLUDED_REGION
	 * @param datum
	 * @return true if no parse erro
	 */
	public boolean p3_set_sa_incl_sl(String datum) {
		String sep = ",";

		String[] pairs = datum.split(sep);
		this.includedRegionStart = Integer.parseInt(pairs[0]);
		this.includedRegionLength = Integer.parseInt(pairs[1]);
		return true;
	}
	/**
	 * SEQUENCE_START_CODON_POSITION
	 * @param datum
	 */
	public void p3_set_sa_start_codon_pos(String datum) {
		this.startCodonPos = Integer.parseInt(datum);

	}

	/**
	 * SEQUENCE_FORCE_LEFT_START
	 * @param datum
	 */
	public void p3_set_sa_force_left_start(String datum) {
		this.forceLeftStart = Integer.parseInt(datum);

	}

	/**
	 * SEQUENCE_FORCE_LEFT_END
	 * @param datum
	 */
	public void p3_set_sa_force_left_end(String datum) {
		this.forceLeftEnd = Integer.parseInt(datum);

	}

	/**
	 * SEQUENCE_FORCE_RIGHT_START
	 * @param datum
	 */
	public void p3_set_sa_force_right_start(String datum) {
		this.forceRightStart = Integer.parseInt(datum);

	}

	/**
	 * SEQUENCE_FORCE_RIGHT_END
	 * @param datum
	 */
	public void p3_set_sa_force_right_end(String datum) {
		this.forceRightEnd = Integer.parseInt(datum);

	}


	public int fake_a_sequence(P3GlobalSettings pa) {
		int i, product_size, space, ns_to_fill;
		char[] rev =null;
		int ns_to_fill_first, ns_to_fill_second;

		/* Determine the product size */
		if ( pa.getProductOptSize() == LibPrimer3.PR_UNDEFINED_INT_OPT){
			product_size = pa.getProductSizeRange(0).getRight() - pa.getProductSizeRange(0).getLeft();
		} else {
			product_size = pa.getProductOptSize();
		} 

		space = product_size + 1;
		ns_to_fill = product_size;


		/* Calculate how many Ns have to be added */
		if (this.leftInput != null ){
			ns_to_fill = ns_to_fill - this.leftInput.length;
		}
		if (this.rightInput != null){
			ns_to_fill = ns_to_fill - this.rightInput.length;
			rev = Sequence.p3_reverse_complement(this.rightInput);
		}
		if (this.internalInput != null){
			ns_to_fill = ns_to_fill - this.internalInput.length;
		}
		/* Return if there are no primers provided */
		if (ns_to_fill == product_size + 1){
			return 0;
		}
		ns_to_fill_first = ns_to_fill / 2;
		ns_to_fill_second = ns_to_fill - ns_to_fill_first;
		/* Allocate the space for the sequence */  
		this.sequence = new char[space];

		/* Copy over the primers */
		int strcatPos = 0;

		if (this.leftInput != null){
			for (i = 0; i < this.leftInput.length; i++ ) {
				this.sequence[strcatPos] =this.leftInput[i];
				strcatPos++;
			}

		}

		/* Add the Ns*/
		for (i = 0; i < ns_to_fill_first; i++ ) {
			this.sequence[strcatPos] = 'N';
			strcatPos++;
		}
		if (this.internalInput != null){
			//			  strcat(this.sequence, this.internal_input);
			for (i = 0; i < this.internalInput.length; i++ ) {
				this.sequence[strcatPos] =this.internalInput[i];
				strcatPos++;
			}

		}
		for (i = 0; i < ns_to_fill_second; i++ ) {
			//		    strcat(this.sequence, "N\0");
			this.sequence[strcatPos] = 'N';
			strcatPos++;
		}
		if (this.rightInput != null ){
			//		    strcat(this.sequence, rev);
			for (i = 0; i < rev.length; i++ ) {
				this.sequence[strcatPos] = rev[i];
				strcatPos++;
			}
		}
		return 0;		
	}



	/**
	 * This function uses the max/min product size info and the max/min
	 * oligo length in order to reduce the ranges of the ok regions. On
	 * some imputs this improves speed dramatically.
	 */
//	public void _optimize_ok_regions_list(P3GlobalSettings pa) {
//		/* We do this only if we enabled the optimization and
//		 * the primers were NOT specified. */
//		if (!LibPrimer3.OPTIMIZE_OK_REGIONS || (this.left_input != null) || (this.right_input!= null)) {
//			return;
//		}
//
//		/* If any pair is allowed, no point in doing this */
//		if (this.ok_regions.any_pair) {
//			return;
//		}
//
//		int pmin = Integer.MAX_VALUE;
//		int pmax = 0;
//		int omin = pa.primersArgs.getMinSize();
//		int omax = pa.primersArgs.getMaxSize();
//
//		/* Determine min/max product size */
//		for (int i=0; i<pa.getProductSizeRangesNumber(); i++) {
//			if (pa.getProductSizeRange(i).getLeft() < pmin) { 
//				pmin = pa.getProductSizeRange(i).getLeft(); 
//			}
//			if (pa.getProductSizeRange(i).getRight() > pmax) { 
//				pmax = pa.getProductSizeRange(i).getRight(); 
//			}
//		}
//
//		/* Update each region */
//		for (int i=0; i<this.ok_regions.count; i++) {
//			int ls = -1, le = -1, rs = -1, re = -1;
//			int new_ls = -1, new_le = -1, new_rs = -1, new_re = -1;
//			if (this.ok_regions.left_pairs[i][0] != -1) {
//				ls = this.ok_regions.left_pairs[i][0];
//				le = this.ok_regions.left_pairs[i][0] + this.ok_regions.left_pairs[i][1] - 1;
//			}
//			if (this.ok_regions.right_pairs[i][0] != -1) {
//				rs = this.ok_regions.right_pairs[i][0];
//				re = this.ok_regions.right_pairs[i][0]	+ this.ok_regions.right_pairs[i][1] - 1;
//			}
//			/* Compute new right region based on left range and min/max values
//		       of product size and oligo length */
//			if (ls != -1) {
//				new_rs = ls + pmin - omax - 1; /* -1 just to be safe */
//				new_re = le - omin + pmax + 1; /* +1 just to be safe */
//				/* Adjust the ranges */
//				if ((rs == -1) || (new_rs > rs)) { 
//					rs = new_rs; 
//				}
//				if ((re == -1) || (new_re < re)) { 
//					re = new_re; 
//				}
//				if (rs < 0) { 
//					rs = 0; 
//				}
//				if (re > (this.sequence.length)) { 
//					re = (this.sequence.length); 
//				}
//			}
//			/* Compute new left region based on right range and min/max values
//		       of product size and oligo length */
//			if (rs != -1) {
//				new_ls = rs + omin - pmax - 1; /* -1 just to be safe */
//				new_le = re - pmin + omax + 1; /* +1 just to be safe */
//				/* Adjust the ranges */
//				if ((ls == -1) || (new_ls > ls)) { 
//					ls = new_ls; 
//				}
//				if ((le == -1) || (new_le < le)) { 
//					le = new_le; 
//				}
//				if (ls < 0) { ls = 0; }
//				if (le > (this.sequence.length)) { 
//					le = (this.sequence.length); 
//				}
//			}
//			/* Temporary testing fSystem.out.format: */
//			/* fSystem.out.format(stderr, "Adjusted range [%d,%d,%d,%d] to [%d,%d,%d,%d],
//			    pmin is %d, pmax is %d, omin is %d, omax is %d\n",
//			    this.ok_regions.left_pairs[i][0],
//			    this.ok_regions.left_pairs[i][0] +
//			    this.ok_regions.left_pairs[i][1] - 1,
//			    this.ok_regions.right_pairs[i][0],
//			    this.ok_regions.right_pairs[i][0] +
//			    this.ok_regions.right_pairs[i][1] - 1, ls, le, rs, re,
//			    pmin, pmax, omin, omax);
//			 */
//			this.ok_regions.left_pairs[i][0]  = ls;
//			this.ok_regions.left_pairs[i][1]  = le - ls + 1;
//			this.ok_regions.right_pairs[i][0] = rs;
//			this.ok_regions.right_pairs[i][1] = re - rs + 1;
//		}
//		/* any_left and any_right not true anymore */
//		this.ok_regions.any_left = false;
//		this.ok_regions.any_right = false;
//	}


	public boolean _check_and_adjust_intervals(int seq_len,
			int first_index, StringBuilder nonfatal_err,
			StringBuilder warning) {

		if (this.targetRegions.checkAndAdjustInterval("TARGET",
				seq_len, first_index, nonfatal_err, this, warning, false)) 
			return true;
		this.startCodonPos -= this.includedRegionStart;
		if ( this.excludedRegions.checkAndAdjustInterval("EXCLUDED_REGION",
				seq_len, first_index, 
				nonfatal_err, this, warning, false)) return true;

		if (this.excludedInternalRegions.checkAndAdjustInterval("PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION",
				seq_len,
				first_index,
				nonfatal_err, this, warning, false)) 
			return true;
		if (IntervalList.checkAndAdjustInterval("PRIMER_PAIR_OK_REGION_LIST",
				this.okRegions.getCount(), 
				this.okRegions.getLeftPairs(),
				seq_len,
				first_index,
				nonfatal_err, this, warning, true)) 
			return true;
		if (IntervalList.checkAndAdjustInterval("PRIMER_PAIR_OK_REGION_LIST",
				this.okRegions.getCount(), 
				this.okRegions.getRightPairs(),
				seq_len,
				first_index,
				nonfatal_err, this, warning, true)) 
			return true;



		return false;
	}


	/**
	 * call _check_and_adjust_overlap_pos with primer_overlap_junctions
	 * @param sa
	 * @param tag
	 * @param seq_len
	 * @param first_index
	 * @param nonfatal_err
	 * @param warning
	 * @return
	 */
	public boolean _check_and_adjust_overlap_pos( String tag, int seq_len, int first_index, StringBuilder  nonfatal_err, StringBuilder warning) {
		return _check_and_adjust_overlap_pos( this, 
				this.primerOverlapJunctionsList, this.primer_overlap_junctions_count,  tag,
				seq_len, first_index, nonfatal_err, warning); 
	}





	/**
	 *  Function to set the included region and fix the start positions 
	 */
	public void _adjust_seq_args(
			P3GlobalSettings pa,
			P3RetVal retval) {



		SeqArgs sa = this;
		StringBuilder nonfatal_err = retval.per_sequence_err;
		StringBuilder warning = retval .warnings;

		int seq_len, inc_len;

		/* Create a seq for check primers if needed */
		if (pa.getPrimerTask() == P3Task.CHECK_PRIMERS) {
			if (null == sa.sequence) {
				sa.fake_a_sequence( pa);
			}
		}
		if (pa.getPrimerTask() == P3Task.PICK_SEQUENCING_PRIMERS && sa.includedRegionLength != -1) {
			nonfatal_err.append("Task pick_sequencing_primers cannot be combined with included region");
			return;
		}

		/*
		     Complain if there is no sequence; We need to check this
		     error here, because this function cannot do its work if
		     sa.sequence == NULL
		 */
		if (null == sa.sequence) {
			if (pa.getPrimerTask() == P3Task.CHECK_PRIMERS) {
				nonfatal_err.append( "No primers provided");
			} else {
				nonfatal_err.append( "Missing SEQUENCE tag");
			}
			return;
		}

		seq_len =sa.sequence.length;

		/* For pick_cloning_primers set the forced positions */
		if (pa.getPrimerTask() == P3Task.PICK_CLONING_PRIMERS) {
			if(sa.includedRegionLength == -1) {
				sa.forceLeftStart = pa.getFirstBaseIndex();
				sa.forceRightStart = seq_len + pa.getFirstBaseIndex() - 1;
			} else {
				sa.forceLeftStart = sa.includedRegionStart;
				sa.forceRightStart = sa.includedRegionStart + sa.includedRegionLength - 1;
			}
		}

		/* For pick_discriminative_primers set the forced positions */
		if (pa.getPrimerTask() == P3Task.PICK_DISCRIMINATIVE_PRIMERS) {
			/* Changed here from incl_s and incl_l to sa.tar2.pairs[0][0/1] */
			if (sa.targetRegions.getCount() != 1) {
				nonfatal_err.append("Task pick_discriminative_primers requires exactly one SEQUENCE_TARGET");
			}
			sa.forceLeftEnd = sa.targetRegions.getInterval(0)[0] - 1;
			sa.forceRightEnd = sa.targetRegions.getInterval(0)[0] + sa.targetRegions.getInterval(0)[1];
		}

		/* If no included region is specified,
		 * use the whole sequence as included region */
		if (sa.includedRegionLength == -1) {
			sa.includedRegionLength = seq_len;
			sa.includedRegionStart = pa.getFirstBaseIndex();
		}

		/* Generate at least one target */
		if (pa.getPrimerTask() == P3Task.PICK_SEQUENCING_PRIMERS && sa.targetRegions.getCount() == 0) {
			
			sa.targetRegions.addInterval(pa.getFirstBaseIndex(),seq_len);
			//sa.tar2.pairs[0][0] = pa.getFirstBaseIndex();
			//sa.tar2.pairs[0][1] = seq_len;
			//sa.tar2.count = 1;
		}

		/* Fix the start of the included region and start codon */
		sa.includedRegionStart -= pa.getFirstBaseIndex();
		sa.startCodonPos -= pa.getFirstBaseIndex();

		/* Fix the start */
		sa.forceLeftStart -= pa.getFirstBaseIndex();
		sa.forceLeftEnd -= pa.getFirstBaseIndex();
		sa.forceRightStart -= pa.getFirstBaseIndex();
		sa.forceRightEnd -= pa.getFirstBaseIndex();

		/* Make it relative to included region */
		sa.forceLeftStart -= sa.includedRegionStart;
		sa.forceLeftEnd -= sa.includedRegionStart;
		sa.forceRightStart -= sa.includedRegionStart;
		sa.forceRightEnd -= sa.includedRegionStart;

		inc_len = sa.includedRegionStart + sa.includedRegionLength - 1;

		if ((sa.includedRegionLength < Integer.MAX_VALUE) && (sa.includedRegionStart > -1)
				&& (sa.includedRegionLength > -1) && (inc_len < seq_len) ) {
			/* Copies inluded region into trimmed_seq */
			sa.trimmedSequence = Sequence._pr_substr(sa.sequence, sa.includedRegionStart, sa.includedRegionLength);

			/* Copies inluded region into trimmed_orig_seq */
			/* edited by T. Koressaar for lowercase masking */
			//		    sa.trimmed_orig_seq = (char *) pr_safe_malloc(sa.incl_l + 1);
			sa.trimmedOrigSequence = Sequence._pr_substr(sa.sequence, sa.includedRegionStart, sa.includedRegionLength);

			/* Masks original trimmed sequence */
			/* edited by M. Lepamets */
			if (pa.isMaskTemplate() && (pa.isPickLeftPrimer()  && pa.isPickRightPrimer() )) {
				input_sequence input_seq = new input_sequence(sa.trimmedOrigSequence);
				output_sequence output_seq =null;

				//		       input_seq = input_sequence.create_input_sequence_from_string (sa.trimmed_orig_seq, nonfatal_err);
				output_seq =  output_sequence.create_output_sequence (sa.includedRegionLength, pa.getMaskingParameters().mdir);

				masker.read_and_mask_sequence(input_seq, output_seq, pa.getMaskingParameters(), false);
				if (output_seq.sequence != null) {
					if (pa.getMaskingParameters().mdir == masking_direction.fwd) {
						sa.trimmedMaskedSeq = output_seq.sequence;
					} else if (pa.getMaskingParameters().mdir == masking_direction.rev) {
						sa.trimmedMaskedSeqRev = output_seq.sequence;
					}
				} else {
					sa.trimmedMaskedSeq = output_seq.sequence_fwd;
					sa.trimmedMaskedSeqRev = output_seq.sequence_rev;
				}
				//		       delete_input_sequence (input_seq);
				//		       delete_output_sequence (output_seq);
			}

			/* Copies the whole sequence into upcased_seq */
			sa.upcasedSeq = Sequence._pr_substr(sa.sequence, 0, sa.sequence.length);
			LibPrimer3.dna_to_upper(sa.upcasedSeq, 1);
			/* We do not need to check for illegal characters in the return
		       from dna_to_upper(), because errors are checked in
		       _pr_data_control sa.trimmed_seq. */

			/* Copies the reverse complement of the whole sequence into upcased_seq_r */
			sa.upcasedSeqRev = Sequence.p3_reverse_complement(sa.upcasedSeq);
		}

		if (sa._check_and_adjust_intervals( seq_len, pa.getFirstBaseIndex(), nonfatal_err, warning)) {
			return;
		}

		if (sa._check_and_adjust_overlap_pos(

				"SEQUENCE_OVERLAP_JUNCTION_LIST",
				seq_len,
				pa.getFirstBaseIndex(),
				nonfatal_err, warning)) {
			return;
		}

		/* Update ok regions, if non empty */
		if (sa.okRegions.getCount() > 0) {
			sa.okRegions.optimizeOkRegionsList(pa, sa);
//			sa._optimize_ok_regions_list(pa);
		}
	}





	static boolean _check_and_adjust_overlap_pos(SeqArgs sa,
			List<Integer> list,
			int count, 
			String tag,
			int seq_len,
			int first_index,
			StringBuilder  nonfatal_err,
			StringBuilder warning) {
		int i;
		boolean outside_warning_issued = false;

		if (count == 0) {
			return false;
		}

		for (i = 0; i < count; i++) {
			/* Subtract first_index from the "must overlap" positions */
			int iThElem = list.get(i);
			iThElem -= first_index;
			list.set(i,iThElem);
			
			if (iThElem >= seq_len) {
				nonfatal_err.append(tag + " beyond end of sequence");
				return true;
			}

			if (iThElem < 0) {
				nonfatal_err.append( "Negative " +tag + " length");
				return true;
			}

			/* Cause the intron positions to be relative to the included region. */
			iThElem -= sa.includedRegionStart;
			list.set(i,iThElem);
			/* Check that intron positions are within the included region. */
			if (iThElem < 0 
					|| iThElem > sa.includedRegionLength) {
				if (!outside_warning_issued) {
					warning.append(tag +  " outside of INCLUDED_REGION");
					outside_warning_issued = true;
				}
			}
		}

		return false;	
	}


	public void p3_print_args() {
		
		
		SeqArgs  s =  this;
		System.out.format("=============\n");
		System.out.format("BEGIN SEQUENCE ARGS\n") ;
		/* TODO: complete the statments for dumping this data
	       System.out.format("interval_array_t2 tar2 %s\n",
	       int pairs[PR_MAX_INTERVAL_ARRAY][2]) ;
	       int count) ;
	       System.out.format("interval_array_t2 excl2 %s\n",
	       int pairs[PR_MAX_INTERVAL_ARRAY][2]) ;
	       int count) ;
	       System.out.format("interval_array_t2 excl_internal2 %s\n",
	       int pairs[PR_MAX_INTERVAL_ARRAY][2]) ;
	       int count) ;
	       System.out.format("ok_regions \n");
		 */

		if (s.primer_overlap_junctions_count > 0) {
			System.out.format("primer_overlap_junctions_count %s\n",
					s.primer_overlap_junctions_count);
			System.out.format("primer_overlap_junctions_list [\n");
			for (int i = 0; i < s.primer_overlap_junctions_count; i++) {
				System.out.format("   %s\n", s.primerOverlapJunctionsList.get(i));
			}
			System.out.format("]\n");
		}

		System.out.format("incl_s %s\n", s.includedRegionStart) ;
		System.out.format("incl_l %s\n", s.includedRegionLength) ;
		System.out.format("start_codon_pos %s\n", s.startCodonPos) ;
		System.out.format("n_quality %s\n", s.sequenceQuality == null ? 0 : s.sequenceQuality.length) ;
		/* TO DO System.out.format("quality%s\", s.quality) ; */
		System.out.format("quality_storage_size %s\n", s.sequenceQuality == null ? 0 : s.sequenceQuality.length) ;
		System.out.format("*sequence %s\n", string(s.sequence)) ;
		System.out.format("*sequence_name %s\n", s.sequenceName) ;
		System.out.format("*sequence_file %s\n",  string(s.sequence_file)) ;
		System.out.format("*trimmed_seq %s\n", string(s.trimmedSequence)) ;
		System.out.format("*trimmed_orig_seq %s\n", string(s.trimmedOrigSequence)) ;
		System.out.format("*trimmed_masked_seq %s\n", string(s.trimmedMaskedSeq)) ;
		System.out.format("*trimmed_masked_seq_r %s\n", string(s.trimmedMaskedSeqRev)) ;
		System.out.format("*upcased_seq %s\n", string(s.upcasedSeq)) ;
		System.out.format("*upcased_seq_r %s\n", string(s.upcasedSeqRev)) ;
		System.out.format("*left_input %s\n", string(s.leftInput)) ;
		System.out.format("*right_input %s\n", string(s.rightInput)) ;
		System.out.format("*internal_input %s\n", s.internalInput) ;
		System.out.format("force_left_start %s\n", s.forceLeftStart) ;
		System.out.format("force_left_end %s\n", s.forceLeftEnd) ;
		System.out.format("force_right_start %s\n", s.forceRightStart) ;
		System.out.format("force_right_end %s\n", s.forceRightEnd) ;
		System.out.format("END SEQUENCE ARGS\n") ;
		System.out.format("=============\n");
		System.out.format("\n");
	}


	
	
	static String string(char[] str) {
		if(str == null)
			return "NULL";
		return new String(str);
	}


	/**
	 * @return the targetRegions
	 */
	public IntervalList getTargetRegions() {
		return targetRegions;
	}


	/**
	 * @param targetRegions the targetRegions to set
	 */
	public void setTargetRegions(IntervalList targetRegions) {
		this.targetRegions = targetRegions;
	}


	/**
	 * @return the excludedRegions
	 */
	public IntervalList getExcludedRegions() {
		return excludedRegions;
	}


	/**
	 * @param excludedRegions the excludedRegions to set
	 */
	public void setExcludedRegions(IntervalList excludedRegions) {
		this.excludedRegions = excludedRegions;
	}


	/**
	 * @return the excludedInternalRegions
	 */
	public IntervalList getExcludedInternalRegions() {
		return excludedInternalRegions;
	}


	/**
	 * @param excludedInternalRegions the excludedInternalRegions to set
	 */
	public void setExcludedInternalRegions(IntervalList excludedInternalRegions) {
		this.excludedInternalRegions = excludedInternalRegions;
	}


	/**
	 * @return the okRegions
	 */
	public PairIntervalList getOkRegions() {
		return okRegions;
	}


	/**
	 * @param okRegions the okRegions to set
	 */
	public void setOkRegions(PairIntervalList okRegions) {
		this.okRegions = okRegions;
	}


	/**
	 * @return the includedRegionStart
	 */
	public int getIncludedRegionStart() {
		return includedRegionStart;
	}


	/**
	 * @param includedRegionStart the includedRegionStart to set
	 */
	public void setIncludedRegionStart(int includedRegionStart) {
		this.includedRegionStart = includedRegionStart;
	}


	/**
	 * @return the includedRegionLength
	 */
	public int getIncludedRegionLength() {
		return includedRegionLength;
	}


	/**
	 * @param includedRegionLength the includedRegionLength to set
	 */
	public void setIncludedRegionLength(int includedRegionLength) {
		this.includedRegionLength = includedRegionLength;
	}


	/**
	 * @return the sequence
	 */
	public char[] getSequence() {
		return sequence;
	}


	/**
	 * SEQUENCE_TEMPLATE
	 * @param sequence to seq
	 */
	public void setSequence(char[] sequence) {
		this.sequence = sequence;		
	}


	/**
	 * @return the startCodonPos
	 */
	public int getStartCodonPos() {
		return startCodonPos;
	}


	/**
	 * @param startCodonPos the startCodonPos to set
	 */
	public void setStartCodonPos(int startCodonPos) {
		this.startCodonPos = startCodonPos;
	}


	/**
	 * @return the sequenceQuality
	 */
	public int[] getSequenceQuality() {
		return sequenceQuality;
	}


	/**
	 * @param sequenceQuality the sequenceQuality to set
	 */
	public void setSequenceQuality(int[] sequenceQuality) {
		this.sequenceQuality = sequenceQuality;
	}


	/**
	 * @return the trimmedSequence
	 */
	public char[] getTrimmedSequence() {
		return trimmedSequence;
	}


	/**
	 * @param trimmedSequence the trimmedSequence to set
	 */
	public void setTrimmedSequence(char[] trimmedSequence) {
		this.trimmedSequence = trimmedSequence;
	}


	/**
	 * @return the trimmedOrigSequence
	 */
	public char[] getTrimmedOrigSequence() {
		return trimmedOrigSequence;
	}


	/**
	 * @param trimmedOrigSequence the trimmedOrigSequence to set
	 */
	public void setTrimmedOrigSequence(char[] trimmedOrigSequence) {
		this.trimmedOrigSequence = trimmedOrigSequence;
	}


	/**
	 * @return the primerOverlapJunctionsList
	 */
	public List<Integer> getPrimerOverlapJunctionsList() {
		return primerOverlapJunctionsList;
	}


	/**
	 * @return the upcasedSeq
	 */
	public char[] getUpcasedSeq() {
		return upcasedSeq;
	}


	/**
	 * @param upcasedSeq the upcasedSeq to set
	 */
	public void setUpcasedSeq(char[] upcasedSeq) {
		this.upcasedSeq = upcasedSeq;
	}


	/**
	 * @return the upcasedSeqRev
	 */
	public char[] getUpcasedSeqRev() {
		return upcasedSeqRev;
	}


	/**
	 * @param upcasedSeqRev the upcasedSeqRev to set
	 */
	public void setUpcasedSeqRev(char[] upcasedSeqRev) {
		this.upcasedSeqRev = upcasedSeqRev;
	}


	/**
	 * @return the leftInput
	 */
	public char[] getLeftInput() {
		return leftInput;
	}


	/**
	 * @param leftInput the leftInput to set
	 */
	public void setLeftInput(char[] leftInput) {
		this.leftInput = leftInput;
	}


	/**
	 * @return the rightInput
	 */
	public char[] getRightInput() {
		return rightInput;
	}


	/**
	 * @param rightInput the rightInput to set
	 */
	public void setRightInput(char[] rightInput) {
		this.rightInput = rightInput;
	}


	/**
	 * @return the internalInput
	 */
	public char[] getInternalInput() {
		return internalInput;
	}


	/**
	 * @param internalInput the internalInput to set
	 */
	public void setInternalInput(char[] internalInput) {
		this.internalInput = internalInput;
	}


	/**
	 * @return the forceLeftStart
	 */
	public int getForceLeftStart() {
		return forceLeftStart;
	}


	/**
	 * @param forceLeftStart the forceLeftStart to set
	 */
	public void setForceLeftStart(int forceLeftStart) {
		this.forceLeftStart = forceLeftStart;
	}


	/**
	 * @return the forceLeftEnd
	 */
	public int getForceLeftEnd() {
		return forceLeftEnd;
	}


	/**
	 * @param forceLeftEnd the forceLeftEnd to set
	 */
	public void setForceLeftEnd(int forceLeftEnd) {
		this.forceLeftEnd = forceLeftEnd;
	}


	/**
	 * @return the forceRightStart
	 */
	public int getForceRightStart() {
		return forceRightStart;
	}


	/**
	 * @param forceRightStart the forceRightStart to set
	 */
	public void setForceRightStart(int forceRightStart) {
		this.forceRightStart = forceRightStart;
	}


	/**
	 * @return the forceRightEnd
	 */
	public int getForceRightEnd() {
		return forceRightEnd;
	}


	/**
	 * @param forceRightEnd the forceRightEnd to set
	 */
	public void setForceRightEnd(int forceRightEnd) {
		this.forceRightEnd = forceRightEnd;
	}


	/**
	 * @return the trimmedMaskedSeq
	 */
	public char[] getTrimmedMaskedSeq() {
		return trimmedMaskedSeq;
	}


	/**
	 * @param trimmedMaskedSeq the trimmedMaskedSeq to set
	 */
	public void setTrimmedMaskedSeq(char[] trimmedMaskedSeq) {
		this.trimmedMaskedSeq = trimmedMaskedSeq;
	}


	/**
	 * @return the trimmedMaskedSeqRev
	 */
	public char[] getTrimmedMaskedSeqRev() {
		return trimmedMaskedSeqRev;
	}


	/**
	 * @param trimmedMaskedSeqRev the trimmedMaskedSeqRev to set
	 */
	public void setTrimmedMaskedSeqRev(char[] trimmedMaskedSeqRev) {
		this.trimmedMaskedSeqRev = trimmedMaskedSeqRev;
	}


	/**
	 * @param primerOverlapJunctionsList the primerOverlapJunctionsList to set
	 */
//	public void setPrimerOverlapJunctionsList(
//			int[] primerOverlapJunctionsList) {
//		this.primerOverlapJunctionsList = primerOverlapJunctionsList;
//	}


}