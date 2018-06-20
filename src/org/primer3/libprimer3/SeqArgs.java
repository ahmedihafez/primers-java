package org.primer3.libprimer3;

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

	IntervalList tar2 = new IntervalList(); /* The targets.  tar2.pairs[i][0] is the start
	 * of the ith target, tar2.pairs[i][1] its length.  */

	IntervalList excl2 = new IntervalList();/* The number of excluded regions. */

	IntervalList excl_internal2 = new IntervalList(); 
	/* Number of excluded regions for internal
	           oligo; similar to excl2.*/

	IntervalArrayT4 ok_regions = new IntervalArrayT4();

	int[] primer_overlap_junctions= new int[LibPrimer3.PR_MAX_INTERVAL_ARRAY]; 
	/* List of overlap junction positions. */

	int primer_overlap_junctions_count;

	public int incl_s;             /* The 0-based start of included region. */

	/**
	 * TRIMMED_SEQ_LEN
	 */
	int incl_l;             /* 
	 * The length of the included region, which is
	 * also the length of the trimmed_seq field.
	 */
	int  start_codon_pos;   /* Index of first base of the start codon. */


	int[] quality;             /* Vector of quality scores. */
	int  n_quality;            /* Number of valid elements in 'quality' */
	int  quality_storage_size; /* Amount of storage quality points to. */

	public char[] sequence;         /* The template sequence itself as input, 
	           not trimmed, not up-cased. */
	public String sequence_name;    /* An identifier for the sequence. */
	char[] sequence_file;    /* Another identifer for the sequence. */
	public char[] trimmed_seq;      /* The included region only, _UPCASED_. */

	/* Element add by T. Koressaar support lowercase masking: */
	char[] trimmed_orig_seq; /* Trimmed version of the original, mixed-case sequence. */

	/* Added by M. Lepamets */
	char[] trimmed_masked_seq; /* Masked version of the trimmed seq */
	char[] trimmed_masked_seq_r; /* Masked version of the other strand
	              of the trimmed seq */                           

	char[] upcased_seq;      /* Upper case version of sequence
	           (_not_ trimmed). */

	char[] upcased_seq_r;    /* Upper case version of sequence, 
	           other strand (_not_ trimmed). */

	public char[] left_input;       /* A left primer to check or design around. */

	public char[] right_input;      /* A right primer to check or design around. */

	public char[] internal_input;   /* An internal oligo to check or design around. */

	int force_left_start;   /* The 0-based forced 5' start left primer. */
	int force_left_end;     /* The 0-based forced 3' end left primer. */
	int force_right_start;  /* The 0-based forced 5' start right primer. */
	int force_right_end;    /* The 0-based forced 3' end right primer. */


	/** 
	 * default init from create_seq_arg()
	 */
	public SeqArgs()
	{
		this.start_codon_pos = LibPrimer3.PR_DEFAULT_START_CODON_POS;
		this.incl_l = -1; /* Indicates logical NULL. */

		this.force_left_start = LibPrimer3.PR_NULL_FORCE_POSITION; /* Indicates logical NULL. */
		this.force_left_end = LibPrimer3.PR_NULL_FORCE_POSITION; /* Indicates logical NULL. */
		this.force_right_start = LibPrimer3.PR_NULL_FORCE_POSITION; /* Indicates logical NULL. */
		this.force_right_end = LibPrimer3.PR_NULL_FORCE_POSITION; /* Indicates logical NULL. */
		this.primer_overlap_junctions_count = 0;

		this.n_quality = 0;
		this.quality = null;
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
		return this.start_codon_pos <= LibPrimer3.PR_NULL_START_CODON_POS;
	}


	/**
	 * SEQUENCE_ID
	 * @param datum
	 */
	public void p3_set_sa_sequence_name(String datum) {
		this.sequence_name = datum;
	}

	/**
	 * SEQUENCE_TEMPLATE
	 * @param datum
	 */
	public void p3_set_sa_sequence(String datum) {
		this.sequence = datum.toCharArray();		
	}

	/**
	 * SEQUENCE_QUALITY
	 * @param datum
	 * @return true if parsed with no errors
	 */
	public boolean set_n_quality(String datum) {

		String[] qsStr =datum.split(" ");
		this.quality  = new int[qsStr.length];
		for(int i = 0 ; i< qsStr.length;i++)
		{
			this.quality[i] = Integer.parseInt(qsStr[i]);
		}
		this.n_quality = quality.length;
		return true;
	}

	/**
	 * SEQUENCE_PRIMER
	 * @param datum
	 */
	public void p3_set_sa_left_input(String datum) {
		this.left_input = datum.toCharArray();		
	}

	/**
	 * SEQUENCE_PRIMER_REVCOMP
	 * @param datum
	 */
	public void p3_set_sa_right_input(String datum) {
		this.right_input = datum.toCharArray();		
	}

	/**
	 * SEQUENCE_INTERNAL_OLIGO
	 * @param datum
	 */
	public void p3_set_sa_internal_input(String datum) {
		this.internal_input = datum.toCharArray();			
	}

	/**
	 * SEQUENCE_PRIMER_PAIR_OK_REGION_LIST
	 * @param datum
	 */
	public void p3_set_sa_ok_regions(String datum) {
		this.ok_regions = new IntervalArrayT4();
		String numSep = ",";
		String intervalSep = ";";

		String[] intervalStrs = datum.split(intervalSep);
		for(String intervalStr : intervalStrs){
			String[] intervalNums = intervalStr.split(numSep);
			int i1,i2,i3,i4;
			i1 = i2 = i3 = i4 = -1;
			if(!intervalNums[0].isEmpty())
				i1 = Integer.parseInt(intervalNums[0]);
			if(!intervalNums[1].isEmpty())
				i2 = Integer.parseInt(intervalNums[1]);
			if(!intervalNums[2].isEmpty())
				i3 = Integer.parseInt(intervalNums[2]);
			if(!intervalNums[3].isEmpty())
				i4 = Integer.parseInt(intervalNums[3]);
			this.ok_regions.p3_add_to_2_interval_array(i1, i2, i3, i4);
		}


	}
	/**
	 * SEQUENCE_TARGET
	 * @param datum
	 */
	public void p3_set_sa_tar2(String datum) {
//		this.tar2 = IntervalArrayT2.append_interval(datum);	
		this.tar2.append_interval(datum);	

	}



	/**
	 * SEQUENCE_EXCLUDED_REGION
	 * @param datum
	 */
	public void p3_set_sa_excl2(String datum) {
//		this.excl2 = IntervalArrayT2.append_interval(datum);	
		this.excl2.append_interval(datum);	

	}
	/**
	 * SEQUENCE_INTERNAL_EXCLUDED_REGION
	 * @param datum
	 */
	public void p3_set_sa_excl_internal2(String datum) {
//		this.excl_internal2 = IntervalArrayT2.append_interval(datum);	
		this.excl_internal2.append_interval(datum);	


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
			this.primer_overlap_junctions[this. primer_overlap_junctions_count] = 
					Integer.parseInt(item);
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
		this.incl_s = Integer.parseInt(pairs[0]);
		this.incl_l = Integer.parseInt(pairs[1]);
		return true;
	}
	/**
	 * SEQUENCE_START_CODON_POSITION
	 * @param datum
	 */
	public void p3_set_sa_start_codon_pos(String datum) {
		this.start_codon_pos = Integer.parseInt(datum);

	}

	/**
	 * SEQUENCE_FORCE_LEFT_START
	 * @param datum
	 */
	public void p3_set_sa_force_left_start(String datum) {
		this.force_left_start = Integer.parseInt(datum);

	}

	/**
	 * SEQUENCE_FORCE_LEFT_END
	 * @param datum
	 */
	public void p3_set_sa_force_left_end(String datum) {
		this.force_left_end = Integer.parseInt(datum);

	}

	/**
	 * SEQUENCE_FORCE_RIGHT_START
	 * @param datum
	 */
	public void p3_set_sa_force_right_start(String datum) {
		this.force_right_start = Integer.parseInt(datum);

	}

	/**
	 * SEQUENCE_FORCE_RIGHT_END
	 * @param datum
	 */
	public void p3_set_sa_force_right_end(String datum) {
		this.force_right_end = Integer.parseInt(datum);

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
		if (this.left_input != null ){
			ns_to_fill = ns_to_fill - this.left_input.length;
		}
		if (this.right_input != null){
			ns_to_fill = ns_to_fill - this.right_input.length;
			rev = Sequence.p3_reverse_complement(this.right_input);
		}
		if (this.internal_input != null){
			ns_to_fill = ns_to_fill - this.internal_input.length;
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

		if (this.left_input != null){
			for (i = 0; i < this.left_input.length; i++ ) {
				this.sequence[strcatPos] =this.left_input[i];
				strcatPos++;
			}

		}

		/* Add the Ns*/
		for (i = 0; i < ns_to_fill_first; i++ ) {
			this.sequence[strcatPos] = 'N';
			strcatPos++;
		}
		if (this.internal_input != null){
			//			  strcat(this.sequence, this.internal_input);
			for (i = 0; i < this.internal_input.length; i++ ) {
				this.sequence[strcatPos] =this.internal_input[i];
				strcatPos++;
			}

		}
		for (i = 0; i < ns_to_fill_second; i++ ) {
			//		    strcat(this.sequence, "N\0");
			this.sequence[strcatPos] = 'N';
			strcatPos++;
		}
		if (this.right_input != null ){
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
	public void _optimize_ok_regions_list(P3GlobalSettings pa) {
		/* We do this only if we enabled the optimization and
		 * the primers were NOT specified. */
		if (!LibPrimer3.OPTIMIZE_OK_REGIONS || (this.left_input != null) || (this.right_input!= null)) {
			return;
		}

		/* If any pair is allowed, no point in doing this */
		if (this.ok_regions.any_pair) {
			return;
		}

		int pmin = Integer.MAX_VALUE;
		int pmax = 0;
		int omin = pa.primersArgs.getMinSize();
		int omax = pa.primersArgs.getMaxSize();

		/* Determine min/max product size */
		for (int i=0; i<pa.getProductSizeRangesNumber(); i++) {
			if (pa.getProductSizeRange(i).getLeft() < pmin) { 
				pmin = pa.getProductSizeRange(i).getLeft(); 
			}
			if (pa.getProductSizeRange(i).getRight() > pmax) { 
				pmax = pa.getProductSizeRange(i).getRight(); 
			}
		}

		/* Update each region */
		for (int i=0; i<this.ok_regions.count; i++) {
			int ls = -1, le = -1, rs = -1, re = -1;
			int new_ls = -1, new_le = -1, new_rs = -1, new_re = -1;
			if (this.ok_regions.left_pairs[i][0] != -1) {
				ls = this.ok_regions.left_pairs[i][0];
				le = this.ok_regions.left_pairs[i][0] + this.ok_regions.left_pairs[i][1] - 1;
			}
			if (this.ok_regions.right_pairs[i][0] != -1) {
				rs = this.ok_regions.right_pairs[i][0];
				re = this.ok_regions.right_pairs[i][0]	+ this.ok_regions.right_pairs[i][1] - 1;
			}
			/* Compute new right region based on left range and min/max values
		       of product size and oligo length */
			if (ls != -1) {
				new_rs = ls + pmin - omax - 1; /* -1 just to be safe */
				new_re = le - omin + pmax + 1; /* +1 just to be safe */
				/* Adjust the ranges */
				if ((rs == -1) || (new_rs > rs)) { 
					rs = new_rs; 
				}
				if ((re == -1) || (new_re < re)) { 
					re = new_re; 
				}
				if (rs < 0) { 
					rs = 0; 
				}
				if (re > (this.sequence.length)) { 
					re = (this.sequence.length); 
				}
			}
			/* Compute new left region based on right range and min/max values
		       of product size and oligo length */
			if (rs != -1) {
				new_ls = rs + omin - pmax - 1; /* -1 just to be safe */
				new_le = re - pmin + omax + 1; /* +1 just to be safe */
				/* Adjust the ranges */
				if ((ls == -1) || (new_ls > ls)) { 
					ls = new_ls; 
				}
				if ((le == -1) || (new_le < le)) { 
					le = new_le; 
				}
				if (ls < 0) { ls = 0; }
				if (le > (this.sequence.length)) { 
					le = (this.sequence.length); 
				}
			}
			/* Temporary testing fSystem.out.format: */
			/* fSystem.out.format(stderr, "Adjusted range [%d,%d,%d,%d] to [%d,%d,%d,%d],
			    pmin is %d, pmax is %d, omin is %d, omax is %d\n",
			    this.ok_regions.left_pairs[i][0],
			    this.ok_regions.left_pairs[i][0] +
			    this.ok_regions.left_pairs[i][1] - 1,
			    this.ok_regions.right_pairs[i][0],
			    this.ok_regions.right_pairs[i][0] +
			    this.ok_regions.right_pairs[i][1] - 1, ls, le, rs, re,
			    pmin, pmax, omin, omax);
			 */
			this.ok_regions.left_pairs[i][0]  = ls;
			this.ok_regions.left_pairs[i][1]  = le - ls + 1;
			this.ok_regions.right_pairs[i][0] = rs;
			this.ok_regions.right_pairs[i][1] = re - rs + 1;
		}
		/* any_left and any_right not true anymore */
		this.ok_regions.any_left = false;
		this.ok_regions.any_right = false;
	}


	public boolean _check_and_adjust_intervals(int seq_len,
			int first_index, StringBuilder nonfatal_err,
			StringBuilder warning) {

		if (this.tar2.checkAndAdjustInterval("TARGET",
				seq_len, first_index, nonfatal_err, this, warning, false)) 
			return true;
		this.start_codon_pos -= this.incl_s;
		if ( this.excl2.checkAndAdjustInterval("EXCLUDED_REGION",
				seq_len, first_index, 
				nonfatal_err, this, warning, false)) return true;

		if (this.excl_internal2.checkAndAdjustInterval("PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION",
				seq_len,
				first_index,
				nonfatal_err, this, warning, false)) 
			return true;
		if (IntervalList.checkAndAdjustInterval("PRIMER_PAIR_OK_REGION_LIST",
				this.ok_regions.count, 
				this.ok_regions.left_pairs,
				seq_len,
				first_index,
				nonfatal_err, this, warning, true)) 
			return true;
		if (IntervalList.checkAndAdjustInterval("PRIMER_PAIR_OK_REGION_LIST",
				this.ok_regions.count, 
				this.ok_regions.right_pairs,
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
				this.primer_overlap_junctions, this.primer_overlap_junctions_count,  tag,
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
		if (pa.getPrimerTask() == P3Task.PICK_SEQUENCING_PRIMERS && sa.incl_l != -1) {
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
			if(sa.incl_l == -1) {
				sa.force_left_start = pa.getFirstBaseIndex();
				sa.force_right_start = seq_len + pa.getFirstBaseIndex() - 1;
			} else {
				sa.force_left_start = sa.incl_s;
				sa.force_right_start = sa.incl_s + sa.incl_l - 1;
			}
		}

		/* For pick_discriminative_primers set the forced positions */
		if (pa.getPrimerTask() == P3Task.PICK_DISCRIMINATIVE_PRIMERS) {
			/* Changed here from incl_s and incl_l to sa.tar2.pairs[0][0/1] */
			if (sa.tar2.getCount() != 1) {
				nonfatal_err.append("Task pick_discriminative_primers requires exactly one SEQUENCE_TARGET");
			}
			sa.force_left_end = sa.tar2.getInterval(0)[0] - 1;
			sa.force_right_end = sa.tar2.getInterval(0)[0] + sa.tar2.getInterval(0)[1];
		}

		/* If no included region is specified,
		 * use the whole sequence as included region */
		if (sa.incl_l == -1) {
			sa.incl_l = seq_len;
			sa.incl_s = pa.getFirstBaseIndex();
		}

		/* Generate at least one target */
		if (pa.getPrimerTask() == P3Task.PICK_SEQUENCING_PRIMERS && sa.tar2.getCount() == 0) {
			
			sa.tar2.addInterval(pa.getFirstBaseIndex(),seq_len);
			//sa.tar2.pairs[0][0] = pa.getFirstBaseIndex();
			//sa.tar2.pairs[0][1] = seq_len;
			//sa.tar2.count = 1;
		}

		/* Fix the start of the included region and start codon */
		sa.incl_s -= pa.getFirstBaseIndex();
		sa.start_codon_pos -= pa.getFirstBaseIndex();

		/* Fix the start */
		sa.force_left_start -= pa.getFirstBaseIndex();
		sa.force_left_end -= pa.getFirstBaseIndex();
		sa.force_right_start -= pa.getFirstBaseIndex();
		sa.force_right_end -= pa.getFirstBaseIndex();

		/* Make it relative to included region */
		sa.force_left_start -= sa.incl_s;
		sa.force_left_end -= sa.incl_s;
		sa.force_right_start -= sa.incl_s;
		sa.force_right_end -= sa.incl_s;

		inc_len = sa.incl_s + sa.incl_l - 1;

		if ((sa.incl_l < Integer.MAX_VALUE) && (sa.incl_s > -1)
				&& (sa.incl_l > -1) && (inc_len < seq_len) ) {
			/* Copies inluded region into trimmed_seq */
			sa.trimmed_seq = Sequence._pr_substr(sa.sequence, sa.incl_s, sa.incl_l);

			/* Copies inluded region into trimmed_orig_seq */
			/* edited by T. Koressaar for lowercase masking */
			//		    sa.trimmed_orig_seq = (char *) pr_safe_malloc(sa.incl_l + 1);
			sa.trimmed_orig_seq = Sequence._pr_substr(sa.sequence, sa.incl_s, sa.incl_l);

			/* Masks original trimmed sequence */
			/* edited by M. Lepamets */
			if (pa.isMaskTemplate() && (pa.isPickLeftPrimer()  && pa.isPickRightPrimer() )) {
				input_sequence input_seq = new input_sequence(sa.trimmed_orig_seq);
				output_sequence output_seq =null;

				//		       input_seq = input_sequence.create_input_sequence_from_string (sa.trimmed_orig_seq, nonfatal_err);
				output_seq =  output_sequence.create_output_sequence (sa.incl_l, pa.getMaskingParameters().mdir);

				masker.read_and_mask_sequence(input_seq, output_seq, pa.getMaskingParameters(), false);
				if (output_seq.sequence != null) {
					if (pa.getMaskingParameters().mdir == masking_direction.fwd) {
						sa.trimmed_masked_seq = output_seq.sequence;
					} else if (pa.getMaskingParameters().mdir == masking_direction.rev) {
						sa.trimmed_masked_seq_r = output_seq.sequence;
					}
				} else {
					sa.trimmed_masked_seq = output_seq.sequence_fwd;
					sa.trimmed_masked_seq_r = output_seq.sequence_rev;
				}
				//		       delete_input_sequence (input_seq);
				//		       delete_output_sequence (output_seq);
			}

			/* Copies the whole sequence into upcased_seq */
			sa.upcased_seq = Sequence._pr_substr(sa.sequence, 0, sa.sequence.length);
			LibPrimer3.dna_to_upper(sa.upcased_seq, 1);
			/* We do not need to check for illegal characters in the return
		       from dna_to_upper(), because errors are checked in
		       _pr_data_control sa.trimmed_seq. */

			/* Copies the reverse complement of the whole sequence into upcased_seq_r */
			sa.upcased_seq_r = Sequence.p3_reverse_complement(sa.upcased_seq);
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
		if (sa.ok_regions.count > 0) {
			sa._optimize_ok_regions_list(pa);
		}
	}





	static boolean _check_and_adjust_overlap_pos(SeqArgs sa,
			int[] list,
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
			list[i] -= first_index;

			if (list[i] >= seq_len) {
				nonfatal_err.append(tag + " beyond end of sequence");
				return true;
			}

			if (list[i] < 0) {
				nonfatal_err.append( "Negative " +tag + " length");
				return true;
			}

			/* Cause the intron positions to be relative to the included region. */
			list[i] -= sa.incl_s;

			/* Check that intron positions are within the included region. */
			if (list[i] < 0 
					|| list[i] > sa.incl_l) {
				if (!outside_warning_issued) {
					warning.append(tag +  " outside of INCLUDED_REGION");
					outside_warning_issued = true;
				}
			}
		}

		return true;	
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
				System.out.format("   %s\n", s.primer_overlap_junctions[i]);
			}
			System.out.format("]\n");
		}

		System.out.format("incl_s %s\n", s.incl_s) ;
		System.out.format("incl_l %s\n", s.incl_l) ;
		System.out.format("start_codon_pos %s\n", s.start_codon_pos) ;
		System.out.format("n_quality %s\n", s.n_quality) ;
		/* TO DO System.out.format("quality%s\", s.quality) ; */
		System.out.format("quality_storage_size %s\n", s.quality_storage_size) ;
		System.out.format("*sequence %s\n", string(s.sequence)) ;
		System.out.format("*sequence_name %s\n", s.sequence_name) ;
		System.out.format("*sequence_file %s\n",  string(s.sequence_file)) ;
		System.out.format("*trimmed_seq %s\n", string(s.trimmed_seq)) ;
		System.out.format("*trimmed_orig_seq %s\n", string(s.trimmed_orig_seq)) ;
		System.out.format("*trimmed_masked_seq %s\n", string(s.trimmed_masked_seq)) ;
		System.out.format("*trimmed_masked_seq_r %s\n", string(s.trimmed_masked_seq_r)) ;
		System.out.format("*upcased_seq %s\n", string(s.upcased_seq)) ;
		System.out.format("*upcased_seq_r %s\n", string(s.upcased_seq_r)) ;
		System.out.format("*left_input %s\n", string(s.left_input)) ;
		System.out.format("*right_input %s\n", string(s.right_input)) ;
		System.out.format("*internal_input %s\n", s.internal_input) ;
		System.out.format("force_left_start %s\n", s.force_left_start) ;
		System.out.format("force_left_end %s\n", s.force_left_end) ;
		System.out.format("force_right_start %s\n", s.force_right_start) ;
		System.out.format("force_right_end %s\n", s.force_right_end) ;
		System.out.format("END SEQUENCE ARGS\n") ;
		System.out.format("=============\n");
		System.out.format("\n");
	}


	
	
	static String string(char[] str) {
		if(str == null)
			return "NULL";
		return new String(str);
	}


}