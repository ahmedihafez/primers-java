package org.primer3.libprimer3;

import org.primer3.primer.PrimerPair;

public class PairIntervalList {


	private int[][] left_pairs = new int[LibPrimer3.PR_MAX_INTERVAL_ARRAY][2];
	private int[][] right_pairs = new int[LibPrimer3.PR_MAX_INTERVAL_ARRAY][2];
	boolean any_left = false ;  /* set to 1 if the empty pair ",," is given for any
			    left interval */
	boolean any_right = false; /* set to 1 if the empty pair ",," is given for any
			    right interval */
	public boolean any_pair = false;  /* set to 1 if both intervals are given as empty */
	private int count = 0;     /* total number of pairs */





	public boolean addIntervalPair(int i1, int i2, int i3, int i4){
		int c = this.getCount();
		if (c >= LibPrimer3.PR_MAX_INTERVAL_ARRAY) 
			throw new  IndexOutOfBoundsException("Too many elements for tag");
		/* for a region either both values are given, or none is given */
		if (((i1 == -1) && (i2 != -1)) || ((i1 != -1) && (i2 == -1)))
			throw new  IndexOutOfBoundsException("Invalid range at tag  ");
		if (((i3 == -1) && (i4 != -1)) || ((i3 != -1) && (i4 == -1)))
			throw new  IndexOutOfBoundsException("Invalid range at tag  ");
		this.left_pairs[c][0] = i1;
		this.left_pairs[c][1] = i2;
		this.right_pairs[c][0] = i3;
		this.right_pairs[c][1] = i4;
		if ((i1 == -1) && (i2 == -1))
			this.any_left = true;
		if ((i3 == -1) && (i4 == -1))
			this.any_right = true;
		if ((i1 == -1) && (i2 == -1) && (i3 == -1) && (i4 == -1))
			this.any_pair = true;
		this.count++ ;
		return  false;
	}
	public boolean checkIncludedInAny(PrimerPair ppair) {

		boolean included = false;
		int l_start = ppair.left.start, l_end = ppair.left.start + ppair.left.length - 1;
		int r_start = ppair.right.start - ppair.right.length + 1, r_end = ppair.right.start;
		for (int i=0; i<this.getCount(); i++) {
			if (this.left_pairs[i][0] == -1) {
				/* any left primer, check right primer */
				if ((r_start >= this.right_pairs[i][0]) &&
						(r_end <= this.right_pairs[i][0] + this.right_pairs[i][1] - 1)) {
					included = true;
					break;
				}
			} else if (this.right_pairs[i][0] == -1) {
				/* any right primer, check the left primer */
				if ((l_start >= this.left_pairs[i][0]) && 
						(l_end <= this.left_pairs[i][0] + this.left_pairs[i][1] - 1)) {
					included = true;
					break;
				}
			} else {
				/* check both primers */
				if ((l_start >= this.left_pairs[i][0]) && 
						(l_end <= this.left_pairs[i][0] + this.left_pairs[i][1] - 1) &&
						(r_start >= this.right_pairs[i][0]) &&
						(r_end <= this.right_pairs[i][0] + this.right_pairs[i][1] - 1)) {
					included = true;
					break;
				}  
			}
		}

		return included;
	}

	public boolean chechIncludedInAny(OligoType otype,int j, int k)
	{
		boolean included = false;

		if(this.getCount() > 0)
		{
			if(otype == OligoType.OT_LEFT)
			{
				if(!this.any_left)
				{
					for (int i = 0; i < this.getCount(); i++) {
						if ((j >= this.left_pairs[i][0])
								&& (k <= this.left_pairs[i][0]
										+ this.left_pairs[i][1] - 1)) {
							return  true;
						}
					}
				}
				else
				{
					return true;
				}
			}
			else if (otype == OligoType.OT_RIGHT)
			{
				if(!this.any_right)
				{
					for (int i = 0; i < this.getCount(); i++) {
						if ((j >= this.right_pairs[i][0])
								&& (k <= this.right_pairs[i][0]
										+ this.right_pairs[i][1] - 1)) {
							return true;
						}
					}
				}
				else
				{
					return true;
				}
			}

		}
		else
		{
			// there is no any restrictions
			return true;
		}
		//		 if ((otype == OligoType.OT_LEFT) && (this.count > 0)
		//					&& (!this.any_left)) {
		//				boolean included = false;
		//				for (int i = 0; i < this.count; i++) {
		//					if ((j >= this.left_pairs[i][0])
		//							&& (k <= this.left_pairs[i][0]
		//									+ this.left_pairs[i][1] - 1)) {
		//						included = true;
		//						break;
		//					}
		//				}
		//			}
		//			if ((otype == OligoType.OT_RIGHT) && (sa.ok_regions.count > 0)
		//					&& (!sa.ok_regions.any_right)) {
		//				boolean included = false;
		//				for (i = 0; i < sa.ok_regions.count; i++) {
		//					if ((j >= sa.ok_regions.right_pairs[i][0])
		//							&& (k <= sa.ok_regions.right_pairs[i][0]
		//									+ sa.ok_regions.right_pairs[i][1] - 1)) {
		//						included = true;
		//						break;
		//					}
		//				}
		//				if (!included) {
		//					op_set_not_in_any_ok_region();
		//					stats.not_in_any_right_ok_region++;
		//					if (!must_use)
		//						return;
		//				}
		//			}

		return included;
	}


	/**
	 * This function uses the max/min product size info and the max/min
	 * oligo length in order to reduce the ranges of the ok regions. On
	 * some imputs this improves speed dramatically.
	 */
	public void optimizeOkRegionsList(P3GlobalSettings pa, SeqArgs sa) {
		/* We do this only if we enabled the optimization and
		 * the primers were NOT specified. */
		if (!LibPrimer3.OPTIMIZE_OK_REGIONS || (sa.getLeftInput() != null) || (sa.getRightInput()!= null)) {
			return;
		}

		/* If any pair is allowed, no point in doing this */
		if (this.any_pair) {
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
		for (int i=0; i<this.getCount(); i++) {
			int ls = -1, le = -1, rs = -1, re = -1;
			int new_ls = -1, new_le = -1, new_rs = -1, new_re = -1;
			if (this.left_pairs[i][0] != -1) {
				ls = this.left_pairs[i][0];
				le = this.left_pairs[i][0] + this.left_pairs[i][1] - 1;
			}
			if (this.right_pairs[i][0] != -1) {
				rs = this.right_pairs[i][0];
				re = this.right_pairs[i][0]	+ this.right_pairs[i][1] - 1;
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
				if (re > (sa.getSequence().length)) { 
					re = (sa.getSequence().length); 
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
				if (le > (sa.getSequence().length)) { 
					le = (sa.getSequence().length); 
				}
			}
			/* Temporary testing fSystem.out.format: */
			/* fSystem.out.format(stderr, "Adjusted range [%d,%d,%d,%d] to [%d,%d,%d,%d],
				    pmin is %d, pmax is %d, omin is %d, omax is %d\n",
				    this.left_pairs[i][0],
				    this.left_pairs[i][0] +
				    this.left_pairs[i][1] - 1,
				    this.right_pairs[i][0],
				    this.right_pairs[i][0] +
				    this.right_pairs[i][1] - 1, ls, le, rs, re,
				    pmin, pmax, omin, omax);
			 */
			this.left_pairs[i][0]  = ls;
			this.left_pairs[i][1]  = le - ls + 1;
			this.right_pairs[i][0] = rs;
			this.right_pairs[i][1] = re - rs + 1;
		}
		/* any_left and any_right not true anymore */
		this.any_left = false;
		this.any_right = false;
	}
	public int[][] getLeftPairs() {
		return left_pairs;
	}
	public int[][] getRightPairs() {
		return right_pairs;
	}
	/**
	 * @return the count
	 */
	public int getCount() {
		return count;
	}
	
}