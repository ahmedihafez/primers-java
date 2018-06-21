package org.primer3.libprimer3;

public class IntervalArrayT4 {
	  
	
	  int[][] left_pairs = new int[LibPrimer3.PR_MAX_INTERVAL_ARRAY][2];
	  int[][] right_pairs= new int[LibPrimer3.PR_MAX_INTERVAL_ARRAY][2];
	  boolean any_left = false ;  /* set to 1 if the empty pair ",," is given for any
			    left interval */
	  boolean any_right = false; /* set to 1 if the empty pair ",," is given for any
			    right interval */
	  boolean any_pair = false;  /* set to 1 if both intervals are given as empty */
	  int count = 0;     /* total number of pairs */
	  
	  
	 public int interval_array_t2_count()
	 {
		 return count;
	 }
	 public int[] interval_array_t2_get_pair(int i)
	 {
		 return null;
	 }
	 
	 public boolean p3_add_to_2_interval_array(int i1, int i2, int i3, int i4){
		 int c = this.count;
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
		 this.count++;
		 return  false;
	 }
	public boolean checkIncludedInAny(PrimerPair ppair) {

		boolean included = false;
		int l_start = ppair.left.start, l_end = ppair.left.start + ppair.left.length - 1;
		int r_start = ppair.right.start - ppair.right.length + 1, r_end = ppair.right.start;
		for (int i=0; i<this.count; i++) {
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

	 

}