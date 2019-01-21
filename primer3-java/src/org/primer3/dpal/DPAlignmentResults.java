package org.primer3.dpal;

public class DPAlignmentResults {
	
	
	//const   char *msg;
	public String msg ;
    int[][]     path = new int[DPAlignment.DPAL_MAX_ALIGN][2];
    int     path_length;
    public int     align_end_1; /* Last alignment position in the 1st sequence. */
    public int     align_end_2; /* Last alignment position in the 2nd sequence. */
    public double  score;
}
