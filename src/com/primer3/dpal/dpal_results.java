package com.primer3.dpal;

public class dpal_results {
	
	
	//const   char *msg;
	public String msg ;
    int[][]     path = new int[dpallib.DPAL_MAX_ALIGN][2];
    int     path_length;
    int     align_end_1; /* Last alignment position in the 1st sequence. */
    int     align_end_2; /* Last alignment position in the 2nd sequence. */
    public double  score;
}
