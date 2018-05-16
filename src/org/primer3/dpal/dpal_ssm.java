package org.primer3.dpal;

/** 
 * The data structure that stores the "scoring system matrix". (The socring
 * system matrix data structure is of size UCHAR_MAX + 1 by UCHAR_MAX + 1.)
 */
public class dpal_ssm {
	static int UCHAR_MAX = 255;
	int[][] ssm = new int[UCHAR_MAX + 1][UCHAR_MAX + 1];
	
	//  alphabet
	String acgtn = "ACGTN";
	public dpal_ssm()
	{
		
		 for (int i = 0; i <= UCHAR_MAX; i++)
		 {
			    for (int j = 0; j <= UCHAR_MAX; j++)
			    {
			    	 if (	('A' == i || 'C' == i || 'G' == i || 'T' == i || 'N' == i) && 
			    			('A' == j || 'C' == j || 'G' == j || 'T' == j || 'N' == j)
			    		)
			    	 {
			    		 if(i == 'N' || j == 'N')
			    		 {
			    			 ssm[i][j] = -25;
			    		 } else if (i == j)
			    		 {
			    			 ssm[i][j] = 100;
			    		 }
			    		 else
			    		 {
			    			 ssm[i][j] = -100;
			    		 }
			    	 }
			    	 else
			    	 {
			    		 ssm[i][j] = Short.MIN_VALUE;
			    	 }
			    }
		 }
		
	}


}
