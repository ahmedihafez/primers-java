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
package org.primer3.dpal;

import java.util.HashMap;

/**
 *  Structure for passing in arguments to the main function, dpal. 
 */
public class DPAlignmentArgs {
	
	/** 
     * If non-0, check for and raise an error on an
     * illegal character in the input strings.
     */ 
	int check_chars;        
	/** 
	  * If non-0, print debugging information to
	  * stderr.
	  */
	int debug;              
	
	 /** Exit with -1 on error. */
	int fail_stop;          
	
	
    /** 
     * One of DPAL_GLOBAL, DPAL_LOCAL,
     * DPAL_GLOBAL_END, DPAL_LOCAL_END
     */
	public int flag;
	
	/**
	 *  Force the use of the generic function. 
	 */
	int force_generic;     
	/** 
	  * Force the use of the long generic no-path
	  * function.
	  */
	int force_long_generic; 
	
	/**
	 *  Force the use of the long maxgap 1 functions. 
	 */
	int force_long_maxgap1; 
	/** The "gap opening" penalty. */
	int gap;               
	
	/**
	 *  The "gap extension" penalty. 
	 */
	int gapl;               
	
	/** 
     * The maximum allowable size for a gap. -1
     * indicates that the gap can be of any size.
     */
	int max_gap;
	/** If greater than 0 stop search as soon as
	  * score > score_max.
	  */
	int score_max;
	 /** 
	   * If non-0, only print the score on
	   * stdout. (Incompatible with debug.)
	   */
	int score_only;
	/**
	 *   
	 */
//	dpal_ssm ssm = new dpal_ssm();      
	/** 
	 * The scoring system matrix.
	 * The data structure that stores the "scoring system matrix". (The socring
	 * system matrix data structure is of size UCHAR_MAX + 1 by UCHAR_MAX + 1.)
	 */
	static int UCHAR_MAX = 255;
	int[][] ssm = new int[UCHAR_MAX + 1][UCHAR_MAX + 1];
	
	
	/**
	 * Initialize the argument to the default matrix for nucleotide matches.
	 */
	public void dpal_set_default_nt_args(){
		set_dpal_args( );
		this.check_chars        = 1;
		this.debug              = 0;
		this.fail_stop          = DPAlignment.DPAL_EXIT_ON_ERROR;
		this.gap                = -100;
		this.gapl               = -100;
		this.max_gap            = 3;
		this.score_only         = 0;
	}
	
	/** 
	 * Routine primarily for testing: sets CC & GG matches to 3, AA & TT 
	 * matches to 2. 
	 */
	void dpal_set_h_nt_matrix(){
		
	}
	
	/** 
	 * The argument, a, must be a DNA scoring matrix.  Modifies a so that it for a
	 * match between any two ambiguity codes (or between ambiguity code and base),
	 * e.g. B and S, the score will be the maximum of score between any base in B
	 * and any base in S, in the example between any pair in {C, G, T} X {C, G}.
	 * This function overwrites any scores already associated with pairs of
	 * ambiguity codes.  Return 0 on error, 1 on success.
	 */
	public int dpal_set_ambiguity_code_matrix(){
		
		char[] amb_codes = "BDHVRYKMSWN".toCharArray();
		char[] all_bases = "ACGT".toCharArray();
		
		
		for(char c1 : amb_codes)
		{	
			char[] bases1 = xlate_ambiguity_code.get(c1);
			if(bases1 == null) return 0;
		    /* Do matches between c1 and all other
		       ambiguity codes. */
			for(char c2 : amb_codes)
			{
				char[] bases2 = xlate_ambiguity_code.get(c2);
				if(bases2 == null) return 0;
				int extreme = Short.MIN_VALUE;
				for(char b1 : bases1)
				{
					for(char b2 : bases2)
					{
						if(ssm[b1][b2] > extreme)
						{
							 extreme = ssm[b1][b2];
						}
					}
				}
				/* extreme is now the maximum score
		         for a match between any 2 bases
		         represented *c1, *c2. */
		      ssm[c1][c2] = extreme;
			}

			
			 /* Do matches between c1 and all bases. */
		    for (char b2 : all_bases) {
		      int extreme =  Short.MIN_VALUE;
		      for (char  b1 : bases1) {
		        if (ssm[b1][b2] > extreme) {
		          extreme = ssm[b1][b2];
		        }
		      }
		      ssm[c1][b2] = extreme;
		      ssm[b2][c1] = extreme;
		    }
			
		}
		
		return 1;
		
	}
	
	public void set_dpal_args()
	{
		init_ssm();
		this.gap                = -200;
		this.gapl               = -200;
		this.flag               = DPAlignment.DPAL_LOCAL;
		this.max_gap            = 1;
		this.fail_stop          = 1;
		this.check_chars        = 0;
		this.debug              = 0;
		this.score_only         = 1;
		this.force_generic      = 0;
		this.force_long_generic = 0;
		this.force_long_maxgap1 = 0;
	}
	
	
	void init_ssm()
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
	
	static HashMap<Character,char[]> xlate_ambiguity_code = new HashMap<Character, char[]>();
	
	
	static {
		
		xlate_ambiguity_code.put('N', "ACGT".toCharArray());
		xlate_ambiguity_code.put('B', "CGT".toCharArray());
		xlate_ambiguity_code.put('D', "AGT".toCharArray());
		xlate_ambiguity_code.put('H', "ACT".toCharArray());
		xlate_ambiguity_code.put('V', "ACG".toCharArray());
		xlate_ambiguity_code.put('R', "AG".toCharArray());
		xlate_ambiguity_code.put('Y', "CT".toCharArray());
		xlate_ambiguity_code.put('K', "GT".toCharArray());
		xlate_ambiguity_code.put('M', "AC".toCharArray());
		xlate_ambiguity_code.put('S', "CG".toCharArray());
		xlate_ambiguity_code.put('W', "AT".toCharArray());

	}
}
