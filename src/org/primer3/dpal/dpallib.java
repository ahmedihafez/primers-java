package org.primer3.dpal;

import java.util.HashMap;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.primer3.sequence.Sequence;

public class dpallib {
	
	public static int DPAL_ERROR_SCORE = Integer.MIN_VALUE;
	
	/* 0 means do not exit on error. */
	public static int DPAL_EXIT_ON_ERROR = 0;

	
	/** 
     * The maximum size of a string that can be
     * aligned with with generic dpal and for which
     * we can return a "path".  Several arrays of
     * size DPAL_MAX_ALIGN X DPAL_MAX_ALIGN are
     * statically allocated in dpal.o
     */
	public static int DPAL_MAX_ALIGN = 1600;
	
	 /**
	  *  Return a local alignment. 
	  */
	public static int DPAL_LOCAL = 0;
	/**
	 * Return a global alignment anchored at the end of the first sequence
	 */
	public static int DPAL_GLOBAL_END = 1;
	
	/** 
     * Return an arbitrary global alignment, that is
     * one anchored at the end of either the first or
     * the second sequence.
     */
	public static int DPAL_GLOBAL = 2;
	
	/** 
     * Return a local alignment that includes the
     * end (but not necessarily the beginning) of
     * the first sequence.
     */
	public static int DPAL_LOCAL_END = 3;

	
	
	
	
	
	
	
	public static dpal_results dpal(Sequence s1 , Sequence s2, dpal_args in	) throws AlignmentException
	{
		dpal_results out = new dpal_results();
		int xlen, ylen;
//		  char msg = "Illegal character in input: ?";

		out.score = DPAL_ERROR_SCORE;
		out.path_length = 0;
		out.msg = null;

		  
		char[] X = s1 .getSequence();
		char[] Y = s2 .getSequence();
		
		  
		if(null == X)
			throw new AlignmentException("NULL first sequence");
		if(null == Y)
			throw new AlignmentException("NULL second sequence");
		if(null == in)
			throw new AlignmentException("NULL 'in' pointer");

//		if (NULL == out) return; /* Leave it to the caller to crash */

		if(in.flag != DPAL_GLOBAL
		              && in.flag != DPAL_GLOBAL_END
		              && in.flag != DPAL_LOCAL_END
		              && in.flag != DPAL_LOCAL)
			throw new AlignmentException("Illegal flag");
//		  if (in.check_chars) {
//		    CHECK_ERROR(illegal_char(X, in.ssm, &msg[28]), msg);
//		    CHECK_ERROR(illegal_char(Y, in.ssm, &msg[28]), msg);
//		  }

		xlen = X.length;
		ylen = Y.length;

		out.align_end_1 = -1;
		out.align_end_2 = -1;

		if (xlen == 0) {
			out.msg = "Empty first sequence";
			out.score = 0;
			return out ;
		}
		if (ylen == 0) {
			out.msg = "Empty second sequence";
			out.score = 0;
			return out ;
		}
		
		if (1 == in.force_generic || in.debug == 1 || 0 == in.score_only) {
		    /* 
		     * A true value of in.debug really means "print alignment on stderr"
		     * and implies 0 == a.score_only.
		     */
		    if(xlen > DPAL_MAX_ALIGN)
		    	throw new AlignmentException("Sequence 1 longer than DPAL_MAX_ALIGN and alignment is requested");
		    if(ylen > DPAL_MAX_ALIGN)
		    	throw new AlignmentException("Sequence 2 longer than DPAL_MAX_ALIGN and alignment is requested");
		    _dpal_generic(X, Y, in, out);
		  } else if (1 == in.force_long_generic) {
		    _dpal_long_nopath_generic(X, Y, in, out);
		  } else if (1 == in.max_gap ) {
		    if (DPAL_LOCAL == in.flag)
		      _dpal_long_nopath_maxgap1_local(X, Y, in, out);
		    else if (DPAL_GLOBAL_END == in.flag)
		      _dpal_long_nopath_maxgap1_global_end(X, Y, in, out);
		    else if (DPAL_LOCAL_END == in.flag)
		      _dpal_long_nopath_maxgap1_local_end(X, Y, in, out);
		    else if (xlen <= DPAL_MAX_ALIGN && ylen <= DPAL_MAX_ALIGN)
		      _dpal_generic(X, Y, in, out);
		    else _dpal_long_nopath_generic(X, Y , in, out);
		  }
		  else if (xlen < DPAL_MAX_ALIGN && ylen < DPAL_MAX_ALIGN)
		    _dpal_generic(X, Y, in, out);
		  else
		    _dpal_long_nopath_generic(X, Y, in, out);

		
		
		
		
		
		return out;
	}
	
	
	public static dpal_results dpal(char[] s1, char[] s2, dpal_args a) throws AlignmentException {
		return dpal(new Sequence(s1),new Sequence(s2),a);		
	}
	
	
	public static void _dpal_generic(char[] X , char[] Y, dpal_args in , dpal_results out	) throws AlignmentException {
	    
		
		int xlen =  X.length,ylen =Y.length;
		/* The "score matrix" (matrix of best scores). */
	    int[][] S = new int[DPAL_MAX_ALIGN][DPAL_MAX_ALIGN];
//#ifndef DPAL_FORGET_PATH
	    /* The matrix of "trace" pointers */
	    int[][][] P = new int[DPAL_MAX_ALIGN][DPAL_MAX_ALIGN][3];
//#endif
	    
	    int i, j, k, mg, c;
	    int gap = in.gap, gapl = in.gapl, max_gap = in.max_gap;

//#ifndef DPAL_FORGET_PATH
	    int i0 = -99, j0 = -99;
	    int saved_k;
//#endif 

	    int I = -99, J = -99; /* Coordinates of the maximum score. */
	    int smax;             /* The optimum score. */
	    int score = -99;      /* Current score. */

	    int a,b,max;

//#ifdef DPAL_PRINT_COVERAGE
//	    fprintf(stderr, "_dpal_generic called\n");
//#endif

	    if(xlen > DPAL_MAX_ALIGN)
	            throw new AlignmentException("First sequence too long for _dpal_generic");
	    if(ylen > DPAL_MAX_ALIGN)
	    		throw new AlignmentException("Second sequence too long for _dpal_generic");

	    /* Initialize the 0th column of the score matrix. */
	    smax = Integer.MIN_VALUE;
	    for(i=0; i < xlen; i++) {
	    	score = in.ssm[X[i]][Y[0]]; 
	    	if (DPAL_LOCAL == in.flag) {
	    		if (score < 0) 
	    			score = 0;
	    		if(score > smax) {
	    			smax = score;
	    			I=i; J=0;
	    		}
	    	}
	        else if (DPAL_LOCAL_END == in.flag) { 
	        	if (score < 0) score = 0;
	        }
	    	S[i][0] = score;
	    }   
	    /* Move code for find global-alignment and end-anchored
	       alignment below? */
	    if (DPAL_LOCAL != in.flag) {
	    	/* 
	         * For a non-local alignment we restrict our search for the maximum
	         * score to the last row.
	         */
	    	smax = S[xlen-1][0]; I=xlen-1; J=0;
	    }
	           
	    /* Initialize the 0th row of the score matrix. */
	    for(j=0; j<ylen; j++) { 
	        score = in.ssm[X[0]][Y[j]]; 

	        if(DPAL_LOCAL == in.flag){
	            if (score < 0) score = 0;
	            if(score > smax){
	                smax = score;
	                I=0; J=j;
	            }
	        }
	        else if (DPAL_LOCAL_END == in.flag) {if (score < 0) score = 0;}
	        S[0][j] = score;
	    }   
	    if(DPAL_GLOBAL == in.flag&&S[0][ylen-1]>smax){
	                smax = S[0][ylen-1];
	                I=0; J=ylen-1;
	    }

	    /* Further is the solution for dynamic programming problem. */
	    for(i=1; i<xlen; i++) {
	        for(j=1; j<ylen; j++) {

	            a=S[i-1][j-1];

	            b = c = Integer.MIN_VALUE;
	            if (1 == max_gap) {
	                if (i > 1) {
	                    b = S[i-2][j-1] + gap;
//#ifndef DPAL_FORGET_PATH
	                    i0 = i - 2;
//#endif
	                }
	                if (j > 1) {
	                    c = S[i-1][j-2] + gap;
//#ifndef DPAL_FORGET_PATH
	                    j0 = j - 2;
//#endif
	                }
	            } else if (max_gap > 1) {
	                max = Integer.MIN_VALUE;
	                mg=(max_gap+1>i||max_gap<0)?i:max_gap+1;
	                for(k=2; k<=mg; k++) {
	                    c = S[i-k][j-1] + gap + gapl*(k-2);
	                    if(c>max){
	                        max=c;
//#ifndef DPAL_FORGET_PATH
	                        i0 = i-k;
//#endif
	                    }
	                }
	                b=max;

	                max=Integer.MIN_VALUE;
	                mg=(max_gap+1>j||max_gap<0)?j:max_gap+1;
	                for(k=2;k<=mg;k++) {
	                    c = S[i-1][j-k] + gap + gapl*(k-2);
	                    if(c>max){
	                        max=c;
//#ifndef DPAL_FORGET_PATH
	                        j0 = j-k;
//#endif
	                    }
	                }
	                c=max;
	            }

	            if(a>=b && a>=c) {
	                score = a + in.ssm[X[i]][Y[j]];
//#ifndef DPAL_FORGET_PATH
	                P[i][j][1] = i-1;
	                P[i][j][2] = j-1;
//#endif
	            } else if (b > a && b >= c) {
	                score = b + in.ssm[X[i]][Y[j]];
//#ifndef DPAL_FORGET_PATH
	                P[i][j][1] = i0;
	                P[i][j][2] = j-1;
//#endif
	            } else if (c > a && c > b) {
	                score = c + in.ssm[X[i]][Y[j]];
//#ifndef DPAL_FORGET_PATH
	                P[i][j][1] = i-1;
	                P[i][j][2] = j0;
//#endif
	            }
	            if (score >= smax)
	                /* 
	                 * Because of comparison '>=' immediately above, dpal reports
	                 * ungapped (i.e. diagonal) alignments if there is a choice
	                 * of more than one optimum alignment.
	                 */
	                /* Move code to get 'g' and 'e' maxima to a separate loop ? */
	                if (DPAL_LOCAL == in.flag 
	                    || (DPAL_GLOBAL_END == in.flag && i == xlen-1)
	                    || (DPAL_LOCAL_END   == in.flag && i == xlen-1)
	                    || (DPAL_GLOBAL == in.flag&& (i==xlen-1||j==ylen-1))) {
	                    /*  
	                     * If in.flag is DPAL_LOCAL, then a cell anywhere within
	                     * S may be the endpoint of the alignment.  If in.flag is
	                     * DPAL_GLOBAL_END, then only cells in the last row may be
	                     * the endpoint.  If in.flag is DPAL_GLOBAL cells in the
	                     * last row _or_ the last column may be the endpoint.
	                     */
	                    smax = score;
	                    I = i;
	                    J = j;
	            } /*  put else here ? */
	            if (score < 0 && (DPAL_LOCAL == in.flag 
	                              || DPAL_LOCAL_END == in.flag))
	                /* 
	                 * For a local alignment, 0 is the lowest score that we record
	                 * in S.
	                 */
	                score = 0;

	            S[i][j]=score;
	        }
	    }
	    /* I and J now specify the last pair of an optimum alignment. */

//#ifndef DPAL_FORGET_PATH    
	        k = (I > J) ? I+1 : J+1;
	        saved_k=k;

	        out.path[k][0]=I; out.path[k][1]=J;
	        while(out.path[k][0]!=0&&out.path[k][1]!=0) {
	            if ((in.flag== DPAL_LOCAL || in.flag == DPAL_LOCAL_END)
	                     &&S[out.path[k][0]][out.path[k][1]]==0) {
	              k++; break;
	            }
	            out.path[k-1][0] = P[out.path[k][0]][out.path[k][1]][1];
	            out.path[k-1][1] = P[out.path[k][0]][out.path[k][1]][2];
	            k--;
	        }
	        if (k>0) {
	            for (i=0;i<=saved_k-k;i++) {
	                out.path[i][0] = out.path[i+k][0];
	                out.path[i][1] = out.path[i+k][1];
	            }
	        }
//#endif

	        if ((DPAL_LOCAL == in.flag 
	             || DPAL_LOCAL_END == in.flag)&& S[I][J] <= 0) {
	            /* There is no alignment at all. */
	            out.score = 0;
	            out.path_length = 0;
	        } else {
	            out.score = smax;
	            out.align_end_1 = I;
	            out.align_end_2 = J;
//#ifndef DPAL_FORGET_PATH        
	            out.path_length = saved_k - k + 1;
//#else
	            // TODO :: see what is best to use
//	            out.path_length = 0;
//#endif
	        }
//#ifndef DPAL_FORGET_PATH
	        if (in.debug == 1) 
	        	print_align(X,Y,P,I,J, in);
//#endif
	        return;
	    
	}
	
	
	/* Linear space, no path, for any value of maxgap and for any alignment. */
	public static void _dpal_long_nopath_generic(char[] X , char[] Y, dpal_args in , dpal_results out	) throws AlignmentException {
		int xlen =  X.length,ylen =Y.length;
	    int i, j, k, mg, mgy, c;
	    int gap = in.gap, gapl = in.gapl, max_gap = in.max_gap;


	    int I = -99, J = -99; /* Coordinates of the maximum score. */
	    int smax;             /* The optimum score. */
	    int score;            /* Current score. */
	    
		/* The "score matrix" (matrix of best scores). */
	    int[][] S  = new int[max_gap+2][xlen];
//	    int[][] P  = new int[max_gap+2][xlen];
	    int[]   SI;
	    
	    out.score = DPAL_ERROR_SCORE;
	    out.path_length = 0;
	    out.msg = null;
	    
	    
	    /* Initialize the 0th column of the score matrix. */
	    smax = Integer.MIN_VALUE;
	    for(i=0; i < xlen; i++) {
	        score = in.ssm[X[i]][Y[0]]; 
	        if (DPAL_LOCAL == in.flag) {
	            if (score < 0) score = 0;
	            if(score > smax) {
	                smax = score;
	                I=i; J=0;
	            }
	        }
	        else if (DPAL_LOCAL_END == in.flag) {if (score < 0) score = 0;}
	        S[0][i] = score;
	        /* Move code for find global-alignment and end-anchored
	        alignment below? */
	     if (DPAL_LOCAL != in.flag) {
	         /* 
	          * For a non-local alignment we restrict our search for the maximum
	          * score to the last row.
	          */
	         smax = S[0][xlen-1]; I=xlen-1; J=0;
	     }
	            
	     /* Initialize the 0th row of the score matrix. 
	     for(j=0; j<ylen; j++) { 
	         score = in.ssm[X[0]][Y[j]]; 

	         if(DPAL_LOCAL == in.flag){
	             if (score < 0) score = 0;
	             if(score > smax){
	                 smax = score;
	                 I=0; J=j;
	             }
	         }
	         S[0][j] = score;
	     }   
	     
	     if(DPAL_GLOBAL == in.flag&&S[0][ylen-1]>smax){
	                 smax = S[0][ylen-1];
	                 I=0; J=ylen-1;
	     }
	     */

	     /* Further is the solution for dynamic programming problem. */
	     for(j=1; j<ylen; j++) {
	         mgy=(max_gap+1>j||max_gap<0)?j:max_gap+1;
	         score = in.ssm[X[0]][Y[j]];
	          if (DPAL_LOCAL == in.flag) {
	              if (score < 0) score = 0;
	              if(score > smax) smax = score;
	          }    
	          else if (DPAL_LOCAL_END == in.flag) { if (score < 0) score = 0;}
	          else if (DPAL_GLOBAL == in.flag && j == ylen-1 && score > smax)
	                         smax = score;
	         S[mgy][0] = score;
	         for(i=1; i<xlen; i++) {

	             score=S[mgy-1][i-1];

	                 mg=(max_gap+1>i||max_gap<0)?i:max_gap+1;
	                 for(k=2; k<=mg; k++) 
	                     if((c = S[mgy-1][i-k] + gap + gapl*(k-2)) > score)score = c;

	                 for(k=2;k<=mgy;k++) 
	                     if((c = S[mgy-k][i-1] + gap + gapl*(k-2)) > score)score=c;

	                 score += in.ssm[X[i]][Y[j]];

	             if (score >= smax)
	                 /* 
	                  * Because of comparison '>=' immediately above, dpal reports
	                  * ungapped (i.e. diagonal) alignments if there is a choice
	                  * of more than one optimum alignment.
	                  */
	                 /* Move code to get 'g' and 'e' maxima to a separate loop ? */
	                 if (DPAL_LOCAL == in.flag 
	                     || ((DPAL_GLOBAL_END == in.flag
	                          || DPAL_LOCAL_END == in.flag) 
	                         && i == xlen-1)
	                     || (DPAL_GLOBAL == in.flag&& (i==xlen-1||j==ylen-1))) {
	                     /*  
	                      * If in.flag is DPAL_LOCAL, then a cell anywhere within
	                      * S may be the endpoint of the alignment.  If in.flag is
	                      * DPAL_GLOBAL_END, then only cells in the last row may be
	                      * the endpoint.  If in.flag is DPAL_GLOBAL cells in the
	                      * last row _or_ the last column may be the endpoint.
	                      */
	                     smax = score;
	                     I = i;
	                     J = j;
	             } /*  put else here ? */
	             if (score < 0 && (DPAL_LOCAL == in.flag
	                               || DPAL_LOCAL_END == in.flag))
	                 /* 
	                  * For a local alignment, 0 is the lowest score that we record
	                  * in S.
	                  */
	                 score = 0;

	             S[mgy][i]=score;
	         }
	         if(mgy == max_gap + 1){
	             SI = S[0];
	             for(i=0; i<mgy; i++)S[i] = S[i+1];
	             S[mgy] = SI;
	         }
	     }
	     /* I and J now specify the last pair of an optimum alignment. */

	     if (DPAL_LOCAL == in.flag && smax <= 0) {
	         /* There is no alignment at all. */
	         out.score = 0;
	         out.path_length = 0;
	     } else {
	         out.score = smax;
	         out.align_end_1 = I;
	         out.align_end_2 = J;
	     }
	    }   
	}
	
	static void _dpal_long_nopath_maxgap1_local(char[] X , char[] Y, dpal_args in , dpal_results out	) throws AlignmentException {
		int xlen =  X.length,ylen =Y.length;
		 
		int i, j;
		int gap = in.gap;
		int smax;           /* The optimum score. */
		int score;          /* Current score. */
		int a;
		
		/* The "score matrix" (matrix of best scores). */
	    int[] S0, S1, S2; 
	    int[] P0, P1, P2;
	    int[] S;
	    
	    P0 = new int[ylen];
	    P1 = new int[ylen];
	    P2 = new int[ylen];
	    
	    S0 = P0; S1 = P1; S2 = P2;

	    smax = 0; /* For local alignment score can never be less than 0. */

	    /* Initialize the 0th row of the score matrix. */
	    for(j=0; j < ylen; j++) { 
	        score = in.ssm[X[0]][Y[j]]; 
	        if (score < 0) score = 0;
	        else if (score > smax) smax = score;
	        /*S[0][j] = score;*/
	        S0[j] = score;
	    }   

	    /* Set the 1st row of the score matrix. */
	    score = in.ssm[X[1]][Y[0]];
	    if(score < 0) score = 0;
	    else if (score > smax) smax = score;
	    S1[0] = score;
	    for(j=1; j < ylen; j++) {
	        score = S0[j-1];
	        if(j>1 && (a=S0[j-2] + gap) > score)score = a;
	        score += in.ssm[X[1]][Y[j]];
	        if (score < 0) score = 0;
	        else if(score > smax) smax = score;
	        S1[j] = score;
	    }

	    for(i=2; i < xlen; i++) {
	        score = in.ssm[X[i]][Y[0]];
	        if (score < 0) score = 0;
	        else if (score > smax) smax = score;
	        S2[0] = score;
	        score = S1[0];
	        if((a=S0[0] + gap) > score) score = a;
	        score += in.ssm[X[i]][Y[1]];
	        if(score < 0) score = 0;
	        else if (score > smax) smax = score;
	        S2[1] = score;
	        for(j=2; j < ylen; j++) {
	            score = S0[j-1];
	            if((a=S1[j-2])>score) score = a;
	            score +=gap;
	            if((a=S1[j-1]) >score) score = a;

	            score += in.ssm[X[i]][Y[j]];       
	            if (score < 0 ) score = 0;
	            else if (score > smax) smax = score;
	            S2[j]=score;
	        }
	        S = S0; S0 = S1; S1 = S2; S2 = S;
	    }
	    out.score = smax;
	    out.path_length=0;
	}
	
	
	static void _dpal_long_nopath_maxgap1_global_end(char[] X , char[] Y, dpal_args in , dpal_results out	) throws AlignmentException {
		int xlen =  X.length,ylen =Y.length;
		int i, j,k;
		int gap = in.gap;
		int smax;           /* The optimum score. */
		int score;          /* Current score. */
		int a,t;
		
		/* The "score matrix" (matrix of best scores). */
	    int[] S0, S1, S2; 
	    int[] P0, P1, P2;
	    int[] S;
	    
	    P0 = new int[xlen];
	    P1 = new int[xlen];
	    P2 = new int[xlen];
	    S0 = P0; S1 = P1; S2 = P2;

	    smax = in.ssm[X[xlen-1]][Y[0]];
	             
	    /* Set the 0th row of the score matrix. */
	    for(j=0; j<xlen; j++) S0[j] = in.ssm[X[j]][Y[0]]; 

	    /* Set the 1st row of the score matrix. */
	    S1[0] = in.ssm[X[0]][Y[1]];
	    for(j=1; j < xlen; j++){
	      score = S0[j-1];
	      if(j>1 && (a=S0[j-2] + gap)> score)score = a;
	      score += in.ssm[X[j]][Y[1]];
	      if(score > smax && j == xlen-1) smax = score;
	      S1[j] = score;
	    }

	    k = ylen - (int)(xlen / 2) + 1;
	    if (k<1) k = 1;

	    /* Set the rectangular part of almost the remainder of the matrix. */
	    for(j=2; j<k+1; j++) {
	      S2[0] = in.ssm[X[0]][Y[j]];
	      score = S1[0];
	      if((a=S0[0]+gap) > score) score = a;
	      score += in.ssm[X[1]][Y[j]];
	      S2[1] = score;
	      for(i=2; i<xlen-1; i++) {
	        score = S1[i-2];
	        if((a=S0[i-1]) > score)score = a;
	        score += gap;
	        if((a=S1[i-1]) > score)score = a;
	        score += in.ssm[X[i]][Y[j]];
	        S2[i] = score;
	      }
	      score = S1[xlen-3];
	      if((a=S0[xlen-2]) > score)score = a;
	      score += gap;
	      if((a=S1[xlen-2]) > score)score = a;
	      score += in.ssm[X[xlen-1]][Y[j]];
	      S2[xlen-1] = score;
	      if(score > smax) smax = score;
	      S = S0; S0 = S1; S1 = S2; S2 = S;
	    }

	    /* Set the triangular part of almost the remainder of the matrix. */
	    t = 2;
	    for(j=k+1; j<ylen; j++) {
	      for(i=t; i<xlen-1; i++) {
	        score = S1[i-2];
	        if((a=S0[i-1]) > score) score = a;
	        score += gap;
	        if((a=S1[i-1]) > score) score = a;
	        score += in.ssm[X[i]][Y[j]];
	        S2[i] = score;
	      }
	      t += 2;
	      score = S1[xlen-3];
	      if((a=S0[xlen-2]) > score)score = a;
	      score += gap;
	      if((a=S1[xlen-2]) > score)score = a;
	      score += in.ssm[X[xlen-1]][Y[j]];
	      S2[xlen-1] = score;
	      if(score > smax) smax = score;
	      S = S0; S0 = S1; S1 = S2; S2 = S;
	    }

	    out.score = smax;
	    out.path_length=0;
	       
	    
	}
	
	static void _dpal_long_nopath_maxgap1_local_end(char[] X , char[] Y, dpal_args in , dpal_results out) throws AlignmentException
	{
		int xlen =  X.length,ylen =Y.length;
		/* The "score matrix" (matrix of best scores). */
		  

		int i, j;
		int gap = in.gap;
		int smax;           /* The optimum score. */
		int score;          /* Current score. */
		int a;

		int[] S0, S1, S2; 
		int[] P0, P1, P2;
		int[] S;
		if(ylen < 3)
			throw new AlignmentException("_dpal_long_nopath_maxgap1_local_end requires ylen >= 3\n");  
	    
		P0 = new int[ylen];
		P1 = new int[ylen];
		P2 = new int[ylen];
		 S0 = P0; S1 = P1; S2 = P2;

		  smax = 0; /* For local alignment score can never be less than 0. */

		  /* Initialize the 0th row of the score matrix. */
		  for(j=0; j < ylen; j++) { 
		    score = in.ssm[X[0]][Y[j]]; 
		    if (score < 0) score = 0;
		    /*S[0][j] = score;*/
		    S0[j] = score;
		  }   

		  /* Set the 1st row of the score matrix. */
		  score = in.ssm[X[1]][Y[0]];
		  if(score < 0) score = 0;
		  S1[0] = score;
		  for(j=1; j < ylen; j++) {
		    score = S0[j-1];
		    if(j>1 && (a=S0[j-2] + gap) > score)score = a;
		    score += in.ssm[X[1]][Y[j]];
		    if (score < 0) score = 0;
		    S1[j] = score;
		  }

		  for(i=2; i < xlen - 1; i++) {
		    score = in.ssm[X[i]][Y[0]];
		    if (score < 0) score = 0;
		    S2[0] = score;
		    score = S1[0];
		    if((a=S0[0] + gap) > score) score = a;
		    score += in.ssm[X[i]][Y[1]];
		    if(score < 0) score = 0;
		    S2[1] = score;
		    for(j=2; j < ylen; j++) {
		      score = S0[j-1];
		      if((a=S1[j-2])>score) score = a;
		      score +=gap;
		      if((a=S1[j-1]) >score) score = a;

		      score += in.ssm[X[i]][Y[j]];       
		      if (score < 0 ) score = 0;
		      S2[j]=score;
		    }
		    S = S0; S0 = S1; S1 = S2; S2 = S;
		  }
		  /* Calculate scores for last row (i = xlen-1) and find smax */
		  i = xlen - 1;
		  score = in.ssm[X[i]][Y[0]];
		  if (score < 0) score = 0;
		  else if (score > smax) smax = score;
		  S2[0] = score;
		  score = S1[0];
		  if((a=S0[0] + gap) > score) score = a;
		  score += in.ssm[X[i]][Y[1]];
		  if(score < 0) score = 0;
		  else if (score > smax) smax = score;
		  S2[1] = score;
		  for(j=2; j < ylen; j++) {
		    score = S0[j-1];
		    if((a=S1[j-2])>score) score = a;
		    score +=gap;
		    if((a=S1[j-1]) >score) score = a;
		    score += in.ssm[X[i]][Y[j]];
		    if (score < 0 ) score = 0;
		    else if (score > smax) smax = score;
		    S2[j]=score;
		  }
		  out.score = smax;
		  out.path_length=0;
		
		
	}
	
	
	private static void print_align(char[] X, char[] Y, int[][][] P, int I,
			int J, dpal_args dargs) {
		 
		
		int xlen =  X.length,ylen =Y.length;
		int[] JX = new int[DPAL_MAX_ALIGN],JY = new int[DPAL_MAX_ALIGN];
		int k,i,j,n,m;
		char[] sx = new char[3*DPAL_MAX_ALIGN],sy = new char [3*DPAL_MAX_ALIGN],sxy = new char[3*DPAL_MAX_ALIGN];

		for (i=0; i < 3*DPAL_MAX_ALIGN; i++) {
			sx[i] = ' '; sy[i] = ' '; sxy[i] = ' ';
		}
		if(I>J)
			k=I+1;
		else 
			k=J+1;

		n=k;
		JX[k] = I;
		JY[k] = J;
		while(JX[k]!=0&&JY[k]!=0){
			JX[k-1] = P[JX[k]][JY[k]][1];
			JY[k-1] = P[JX[k]][JY[k]][2];
			k--;
		}
		if(JX[k]>JY[k]){
			for(i=0;i<JX[k];i++)sx[i] = X[i];
			for(i=0;i<JX[k]-JY[k];i++)sy[i] = ' ';
			j = JX[k]-JY[k];
			for(i=JX[k]-JY[k];i<JX[k];i++)sy[i] = Y[i-j];
			m = JX[k];
		}
		else{
			for(i=0;i<JY[k];i++)sy[i] = Y[i];
			for(i=0;i<JY[k]-JX[k];i++)sx[i] = ' ';
		    j= JY[k]-JX[k];
		    for(i=j;i<JY[k];i++)sx[i] = X[i-j];
		    m = JY[k];
		}
		for(i=0;i<m;i++)sxy[i] = ' ';
		for(i=k;i<n;i++){
			sx[m] = X[JX[i]];
			sy[m] = Y[JY[i]];
			/* if(sx[m]==sy[m]&&sx[m]!='N') sxy[m] = '|'; */
			if (dargs.ssm[sx[m]][sy[m]] > 0)
				sxy[m] = '|';
			else sxy[m]=' ';
			if(JX[i+1]-JX[i]>JY[i+1]-JY[i]){
				for(j=1;j<JX[i+1]-JX[i];j++){
					sy[m+j] = '-';
					sx[m+j] = X[JX[i]+j];
					sxy[m+j] = ' ';
				}
				m += JX[i+1]-JX[i]-1;
			}
			if(JY[i+1]-JY[i]>JX[i+1]-JX[i]){
				for(j=1;j<JY[i+1]-JY[i];j++){
					sx[m+j] = '-';
					sy[m+j] = Y[JY[i]+j];
					sxy[m+j] = ' ';
				}
				m += JY[i+1]-JY[i]-1;
			}
			m++;
		}
		sx[m] = X[I];
		sy[m] = Y[J];
		for (i=m+1; i < (m + xlen - I); i++) 
			sx[i]=X[i-m+I];
		for (i=m+1; i < (m + ylen - J); i++) 
			sy[i]=Y[i-m+J];

		if (dargs.ssm[sx[m]][sy[m]] > 0)
			sxy[m] = '|';
		else sxy[m]=' ';
		m++;
		if (xlen - I > ylen -J) {
			k = m + xlen - I;
		} else {
			k = m + ylen - J;
		}
		        
		j=0;
		while(j<k){
			for(i=j;i<j+70;i++) System.err.print(sx[i]);
			System.err.print("\n");
			for(i=j;i<j+70;i++) System.err.print(sxy[i]);
			System.err.print("\n");
		    for(i=j;i<j+70;i++) System.err.print(sy[i]); 
		    System.err.print("\n");
		    for(i=0;i<70;i++)   System.err.print("_");
		    System.err.print( "\n");
		    j +=70;
		}
		
	}











	public static void main(String[] args)
	{
		// AAGATCAGCTAGCTAGCTAGCTTTTT TTTTTCCGATCAGCTAGCTAGCTAGCTTTTTTT l
		String[] testCases = new String[]{
				"-s1 GTGAAGCCTCAGGTAGTGCA -s2 CTCTGTCGACTTTGCCACCA -mode g  "
		}; 
		
		
		args = testCases[0].split(" ");
		
		String s1 = null, s2 = null;
		
		dpal_args a = new dpal_args();
	    dpal_results r;
	    int tmp_ret;
	    int i;
	    int print_align_end = 0; /* 
	                              * Print align_end_1 and align_end_2 from
	                              * dpal_results.
				      */
	    int use_ambiguity_codes = 0;
	    int use_h_matrix = 0;
	    char mode;
		
		
		
		CommandLine commandLine = setArgs(args);
		
		
		if(commandLine != null)
		{
			if(commandLine.hasOption("s1")) //s1
			{
				s1 = commandLine.getOptionValue("s1");	
			}
			if(commandLine.hasOption("s2")) //s1
			{
				s2 = commandLine.getOptionValue("s2");	
			}
			
			
			
			
			
			a.dpal_set_default_nt_args();
			
			
			mode = 'G';
			if ('l' == mode)
				a.flag = DPAL_LOCAL;
			else if ('e' == mode || 'G' == mode)
				a.flag = DPAL_GLOBAL_END;
			else if ('g' == mode)
				a.flag = DPAL_GLOBAL;
			else if ('L' == mode)
				a.flag = DPAL_LOCAL_END;
		    if(print_align_end == 1) a.force_long_generic = 1;

			try {
				r = dpallib.dpal(new Sequence(s1.toCharArray()),new Sequence(s2.toCharArray()), a);
				if (r.score == DPAL_ERROR_SCORE) {
				      System.err.format("Error: %s\n", r.msg);
				      System.exit(-1);
				    }
				if (a.score_only == 1) {
					System.out.format("%.2f\n", 0.01 * r.score);
					if (print_align_end == 1) {
						if(r.align_end_1 >= 0) 
							System.out.format("align_end_1=%d ",r.align_end_1);
						if(r.align_end_2 >= 0) 
							System.out.format("align_end_2=%d\n ",r.align_end_2);
					}
				} else {
					System.out.format("|%s|  |%s| %c ", s1, s2, mode); 
					System.out.format("score=%.2f len=%d ", (0.01 * r.score), r.path_length);
					if (print_align_end == 1) {
						if(r.align_end_1 >= 0) System.out.format("align_end_1=%d ",r.align_end_1);
						if(r.align_end_2 >= 0) System.out.format("align_end_2=%d ",r.align_end_2);
					}
					for (i=0; i<r.path_length; i++)
						System.out.format("|%d,%d", r.path[i][0],r.path[i][1]);
					System.out.format("|\n");
				}
				
				
			} catch (AlignmentException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}


	static CommandLine setArgs(String[] args)
	{
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		Option  argOption   = Option.builder("g").argName( "gval" )
                .hasArg()
                .desc("gval> are (positive) float (.01 precision) specifying penalties for creating a gap." + 
                " (the penalties are subtracted from the output score")
                .build();
		options.addOption(argOption);
		
		argOption   = Option.builder("l").argName( "lval" )
                .hasArg()
                .desc("lval are (positive) float (.01 precision) specifying penalties for lengthening a gap." + 
                " (the penalties are subtracted from the output score")
                .build();
		options.addOption(argOption);
		
		argOption   = Option.builder("a")
                .desc("causes the scoring matrix to be modified by dpal_set_ambiguity_codes." )
                .build();
		options.addOption(argOption);
		
		argOption   = Option.builder("e")
                .hasArg()
                .desc("causes the end postion of the alignment in both sequences to be printed")
                .build();
		options.addOption(argOption);
		
		argOption   = Option.builder("mode").argName( "mode" )
                .hasArg()
                .required()
                .desc("is one of g, G, l, or L")
                .build();
		options.addOption(argOption);
		
	
		
		argOption = Option.builder("s1").argName( "seq1" )
                .hasArg()
                
                .required()
                .desc("seq1")
                .build();
		options.addOption(argOption);
		argOption   = Option.builder("s2").argName( "seq2" )
                .hasArg()
                .required()
                .desc("seq2" )
                .build();
		
		options.addOption(argOption);
		
		Option helpOption = new Option( "help", "Print Usage" );

		options.addOption(helpOption);

		// automatically generate the help statement
		HelpFormatter formatter = new HelpFormatter();
		formatter.setOptionComparator(null);
		try {
			// parse the command line arguments
		    CommandLine line = parser.parse( options, args );
		    
		    if(line.getOptions().length == 0 )
		    {
				formatter.printHelp( "ntdpal", options );
		    }
		    
		    return line;
		}
		catch( ParseException exp ) {
			// oops, something went wrong
		    System.err.println( "Parsing failed.  Reason: " + exp.getMessage() );
			formatter.printHelp( "ntdpal", options );
		}
		
		return null;
	}


	
}















