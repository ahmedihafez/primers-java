package org.pooler;

import java.io.PrintStream;

public final class HelperMethods {


	static final public long rm_unstable_bonds64(long bonds,long overlap) {
		/* MPprimer_dimer_check.pl's "primer_dimer" function
		     does this before its deltaG calculation.  Says
		     single matched bases surrounded by mismatches are
		     unstable, removes 01 when not followed by valid 1 */
		long lone_1s = bonds & ((~bonds)<<1) & ((~bonds)>>>1);
		lone_1s &= (overlap>>>1); /* left-hand bit is never got rid of */
		return bonds & ~lone_1s;
	}

	static final public float minDGdegenerate(
			long MaybeA12,long MaybeC12,long MaybeG12,long MaybeT12,long MaybeA34,long MaybeC34,long MaybeG34,long MaybeT34,
			final float[] table) {
		int[] inBits = new int[16]; // [a,c,g,t, a,c,g,t, a,c,g,t, a,c,g,t]
		inBits[0] = (int) (MaybeA12 >>> 1); 
		inBits[1] = (int) (MaybeC12 >>> 1);
		inBits[2] = (int) (MaybeG12 >>> 1); 
		inBits[3] = (int) (MaybeT12 >>> 1);
		inBits[4] = (int) (MaybeA12 & 1); 
		inBits[5] = (int) (MaybeC12 & 1);
		inBits[6] = (int) (MaybeG12 & 1); 
		inBits[7] = (int) (MaybeT12 & 1);
		inBits[8] = (int) (MaybeA34 >>> 1); 
		inBits[9] = (int) (MaybeC34 >>> 1);
		inBits[10] = (int) (MaybeG34 >>> 1); 
		inBits[11] = (int) (MaybeT34 >>> 1);
		inBits[12] = (int) (MaybeA34 & 1); 
		inBits[13] = (int) (MaybeC34 & 1);
		inBits[14] = (int) (MaybeG34 & 1); 
		inBits[15] = (int) (MaybeT34 & 1);
		return minDGdegenerate(0,5,inBits,0,table);
	}

	static final float minDGdegenerate(int valSoFar,int shift,final int[] inBits,int inBitsOffset, final float[] table) {
		/* AorT 0 GorT */ final int A=4,C=0,G=1,T=5 ; /* shift 5,4,1,0 */
		if(shift == 3 ) 
			shift=1; 
		else if (shift<0) 
			return table[valSoFar]; 
		float m=0;
		if(inBits[inBitsOffset + 0] != 0) 
			m= Float.min(m,minDGdegenerate(valSoFar|(A<<shift),shift-1,inBits, inBitsOffset + 4 , table));
		if(inBits[inBitsOffset + 1] != 0) 
			m= Float.min(m,minDGdegenerate(valSoFar|(C<<shift),shift-1,inBits, inBitsOffset + 4, table));
		if(inBits[inBitsOffset + 2] != 0) 
			m= Float.min(m,minDGdegenerate(valSoFar|(G<<shift),shift-1,inBits, inBitsOffset + 4,table));
		if(inBits[inBitsOffset + 3] != 0) 
			m= Float.min(m,minDGdegenerate(valSoFar|(T<<shift),shift-1,inBits, inBitsOffset + 4 ,table));
		return m;
	}
	static final int dGbucket(float dG,int max) {
		  /* Added in v1.32.  Previously had min2(-min2(0,2*dG),max)
		     with min2 declared as int, but it segfaulted dGsCounts64
		     in the case of dG calculation going to infinity
		     (e.g. user had set all parameters to 0, causing a
		     log of 0 to be attempted in deltaG_table) and
		     the cast being undefined: could get 0x80000000 which doesn't
		     play nicely with sign flipping etc */
		  if (dG >= 0) return 0;
		  dG *= -2.0;
		  if (dG >= max) return max;
		  /* OK, safe to cast */
		  int dGi = (int) dG;
		  return dGi > max ? max : dGi; /* repeat just in case */
		}

	static final void prnSeconds(long numSecs) {
		if(numSecs>=3600) 
			System.err.format(" (%d:%02d:%02d)",(int)(numSecs/3600L),(int)(numSecs%3600L)/60,(int)(numSecs%60L));
		else {
			int secs = (int)numSecs;
			if(secs>=60) 
				System.err.format(" (%dmin %dsec)",secs/60,secs%60);
			else if(secs>1) 
				System.err.format(" (%d seconds)",secs);
			else System.err.format(" (1 second)");
		}
	}
	
	public static final boolean nearlyEqual(float a,float b) {
		  /* After having calculated minDG, the output code
		     ought to be able to say if dG==minDG to check
		     it's seeing the minimum one again.  In theory,
		     rounding error shouldn't matter as long as we
		     take the EXACT same steps for both calculations
		     (rouding errors should be exactly identical).
		     In practice, when compiling with optimisation,
		     some of GCC's optimisations result in floating-
		     point calculations becoming slightly different
		     in different places (reordering etc), and it's
		     NOT guaranteed that you'll get the exact same
		     bits from two different parts of the program
		     that seemingly do the same thing!  So we have
		     to do the "nearly equal" thing.
		     
		     TODO: this won't work at threshold 0 if there
		     are very small <0 deltaG's (but those are not
		     likely to be accurate anyway: take them out?
		     insist on a threshold of at least -1 ?) could
		     adapt the 0.001 but beware of slowing it down
		     (if b is a min, try b += fabs(b)*0.0001 before
		     the loop)
		  */
		  return Math.abs(a-b) <= 0.001;
		}
	
	
	static public final void indent(int n,PrintStream out){
		while(n-- != 0)
			out.print(' ');
	}
	
}
