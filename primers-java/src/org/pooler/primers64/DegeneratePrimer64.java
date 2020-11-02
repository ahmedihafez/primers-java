/*
  This file is part of Java porting of Primer Pooler (https://github.com/ssb22/PrimerPooler)
  Primer Pooler (c) Silas S. Brown.  For Wen.
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
    http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
package org.pooler.primers64;

import java.io.PrintStream;

import org.pooler.Amplicons.PrimerToFind;
import org.pooler.HelperMethods;
import org.pooler.primers.CountResult;
import org.pooler.primers.DegeneratePrimer;
import org.pooler.primers.IPrimer;

class DegeneratePrimer64  extends DegeneratePrimer {
	static String degenerateCombos="TGKCYSBAWRDMHVN";


	long 	MaybeA = 0,
			MaybeC = 0,
			MaybeG = 0,
			MaybeT = 0;

	public DegeneratePrimer64(String primerSeq) {

		for(int i =0;i<primerSeq.length() ; i++)
		{
			char l = Character.toUpperCase(primerSeq.charAt(i));
			int c = degenerateCombos.indexOf(l) + 1;
			this.MaybeA = (this.MaybeA <<1) | ((c&8)==0? 0L:01L);
			this.MaybeC = (this.MaybeC <<1) | ((c&4)==0? 0L:01L); 
			this.MaybeG = (this.MaybeG <<1) | ((c&2)==0? 0L:01L); 
			this.MaybeT = (this.MaybeT <<1) | ((c&1)==0? 0L:01L); 
		}
	}



	protected DegeneratePrimer64(long MaybeA ,
			long MaybeC ,
			long MaybeG,
			long MaybeT ) {
		this.MaybeA = MaybeA;
		this.MaybeC = MaybeC;
		this.MaybeG = MaybeG;
		this.MaybeT = MaybeT;
	}







	public DegeneratePrimer64(long AorT, long GorT, long valid ) {
		this.MaybeA = valid & AorT & ~GorT;
		this.MaybeC = valid & ~(AorT | GorT);
		this.MaybeG = valid & ~AorT & GorT;
		this.MaybeT = valid & AorT & GorT;
	}



	@Override
	public String toString() {
		return 
				"MaybeA: "  + Long.toBinaryString(MaybeA) + "\n" + 
				"MaybeC: "  + Long.toBinaryString(MaybeC) + "\n" +
				"MaybeG: "  + Long.toBinaryString(MaybeG) + "\n" +
				"MaybeT: "  + Long.toBinaryString(MaybeT) + "\n" ;
	}

	@Override
	public IPrimer getReverse() {
		/* assumes right-aligned */
		long    i1= this.MaybeA,
				i2=this.MaybeC,
				i3=this.MaybeG, 
				i4=this.MaybeT, 
				o1=0, o2=0, o3=0, o4=0;
		while(i1 != 0 || i2 != 0 || i3 != 0 || i4 != 0) {
			o1=(o1<<1)|(i1&1); o2=(o2<<1)|(i2&1);
			o3=(o3<<1)|(i3&1); o4=(o4<<1)|(i4&1);
			i1>>>=1; i2>>>=1; i3>>>=1; i4>>>=1;
		}		
		return new DegeneratePrimer64(o1,o2,o3,o4);
	}



	@Override
	public void addTag(IPrimer tag) {
		DegeneratePrimer64 degTag = (DegeneratePrimer64) tag;	
		int sL = Long.bitCount(this.MaybeA | this.MaybeC | this.MaybeG | this.MaybeT);
		this.MaybeA |= (degTag.MaybeA << sL);
		this.MaybeC |= (degTag.MaybeC << sL);
		this.MaybeG |= (degTag.MaybeG << sL);
		this.MaybeT |= (degTag.MaybeT << sL);
	}



	@Override
	public void addTagBackward(IPrimer tag) {
		DegeneratePrimer64 degTag = (DegeneratePrimer64) tag.getReverse();	
		int sL = Long.bitCount(degTag.MaybeA | degTag.MaybeC | degTag.MaybeG | degTag.MaybeT);
		this.MaybeA = ((this.MaybeA) << sL) | degTag.MaybeA;
		this.MaybeC = ((this.MaybeC) << sL) | degTag.MaybeC;
		this.MaybeG = ((this.MaybeG) << sL) | degTag.MaybeG;
		this.MaybeT = ((this.MaybeT) << sL) | degTag.MaybeT;		
	}



	@Override
	public void removeTag(IPrimer tag) {
		DegeneratePrimer64 degTag = (DegeneratePrimer64) tag;	
		long tValid = degTag.MaybeA | degTag.MaybeC | degTag.MaybeG | degTag.MaybeT,
				pValid = this.MaybeA | this.MaybeC | this.MaybeG | this.MaybeT;
		int toRM = Long.bitCount(tValid);
		long mask = ~(tValid << (Long.bitCount(pValid) - toRM));
		this.MaybeA &= mask; 
		this.MaybeC &= mask; 
		this.MaybeG &= mask;
		this.MaybeT &= mask;		
	}



	@Override
	public void removeTagBackward(IPrimer tag) {
		DegeneratePrimer64 degTag = (DegeneratePrimer64) tag;	
		int sR = Long.bitCount(degTag.MaybeA | degTag.MaybeC | degTag.MaybeG | degTag.MaybeT);
		this.MaybeA >>>= sR; 
		this.MaybeC >>>= sR; 
		this.MaybeG >>>= sR;
		this.MaybeT >>>= sR;		
	}



	@Override
	public int calcScore(IPrimer primer) {
		DegeneratePrimer64 degPrimer64 = (DegeneratePrimer64) primer;
		long p1Valid = this.MaybeA | this.MaybeC | this.MaybeG | this.MaybeT;
		long p2Valid = degPrimer64.MaybeA | degPrimer64.MaybeC | degPrimer64.MaybeG | degPrimer64.MaybeT;
		int sL=(64 - 1/*threshold*/) - Long.numberOfLeadingZeros(p2Valid);
		int maxScore = 0;

		long maybeA =  this.MaybeA << sL;
		long maybeC =  this.MaybeC << sL;
		long maybeG =  this.MaybeG << sL;
		long maybeT =  this.MaybeT << sL;
	 	long p1Bvalid = p1Valid << sL;
		int reload = sL - Long.numberOfLeadingZeros(p1Valid);
		if(reload< 0 )  
			reload=0;
//		int iterIndex = 0;
		while(true) {
			long overlap = p1Bvalid & p2Valid;
			if( overlap ==  0) 
				return maxScore;
			long bonds = (maybeA & degPrimer64.MaybeT) |
					(maybeC & degPrimer64.MaybeG) |
					(maybeG & degPrimer64.MaybeC) |
					(maybeT & degPrimer64.MaybeA);
			int score = 2 * Long.bitCount(bonds) - Long.bitCount(overlap);
			maxScore = (score > maxScore ? score : maxScore);
			if(reload != 0) {
				--sL;
				maybeA   = this.MaybeA << sL;
				maybeC   = this.MaybeC << sL;
				maybeG   = this.MaybeG << sL;
				maybeT   = this.MaybeT << sL;
				p1Bvalid = p1Valid << sL; 
				--reload;
			} else {
				maybeA   >>>=1; 
				maybeC   >>>=1; 
				maybeG   >>>=1;
				maybeT   >>>=1; 
				p1Bvalid >>>=1;
			}
//			iterIndex++;
		}
	}



	@Override
	public float calcDeltaG(IPrimer primer, float[] table) {
		DegeneratePrimer64 degPrimer64 = (DegeneratePrimer64) primer;
		long p1Valid = this.MaybeA | this.MaybeC | this.MaybeG | this.MaybeT;
		long p2Valid = degPrimer64.MaybeA | degPrimer64.MaybeC | degPrimer64.MaybeG | degPrimer64.MaybeT;
		int sL=(64 - 1/*threshold*/) - Long.numberOfLeadingZeros(p2Valid);
		float minDG = Float.POSITIVE_INFINITY;
		long maybeA = this.MaybeA << sL;
		long maybeC = this.MaybeC << sL;
		long maybeG = this.MaybeG << sL;
		long maybeT = this.MaybeT << sL;
		long p1Bvalid = p1Valid << sL;
		int reload = sL - Long.numberOfLeadingZeros(p1Valid);
		if(reload<0) reload=0;
		while(true) {
			long overlap = p1Bvalid & p2Valid;
		    if(overlap == 0) 
		    	return minDG;
		    long bonds = HelperMethods.rm_unstable_bonds64( (maybeA & degPrimer64.MaybeT) |
		                                      (maybeC & degPrimer64.MaybeG) |
		                                      (maybeG & degPrimer64.MaybeC) |
		                                      (maybeT & degPrimer64.MaybeA),overlap);
		    
		    // TODO :: this (bonds == 0 ) ? 0 : bonds is not need for now,
		    // the difference is the operation in C is undefinde and give 1 for 0 ?? 
		    // here it return 0 for 0
		    int shift = 64-2-Long.numberOfLeadingZeros((bonds == 0 ) ? 0 : bonds); 
		    long mask= ((long) 3) << shift;
		    long maskEnd = ((long) 3) << Long.numberOfTrailingZeros(bonds);
		    float dG = table[256 + ((((maybeC|maybeG) & mask)>>>(shift+1) ) == 0 ? 1 : 0 )]; // init (worst-case scenario is C or G)
		    for(; Long.compareUnsigned(mask,maskEnd) >= 0   ; mask>>>=1,shift--)
		      dG += HelperMethods.minDGdegenerate( 	
		    		  					(maybeA & mask) >>> shift,
				 						(maybeC & mask) >>> shift,
										(maybeG & mask) >>> shift,
										(maybeT & mask) >>> shift,
										(degPrimer64.MaybeA & mask) >>> shift,
										(degPrimer64.MaybeC & mask) >>> shift,
										(degPrimer64.MaybeG & mask) >>> shift, 
										(degPrimer64.MaybeT & mask) >>> shift, table);
		    dG += table[256 +  ((((maybeC|maybeG) & mask) >>> (shift+1)) == 0 ? 1 : 0)  ];
		    minDG = (dG < minDG ? dG : minDG);
		    if(reload != 0) {
		      --sL;
		      maybeA = this.MaybeA << sL;
		      maybeC = this.MaybeC << sL;
		      maybeG = this.MaybeG << sL;
		      maybeT = this.MaybeT << sL;
		      p1Bvalid = p1Valid << sL; --reload;
		    } else {
		    	maybeA   >>>=1; 
		      	maybeC   >>>=1; 
		      	maybeG   >>>=1;
		      	maybeT   >>>=1; 
		    	p1Bvalid >>>=1;
		    }
		  }
	}



	@Override
	public CountResult calcCount(IPrimer backwardPrimer) {
		DegeneratePrimer64 degPrimer64 = (DegeneratePrimer64) backwardPrimer;
		long p1Valid = this.MaybeA | this.MaybeC | this.MaybeG | this.MaybeT;
		long p2Valid = degPrimer64.MaybeA | degPrimer64.MaybeC | degPrimer64.MaybeG | degPrimer64.MaybeT;
		int sL=(64 - 1/*threshold*/) - Long.numberOfLeadingZeros(p2Valid);
		CountResult  count = new CountResult();
	
		long maybeA = this.MaybeA << sL;
		long maybeC = this.MaybeC << sL;
		long maybeG = this.MaybeG << sL;
		long maybeT = this.MaybeT << sL;
		long p1Bvalid = p1Valid << sL;
		int reload = sL - Long.numberOfLeadingZeros(p1Valid);
		if(reload<0) 
			reload=0;
		while(true) {
			long overlap = p1Bvalid & p2Valid;
			if(overlap == 0 ) return count;
			count.tried++;
			if (   ((maybeA & degPrimer64.MaybeT) |
					(maybeC & degPrimer64.MaybeG) |
					(maybeG & degPrimer64.MaybeC) |
					(maybeT & degPrimer64.MaybeA) ) != 0) count.count++;
			if(reload != 0) {
				--sL;
				maybeA = this.MaybeA << sL;
				maybeC = this.MaybeC << sL;
				maybeG = this.MaybeG << sL;
				maybeT = this.MaybeT << sL;
				p1Bvalid = p1Valid << sL; 
				--reload;
			} else {
				maybeA >>>=1; 
				maybeC >>>=1; 
				maybeG >>>=1;
				maybeT >>>=1; 
				p1Bvalid >>>=1;
			}
		}
	}



	@Override
	public void dGprint(IPrimer backwardPrimer, float minDG, PrintStream out, float[] table) {
		DegeneratePrimer64 degPrimer64 = (DegeneratePrimer64) backwardPrimer;
		long p1Valid = this.MaybeA | this.MaybeC | this.MaybeG | this.MaybeT;
		long p2Valid = degPrimer64.MaybeA | degPrimer64.MaybeC | degPrimer64.MaybeG | degPrimer64.MaybeT;
  
		int sL=(64 - 1) - Long.numberOfLeadingZeros(p2Valid);
		  int sR = 0; 
		  long p1B_MaybeA , p1B_MaybeC , p1B_MaybeG , p1B_MaybeT;
		  while(true) {
		    if(sL != 0) {
		      p1B_MaybeA = this.MaybeA << sL;
		      p1B_MaybeC = this.MaybeC << sL;
		      p1B_MaybeG = this.MaybeG << sL;
		      p1B_MaybeT = this.MaybeT << sL;
		    } else {
		      p1B_MaybeA = this.MaybeA >>> sR;
		      p1B_MaybeC = this.MaybeC >>> sR;
		      p1B_MaybeG = this.MaybeG >>> sR;
		      p1B_MaybeT = this.MaybeT >>> sR;
		      assert(sR < 64); /* if this breaks, check the range of nearlyEqual */
		    }
		    long overlap = (p1B_MaybeA|p1B_MaybeC|p1B_MaybeG|p1B_MaybeT) & p2Valid;
		    long bonds0 = (p1B_MaybeA & degPrimer64.MaybeT) |
		      (p1B_MaybeC & degPrimer64.MaybeG) |
		      (p1B_MaybeG & degPrimer64.MaybeC) |
		      (p1B_MaybeT & degPrimer64.MaybeA),
		      bonds = HelperMethods.rm_unstable_bonds64(bonds0,overlap);
		    int shift = 64-2-Long.numberOfLeadingZeros(bonds); 
		    long mask=(long)3 << shift,
		      maskEnd = (long)3 << Long.numberOfTrailingZeros(bonds);
		    float dG = table[256+ ( (((p1B_MaybeC|p1B_MaybeG) & mask)>>>(shift+1)) == 0  ? 1 : 0  )  ]; // init (worst-case scenario is C or G)
		    for(; Long.compareUnsigned(mask,maskEnd) >= 0; mask>>>=1,shift--)
		      dG += HelperMethods.minDGdegenerate(
		    		  	(p1B_MaybeA & mask)>>>shift,
						(p1B_MaybeC & mask)>>>shift,
						(p1B_MaybeG & mask)>>>shift,
						(p1B_MaybeT & mask)>>>shift,
						(degPrimer64.MaybeA & mask)>>>shift,
						(degPrimer64.MaybeC & mask)>>>shift,
						(degPrimer64.MaybeG & mask)>>>shift,
						(degPrimer64.MaybeT & mask)>>>shift,table);
		    dG += table[256+  ( (((p1B_MaybeC|p1B_MaybeG) & mask)>>>(shift+1)) == 0 ? 1 : 0  )];
		    if(HelperMethods.nearlyEqual(dG,minDG)) {
		      out.format("dG = %.3g\n",dG);
		      print64D_inner(sL-sR,degPrimer64,overlap,bonds0,out);
		      return;
		    }
		    if(sL != 0) sL--; else sR++;		
		  }
	}
	
	@Override
	public void print(IPrimer backwardPrimer, int maxScore, PrintStream out) {
		DegeneratePrimer64 degPrimer64 = (DegeneratePrimer64) backwardPrimer;
//		long p1Valid = this.MaybeA | this.MaybeC | this.MaybeG | this.MaybeT;
		long p2Valid = degPrimer64.MaybeA | degPrimer64.MaybeC | degPrimer64.MaybeG | degPrimer64.MaybeT;
		
		int sL=(64 - 1) - Long.numberOfLeadingZeros(p2Valid);
		int sR = 0; 
		long p1B_MaybeA , p1B_MaybeC , p1B_MaybeG , p1B_MaybeT;

		while(true) {
			if(sL != 0) {
				p1B_MaybeA = this.MaybeA << sL;
				p1B_MaybeC = this.MaybeC << sL;
				p1B_MaybeG = this.MaybeG << sL;
				p1B_MaybeT = this.MaybeT << sL;
		    } else {
		      p1B_MaybeA = this.MaybeA >>> sR;
		      p1B_MaybeC = this.MaybeC >>> sR;
		      p1B_MaybeG = this.MaybeG >>> sR;
		      p1B_MaybeT = this.MaybeT >>> sR;
		    }
		    long overlap = (p1B_MaybeA|p1B_MaybeC|p1B_MaybeG|p1B_MaybeT) & p2Valid;
		    // bit64 overlap = DegenerateValid64(p1B) & DegenerateValid64(p2);
		    long bonds = (p1B_MaybeA & degPrimer64.MaybeT) |
		      (p1B_MaybeC & degPrimer64.MaybeG) |
		      (p1B_MaybeG & degPrimer64.MaybeC) |
		      (p1B_MaybeT & degPrimer64.MaybeA);
		    
		    int score = 2*Long.bitCount(bonds) - Long.bitCount(overlap);
		    if(score == maxScore) {
		      /* TODO: if more than one ==maxScore, how to
		         prioritise the links in the degenerate case? */
		      out.format("Matches = %d\n", Long.bitCount(bonds));
		      out.format("Score = %d\n",maxScore);
		      print64D_inner(sL-sR,degPrimer64,overlap,bonds,out);
//		      print64D_inner(sL-sR,p2,overlap,bonds,f);
		      //return; /* comment out to print ALL maxScore matches */
		    }
		    if(!   (overlap != 0)) return; /* needed if not returning above */
		    if(sL != 0) sL--; else sR++;
		  }
		
	}
	
	
	 void print64D_inner(int sL,DegeneratePrimer64 p2,long overlap,long bonds,PrintStream out) {
		 int i1 = Long.numberOfLeadingZeros((this.MaybeA | this.MaybeC | this.MaybeG | this.MaybeT))-sL, 
				 i2 = Long.numberOfLeadingZeros(p2.MaybeA | p2.MaybeC | p2.MaybeG | p2.MaybeT), 
				 iMid = Long.numberOfLeadingZeros(overlap); 
		if(i1<i2){ i2-=i1; iMid-=i1; i1=0; } else { i1-=i2; iMid-=i2; i2=0; }
		HelperMethods.indent(i1,out); 
		out.print("5'-"); 
		this.printBases64D(out); 
		out.print("-3'\n");
//		HelperMethods.indent(iMid+(sizeof("5'-")-1), out);
		HelperMethods.indent(iMid+(3), out);

		long bond = (long)1 << (64-1-Long.numberOfLeadingZeros(overlap));
		for(; (bond&overlap) != 0; bond>>>=1) 
			out.print((bond&bonds) != 0 ?'|':'x');
		out.print('\n'); 
		HelperMethods.indent(i2,out);
		out.print("3'-"); 
		p2.printBases64D(out); 
		out.print("-5'\n");
	 }



	private void printBases64D(PrintStream out) {
		long valid = this.MaybeA | this.MaybeC | this.MaybeG | this.MaybeT;
		long i = (long)1 << (64-1-Long.numberOfLeadingZeros(valid));
		  for(; (i&valid)  != 0; i>>>=1) {
		    int j =  (((this.MaybeA & i)!=0 ? 1 : 0   )  << 3 ) |
		    		 (((this.MaybeC & i)!=0 ? 1 : 0   )  << 2 ) |
		    		 (((this.MaybeG & i)!=0 ? 1 : 0   )  << 1 ) |
		    		  ((this.MaybeT & i)!=0  ? 1 : 0  );
		    out.print(degenerateCombos.charAt(j-1));
		  }		
	}



	@Override
	public int NumPossibilities_32bases() {
		  /* number of possible primers this degenerate primer might be equal to (assumes aligned right) */
		  int r = 1, sR, poss;
		  for(sR=0; sR<32; sR++) { /* 32, not sixty four etc, because we want the return value to work with Make2bitFrom64D  */
		    poss = (int) (((this.MaybeA >>> sR) & 1) + ((this.MaybeC >>> sR) & 1)
		      + ((this.MaybeG >>> sR) & 1) + ((this.MaybeT >>> sR) & 1));
		    if (poss != 0) 
		    	r *= poss; 
		    else 
		    	break;
		  }
		  return r;		
		  
	}



	@Override
	public IPrimer getComplement() {
		DegeneratePrimer64 complment = new DegeneratePrimer64(MaybeA, MaybeC, MaybeG, MaybeT);
		  long tmp = complment.MaybeA; 
		  complment.MaybeA = complment.MaybeT; 
		  complment.MaybeT=tmp;
		  tmp = complment.MaybeG; 
		  complment.MaybeG = complment.MaybeC; 
		  complment.MaybeC = tmp;
		  return complment;
	}



	@Override
	public boolean make2bit(PrimerToFind newPrimerPoss, int possNo, int nPoss) {
		 /* note: ULL across all 3 .h variants
	     - more than 32 bases here will be truncated
	     (return value is 1 if it got truncated).
	     TODO: we could make a "not degenerate" version of
	     this which just shifts bits around, but I'd be very
	     surprised if this function is anywhere near the top
	     of a profile trace, so for now I'll leave it as you
	     have to call upgradeToDegenerate64 from the Maybe.

	     Note: for ease of binary search (see amplicons.c),
	     bases are shifted in from the LEFT (not from the
	     right as in the other functions), and the result
	     reads backwards from the 'end' cursor at left.
	  */
	  int sR, poss;
	  long toOut=0,toOutValid=0;
	  for(sR=0; sR<32; sR++) { /* 32, not sixty four etc */
	    int MaybeA = (int) ((this.MaybeA >>> sR) & 1),
	      MaybeC = (int) ((this.MaybeC >>> sR) & 1),
	      MaybeG = (int) ((this.MaybeG >>> sR) & 1),
	      MaybeT = (int) ((this.MaybeT >>> sR) & 1);
	    poss = MaybeA + MaybeC + MaybeG + MaybeT;
	    if (poss != 0) {
	      int partitionSize = nPoss / poss;
	      int possToTake = possNo / partitionSize;
	      possNo %= partitionSize; nPoss /= poss;
	      long bits = 0; /* we set it to a value to stop the "might be uninitialised" warning */
	      if (MaybeT != 0 ) possToTake--; /* if(!possToTake--) bits=0; but it's at 0 anyway */
	      if (MaybeC != 0  && possToTake-- == 0) bits=1;
	      if (MaybeA != 0  && possToTake-- == 0) bits=2;
	      if (MaybeG != 0  && possToTake   == 0) bits=3;
	      int sL = 62-2*sR; /* IMPORTANT: don't write (64-2) as it'll be changed to 32-2 in 32.h; this is ULL */
	      toOut |= (bits << sL);
	      toOutValid |= ((long)3 << sL);
	    } else break;
	  } 
	  newPrimerPoss.p = toOut; 
	  newPrimerPoss.valid = toOutValid;
//	#define Is_64bit /* will change to Is_32bit in 32.h */
//	#ifdef Is_32bit
//	  return 0; // avoid compiler warnings
//	#else
	  return ((this.MaybeA >>> 32) | (this.MaybeC >>> 32) | (this.MaybeG >>> 32) | (this.MaybeT >>> 32)) != 0;
//	#endif
	}



	@Override
	public void printBases(PrintStream outstream) {
		printBases64D(outstream);
		
	}



	
}