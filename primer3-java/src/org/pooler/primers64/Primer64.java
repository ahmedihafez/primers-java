package org.pooler.primers64;

import org.pooler.primers.NonDegeneratePrimer;

import java.io.PrintStream;

import org.pooler.HelperMethods;
import org.pooler.primers.CountResult;
import org.pooler.primers.IPrimer;

class Primer64  extends NonDegeneratePrimer {
	public Primer64(String primerSeq) {
		for(int i =0;i<primerSeq.length() ; i++)
		{
			char l = Character.toUpperCase(primerSeq.charAt(i));
			// l should be AGCT only
			this.AorT = (this.AorT<<1) | ( (l=='A'||l=='T')? 1:0 );
			this.GorT = (this.GorT<<1) | ( (l=='G'||l=='T') ? 1 :0 );
			this.valid = (this.valid<<1) | 1;
		}
	}
	protected Primer64(long AorT, long GorT, long valid) {
		this.AorT = AorT;
		this.GorT = GorT;
		this.valid= valid;
	}
	long AorT = 0;
	long GorT = 0;
	long valid = 0;

	@Override
	public String toString() {
		return 
				"AorT: "  + Long.toBinaryString(AorT) + "\n" + 
				"GorT: "  + Long.toBinaryString(GorT) + "\n" +
				"valid: "  + Long.toBinaryString(valid) + "\n" ;
	}

	@Override
	public IPrimer getReverse() {
		/* assumes 'valid' is right-aligned */
		long 
		v=this.valid,
		i1=this.AorT,
		i2=this.GorT,
		o1=0,
		o2=0;
		while(v != 0) {
			o1=(o1<<1)|(i1&1); o2=(o2<<1)|(i2&1);
			i1>>>=1; i2>>>=1; v>>>=1;
		}
		return new Primer64(o1,o2,this.valid);
	}
	@Override
	public void addTag(IPrimer tag) {
		Primer64 tag64 =  (Primer64) tag;
		/* assumes 'valid' is right-aligned in both p and tag */
		int sL = Long.bitCount(this.valid); /* = 64-leading0_64 */
		this.valid |= (tag64.valid << sL);
		this.AorT  |= (tag64.AorT  << sL); 
		this.GorT  |= (tag64.GorT  << sL);

	}
	@Override
	public void addTagBackward(IPrimer tag) {
		/* for backward primers: reverse the tag first, and add it to the lsb of the primer rather than the msb */
		Primer64 tag64 =  (Primer64) tag.getReverse();
		int sL = Long.bitCount(tag64.valid);
		this.valid = ((this.valid) << sL) | tag64.valid;
		this.AorT  = ((this.AorT)  << sL) | tag64.AorT;
		this.GorT  = ((this.GorT)  << sL) | tag64.GorT;
	}
	@Override
	public void removeTag(IPrimer tag) {
		Primer64 tag64 =  (Primer64) tag;
		int toRM   =  Long.bitCount(tag64.valid);
		long mask  =  ~(tag64.valid << (Long.bitCount(this.valid) - toRM));
		this.valid &= mask;
		this.AorT  &= mask; 
		this.GorT  &= mask;
	}
	@Override
	public void removeTagBackward(IPrimer tag) {
		Primer64 tag64 =  (Primer64) tag;
		int sR = Long.bitCount(tag64.valid);
		this.valid >>>= sR; 
		this.AorT  >>>= sR; 
		this.GorT  >>>= sR;
	}
	@Override
	public int calcScore(IPrimer primer) {
		Primer64 primer64 = (Primer64) primer;
		/* score the interaction between p1 and p2, fast.
	     p2 must have been PrimerReverse64'd by the caller. */
		int sL=(64 - 1/*threshold*/) - Long.numberOfLeadingZeros(primer64.valid);

		int maxScore = 0; /* this initial value of maxScore is also the minimum score that can be returned.  Do not make it negative without reviewing code that assumes it's >=0 */
		long AorT   = this.AorT  << sL; /* we start with p1 shifted left */
		long GorT   = this.GorT  << sL;
		long valid  = this.valid << sL;
		int reload = sL - Long.numberOfLeadingZeros(this.valid);
		if(reload < 0) 
			reload = 0;
		/* TODO: if rewritten into 2 loops, can reload only once: load in the middle, 
		 * shift to the right until gone, reload in the middle <<1, while overlap test + shift left.  
		 * Of course if the above reload <= 0 then do just the 1 loop as below because it'll be faster in that case. */
		while(true) {
			long overlap = valid & primer64.valid;
			if(overlap == 0) 
				return maxScore; /* all done */
			long bonds = (~(AorT ^ primer64.AorT)) & (GorT ^ primer64.GorT) & overlap;
			int score = 2 * Long.bitCount(bonds) - Long.bitCount(overlap);
			maxScore = (score > maxScore ? score : maxScore);
			if(reload != 0) {
				--sL;
				AorT =  this.AorT  << sL; 
				GorT =  this.GorT  << sL;
				valid = this.valid << sL; 
				--reload;
			} else {
				AorT  >>>= 1; 
				GorT  >>>= 1; 
				valid >>>= 1;
			}
		}
	}
	@Override
	public float calcDeltaG(IPrimer primer, float[] table) {
		Primer64 primer64 = (Primer64) primer;

		  /* like score64 but does deltaG instead */
		int sL=(64 - 1/*threshold*/) - Long.numberOfLeadingZeros(primer64.valid);
		  float minDG = Float.POSITIVE_INFINITY;
		  long AorT  = this.AorT  << sL; 
		  long GorT  = this.GorT  << sL;
		  long valid = this.valid << sL;
		  int reload = sL - Long.numberOfLeadingZeros(this.valid);
		  if(reload<0) reload=0;
		  while(true) {
		    long overlap =  valid & primer64.valid;
		    if( overlap == 0) 
		    	return minDG;
		    long bonds = HelperMethods.rm_unstable_bonds64((~(AorT ^ primer64.AorT)) & (GorT ^ primer64.GorT) & overlap, overlap);
		    int shift = 64-2- Long.numberOfLeadingZeros(bonds); 
		    long mask= ((long)3) << shift; 
		    long maskEnd = ((long)3) << Long.numberOfTrailingZeros(bonds);
		    float dG = table[256 + (((AorT & mask) >>> (shift+1)) != 0 ? 1 : 0 ) ]; // init
		    for(; Long.compareUnsigned(mask,maskEnd) >= 0; mask>>>=1,shift--)
		      dG += table[ (int) 
		                   	( 	(((AorT & mask) >>> shift) << 6) | 
		                		(((GorT & mask) >>> shift) << 4) | 
		                		(((primer64.AorT & mask) >>> shift) << 2) | 
		                		 ((primer64.GorT & mask) >>> shift)
		                	)];
		    dG += table[256 + (((AorT & mask) >>> (shift+1)) !=0 ? 1 : 0  )]; // init at end
		    minDG = (dG < minDG ? dG : minDG);
		    if(reload != 0) {
		      --sL;
		      AorT = this.AorT << sL; 
		      GorT = this.GorT << sL;
		      valid = this.valid << sL; 
		      --reload;
		    } else {
		      AorT >>>=1; 
		      GorT >>>=1; 
		      valid >>>=1;
		    }
		  }
	}
	@Override
	public CountResult calcCount(IPrimer primer) {
		/* count the number of alignments of >0 bonds,
	     for information only.  Similar to score64, but called
	     only when outputting interaction data.
	     (TODO: could make this return maxScore or minDG as
	     well, to save having to call that func separately,
	     but low priority because this is called only when
	     printing out bonds in excess of threshold) */
		Primer64 primer64 = (Primer64) primer;
		/* score the interaction between p1 and p2, fast.
	     p2 must have been PrimerReverse64'd by the caller. */
		int sL=(64 - 1/*threshold*/) - Long.numberOfLeadingZeros(primer64.valid);
		
		CountResult count =  new CountResult();
	  long AorT = this.AorT << sL; 
	  long GorT = this.GorT << sL;
	  long valid = this.valid << sL;
	  int reload = sL - Long.numberOfLeadingZeros(this.valid);
	  if(reload<0) 
		  reload=0;
	  while(true) {
	    long overlap = valid & primer64.valid;
	    if(overlap == 0) return count;
	    count.tried++;
	    if((( ~( AorT ^ primer64.AorT)) & ( GorT ^primer64.GorT) & overlap ) != 0 ) count.count++;
	    if(reload != 0) {
	      --sL;
	      AorT = this.AorT << sL; 
	      GorT = this.GorT << sL;
	      valid = this.valid << sL; 
	      --reload;
	    } else {
	      AorT  >>>=1; 
	      GorT  >>>=1; 
	      valid >>>=1;
	    }
	  }
	}
	@Override
	public void dGprint(IPrimer primer, float minDG, PrintStream out, float[] table) {
		Primer64 primer64 = (Primer64) primer;
		 int sL= (64 - 1) - Long.numberOfLeadingZeros(primer64.valid);
		  int sR = 0; 
//		  Primer64 p1B;
		  long AorT, GorT,valid;
		  while(true) {
		    if(sL != 0) {
		      AorT = this.AorT << sL; 
		      GorT = this.GorT << sL;
		      valid = this.valid << sL;
		    } else {
		      AorT = this.AorT >>> sR; 
		      GorT = this.GorT >>> sR;
		      valid = this.valid >>> sR;
		      assert(valid != 0); /* if this breaks, check the range of nearlyEqual */
		    }
		    long overlap = valid & primer64.valid;
		    long bonds0 = (~(AorT ^ primer64.AorT)) & (GorT ^primer64.GorT) & overlap;
		    long bonds = HelperMethods.rm_unstable_bonds64(bonds0, overlap);
		    // see if it affect the result bonds == 0 ? 0 : bonds
		    int shift = 64-2-Long.numberOfLeadingZeros(bonds == 0 ? 0 : bonds); 
		    long mask= (long)3 << shift,
		      maskEnd = (long)3 << Long.numberOfTrailingZeros(bonds);
		    float dG = table[256+  (((AorT & mask)>>>(shift+1)) !=0 ? 1 : 0 )   ]; // init
		    for(; Long.compareUnsigned(mask,maskEnd) >= 0; mask>>>=1,shift--)
		      dG += table[ (int) ((((AorT & mask) >>> shift)<<6) | 
		                   (((GorT & mask) >>> shift)<<4) | 
		                   (((primer64.AorT & mask) >>> shift)<<2) |
		                    ((primer64.GorT & mask) >>> shift))];
		    dG += table[256+  (((AorT & mask) >>> (shift+1)) != 0 ? 1 : 0   )]; // init at end
		    if(HelperMethods.nearlyEqual(dG,minDG)) {
		      out.printf("dG = %.3g\n",dG);
		      print64_inner(sL-sR,primer64,overlap,bonds0,out);
		      return;
		    }
		    if(sL != 0) 
		    	sL--; 
		    else sR++;
		  }		
	}
	
	@Override
	public void print(IPrimer primer, int maxScore, PrintStream out) {
		/* maxScore has been found by score64; print a representation
	     of the interaction, along with the score */
		Primer64 primer64 = (Primer64) primer;
		int sL= (64 - 1) - Long.numberOfLeadingZeros(primer64.valid);

//		int sL=(64 - 1) - leading0_64(p2.valid);
	  int sR = 0; 
	  long AorT, GorT,valid;
//	  Primer64 p1B;
	  	while(true) {
	    if(sL != 0) {
	      AorT = this.AorT << sL; 
	      GorT = this.GorT << sL;
	      valid = this.valid << sL;
	    } else {
	      /* this function is allowed to be a bit slower than
	         score64, and we need to keep all bits */
	      AorT = this.AorT >>> sR; 
	      GorT = this.GorT >>> sR;
	      valid = this.valid >>> sR;
	    }
	    long overlap = valid & primer64.valid;
	    long bonds = (~(AorT ^ primer64.AorT)) & (GorT ^ primer64.GorT) & overlap;
	    int score = 2*Long.bitCount(bonds) - Long.bitCount(overlap);
	    if(score == maxScore) {
	      /* TODO: if more than one ==maxScore, prioritise
	         any that has more C-G links, stronger than A-T */
	      out.format("Matches = %d\n",Long.bitCount(bonds));
	      out.format("Score = %d\n",maxScore);
	      print64_inner(sL-sR,primer64,overlap,bonds,out);
	      //return; /* comment out to print ALL maxScore matches */
	    }
	    if(! (overlap !=0)) return; /* needed if not returning above */
	    if(sL != 0 ) sL--; else sR++;
	  }
		
	}
	
	
	void print64_inner(int sL,Primer64 p2,long overlap,long bonds,PrintStream out) {
		  /* code common to print64 and dGprint64 */
		  int i1 = Long.numberOfLeadingZeros(this.valid)-sL, 
				  i2 = Long.numberOfLeadingZeros(p2.valid), 
				  iMid = Long.numberOfLeadingZeros(overlap); 
		  if(i1<i2) {
			  i2-=i1; 
			  iMid-=i1; 
			  i1=0; 
		  } else {
			i1-=i2;
			iMid-=i2; 
			i2=0; 
		  }
		  HelperMethods.indent(i1,out); 
		  out.print("5'-"); 
		  this.printBases64(out); 
		  out.print("-3'\n");
//		  indent(iMid+(sizeof("5'-")-1), out);
		  HelperMethods.indent(iMid+(3), out);

		  long bond = (long)1 << (64-1-Long.numberOfLeadingZeros(overlap));
		  for(; (bond&overlap) != 0; bond>>>=1) 
			  out.print((bond&bonds) != 0 ?'|':'x');
		  out.print('\n'); 
		  HelperMethods.indent(i2,out);
		  out.print("3'-"); 
		  p2.printBases64(out);
		  out.print("-5'\n");
		}
	
	void printBases64(PrintStream out) {
		  long i = (long)1 << (64-1-Long.numberOfLeadingZeros(this.valid));
		  for(; (i&this.valid) != 0 ; i >>>=1)
		    out.print(
		         (this.AorT & i) != 0 ?
		         ((this.GorT & i) != 0 ?'T':'A') :
		         ((this.GorT & i) !=  0 ?'G':'C'));
		}
	@Override
	public IPrimer getComplement() {
		Primer64 complement = new Primer64(AorT, GorT, valid);
		complement.GorT = ~(complement.GorT);
		return complement;
	}
	@Override
	public void printBases(PrintStream outstream) {
		printBases64(outstream);	
	}

	
}