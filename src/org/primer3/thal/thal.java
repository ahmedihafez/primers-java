package org.primer3.thal;

import java.util.Stack;

import org.primer3.sequence.Sequence;

public class thal {

	// statics 

	static public double THAL_ERROR_SCORE = Double.NEGATIVE_INFINITY;
	
	/* The maximum length of _one_ of the two sequences being aligned in a
	   thermodynamic alignment. In other words, the length of one sequence
	   must be <= THAL_MAX_ALIGN, but the other sequence can be longer.
	   The rationale behind this value (60) is that this is the maxium
	   reasonable length for nearest neighbor models. It is the maxium
	   length at which we can restrict our model to only two states of
	   melting: fully intact duplex or completely dissociated single
	   strands. */
	static public int THAL_MAX_ALIGN = 60;
	/* The maxium length of the other sequence in a thermodynamic
	   alignment. This value can be increased, though alignments against
	   very long sequences will be quite slow. As of 2012-05-18, we only
	   potentially see sequences longer this when checking for mispriming
	   in the template ('max_template_mispriming') in libprimer3.c, which
	   is really designed to find sites of ectopic primer very close (a
	   few kilobases) from the location of the cadidate primer. */
	static public int THAL_MAX_SEQ  = 10000;
	
	

	
	/**
	 *   minimum size of hairpin loop 
	 */
	static public int MIN_HRPN_LOOP = 3;


	static final double R = 1.9872; /* cal/Kmol */
	static final double ILAS = (-300 / 310.15); /* Internal Loop Entropy ASymmetry correction -0.3kcal/mol*/
	static final double ILAH = 0.0; /* Internal Loop EntHalpy Asymmetry correction */
	static final double AT_H = 2200.0; /* AT penalty */
	static final double AT_S = 6.9; /* AT penalty */
	static final double MinEntropyCutoff = -2500.0; /* to filter out non-existing entropies */
	static final double MinEntropy = -3224.0; /* initiation */
	static final double G2 = 0.0; /* structures w higher G are considered to be unstabile */
	static final double ABSOLUTE_ZERO = 273.15;
	static final double TEMP_KELVIN = 310.15;
	static final int MAX_LOOP = 30; /* the maximum size of loop that can be calculated; for larger loops formula must be implemented */
	static final int MIN_LOOP = 0;
	
	/* matrix for allowed; bp 0 - no bp, watson crick bp - 1 */
	// int BPI[5][5]
	static final int[][] BPI = new int[][]{
	     {0, 0, 0, 1, 0}, /* A, C, G, T, N; */
	     {0, 0, 1, 0, 0},
	     {0, 1, 0, 0, 0},
	     {1, 0, 0, 0, 0},
	     {0, 0, 0, 0, 0}};
	
	// ################################################################################## 

	
	
	
	
	
	Sequence oligo_f;
	Sequence oligo_r;
	thal_args a;
	
	
	
	
	
	
	// for calculation
	int len1, len2, len3; /* length of sequense 1 and 2 *//* 17.02.2009 int temponly;*/ /* print only temperature of the predicted structure */

	double[] send5;
	double[] hend5;
	int[] numSeq1 ;
	int[] numSeq2 ;  /* same as oligo1 and oligo2 but converted to numbers */
	double[] enthalpyDPT; /* matrix for values of enthalpy */
	double[] entropyDPT; /* matrix for values of entropy */
	double SHleft; /* var that helps to find str w highest melting temperature */
	double RC = 0,dplx_init_H = 0 ,dplx_init_S = 0;
	double saltCorrection ;
	char[] oligo1 = null;
	char[] oligo2 = null;
	
	public thal(Sequence oligo_f,Sequence oligo_r, thal_args a){
		this.a = a;
		this.oligo_f =oligo_f;
		this.oligo_r = oligo_r;
	}


	public thal(char[] s1, char[] s2, thal_args a) {
		this.a = a;			
		this.oligo_f = new Sequence(s1);
		if(s1 != s2)
			this.oligo_r = new Sequence(s2);
		else
			this.oligo_r = this.oligo_f;

	}

	public thal(char[] s, thal_args a) {
		this.a = a;
		this.oligo_f = new Sequence(s);
		this.oligo_r = this.oligo_f;

	}
	public thal_results calc_thal() throws ThermodynamicAlignmentException
	{
		thal_results o = new thal_results();
		
		// TODO :: check inputs 
		
		
		int len_f = oligo_f.length();
		int len_r = oligo_r.length();
		 
		if( (len_f > THAL_MAX_ALIGN) && (len_r > THAL_MAX_ALIGN))
		{
			throw new ThermodynamicAlignmentException("Both sequences longer than " + THAL_MAX_ALIGN + " for thermodynamic alignment");
		}
		if(len_f > THAL_MAX_SEQ)
		{
			throw new ThermodynamicAlignmentException("Target sequence length > maximum allowed (" +THAL_MAX_SEQ +") in thermodynamic alignment (f)");
		}
		if(len_r > THAL_MAX_SEQ)
		{
			throw new ThermodynamicAlignmentException("Target sequence length > maximum allowed (" +THAL_MAX_SEQ +") in thermodynamic alignment (r)");
		}
		
		if (len_f == 0) {
			o.msg = "Empty first sequence";
			o.temp = 0.0;
			return o;
		}
		if (len_r == 0) {
			o.msg="Empty second sequence";
		    o.temp = 0.0;
		    return o;
		}
		
		
//		char[] oligo1 = null;
//		char[] oligo2 = null;
		char[] oligo2_rev = null;

		if(a.type != thal_alignment_type.thal_end2) // !=3 
		{
			oligo1 = oligo_f.getSequence();
			oligo2 = oligo_r.getSequence();
		}
		else
		{
			oligo1 = oligo_r.getSequence();
			oligo2 = oligo_f.getSequence();
		}
		

		
		
		/*** INIT values for unimolecular and bimolecular structures ***/
		if(a.type == thal_alignment_type.thal_hairpin) // == 4
		{ /* unimolecular folding */
			len2 = oligo2.length;
		    len3 = len2 -1;
		    dplx_init_H = 0.0;
		    dplx_init_S = -0.00000000001;
		    RC = 0;
		} 
		else if (a.type != thal_alignment_type.thal_hairpin)  // != 4
		{
			/* hybridization of two oligos */
			dplx_init_H = 200;
		    dplx_init_S = -5.7;
		    // (symmetry_thermo(oligo1) && symmetry_thermo(oligo2)
		    if(oligo_f.symmetry() && oligo_r.symmetry()) {
		    	RC = R  * Math.log(a.dna_conc/1000000000.0);
		    } else {
		    	RC = R  * Math.log(a.dna_conc/4000000000.0);
		    }
		    
		    if(a.type != thal_alignment_type.thal_end2) // != 3
		    {
		    	oligo2_rev = oligo_r.getReverse().getSequence(); 
		    } else {
		    	oligo2_rev = oligo_f.getReverse().getSequence();
		    }
//		    reverse(oligo2_rev); /* REVERSE oligo2, so it goes to dpt 3'->5' direction */
//		    free(oligo2);
//		    oligo2=NULL;
		    oligo2=oligo2_rev;	
		}
		else
		{
			// Does not make any sense here
//			strcpy(o->msg, "Wrong alignment type!");
//		    o->temp = THAL_ERROR_SCORE;
//		      errno=0;
		}
		
		
		
		len1 = oligo1.length;
		len2 = oligo2.length;
		/* convert nucleotides to numbers */
		numSeq1 = seqToInt(oligo1,1);
		numSeq2 = seqToInt(oligo2,1);// new int[len2+2];
		/*** Calc part of the salt correction ***/
		saltCorrection = saltCorrectS(a.mv,a.dv,a.dntp); /* salt correction for entropy, must be multiplied with N, which is
								   the total number of phosphates in the duplex divided by 2; 8bp dplx N=7 */
		

		
		if(a.type == thal_alignment_type.thal_hairpin) // == 4 /* monomer */
		{
			/* terminal basepairs */
			// init
			send5 = new double[len1+1];
			hend5 = new double[len1+1];
		}
		// make sure that oligo1,oligo2 is uppercase
		// covert seq char into index
//		for(int i = 1; i <= len1; ++i) numSeq1[i] = thallib.str2int(oligo1[i - 1]);
//		for(int i = 1; i <= len2; ++i) numSeq2[i] = thallib.str2int(oligo2[i - 1]);
		numSeq1[0] = numSeq1[len1 + 1] = numSeq2[0] = numSeq2[len2 + 1] = 4; /* mark as N-s */


		
		double mh,ms;
		if(a.type == thal_alignment_type.thal_hairpin) /* calculate structure of monomer */
		{
			enthalpyDPT = new double[len1*len2];
			entropyDPT  = new double[len1*len2];
			initMatrix2();
		    fillMatrix2();
		    calc_terminal_bp(a.temp);
		    mh = hend5[len1];
		    ms = send5[len1];
		    o.align_end_1 = (int) mh;
		    o.align_end_2 = (int) ms;
		    int[] bp = new int[len1];
		    if(Double.isFinite(mh)) {
		    	tracebacku(bp);
		    	/* traceback for unimolecular structure */
		    	drawHairpin(bp, mh, ms, o); /* if temponly=1 then return after printing basic therm data */
		    } else if(a.temponly == 0) {
		    	System.err.println("No secondary structure could be calculated\n");
		    }
		    if(o.temp == Double.NEGATIVE_INFINITY && o.msg.isEmpty() ) 
		    	o.temp=0.0;
		    return o;
		}
		else
		{
			len3 = len2;
			enthalpyDPT = new double[len1*len2];
			entropyDPT  = new double[len1*len2];
			initMatrix();
			fillMatrix();
			double[] SH = new double[2];
			int bestI = 0, bestJ = 0; 
			double G1 = Double.POSITIVE_INFINITY, bestG = Double.POSITIVE_INFINITY;
			if(a.type == thal_alignment_type.thal_any) // == 1
			{
				for (int i = 1; i <= len1; i++) {
					for (int j = 1; j <= len2; j++) {
						RSH(i, j, SH);
						SH[0] = SH[0]+ Double.MIN_VALUE ;// SMALL_NON_ZERO; /* this adding is done for compiler, optimization -O2 vs -O0 */
						SH[1] = SH[1]+ Double.MIN_VALUE ;//SMALL_NON_ZERO;
						G1 = (EnthalpyDPT(i,j)+ SH[1] + dplx_init_H) - TEMP_KELVIN*(EntropyDPT(i,j) + SH[0] + dplx_init_S);  
						if(G1 < bestG){
							bestG = G1;
							bestI = i;
							bestJ = j;
						}
					}
				}
			}
			
			
			int[] ps1 = new int[len1];
			int[] ps2 = new int[len2];
			
			if(a.type == thal_alignment_type.thal_end1 || a.type ==thal_alignment_type.thal_end2) // 2 or 3
			{
				/* THAL_END1 */
				bestI = bestJ = 0;
				bestI = len1;
				int i = len1;
				G1 = bestG = Double.POSITIVE_INFINITY;
				for (int j = 1; j <= len2; ++j) {
					RSH(i, j, SH);
					SH[0] = SH[0]+Double.MIN_VALUE; /* this adding is done for compiler, optimization -O2 vs -O0,
								   that compiler could understand that SH is changed in this cycle */
					SH[1] = SH[1]+Double.MIN_VALUE;
					G1 = (EnthalpyDPT(i,j)+ SH[1] + dplx_init_H) - TEMP_KELVIN*(EntropyDPT(i,j) + SH[0] + dplx_init_S);  
					if(G1<bestG){
						bestG = G1;
						bestJ = j;
					}
				 }
			}
			if (!Double.isFinite(bestG)) 
				bestI = bestJ = 1;
			double dH, dS;
			RSH(bestI, bestJ, SH);
			dH =  EnthalpyDPT(bestI,bestJ)+ SH[1] + dplx_init_H;
			dS = (EntropyDPT(bestI,bestJ) + SH[0] + dplx_init_S);
			/* tracebacking */
			if(Double.isFinite(EnthalpyDPT(bestI,bestJ))){
				traceback(bestI, bestJ, RC, ps1, ps2);
				drawDimer(ps1, ps2, SHleft, dH, dS, o);
				o.align_end_1=bestI;
				o.align_end_2=bestJ;
			} else  {
				o.temp = 0.0;
				/* fputs("No secondary structure could be calculated\n",stderr); */
			}
		}
		
		return o;
	}


	/**
	 * 
	 * @param charSeq
	 * @param padSize -- number of elements to add left and right 1 mean add 1 to the left and one to the right
	 * @return
	 */
	private static int[] seqToInt(char[] charSeq, int padSize) {

		int[] seq = new int[charSeq.length + padSize*2];
		for(int i = padSize,j=0 ; j < charSeq.length  ;j++, i++)
		{
			seq[i] =  thallib.str2int(Character.toUpperCase(charSeq[j]));
		}
		
		return seq;
	}


	private void drawDimer(int[] ps1, int[] ps2, double temp, double H, double S, thal_results o) {

		int temponly = a.temponly;
		double t37 = a.temp;
		int i, j, k, numSS1, numSS2, N;
//		char* duplex[4];
		double G, t;
		t = G = 0;
		  if (!Double.isFinite(temp)){
			  if(temponly==0) {
				  System.out.println("No predicted secondary structures for given sequences\n");
			  }
			  o.temp = 0.0; /* lets use generalization here; this should rather be very negative value */
			  o.msg= "No predicted sec struc for given seq";
			  return;
		  } else {
			  N=0;
			  for(i=0;i<len1;i++){
				  if(ps1[i]>0) ++N;
			  }
			  for(i=0;i<len2;i++) {
				  if(ps2[i]>0) ++N;
			  }
			  N = (N/2) -1;
			  t = ((H) / (S + (N * saltCorrection) + RC)) - ABSOLUTE_ZERO;
			  if(temponly==0) {
				  G = (H) - (t37 * (S + (N * saltCorrection)));
				  S = S + (N * saltCorrection);
				  o.temp = (double) t;
				  /* maybe user does not need as precise as that */
				  /* printf("Thermodynamical values:\t%d\tdS = %g\tdH = %g\tdG = %g\tt = %g\tN = %d, SaltC=%f, RC=%f\n",
					len1, (double) S, (double) H, (double) G, (double) t, (int) N, saltCorrection, RC); */
				  System.out.format("Calculated thermodynamical parameters for dimer:\tdS = %g\tdH = %g\tdG = %g\tt = %g\n",
						  (double) S, (double) H, (double) G, (double) t);
			  } else {
		    	  o.temp = (double) t;
		    	  return;
		      }
		  }
		  
		  StringBuilder[] duplex = new StringBuilder[4];
		  duplex[0] = new StringBuilder();
		  duplex[1] = new StringBuilder();
		  duplex[2] = new StringBuilder();
		  duplex[3] = new StringBuilder();
		  
		  
		  i = 0;
		  numSS1 = 0;
		  while (ps1[i++] == 0) ++numSS1;
		  j = 0;
		  numSS2 = 0;
		  while (ps2[j++] == 0) ++numSS2;

		  if (numSS1 >= numSS2){
			  for (i = 0; i < numSS1; ++i) {
				  duplex[0].append(oligo1[i]) ;// strcatc(duplex[0], oligo1[i]);
				  duplex[1].append(' ');//strcatc(duplex[1], ' ');
				  duplex[2].append(' ');//strcatc(duplex[2], ' ');
			  }
			  for (j = 0; j < numSS1 - numSS2; ++j) duplex[3].append(' '); //strcatc(duplex[3], ' ');
			  for (j = 0; j < numSS2; ++j) duplex[3].append(oligo2[j]); //strcatc(duplex[3], oligo2[j]);
		  } else {
			  for (j = 0; j < numSS2; ++j) {
				  duplex[3].append(oligo2[j]);
				  duplex[1].append(' ');
				  duplex[2].append(' ');
			  }
			  for (i = 0; i < numSS2 - numSS1; ++i)
				  duplex[0].append(' ');
			  for (i = 0; i < numSS1; ++i)
				  duplex[0].append(oligo1[i]);
		  }
		  i = numSS1 + 1;
		  j = numSS2 + 1;

		  while (i <= len1) {
			  while (i <= len1 && ps1[i - 1] != 0 && j <= len2 && ps2[j - 1] != 0) {
				  duplex[0].append(' ');
				  duplex[1].append(oligo1[i - 1]);
				  duplex[2].append(oligo2[j - 1]);
				  duplex[3].append(' ');
				  ++i;
				  ++j;
			  }
			  numSS1 = 0;
			  while (i <= len1 && ps1[i - 1] == 0) {
				  duplex[0].append(oligo1[i - 1]);
				  duplex[1].append(' ');
				  ++numSS1;
				  ++i;
			  }
			  numSS2 = 0;
			  while (j <= len2 && ps2[j - 1] == 0) {
				  duplex[2].append(' ');
				  duplex[3].append(oligo2[j - 1]);
				  ++numSS2;
				  ++j;
			  }
		      if (numSS1 < numSS2)
		    	  for (k = 0; k < numSS2 - numSS1; ++k) {
		    		  duplex[0].append('-');
		    		  duplex[1].append(' ');
		    	  }
		      else if (numSS1 > numSS2)
		    	  for (k = 0; k < numSS1 - numSS2; ++k) {
		    		  duplex[2].append(' ');
		    		  duplex[3].append('-');
		    	  }
		  }
		   System.out.println("SEQ\t" + duplex[0].toString() );
		   System.out.println("SEQ\t" + duplex[1].toString() );
		   System.out.println("STR\t" + duplex[2].toString() );
		   System.out.println("STR\t" + duplex[3].toString() );
		
	}

	/**
	 *  traceback for dimers 
	 */
	private void traceback(int i, int j, double RT, int[] ps1,
			int[] ps2) {
		int maxLoop = a.maxLoop;
		int d, ii, jj;
		boolean done;
		double[] SH = new double[2];
		ps1[i - 1] = j;
		ps2[j - 1] = i;
		while(true) {
			SH[0] = -1.0;
			SH[1] = Double.POSITIVE_INFINITY;
			LSH(i,j,SH);
			if(equal(EntropyDPT(i,j),SH[0]) && equal(EnthalpyDPT(i,j),SH[1])) {
				break;
			}
			done = false;
			if (i > 1 && j > 1 && equal(EntropyDPT(i,j), Ss(i - 1, j - 1, 1) + EntropyDPT(i - 1,j - 1)) && equal(EnthalpyDPT(i,j), Hs(i - 1, j - 1, 1) + EnthalpyDPT(i - 1,j - 1))) {
				i = i - 1;
				j = j - 1;
				ps1[i - 1] = j;
				ps2[j - 1] = i;
				done = true;
			}
			for (d = 3; !done && d <= maxLoop + 2; ++d) {
				ii = i - 1;
				jj = -ii - d + (j + i);
				if (jj < 1) {
					ii -= Math.abs(jj-1);
					jj = 1;
				}
				for (; !done && ii > 0 && jj < j; --ii, ++jj) {
					SH[0] = -1.0;
					SH[1] = Double.POSITIVE_INFINITY;
					calc_bulge_internal(ii, jj, i, j, SH,true,maxLoop);
					if (equal(EntropyDPT(i,j), SH[0]) && equal(EnthalpyDPT(i,j), SH[1])) {
						i = ii;
						j = jj;
						ps1[i - 1] = j;
						ps2[j - 1] = i;
						done = true;
						break;
					}
				}
			}
		}
	}

	

//	private boolean DBL_EQ(double a, double b) {
//		// TODO Auto-generated method stub
//		return Math.abs(a - b ) < 1e-5;;
//	}
	private boolean equal(double a, double b) {
		
		if(Double.isInfinite(a) || Double.isInfinite(b))
			return false;
 		return Math.abs(a - b ) < 1e-5;
	}
	
	private boolean isPositive(double v) {
		
		return v > 0;
	}

	/**
	 *   prints ascii output of hairpin structure 
	 */
	private void drawHairpin(int[] bp,  double mh, double ms, thal_results o) {
		int temponly = a.temponly;
		double temp = a.temp;
		int i, N;
		N = 0;
		double mg, t;
		if (!Double.isFinite(ms) || !Double.isFinite(mh)) {
			if(temponly == 0) {
				System.out.format("0\tdS = %g\tdH = %g\tinf\tinf\n", (double) ms,(double) mh);
			} else {
				o.temp = 0.0; /* lets use generalization here */
				o.msg= "No predicted sec struc for given seq\n";
			}
		} else {
			if(temponly == 0) {
				for (i = 1; i < len1; ++i) {
					if(bp[i-1] > 0) N++;
				}
			} else {
				for (i = 1; i < len1; ++i) {
					if(bp[i-1] > 0) N++;
				}
			}
			t = (mh / (ms + (((N/2)-1) * saltCorrection))) - ABSOLUTE_ZERO;
			if(temponly == 0) {
				mg = mh - (temp * (ms + (((N/2)-1) * saltCorrection)));
				ms = ms + (((N/2)-1) * saltCorrection);
				o.temp = (double) t;
				System.out.format("Calculated thermodynamical parameters for dimer:\t%d\tdS = %g\tdH = %g\tdG = %g\tt = %g\n",
						len1, (double) ms, (double) mh, (double) mg, (double) t);
			} else {
				o.temp = (double) t;
				return;
			}
		}
		
		 /* plain-text output */
		char[] asciiRow = new char[len1];
		for(i = 0; i < len1; ++i) asciiRow[i] = '0';
		for(i = 1; i < len1+1; ++i) {
			if(bp[i-1] == 0) {
				asciiRow[(i-1)] = '-';
			} else {
				if(bp[i-1] > (i-1)) {
					asciiRow[(bp[i-1]-1)]='\\';
				} else  {
					asciiRow[(bp[i-1]-1)]='/';
				}
			}
		}
		System.out.println("SEQ\t" + new String(asciiRow));		 
		System.out.println("STR\t" + new String(oligo1));
	}

	/** 
	 * traceback for hairpins 
	 */
	private void tracebacku(int[] bp) {
		int maxLoop = a.maxLoop;
		int i, j;
		int ii, jj, k;
		Stack<thal.tracer> stack = new Stack<thal.tracer>();//
		double[] SH1 = new double[2];
		double[] SH2 = new double[2];
		double[] EntropyEnthalpy = new double[2];
		stack.push(new tracer(len1, 0, 1));
		while(!stack.isEmpty())
		{
			tracer top = stack.pop();
			i = top.i;
			j = top.j;
			if(top.mtrx==1) {
				while (equal(send5[i], send5[i - 1]) && equal(hend5[i], hend5[i - 1])) /* if previous structure is the same as this one */
		    	   --i;
				if (i == 0)
					continue;
				if (equal(send5[i], END5_1(i,2)) && equal(hend5[i], END5_1(i,1))) {
					for (k = 0; k <= i - MIN_HRPN_LOOP - 2; ++k)
						if (equal(send5[i], thallib.atpS[numSeq1[k + 1]][numSeq1[i]] + EntropyDPT(k + 1,i)) &&
								equal(hend5[i], thallib.atpH[numSeq1[k + 1]][numSeq1[i]] + EnthalpyDPT(k + 1,i))) {
							stack.push(new tracer(k + 1, i,0));
							break;
						}
						else if (equal(send5[i], send5[k] + thallib.atpS[numSeq1[k + 1]][numSeq1[i]] + EntropyDPT(k + 1,i)) &&
								equal(hend5[i], hend5[k] + thallib.atpH[numSeq1[k + 1]][numSeq1[i]] + EnthalpyDPT(k + 1,i))) {
							stack.push(new tracer(k + 1, i, 0));
							stack.push(new tracer(k, 0, 1));
							break;
						}
				}
				else if (equal(send5[i], END5_2(i,2)) && equal(hend5[i], END5_2(i,1))) {
					for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k)
						if (equal(send5[i], thallib.atpS[numSeq1[k + 2]][numSeq1[i]] + Sd5(i, k + 2) + EntropyDPT(k + 2,i)) &&
								equal(hend5[i], thallib.atpH[numSeq1[k + 2]][numSeq1[i]] + Hd5(i, k + 2) + EnthalpyDPT(k + 2,i))) {
							stack.push(new tracer (k + 2, i, 0));
							break;
						}
						else if (equal(send5[i], send5[k] + thallib.atpS[numSeq1[k + 2]][numSeq1[i]] + Sd5(i, k + 2) + EntropyDPT(k + 2,i)) &&
								equal(hend5[i], hend5[k] + thallib.atpH[numSeq1[k + 2]][numSeq1[i]] + Hd5(i, k + 2) + EnthalpyDPT(k + 2,i))) {
							stack.push(new tracer(k + 2, i, 0));
							stack.push(new tracer(k, 0, 1));
							break;
						}
		    	 }
				else if (equal(send5[i], END5_3(i,2)) && equal(hend5[i], END5_3(i,1))) {
					for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k)
						if (equal(SEND5(i), atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) + Sd3(i - 1, k + 1) + EntropyDPT(k + 1, i - 1))
								&& equal(HEND5(i), atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) + Hd3(i - 1, k + 1) + EnthalpyDPT(k + 1, i - 1))) {
							stack.push(new tracer( k + 1, i - 1, 0));
							break;
						}
						else if (equal(SEND5(i), SEND5(k) + atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) + Sd3(i - 1, k + 1) + EntropyDPT(k + 1, i - 1)) &&
								equal(HEND5(i), HEND5(k) + atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) + Hd3(i - 1, k + 1) + EnthalpyDPT(k + 1, i - 1))) {
							stack.push(new tracer( k + 1, i - 1, 0)); /* matrix 0  */
							stack.push(new tracer( k, 0, 1)); /* matrix 3 */
							break;
						}
				}
				else if(equal(SEND5(i), END5_4(i,2)) && equal(HEND5(i), END5_4(i,1))) {
					for (k = 0; k <= i - MIN_HRPN_LOOP - 4; ++k)
						if (equal(SEND5(i), atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) + Ststack(i - 1, k + 2) + EntropyDPT(k + 2, i - 1)) &&
								equal(HEND5(i), atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) + Htstack(i - 1, k + 2) + EnthalpyDPT(k + 2, i - 1))) {
							stack.push(new tracer( k + 2, i - 1, 0));
							break;
						}
						else if (equal(SEND5(i), SEND5(k) + atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) + Ststack(i - 1, k + 2) + EntropyDPT(k + 2, i - 1)) &&
								equal(HEND5(i), HEND5(k) + atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) + Htstack(i - 1, k + 2) + EnthalpyDPT(k + 2, i - 1)) ) {
							stack.push(new tracer( k + 2, i - 1, 0));
							stack.push(new tracer( k, 0, 1));
							break;
						}
				}
			}
			else if(top.mtrx==0) {
				 bp[i - 1] = j;
				 bp[j - 1] = i;
				 SH1[0] = -1.0;
				 SH1[1] = Double.POSITIVE_INFINITY;
				 calc_hairpin(i, j, SH1, 1); /* 1 means that we use this method in traceback */
				 SH2[0] = -1.0;
				 SH2[1] = Double.POSITIVE_INFINITY;
				 CBI(i,j,SH2,true,maxLoop);
				 if (equal(EntropyDPT(i, j), Ss(i, j, 2) + EntropyDPT(i + 1, j - 1)) &&
						 equal(EnthalpyDPT(i, j), Hs(i, j, 2) + EnthalpyDPT(i + 1, j - 1))) {
					 stack.push(new tracer(i + 1, j - 1, 0));
				 }
				 else if (equal(EntropyDPT(i, j), SH1[0]) && equal(EnthalpyDPT(i,j), SH1[1]))
					 ;
				 else if (equal(EntropyDPT(i, j), SH2[0]) && equal(EnthalpyDPT(i, j), SH2[1])) {
					 int d;
					 boolean done;
					 for (done = false, d = j - i - 3; d >= MIN_HRPN_LOOP + 1 && d >= j - i - 2 - maxLoop && !done; --d)
						 for (ii = i + 1; ii < j - d; ++ii) {
							 jj = d + ii;
							 EntropyEnthalpy[0] = -1.0;
							 EntropyEnthalpy[1] = Double.POSITIVE_INFINITY;
							 calc_bulge_internal2(i, j, ii, jj,EntropyEnthalpy,true,maxLoop);
							 if (equal(EntropyDPT(i, j), EntropyEnthalpy[0] + EntropyDPT(ii, jj)) &&
									 equal(EnthalpyDPT(i, j), EntropyEnthalpy[1] + EnthalpyDPT(ii, jj))) {
								 stack.push(new tracer(ii, jj, 0));
								 done=true;
								 break;
							 }
						 }
				 } else {
				 }
			}
		}
		
	}


	
	private double SEND5(int i) {
		return send5[i] ;
	}

	private double atPenaltyS(int i, int j) {
		return thallib.atpS[i][j];
	}

	private double atPenaltyH(int i, int j) {
		return thallib.atpH[i][j];
	}

	private double EntropyDPT(int i, int j) {
//		return entropyDPT[i][j];
		return entropyDPT[(j) + ((i-1)*len3) - (1)];
	}
	private void setEntropyDPT(int i, int j,double value) {
//		return entropyDPT[i][j];
		entropyDPT[(j) + ((i-1)*len3) - (1)] = value;
	}

	private double EnthalpyDPT(int i, int j) {
		// return enthalpyDPT[i][j];
		return enthalpyDPT[(j) + ((i-1)*len3) - (1)];
	}
	private void setEnthalpyDPT(int i, int j,double value) {
		// return enthalpyDPT[i][j];
		enthalpyDPT[(j) + ((i-1)*len3) - (1)] = value;
	}
	private double HEND5(int i) {
		return hend5[i];
	}






	/* executed in calc_terminal_bp; to find structure that corresponds to max Tm for terminal bp */
	private double END5_1(int i,int hs) {
		int k;
		double max_tm; /* energy min */
		double T1, T2;
		double H, S;
		double H_max, S_max;
		H_max = H = Double.POSITIVE_INFINITY;
		S_max = S = -1.0;
		T1 = T2 = Double.NEGATIVE_INFINITY;
		max_tm = Double.NEGATIVE_INFINITY;
		for(k = 0; k <= i - MIN_HRPN_LOOP - 2; ++k) {
			T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
			T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC);
			if(T1 >= T2) {
				H = HEND5(k) + atPenaltyH(numSeq1[k + 1], numSeq1[i]) + EnthalpyDPT(k + 1, i);
				S = SEND5(k) + atPenaltyS(numSeq1[k + 1], numSeq1[i]) + EntropyDPT(k + 1, i);
				if(!Double.isFinite(H) || H > 0 || S > 0) { /* H and S must be greater than 0 to avoid BS */
					H = Double.POSITIVE_INFINITY;
					S = -1.0;
				}
				T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
			} else {
				H = 0 + atPenaltyH(numSeq1[k + 1], numSeq1[i]) + EnthalpyDPT(k + 1, i);
				S = 0 + atPenaltyS(numSeq1[k + 1], numSeq1[i]) + EntropyDPT(k + 1, i);
				if(!Double.isFinite(H) || H > 0 || S > 0) {
					H = Double.POSITIVE_INFINITY;
					S = -1.0;
				}
				T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
			}
			if(max_tm < T1) {
				if(S > MinEntropyCutoff) {
					H_max = H;
					S_max = S;
					max_tm = T1;
				}
			}
		}
		if (hs == 1) 
			return H_max;
		return S_max;
	}
	private double END5_2(int i, int hs) {
		int k;
		double max_tm;
		double T1, T2;
		double H, S;
		double H_max, S_max;
		H_max = H = Double.POSITIVE_INFINITY;
		T1 = T2 = max_tm = Double.NEGATIVE_INFINITY;
		S_max = S = -1.0;
		for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k) {
			T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
			T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC);
			if(T1 >= T2) {
				H = HEND5(k) + atPenaltyH(numSeq1[k + 2], numSeq1[i]) + Hd5(i, k + 2) + EnthalpyDPT(k + 2, i);
				S = SEND5(k) + atPenaltyS(numSeq1[k + 2], numSeq1[i]) + Sd5(i, k + 2) + EntropyDPT(k + 2, i);
				if(!Double.isFinite(H) || H > 0 || S > 0) {
					H = Double.POSITIVE_INFINITY;
					S = -1.0;
				}
				T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
			} else {
				H = 0 + atPenaltyH(numSeq1[k + 2], numSeq1[i]) + Hd5(i, k + 2) + EnthalpyDPT(k + 2, i);
				S = 0 + atPenaltyS(numSeq1[k + 2], numSeq1[i]) + Sd5(i, k + 2) + EntropyDPT(k + 2, i);
				if(!Double.isFinite(H) || H > 0 || S > 0) {
					H = Double.POSITIVE_INFINITY;
					S = -1.0;
				}
				T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
			}
			if(max_tm < T1) {
				if(S > MinEntropyCutoff) {
					H_max = H;
					S_max = S;
					max_tm = T1;
				}
			}
		}
		if (hs == 1) return H_max;
		return S_max;
	}

	private double END5_3(int i, int hs) {
		int k;
		double max_tm;
		double T1, T2;
		double H, S;
		double H_max, S_max;
		H_max = H = Double.POSITIVE_INFINITY;;
		T1 = T2 = max_tm = Double.NEGATIVE_INFINITY;
		S_max = S = -1.0;
		for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k) {
			T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
			T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC);
			if(T1 >= T2) {
				H = HEND5(k) + atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) + Hd3(i - 1, k + 1) + EnthalpyDPT(k + 1, i - 1);
				S = SEND5(k) + atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) + Sd3(i - 1, k + 1) + EntropyDPT(k + 1, i - 1);
				if(!Double.isFinite(H) || H > 0 || S > 0) {
					H = Double.POSITIVE_INFINITY;
					S = -1.0;
				}
				T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
			} else {
				H = 0 + atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) + Hd3(i - 1, k + 1) + EnthalpyDPT(k + 1, i - 1);
				S = 0 + atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) + Sd3(i - 1, k + 1) + EntropyDPT(k + 1, i - 1);
				if(!Double.isFinite(H) || H > 0 || S > 0) {
					H = Double.POSITIVE_INFINITY;
					S = -1.0;
				}
				T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
			}
			if(max_tm < T1) {
				if(S > MinEntropyCutoff) {
					H_max = H;
					S_max = S;
					max_tm = T1;
				}
			}
		}
		if (hs == 1) return H_max;
		return S_max;
	}
	private double END5_4(int i, int hs) {
		int k;
		double max_tm;
		double T1, T2;
		double H, S;
		double H_max, S_max;
		H_max = H = Double.POSITIVE_INFINITY;
		T1 = T2 = max_tm = Double.NEGATIVE_INFINITY;
		S_max = S = -1.0;
		for(k = 0; k <= i - MIN_HRPN_LOOP - 4; ++k) {
			T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
			T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC);
			if(T1 >= T2) {
				H = HEND5(k) + atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) + Htstack(i - 1, k + 2) + EnthalpyDPT(k + 2, i - 1);
				S = SEND5(k) + atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) + Ststack(i - 1, k + 2) + EntropyDPT(k + 2, i - 1);
				if(!Double.isFinite(H) || H > 0 || S > 0) {
					H = Double.POSITIVE_INFINITY;
					S = -1.0;
				}
				T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
			} else {
				H = 0 + atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) + Htstack(i - 1, k + 2) + EnthalpyDPT(k + 2, i - 1);
				S = 0 + atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) + Ststack(i - 1, k + 2) + EntropyDPT(k + 2, i - 1);
				if(!Double.isFinite(H) || H > 0 || S > 0) {
					H = Double.POSITIVE_INFINITY;
					S = -1.0;
				}
				T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
			}
			if(max_tm < T1) {
				if(S > MinEntropyCutoff) {
					H_max = H;
					S_max = S;
					max_tm = T1;
				}
			}
		}
		if (hs == 1) return H_max;
		return S_max;
	}






	private double Sd5(int i, int j) {
		return thallib.dangleEntropies5[numSeq1[i]][numSeq1[j]][numSeq1[j - 1]];
	}
	
	private double Hd5(int i, int j) {
		return thallib.dangleEnthalpies5[numSeq1[i]][numSeq1[j]][numSeq1[j - 1]];
	}

	private double Sd3(int i, int j) {
		return thallib.dangleEntropies3[numSeq1[i]][numSeq1[i+1]][numSeq1[j]];

	}
	private double Hd3(int i, int j) {
		return thallib.dangleEnthalpies3[numSeq1[i]][numSeq1[i+1]][numSeq1[j]];
	}
	
	

	private double Ststack(int i, int j) {
		return thallib.tstack2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq1[j]][numSeq1[j-1]];
	}
	private double Htstack(int i, int j) {
		return thallib.tstack2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq1[j]][numSeq1[j-1]];
	}
	


	

	/**
	 *  finds monomer structure that has maximum Tm
	 */
	private void calc_hairpin(int i, int j, double[] EntropyEnthalpy, int traceback) {
		int loopSize = j - i - 1;
		double G1, G2;
		G1 = G2 = Double.NEGATIVE_INFINITY;
		double[] SH=new double[2];
		SH[0] = -1.0;
		SH[1] = Double.POSITIVE_INFINITY;
		if(loopSize < MIN_HRPN_LOOP) {
			EntropyEnthalpy[0] = -1.0;
			EntropyEnthalpy[1] = Double.POSITIVE_INFINITY;
			return;
		}
		if (i <= len1 && len2 < j) {
			EntropyEnthalpy[0] = -1.0;
			EntropyEnthalpy[1] = Double.POSITIVE_INFINITY;
			return;
		} else if (i > len2) {
			i -= len1;
			j -= len2;
		}
		if(loopSize <= 30) {
			EntropyEnthalpy[1] = thallib.hairpinLoopEnthalpies[loopSize - 1];
			EntropyEnthalpy[0] = thallib.hairpinLoopEntropies[loopSize - 1];
		} else {
			EntropyEnthalpy[1] = thallib.hairpinLoopEnthalpies[29];
			EntropyEnthalpy[0] = thallib.hairpinLoopEntropies[29];
		}
		if (loopSize > 3) { /* for loops 4 bp and more in length, terminal mm are accounted */
			EntropyEnthalpy[1] += thallib.tstack2Enthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq1[j]][numSeq1[j - 1]];
			EntropyEnthalpy[0] += thallib.tstack2Entropies[numSeq1[i]][numSeq1[i + 1]][numSeq1[j]][numSeq1[j - 1]];
		} else if(loopSize == 3){ /* for loops 3 bp in length at-penalty is considered */
			EntropyEnthalpy[1] += atPenaltyH(numSeq1[i], numSeq1[j]);
			EntropyEnthalpy[0] += atPenaltyS(numSeq1[i], numSeq1[j]);
		}

		if (loopSize == 3) {	 /* closing AT-penalty (+), triloop bonus, hairpin of 3 (+) */
//			triloop loop = null;
//			if (numTriloops) 
			{
				// (loop = (struct triloop*) bsearch(numSeq1 + i, triloopEnthalpies, numTriloops, sizeof(struct triloop), comp3loop))
				int loopKey = getHashkey(numSeq1,i,5);
				if (thallib.triloopEnthalpies.containsKey(loopKey))
					EntropyEnthalpy[1] += thallib.triloopEnthalpies.get(loopKey);
				// (loop = (struct triloop*) bsearch(numSeq1 + i, triloopEntropies, numTriloops, sizeof(struct triloop), comp3loop))
				if (thallib.triloopEntropies.containsKey(loopKey))
					EntropyEnthalpy[0] += thallib.triloopEntropies.get(loopKey);
			}
		} else if (loopSize == 4) { /* terminal mismatch, tetraloop bonus, hairpin of 4 */
//			tetraloop loop = null;
//			if (numTetraloops) 
			{
				//(loop = (struct tetraloop*) bsearch(numSeq1 + i, tetraloopEnthalpies, numTetraloops, sizeof(struct tetraloop), comp4loop))
				int loopKey = getHashkey(numSeq1,i,6);
				if (thallib.tetraloopEnthalpies.containsKey(loopKey)) {
					EntropyEnthalpy[1] += thallib.tetraloopEnthalpies.get(loopKey);
				}
				// (loop = (struct tetraloop*) bsearch(numSeq1 + i, tetraloopEntropies, numTetraloops, sizeof(struct tetraloop), comp4loop))
				if (thallib.tetraloopEntropies.containsKey(loopKey)) {
					EntropyEnthalpy[0] += thallib.tetraloopEntropies.get(loopKey);
				}
			}
		}
		if(!Double.isFinite(EntropyEnthalpy[1])) {
			EntropyEnthalpy[1] = Double.POSITIVE_INFINITY;
			EntropyEnthalpy[0] = -1.0;
		}
		if(isPositive(EntropyEnthalpy[1]) && isPositive(EntropyEnthalpy[0]) && (!isPositive(EnthalpyDPT(i, j)) || !isPositive(EntropyDPT(i, j)))) { /* if both, S and H are positive */
			EntropyEnthalpy[1] = Double.POSITIVE_INFINITY;
			EntropyEnthalpy[0] = -1.0;
		}
		RSH(i,j,SH);
		G1 = EntropyEnthalpy[1]+SH[1] -TEMP_KELVIN*(EntropyEnthalpy[0]+SH[0]);
		G2 = EnthalpyDPT(i, j)+SH[1] -TEMP_KELVIN*(EntropyDPT(i, j)+SH[0]);
		if(G2 < G1 && traceback == 0) {
			EntropyEnthalpy[0] = EntropyDPT(i, j);
			EntropyEnthalpy[1] = EnthalpyDPT(i, j);
		}		
	}





	private void CBI(int i, int j, double[] EntropyEnthalpy, boolean traceback, int maxLoop) {
		for (int d = j - i - 3; d >= MIN_HRPN_LOOP + 1 && d >= j - i - 2 - maxLoop; --d)
			for (int ii = i + 1; ii < j - d && ii <= len1; ++ii) {
				int jj = d + ii;
				if(traceback==false) {
					EntropyEnthalpy[0] = -1.0;
					EntropyEnthalpy[1] = Double.POSITIVE_INFINITY;
				}
				if (Double.isFinite(EnthalpyDPT(ii,jj)) && Double.isFinite(EnthalpyDPT(i,j))) {
					calc_bulge_internal2(i, j, ii, jj, EntropyEnthalpy, traceback,maxLoop);
					if(Double.isFinite(EntropyEnthalpy[1])) {
						if(EntropyEnthalpy[0] < MinEntropyCutoff) {
							EntropyEnthalpy[0] = MinEntropy;
							EntropyEnthalpy[1] = 0.0;
						}
						if(traceback==false) {
							setEnthalpyDPT(i,j,EntropyEnthalpy[1]) ;
							setEntropyDPT(i,j, EntropyEnthalpy[0]) ;
						}
					}
				}
			}		
	}

	private void calc_bulge_internal2(int i, int j, int ii, int jj, 
			double[] EntropyEnthalpy, boolean traceback, int maxLoop) {
		 int loopSize1, loopSize2, loopSize;
		 double T1, T2;
		 double S,H;
		 /* int N, N_loop; Triinu, please review */
		 T1 = T2 = Double.NEGATIVE_INFINITY;
		 S = MinEntropy;
		 H = 0.0;
		 loopSize1 = ii - i - 1;
		 loopSize2 = j - jj - 1;
		 if (loopSize1 + loopSize2 > maxLoop) {
			 EntropyEnthalpy[0] = -1.0;
			 EntropyEnthalpy[1] = Double.POSITIVE_INFINITY;
			 return;
		 }
//		 #ifdef DEBUG
//		   if (ii <= i)
//		     fputs("Error in calc_bulge_internal(): ii isn't greater than i\n", stderr);
//		   if (jj >= j)
//		     fputs("Error in calc_bulge_internal(): jj isn't less than j\n", stderr);
//		   if (ii >= jj)
//		     fputs("Error in calc_bulge_internal(): jj isn't greater than ii\n", stderr);
//
//		   if ((i <= len1 && len1 < ii) || (jj <= len2 && len2 < j))  {
//		      EntropyEnthalpy[0] = -1.0;
//		      EntropyEnthalpy[1] = _INFINITY;
//		      return;
//		   }
//		#endif
//		 #ifdef DEBUG
//		   if (loopSize1 + loopSize2 > maxLoop) {
//		      fputs("Error: calc_bulge_internal() called with loopSize1 + loopSize2 > maxLoop\n", stderr);
//		      return;
//		   }
//		#endif
		 
		 
//		#ifdef DEBUG
//		   if (loopSize1 == 0 && loopSize2 == 0) {
//		      fputs("Error: calc_bulge_internal() called with nonsense\n", stderr);
//		      return;
//		   }
//		#endif

//		#ifdef DEBUG
//		   if (i > len1)
//		     i -= len1;
//		   if (ii > len1)
//		     ii -= len1;
//		   if (j > len2)
//		     j -= len2;
//		   if (jj > len2)
//		     jj -= len2;
//		#endif 
		 loopSize = loopSize1 + loopSize2 -1; /* for indx only */
		 if((loopSize1 == 0 && loopSize2 > 0) || (loopSize2 == 0 && loopSize1 > 0)) { /* only bulges have to be considered */
			 if(loopSize2 == 1 || loopSize1 == 1) { 
				 /* bulge loop of size one is treated differently
					the intervening nn-pair must be added */
				 if((loopSize2 == 1 && loopSize1 == 0) || (loopSize2 == 0 && loopSize1 == 1)) {
					 H = thallib.bulgeLoopEnthalpies[loopSize] +
							 thallib.stackEnthalpies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
					 S = thallib.bulgeLoopEntropies[loopSize] +
							 thallib.stackEntropies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
				 }
				 if(!traceback) {
					 H += EnthalpyDPT(ii, jj); /* bulge koos otsaga, st bulge i,j-ni */
					 S += EntropyDPT(ii, jj);
				 }

				 if(!Double.isFinite(H)) {
					 H = Double.POSITIVE_INFINITY;
					 S = -1.0;
				 }
			
				 T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
				 T2 = (EnthalpyDPT(i, j) + dplx_init_H) / ((EntropyDPT(i, j)) + dplx_init_S + RC);

				 if((T1 > T2) || ((traceback && T1 >= T2) || traceback)) {
					 EntropyEnthalpy[0] = S;
					 EntropyEnthalpy[1] = H;
				 }

			 } else { /* we have _not_ implemented Jacobson-Stockaymayer equation; the maximum bulgeloop size is 30 */

				 H = thallib.bulgeLoopEnthalpies[loopSize] + atPenaltyH(numSeq1[i], numSeq2[j]) + atPenaltyH(numSeq1[ii], numSeq2[jj]);
				 if(!traceback)
					 H += EnthalpyDPT(ii, jj);

				 S = thallib.bulgeLoopEntropies[loopSize] + atPenaltyS(numSeq1[i], numSeq2[j]) + atPenaltyS(numSeq1[ii], numSeq2[jj]);
				 if(!traceback)
					 S += EntropyDPT(ii, jj);
				 if(!Double.isFinite(H)) {
					 H = Double.POSITIVE_INFINITY;
					 S = -1.0;
				 }
			  
				 T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
				 T2 = (EnthalpyDPT(i, j) + dplx_init_H) / (EntropyDPT(i, j) + dplx_init_S + RC);

				 if((T1 > T2) || ((traceback && T1 >= T2) || (traceback))) {
					 EntropyEnthalpy[0] = S;
					 EntropyEnthalpy[1] = H;
				 }
			 }
		 } /* end of calculating bulges */
		 else if (loopSize1 == 1 && loopSize2 == 1) {
			 /* mismatch nearest neighbor parameters */

			 S = thallib.stackint2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]] +
		    		  thallib.stackint2Entropies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]];
			 if(!traceback)
				 S += EntropyDPT(ii, jj);

			 H = thallib.stackint2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]] +
		    		  thallib.stackint2Enthalpies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]];
			 if(!traceback)
				 H += EnthalpyDPT(ii, jj);
			 if(!Double.isFinite(H)) {
				 H = Double.POSITIVE_INFINITY;
				 S = -1.0;
			 }
		    
			 T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
			 T2 = (EnthalpyDPT(i, j) + dplx_init_H) / (EntropyDPT(i, j) + dplx_init_S + RC);
			 // *** DBL_EQ(T1,T2) == 2
			 if((!equal (T1,T2) ) || traceback) {
				 if((T1 > T2) || ((traceback && T1 >= T2) || traceback)) {
					 EntropyEnthalpy[0] = S;
					 EntropyEnthalpy[1] = H;
				 }
			 }
			 return;
		 } else { /* only internal loops */

			 H = thallib.interiorLoopEnthalpies[loopSize] + thallib.tstackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]] +
					 thallib.tstackEnthalpies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]]
							 + (ILAH * Math.abs(loopSize1 - loopSize2));
			 if(!traceback)
				 H += EnthalpyDPT(ii, jj);

			 S = thallib.interiorLoopEntropies[loopSize] + thallib.tstackEntropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]] +
		    		  thallib.tstackEntropies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]] + (ILAS * Math.abs(loopSize1 - loopSize2));
			 if(!traceback)
				 S += EntropyDPT(ii, jj);
			 if(!Double.isFinite(H)) {
				 H = Double.POSITIVE_INFINITY;
				 S = -1.0;
			 }
		    
			 T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
			 T2 = (EnthalpyDPT(i, j) + dplx_init_H) / ((EntropyDPT(i, j)) + dplx_init_S + RC);
			 if((T1 > T2) || ((traceback && T1 >= T2) || (traceback))) {
				 EntropyEnthalpy[0] = S;
				 EntropyEnthalpy[1] = H;
			 }
		 }
		 
	}





	/**
	 *  finds max Tm while filling the dyn progr table using stacking S and stacking H (monomer) 
	 */
	private void maxTM2(int i, int j) {
		double T0, T1;
		double S0, S1;
		double H0, H1;
		T0 = T1 = Double.NEGATIVE_INFINITY;
		S0 = EntropyDPT(i,j);
		H0 = EnthalpyDPT(i,j);
		T0 = (H0 + dplx_init_H) /(S0 + dplx_init_S + RC);
		if(Double.isFinite(EnthalpyDPT(i,j))) {
			S1 = (EntropyDPT(i + 1, j - 1) + Ss(i, j, 2));
			H1 = (EnthalpyDPT(i + 1,j - 1) + Hs(i, j, 2));
		} else {
			S1 = -1.0;
			H1 = Double.POSITIVE_INFINITY;
		}
		T1 = (H1 + dplx_init_H) /(S1 + dplx_init_S + RC);
		if(S1 < MinEntropyCutoff) {
			S1 = MinEntropy;
			H1 = 0.0;
		}
		if(S0 < MinEntropyCutoff) {
			S0 = MinEntropy;
			H0 = 0.0;
		}

		if(T1 > T0) {
			setEntropyDPT(i,j,S1 );
			setEnthalpyDPT(i,j, H1);
		} else {
			setEntropyDPT(i,j,S0);
			setEnthalpyDPT(i,j,  H0);
		}
		
	}


	

	/**
	 * terminal bp for monomer structure
	 * @throws ThermodynamicAlignmentException 
	 */
	private void calc_terminal_bp(double temp) throws ThermodynamicAlignmentException {
		int i;
		int max;
		send5[0] = send5[1] = -1.0;
		hend5[0] = hend5[1] = Double.POSITIVE_INFINITY;
		for(i = 2; i<=(len1); i++) {
			send5[i] = MinEntropy;
			hend5[i] = 0;
		}

		double T1, T2, T3, T4, T5;
		T1 = T2 = T3 = T4 = T5 = Double.NEGATIVE_INFINITY;
		double G;
		/* adding terminal penalties to 3' end and to 5' end */
		for(i = 2; i <= len1; ++i) {
			max = 0;
			T1 = T2 = T3 = T4 = T5 = Double.NEGATIVE_INFINITY;
			T1 = (HEND5(i - 1) + dplx_init_H) / (SEND5(i - 1) + dplx_init_S + RC);
			T2 = (END5_1(i,1) + dplx_init_H) / (END5_1(i,2) + dplx_init_S + RC);
			T3 = (END5_2(i,1) + dplx_init_H) / (END5_2(i,2) + dplx_init_S + RC);
			T4 = (END5_3(i,1) + dplx_init_H) / (END5_3(i,2) + dplx_init_S + RC);
			T5 = (END5_4(i,1) + dplx_init_H) / (END5_4(i,2) + dplx_init_S + RC);
			max = max5(T1,T2,T3,T4,T5);
			switch (max) {
			case 1:
				send5[i] = SEND5(i - 1);
				hend5[i] = HEND5(i - 1);
				break;
			case 2:
				G = END5_1(i,1) - (temp * (END5_1(i,2)));
				if(G < G2) {
					send5[i] = END5_1(i,2);
					hend5[i] = END5_1(i,1);
				} else {
					send5[i] = SEND5(i - 1);
					hend5[i] = HEND5(i - 1);
				}
				break;
			case 3:
				G = END5_2(i,1) - (temp * (END5_2(i,2)));
				if(G < G2) {
					send5[i] = END5_2(i,2);
					hend5[i] = END5_2(i,1);
				} else {
					send5[i] = SEND5(i - 1);
					hend5[i] = HEND5(i - 1);
				}
				break;
			case 4:
				G = END5_3(i,1) - (temp * (END5_3(i,2)));
				if(G < G2) {
					send5[i] = END5_3(i,2);
					hend5[i] = END5_3(i,1);
				} else {
					send5[i] = SEND5(i - 1);
					hend5[i] = HEND5(i - 1);
				}
				break;
			case 5:
				G = END5_4(i,1) - (temp * (END5_4(i,2)));
				if(G < G2) {
					send5[i] = END5_4(i,2);
					hend5[i] = END5_4(i,1);
				} else {
					send5[i] = SEND5(i - 1);
					hend5[i] = HEND5(i - 1);
				}
				break;
			default:
//		#ifdef DEBUG
//			 printf ("WARNING: max5 returned character code %d ??\n", max);
//		#endif
				throw new ThermodynamicAlignmentException("returned character code is not in range 1-5");
			}
		}		
	}


	
	private int max5(double a, double b, double c, double d, double e) {
		   if(a > b && a > c && a > d && a > e) return 1;
		   else if(b > c && b > d && b > e) return 2;
		   else if(c > d && c > e) return 3;
		   else if(d > e) return 4;
		   else return 5;
	}


	/**
	 * calculates bulges and internal loops for dimer structures
	 * @param ii
	 * @param jj
	 * @param i
	 * @param j
	 * @param EntropyEnthalpy
	 * @param k
	 * @param maxLoop
	 */
	private void calc_bulge_internal(int i, int j, int ii, int jj, double[] EntropyEnthalpy,
			boolean traceback, int maxLoop) {
		int loopSize1, loopSize2, loopSize;
		double S,H,G1,G2;
		int N, N_loop;
		double[] SH = new double[2];
		SH[0] = -1.0;
		SH[1] = Double.POSITIVE_INFINITY;
		S = -1.0;
		H = Double.POSITIVE_INFINITY;
		loopSize1 = ii - i - 1;
		loopSize2 = jj - j - 1;
		if(ii < jj) {
			N = ((2 * i)/2);
			N_loop = N;
			if(loopSize1 > 2) N_loop -= (loopSize1 - 2);
			if(loopSize2 > 2) N_loop -= (loopSize2 - 2);
		} else {
			N = ((2 * j)/2);
			N_loop = 2 * jj;
			if(loopSize1 > 2) N_loop -= (loopSize1 - 2);
			if(loopSize2 > 2) N_loop -= (loopSize2 - 2);
			N_loop = (N_loop/2) - 1;
		}
//		#ifdef DEBUG
//		   if (ii <= i){
//		      fputs("Error in calc_bulge_internal(): ii is not greater than i\n", stderr);
//		   }
//		   if (jj <= j)
//		     fputs("Error in calc_bulge_internal(): jj is not greater than j\n", stderr);
//		#endif

//		#ifdef DEBUG
//		   if (loopSize1 + loopSize2 > maxLoop) {
//		      fputs("Error: calc_bulge_internal() called with loopSize1 + loopSize2 > maxLoop\n", stderr);
//		      free(SH);
//		      return;
//		   }
//		#endif
//		#ifdef DEBUG
//		   if (loopSize1 == 0 && loopSize2 == 0) {
//		      fputs("Error: calc_bulge_internal() called with nonsense\n", stderr);
//		      free(SH);
//		      return;
//		   }
//		#endif
		loopSize = loopSize1 + loopSize2-1;
		if((loopSize1 == 0 && loopSize2 > 0) || (loopSize2 == 0 && loopSize1 > 0)) { /* only bulges have to be considered */
			if(loopSize2 == 1 || loopSize1 == 1) { /* bulge loop of size one is treated differently
							      the intervening nn-pair must be added */

				if((loopSize2 == 1 && loopSize1 == 0) || (loopSize2 == 0 && loopSize1 == 1)) {
					H = thallib.bulgeLoopEnthalpies[loopSize] +
							thallib.stackEnthalpies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
					S = thallib.bulgeLoopEntropies[loopSize] +
							thallib.stackEntropies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
				}
				if(isPositive(H) || isPositive(S)){
					H = Double.POSITIVE_INFINITY;
					S = -1.0;
				}
				H += EnthalpyDPT(i, j);
			 S += EntropyDPT(i, j);
			 if(!Double.isFinite(H)) {
			    H = Double.POSITIVE_INFINITY;
			    S = -1.0;
			 }
			 RSH(ii,jj,SH);
			 G1 = H+SH[1] -TEMP_KELVIN*(S+SH[0]);
			 G2 = EnthalpyDPT(ii, jj)+SH[1]-TEMP_KELVIN*((EntropyDPT(ii, jj)+SH[0]));
			 if((G1< G2) || (traceback)) {
			    EntropyEnthalpy[0] = S;
			    EntropyEnthalpy[1] = H;
			 }
			} else { /* we have _not_ implemented Jacobson-Stockaymayer equation; the maximum bulgeloop size is 30 */

				H = thallib.bulgeLoopEnthalpies[loopSize] + atPenaltyH(numSeq1[i], numSeq2[j]) + atPenaltyH(numSeq1[ii], numSeq2[jj]);
				H += EnthalpyDPT(i, j);

				S = thallib.bulgeLoopEntropies[loopSize] + atPenaltyS(numSeq1[i], numSeq2[j]) + atPenaltyS(numSeq1[ii], numSeq2[jj]);
				S += EntropyDPT(i, j);
				if(!Double.isFinite(H)) {
					H = Double.POSITIVE_INFINITY;
					S = -1.0;
				}
				if(isPositive(H) && isPositive(S)){ 
					H = Double.POSITIVE_INFINITY;
					S = -1.0;
				}
			
				RSH(ii,jj,SH);
				G1 = H+SH[1] -TEMP_KELVIN*(S+SH[0]);
				G2 = EnthalpyDPT(ii, jj)+SH[1]-TEMP_KELVIN*(EntropyDPT(ii, jj)+SH[0]);
				if(G1< G2 || (traceback)){
					EntropyEnthalpy[0] = S;
					EntropyEnthalpy[1] = H;
				}
			 
			}
		} else if (loopSize1 == 1 && loopSize2 == 1) {
			S = thallib.stackint2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
					thallib.stackint2Entropies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]];
			S += EntropyDPT(i, j);

			H = thallib.stackint2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
					thallib.stackint2Enthalpies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]];
			H += EnthalpyDPT(i, j);
			if(!Double.isFinite(H)) {
				H = Double.POSITIVE_INFINITY;
				S = -1.0;
			}
			if(isPositive(H) && isPositive(S)){
				H = Double.POSITIVE_INFINITY;
				S = -1.0;
			}    
			RSH(ii,jj,SH);
			G1 = H+SH[1] -TEMP_KELVIN*(S+SH[0]);
			G2 = EnthalpyDPT(ii, jj)+SH[1]-TEMP_KELVIN*(EntropyDPT(ii, jj)+SH[0]);
			if((G1< G2) || traceback) {
				EntropyEnthalpy[0] = S;
				EntropyEnthalpy[1] = H;
			}
			return;
		} else { /* only internal loops */
			H = thallib.interiorLoopEnthalpies[loopSize] + thallib.tstackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
					thallib.tstackEnthalpies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]]
							+ (ILAH * Math.abs(loopSize1 - loopSize2));
			H += EnthalpyDPT(i, j);

			S = thallib.interiorLoopEntropies[loopSize] + thallib.tstackEntropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
		    		  thallib.tstackEntropies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]] + (ILAS * Math.abs(loopSize1 - loopSize2));
			S += EntropyDPT(i, j);
			if(!Double.isFinite(H)) {
				H = Double.POSITIVE_INFINITY;
				S = -1.0;
			}
			if(isPositive(H) && isPositive(S)){ 
				H = Double.POSITIVE_INFINITY;
				S = -1.0;
			}
			RSH(ii,jj,SH);
			G1 = H+SH[1] -TEMP_KELVIN*(S+SH[0]);
			G2 = EnthalpyDPT(ii, jj)+SH[1]-TEMP_KELVIN*(EntropyDPT(ii, jj)+SH[0]);
			if((G1< G2) || (traceback)){
				EntropyEnthalpy[0] = S;
				EntropyEnthalpy[1] = H;
			}
		}		
	}


	/**
	 * finds max Tm while filling the dyn progr table using stacking S and stacking H (dimer) 
	 */
	private void maxTM(int i, int j) {
		double T0, T1;
		double S0, S1;
		double H0, H1;
		double[] SH = new double[2];
		T0 = T1 = Double.NEGATIVE_INFINITY;
		S0 = EntropyDPT(i,j);
		H0 = EnthalpyDPT(i,j);
		RSH(i,j,SH);
		T0 = (H0 + dplx_init_H + SH[1]) /(S0 + dplx_init_S + SH[0] + RC); /* at current position */
		if(Double.isFinite(EnthalpyDPT(i - 1,j - 1)) && Double.isFinite(Hs(i - 1, j - 1, 1))) {
			S1 = (EntropyDPT(i - 1,j - 1) + Ss(i - 1, j - 1, 1));
			H1 = (EnthalpyDPT(i - 1,j - 1) + Hs(i - 1, j - 1, 1));
			T1 = (H1 + dplx_init_H + SH[1]) /(S1 + dplx_init_S + SH[0] + RC);
		} else {
			S1 = -1.0;
			H1 = Double.POSITIVE_INFINITY;
			T1 = (H1 + dplx_init_H) /(S1 + dplx_init_S + RC);
		}
		   
		if(S1 < MinEntropyCutoff) {
		      /* to not give dH any value if dS is unreasonable */
		      S1 = MinEntropy;
		      H1 = 0.0;
		}
		if(S0 < MinEntropyCutoff) {
			/* to not give dH any value if dS is unreasonable */
			S0 = MinEntropy;
			H0 = 0.0;
		}
		if(T1 > T0) { 
			setEntropyDPT(i,j, S1);
			setEnthalpyDPT(i,j, H1);
		} else if(T0 >= T1) {
			setEntropyDPT(i,j, S0);
			setEnthalpyDPT(i,j, H0);
		}
	}

	/**
	 *  returns stack entropy 
	 */
	private double Ss(int i, int j, int k) {
		if(k==2) {
			if (i >= j)
				return -1.0;
			if (i == len1 || j == len2 + 1)
				return -1.0;

			if (i > len1)
				i -= len1;
			if (j > len2)
				j -= len2;
			return thallib.stackEntropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]];
		} else {
			return thallib.stackEntropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
		}
	}

	/**
	 * returns stack enthalpy
	 * @param i
	 * @param j
	 * @param k
	 * @return
	 */
	private double Hs(int i, int j, int k) {
		if(k==2) {
			if (i >= j)
				return Double.POSITIVE_INFINITY;
			if (i == len1 || j == len2 + 1)
				return Double.POSITIVE_INFINITY;

			if (i > len1)
				i -= len1;
			if (j > len2)
				j -= len2;
			if(Double.isFinite(thallib.stackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]])) {
				return thallib.stackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]];
			} else {
				return Double.POSITIVE_INFINITY;
			}
		} else {
			return thallib.stackEnthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
		}
	}

	/* calculate terminal entropy S and terminal enthalpy H starting reading from 5'end (Left hand/3' end - Right end) */
	private void LSH(int i, int j, double[] EntropyEnthalpy) {
		double S1, H1, T1, G1;
		double S2, H2, T2, G2;
		S1 = S2 = -1.0;
		H1 = H2 = Double.NEGATIVE_INFINITY;
		T1 = T2 = Double.NEGATIVE_INFINITY;
		if (BPI[numSeq1[i]][numSeq2[j]] == 0) {
			setEntropyDPT(i,j,-1.0);
			setEnthalpyDPT(i,j, Double.POSITIVE_INFINITY);
			return;
		}
		S1 = thallib.atpS[numSeq1[i]][numSeq2[j]] + thallib.tstack2Entropies[numSeq2[j]][numSeq2[j-1]][numSeq1[i]][numSeq1[i-1]];
		H1 = thallib.atpH[numSeq1[i]][numSeq2[j]] + thallib.tstack2Enthalpies[numSeq2[j]][numSeq2[j-1]][numSeq1[i]][numSeq1[i-1]];
		G1 = H1 - TEMP_KELVIN*S1;
		if(!Double.isFinite(H1) || G1>0) {
			H1 = Double.POSITIVE_INFINITY;
			S1 = -1.0;
			G1 = 1.0;
		}
		/** If there is two dangling ends at the same end of duplex **/
		if((BPI[numSeq1[i-1]][numSeq2[j-1]] != 1 ) && Double.isFinite(thallib.dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]]) && Double.isFinite(thallib.dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]])) {
			S2 = thallib.atpS[numSeq1[i]][numSeq2[j]] + thallib.dangleEntropies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]] +
					thallib.dangleEntropies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
			H2 = thallib.atpH[numSeq1[i]][numSeq2[j]] + thallib.dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]] +
					thallib.dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
			G2 = H2 - TEMP_KELVIN*S2;
			if(!Double.isFinite(H2) || G2>0) {
				H2 = Double.POSITIVE_INFINITY;
				S2 = -1.0;
				G2 = 1.0;
			}
			T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
			if(Double.isFinite(H1) && G1<0) {
				T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
				if(T1<T2 && G2<0) {
					S1 = S2;
					H1 = H2;
					T1 = T2;
				}
			} else if(G2<0) {
				S1 = S2;
				H1 = H2;
				T1 = T2;
			}
		} else if ((BPI[numSeq1[i-1]][numSeq2[j-1]] != 1) && Double.isFinite(thallib.dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]])) {
			S2 = thallib.atpS[numSeq1[i]][numSeq2[j]] + thallib.dangleEntropies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]];
			H2 = thallib.atpH[numSeq1[i]][numSeq2[j]] + thallib.dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]];
			G2 = H2 - TEMP_KELVIN*S2;
			if(!Double.isFinite(H2) || G2>0) {
				H2 = Double.POSITIVE_INFINITY;
				S2 = -1.0;
				G2 = 1.0;
			}
			T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
			if(Double.isFinite(H1) && G1<0) {
				T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
				if(T1<T2 && G2<0) {
					S1 = S2;
					H1 = H2;
					T1 = T2;
				}
			} else if (G2<0){
				S1 = S2;
				H1 = H2;
				T1 = T2;
			}
		} else if ((BPI[numSeq1[i-1]][numSeq2[j-1]] != 1) && Double.isFinite(thallib.dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]])) {
			S2 = thallib.atpS[numSeq1[i]][numSeq2[j]] + thallib.dangleEntropies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
			H2 = thallib.atpH[numSeq1[i]][numSeq2[j]] + thallib.dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
			G2 = H2 - TEMP_KELVIN*S2;
			if(!Double.isFinite(H2) || G2>0) {
				H2 = Double.POSITIVE_INFINITY;
				S2 = -1.0;
				G2 = 1.0;     
			}
			T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
			if(Double.isFinite(H1) && G1<0) {
				T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
				if(T1 < T2  && G2<0) {
					S1 = S2;
					H1 = H2;
					T1 = T2;
				}
			} else if(G2<0) {
				S1 = S2;
				H1 = H2;
				T1 = T2;
			}
		}
		S2 = thallib.atpS[numSeq1[i]][numSeq2[j]];
		H2 = thallib.atpH[numSeq1[i]][numSeq2[j]];
		T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
		G1 = H1 -TEMP_KELVIN*S1;   
		G2 = H2 -TEMP_KELVIN*S2;
		if(Double.isFinite(H1)) {
			if(T1 < T2) {
				EntropyEnthalpy[0] = S2;
				EntropyEnthalpy[1] = H2;
			} else {
				EntropyEnthalpy[0] = S1;
				EntropyEnthalpy[1] = H1;
			}
		} else {
			EntropyEnthalpy[0] = S2;
			EntropyEnthalpy[1] = H2;
		}
		
	}

	/**
	 *  calculate terminal entropy S and terminal enthalpy H starting reading from 5'end
	 *  (Left hand/3' end - Right end) 
	 */
	private void RSH(int i, int j, double[] EntropyEnthalpy) {
		double G1, G2;
		double S1, S2;
		double H1, H2;
		double T1, T2;
		S1 = S2 = -1.0;
		H1 = H2 = Double.POSITIVE_INFINITY;
		T1 = T2 = Double.NEGATIVE_INFINITY;
		if (BPI[numSeq1[i]][ numSeq2[j]] == 0) {
			EntropyEnthalpy[0] = -1.0;
			EntropyEnthalpy[1] = Double.POSITIVE_INFINITY;
			return;
		}
		S1 = thallib.atpS[numSeq1[i]][ numSeq2[j]] + thallib.tstack2Entropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
		H1 = thallib.atpH[numSeq1[i]][ numSeq2[j]] + thallib.tstack2Enthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
		G1 = H1 - TEMP_KELVIN*S1;
		if(!Double.isFinite(H1) || G1>0) {
			H1 = Double.POSITIVE_INFINITY;
			S1 = -1.0;
			G1 = 1.0;
		}
		   
		if(BPI[numSeq1[i+1]][numSeq2[j+1]] == 0 && Double.isFinite(thallib.dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]]) && Double.isFinite(thallib.dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]])) {
			S2 = thallib.atpS[numSeq1[i]][numSeq2[j]] + thallib.dangleEntropies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]] +
					thallib.dangleEntropies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
			H2 = thallib.atpH[numSeq1[i]][numSeq2[j]] + thallib.dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]] +
					thallib.dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
			G2 = H2 - TEMP_KELVIN*S2;
			if(!Double.isFinite(H2) || G2>0) {
				H2 = Double.POSITIVE_INFINITY;
				S2 = -1.0;
				G2 = 1.0;
			}
		      
			T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
			if(Double.isFinite(H1) && G1<0) {
				T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
				if(T1 < T2 && G2<0) {
					S1 = S2;
					H1 = H2;
					T1 = T2;
				}
			} else if(G2<0){
				S1 = S2;
				H1 = H2;
				T1 = T2;
			}
		}

		else if(BPI[numSeq1[i+1]][numSeq2[j+1]] == 0 && Double.isFinite(thallib.dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]])) {
			S2 = thallib.atpS[numSeq1[i]][numSeq2[j]] + thallib.dangleEntropies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]];
			H2 = thallib.atpH[numSeq1[i]][numSeq2[j]] + thallib.dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]];
			G2 = H2 - TEMP_KELVIN*S2;
			if(!Double.isFinite(H2) || G2 >0) {
				H2 = Double.POSITIVE_INFINITY;
				S2 = -1.0;
				G2 = 1.0;
			}
			T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
			if(Double.isFinite(H1) && G1<0) {
				T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
				if(T1 < T2 && G2<0) {
					S1 = S2;
					H1 = H2;
					T1 = T2;
				}
			} else if(G2<0){
				S1 = S2;
				H1 = H2;
				T1 = T2;
			}
		} else if(BPI[numSeq1[i+1]][numSeq2[j+1]] == 0 && Double.isFinite(thallib.dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]])) {
			S2 = thallib.atpS[numSeq1[i]][numSeq2[j]] + thallib.dangleEntropies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
			H2 = thallib.atpH[numSeq1[i]][numSeq2[j]] + thallib.dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
			G2 = H2 - TEMP_KELVIN*S2;
			if(!Double.isFinite(H2) || G2>0) {
				H2 = Double.POSITIVE_INFINITY;
				S2 = -1.0;
				G2 = 1.0;
			}
			T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
			if(Double.isFinite(H1) && G1<0) {
				T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
				if(T1 < T2 && G2<0) {
					S1 = S2;
					H1 = H2;
					T1 = T2;
				}
			} else if (G2<0){
				S1 = S2;
				H1 = H2;
				T1 = T2;
			}
		}
		S2 = thallib.atpS[numSeq1[i]][numSeq2[j]];
		H2 = thallib.atpH[numSeq1[i]][numSeq2[j]];
		T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
		G1 = H1 -TEMP_KELVIN*S1;
		G2 =  H2 -TEMP_KELVIN*S2;
		if(Double.isFinite(H1)) {
			if(T1 < T2) {
				EntropyEnthalpy[0] = S2;
				EntropyEnthalpy[1] = H2;
			} else {
				EntropyEnthalpy[0] = S1;
				EntropyEnthalpy[1] = H1;
			}
		} else {
			EntropyEnthalpy[0] = S2;
			EntropyEnthalpy[1] = H2;
		}
	}
	
	
	 /**
	  *  initiates thermodynamic parameter tables of entropy and enthalpy for dimer 
	  */
	private void initMatrix() {
		   
		for (int i = 1; i <= len1; ++i) {
			for (int j = 1; j <= len2; ++j) {
				if (BPI[numSeq1[i]][numSeq2[j]] == 0)  {
					setEnthalpyDPT(i,j, Double.POSITIVE_INFINITY);
					setEntropyDPT(i,j, -1.0);
				} else {
					setEnthalpyDPT(i,j,0.0);
					setEntropyDPT(i,j,MinEntropy);
				}
			}
		}
		     		
	}
	/**
	 * initiates thermodynamic parameter tables of entropy and enthalpy for monomer 
	 */
	private void initMatrix2() {
		   
		for (int i = 1; i <= len1; ++i)
		{
			for (int j = i; j <= len2; ++j)
			{
				if (j - i < MIN_HRPN_LOOP + 1 || (BPI[numSeq1[i]][numSeq1[j]] == 0)) {
					setEnthalpyDPT(i,j,Double.POSITIVE_INFINITY);
					setEntropyDPT(i,j,-1.0);
				} else {
					setEnthalpyDPT(i,j, 0.0);
					setEntropyDPT(i,j,MinEntropy);
				}
			}
		}
		     		
	}
	
	/**
	 * calc-s thermod values into dynamic progr table (dimer) 
	 */
	private void fillMatrix() {
		int maxLoop	= a.maxLoop;
		double[] SH = new double[2];
		for (int i = 1; i <= len1; ++i) {
			for (int j = 1; j <= len2; ++j) {
				if(Double.isFinite(EnthalpyDPT(i,j))) { /* if finite */
					SH[0] = -1.0;
					SH[1] = Double.POSITIVE_INFINITY;
					LSH(i,j,SH);
					if(Double.isFinite(SH[1])) {
						setEntropyDPT(i,j, SH[0]);
						setEnthalpyDPT(i,j, SH[1]);
					}
					if (i > 1 && j > 1) {
						maxTM(i, j); /* stack: sets EntropyDPT(i, j) and EnthalpyDPT(i, j) */
						for(int d = 3; d <= maxLoop + 2; d++) {  /* max=30, length over 30 is not allowed */
							int ii = i - 1;
							int jj = - ii - d + (j + i);
							if (jj < 1) {
							     ii -= Math.abs(jj-1);
							     jj = 1;
							}
							for (; ii > 0 && jj < j; --ii, ++jj) {
								if(Double.isFinite(EnthalpyDPT(ii,jj))){
									SH[0] = -1.0;
									SH[1] = Double.POSITIVE_INFINITY;
									calc_bulge_internal(ii, jj, i, j, SH,false,maxLoop);
									if(SH[0] < MinEntropyCutoff) {
										   /* to not give dH any value if dS is unreasonable */
										   SH[0] = MinEntropy;
										   SH[1] = 0.0;
										}
									if(Double.isFinite(SH[1])) {
										   setEnthalpyDPT(i,j,SH[1]);
										   setEntropyDPT(i,j,SH[0]);
										}
								}
							}
						}
					}
				}
			}
		}
	}

	/**
	 * calc-s thermod values into dynamic progr table (monomer)
	 */
	private void fillMatrix2() {
		int maxLoop	= a.maxLoop;
		double[] SH = new double[2];
		for (int j = 2; j <= len2; ++j) {
			for (int i = j - MIN_HRPN_LOOP - 1; i >= 1; --i) {
				if(Double.isFinite(EnthalpyDPT(i,j))){
					SH[0] = -1.0;
					SH[1] = Double.POSITIVE_INFINITY;
					maxTM2(i,j); /* calculate stack */
					CBI(i, j, SH, false,maxLoop); /* calculate Bulge and Internal loop and stack */
					SH[0] = -1.0;
					SH[1] = Double.POSITIVE_INFINITY;
					calc_hairpin(i, j, SH, 0);
					if(Double.isFinite(SH[1])) {
						if(SH[0] < MinEntropyCutoff){ /* to not give dH any value if dS is unreasonable */
							SH[0] = MinEntropy;
							SH[1] = 0.0;
						}
						setEntropyDPT(i,j,SH[0]);
						setEnthalpyDPT(i,j, SH[1]);
					}
				}
			}
		}

	}
	private double saltCorrectS(double mv, double dv, double dntp) {
		 if(dv<=0) 
			 dntp=dv;
		 return 0.368*((Math.log((mv+120*(Math.sqrt(Math.max(0.0, dv-dntp))))/1000.0)));
	}


//	private boolean symmetry_thermo(char[] oligo1) {
//		// TODO Auto-generated method stub
//		return false;
//	}
	
	
	
	class tracer /* structure for tracebacku - unimolecular str */ {
		 
		public tracer(int i , int j, int mtrx)
		{
			this.i = i;
			this.j = j ;
			this.mtrx = mtrx;
		}
		
		int i;
		int j;
		int mtrx; /* [0 1] EntropyDPT/EnthalpyDPT*/
	};
	
	public static int getHashkey(int[]  seq,int start , int len)
	{
		int[] seqLoop = new int[len];
		for(int i =  start,j=0;j<len;i++,j++)
		{
			seqLoop[j] = seq[i];
		}
		return thallib.hashLoop(seqLoop);
	}
}

