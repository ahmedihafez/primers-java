package org.primer3.oligotm;

import java.util.HashMap;

import org.primer3.sequence.Sequence;
import org.primer3.sequence.Sequence.SequenceIterator;

public class OligoTMCalculator {

	
// 	#define OLIGOTM_ERROR -999999.9999
	static public double OLIGOTM_ERROR = -999999.9999;
	
	static final double T_KELVIN = 273.15;
	static char A_CHAR = 'A';
	static char G_CHAR = 'G';
	static char T_CHAR = 'T';
	static char C_CHAR = 'C';
	static char N_CHAR = 'N';
	
	static final double crossover_point = 0.22;

	
	
	
	
	/** 
	 * Return the delta G of the last len bases of oligo if oligo is at least len
	 * bases long; otherwise return the delta G of oligo. 
	 */
	static double end_oligodg(Sequence s, int len, MeltingTemperatureMethod tm_method) {
		
		if (tm_method != MeltingTemperatureMethod.breslauer_auto
			      && tm_method != MeltingTemperatureMethod.santalucia_auto)
			    return OLIGOTM_ERROR;
		int x = s.length();
		
		return x < len ? oligoDeltaG(s,tm_method) : oligoDeltaG(s.subSeqRange(x-len,x-1),tm_method);
	}

	public static double end_oligodg(char[] seq, int len, MeltingTemperatureMethod tm_method) {
		return end_oligodg(new Sequence(seq),len,tm_method);
	}
	/** Calculate the melting temperature of substr(seq, start, length) using the
	   formula from Bolton and McCarthy, PNAS 84:1390 (1962) as presented in
	   Sambrook, Fritsch and Maniatis, Molecular Cloning, p 11.46 (1989, CSHL
	   Press).

	   Tm = 81.5 + 16.6(log10([Na+])) + .41*(%GC) - 600/length

	   Where [Na+] is the molar sodium concentration, (%GC) is the percent of Gs
	   and Cs in the sequence, and length is the length of the sequence.

	   A similar formula is used by the prime primer selection program in GCG
	   (http://www.gcg.com), which instead uses 675.0 / length in the last term
	   (after F. Baldino, Jr, M.-F. Chesselet, and M.E.  Lewis, Methods in
	   Enzymology 168:766 (1989) eqn (1) on page 766 without the mismatch and
	   formamide terms).  The formulas here and in Baldino et al. assume Na+ rather
	   than K+.  According to J.G. Wetmur, Critical Reviews in BioChem. and
	   Mol. Bio. 26:227 (1991) 50 mM K+ should be equivalent in these formulae to .2
	   M Na+.

	   This function takes salt_conc to be the millimolar (mM) concentration,
	   since mM is the usual units in PCR applications.

	 */
	static double longSeqTM(Sequence seq, 
	                   int start, 
	                   int length, 
	                   double salt_conc, 
	                   double divalent_conc, 
	                   double dntp_conc) {
		
		
		if (divalentToMonovalent(divalent_conc, dntp_conc) == OLIGOTM_ERROR)
			return OLIGOTM_ERROR;
		int GC_count = 0;
		salt_conc = salt_conc + divalentToMonovalent(divalent_conc, dntp_conc);

		if ((start + length) > seq.length() || start < 0 || length <= 0)
			    return OLIGOTM_ERROR;
		int end = start + length;
		for(int i = start ; i < end ; i ++)
		{
			char c = seq.get(i);
			if ('G' == c || 'C' == c)
				GC_count++;
		}
		return  81.5 + (16.6 * Math.log10(salt_conc / 1000.0)) + (41.0 * (((double) GC_count) / length)) - (600.0 / length);
	}
	
	
	
	public static double longSeqTM(char[] seq, 
            int start, 
            int length, 
            double salt_conc, 
            double divalent_conc, 
            double dntp_conc) {
	
		if (divalentToMonovalent(divalent_conc, dntp_conc) == OLIGOTM_ERROR)
			return OLIGOTM_ERROR;
		int GC_count = 0;
		salt_conc = salt_conc + divalentToMonovalent(divalent_conc, dntp_conc);

		if ((start + length) > seq.length || start < 0 || length <= 0)
			    return OLIGOTM_ERROR;
		int end = start + length;
		for(int i = start ; i < end ; i ++)
		{
			char c = seq[i];
			if ('G' == c || 'C' == c)
				GC_count++;
		}
		return  81.5
			    + (16.6 * Math.log10(salt_conc / 1000.0))
			    + (41.0 * (((double) GC_count) / length))
			    - (600.0 / length);
	}
	
	/* 
	   For olgigotm() and seqtm()

	   Both functions return the melting temperature of the given oligo
	   calculated as specified by user, but oligotm _should_ only be used on
	   DNA sequences of length <= MAX_PRIMER_LENGTH (which is defined
	   elsewhere).  seqtm uses oligotm for sequences of length <=
	   MAX_PRIMER_LENGTH, and a different, G+C% based formula for longer
	   sequences.  For oligotm(), no error is generated on sequences
	   longer than MAX_PRIMER_LENGTH, but the formula becomes less
	   accurate as the sequence grows longer.  Caveat emptor.

	   We use the folowing typedefs:
	*/
	
	/**
	 * 
	 * @param seq The sequence.
	 * @param DNA_nM DNA concentration (nanomolar).
	 * @param K_mM Salt concentration (millimolar).
	 * @param divalent_conc Concentration of divalent cations (millimolar)
	 * @param dntp_conc Concentration of dNTPs (millimolar)
	 * @param tm_method
	 * @param salt_corrections
	 * @return
	 */
	static double oligoTM(Sequence seq,    
            double DNA_nM,     
            double K_mM,     
            double divalent_conc,
            double dntp_conc,     
            MeltingTemperatureMethod tm_method,    
            SaltCorrectionMethod salt_corrections  
            ) {
		
		int dh = 0, ds = 0 ;
		char c;
		double delta_H,delta_S;
		double Tm; /* Melting temperature */
		double correction;
		int len ;
		boolean sym;
		
		
		 if(divalentToMonovalent(divalent_conc, dntp_conc) == OLIGOTM_ERROR)
			     return OLIGOTM_ERROR;
		  /** K_mM = K_mM + divalent_to_monovalent(divalent_conc, dntp_conc); **/
		 if (tm_method != MeltingTemperatureMethod.breslauer_auto
			      && tm_method != MeltingTemperatureMethod.santalucia_auto)
			    return OLIGOTM_ERROR;
		 if (salt_corrections != SaltCorrectionMethod.schildkraut
			      && salt_corrections != SaltCorrectionMethod.santalucia
			      && salt_corrections != SaltCorrectionMethod.owczarzy)
			    return OLIGOTM_ERROR;
		// len = (strlen(s)-1); no need for this 
//		len = seq.length()-1;
		sym = seq.symmetry();  /*Add symmetry correction if seq is symmetrical*/
		
		SequenceIterator s = (SequenceIterator) seq.iterator();
		c = s.current(); s.next(); // or just c = s.next();
		HashMap<String, Integer> dh_Table = DH;
		HashMap<String, Integer> ds_Table = DS;

		
		if( tm_method == MeltingTemperatureMethod.breslauer_auto ) {
		    ds=108;
		    dh_Table = H;
		    ds_Table = S;
		}
		else {
		  if(sym) {
			  ds+=14;
		  }
		  /** Terminal AT penalty **/
		  if(c == 'A' || c =='T')
		  {
		      ds += -41;
		      dh += -23;
		  } else if (c == 'C' || c =='G') {
		      ds += 28;
		      dh += -1;
		  }
		  
		}
		
		char last = c;
		// TODO :: add try catch here to return OLIGOTM_ERROR in case if error
		while( s.hasNext() )
		{
			last = c;
			c = s.current(); s.next();
			ds += ds_Table.get(STATE(last, c));
			dh += dh_Table.get(STATE(last, c));
		}
		// get last one again
		if( tm_method != MeltingTemperatureMethod.breslauer_auto )
		{
			c = s.previous();
			/** Terminal AT penalty **/
			if(c == 'A' || c =='T')
			{
				ds += -41;
				dh += -23;
			} else if (c == 'C' || c =='G') {
				ds += 28;
			    dh += -1;
			}
		}
		
		
		
		/* dh and ds are now computed for the given sequence. */
		
		/* 
		 * Nearest-neighbor thermodynamic values for dh
		 * are given in 100 cal/mol of interaction.
		 */
		delta_H = dh * -100.0;
		
		 /*
		  * Nearest-neighbor thermodynamic values for ds
		  * are in in .1 cal/K per mol of interaction.
		  */
		delta_S = ds * -0.1;    
		Tm=0;  /* Melting temperature */
		len = seq.length();
		
		
		if (salt_corrections == SaltCorrectionMethod.schildkraut) {
		     K_mM = K_mM + divalentToMonovalent(divalent_conc, dntp_conc);
		     correction = 16.6 * Math.log10(K_mM/1000.0) - T_KELVIN;
		     Tm = delta_H / (delta_S + 1.987 * Math.log(DNA_nM/4000000000.0)) + correction;
		} else if (salt_corrections== SaltCorrectionMethod.santalucia) {
		    K_mM = K_mM + divalentToMonovalent(divalent_conc, dntp_conc);
		    delta_S = delta_S + 0.368 * (len - 1) * Math.log(K_mM / 1000.0 );
		    if(sym ) { /* primer is symmetrical */
		      /* Equation A */
		      Tm = delta_H / (delta_S + 1.987 * Math.log(DNA_nM/1000000000.0)) - T_KELVIN;
		    } else {
		      /* Equation B */
		      Tm = delta_H / (delta_S + 1.987 * Math.log(DNA_nM/4000000000.0)) - T_KELVIN;
		    }      
		 } else if (salt_corrections == SaltCorrectionMethod.owczarzy) {
		     double gcPercent= seq.getGCPrercent();
		     double free_divalent; /* conc of divalent cations minus dNTP conc */
		     /**** BEGIN: UPDATED SALT BY OWCZARZY *****/
		      /* different salt corrections for monovalent (Owczarzy et al.,2004) and divalent cations (Owczarzy et al.,2008) */
		     /* competition bw magnesium and monovalent cations, see Owczarzy et al., 2008 Figure 9 and Equation 16 */
		     double div_monov_ratio;
		     if(dntp_conc >= divalent_conc) {
		    	free_divalent = 0.00000000001; /* to not to get log(0) */
		     	} 
		     else {
		    	free_divalent = (divalent_conc - dntp_conc)/1000.0;
		     }
		     double A = 0,B = 0,C = 0,D = 0,E = 0,F = 0,G = 0;
		      if(K_mM==0) {
		    	  div_monov_ratio = 6.0;
		      } else {
		    	  div_monov_ratio = (Math.sqrt(free_divalent))/(K_mM/1000); /* if conc of monov cations is provided
		    									    a ratio is calculated to further calculate
		    									    the _correct_ correction */
		      }
		      
		      if (div_monov_ratio < crossover_point) {
		    	 /* use only monovalent salt correction, Eq 22 (Owczarzy et al., 2004) */
		    	 correction = (((4.29 * gcPercent) - 3.95) * Math.pow(10,-5) * Math.log(K_mM / 1000.0))
		    		    + (9.40 * Math.pow(10,-6) * (Math.pow(Math.log(K_mM / 1000.0),2)));
		      } 
		      else 
		      {
		    	  /* magnesium effects are dominant, Eq 16 (Owczarzy et al., 2008) is used */
		    	  B =- 9.11 * Math.pow(10,-6);
		    	  C = 6.26 * Math.pow(10,-5);
		    	  E =- 4.82 * Math.pow(10,-4);
		    	  F = 5.25 * Math.pow(10,-4);
		    	  A = 3.92 * Math.pow(10,-5);
		    	  D = 1.42 * Math.pow(10,-5);
		    	  G = 8.31 * Math.pow(10,-5);
		    	  if(div_monov_ratio < 6.0) {
		    		   /* in particular ratio of conc of monov and div cations
		    		    *             some parameters of Eq 16 must be corrected (a,d,g) */
		    		   A = 3.92 * Math.pow(10,-5) * (0.843 - (0.352 * Math.sqrt(K_mM/1000.0) * Math.log(K_mM/1000.0)));
		    		   D = 1.42 * Math.pow(10,-5) * (1.279 - 4.03 * Math.pow(10,-3) * Math.log(K_mM/1000.0) - 8.03 * Math.pow(10,-3) * Math.pow(Math.log(K_mM/1000.0),2));
		    		   G = 8.31 * Math.pow(10,-5) * (0.486 - 0.258 * Math.log(K_mM/1000.0) + 5.25 * Math.pow(10,-3) * Math.pow(Math.log(K_mM/1000.0),3));
		    	  }
		    		
		    	  correction = A + (B * Math.log(free_divalent))
		    		  + gcPercent * (C + (D * Math.log(free_divalent)))
		    		    + (1/(2 * (len - 1))) * (E + (F * Math.log(free_divalent))
		    					     + G * (Math.pow((Math.log(free_divalent)),2)));
		    }
		    /**** END: UPDATED SALT BY OWCZARZY *****/
		    if (sym) {
		   	/* primer is symmetrical */
		   	/* Equation A */
		    	Tm = 1/((1/(delta_H 
		   		    /
		   		    (delta_S + 1.9872 * Math.log(DNA_nM/1000000000.0)))) + correction) - T_KELVIN;
		    } 
		    else 
		    {
		    	/* Equation B */
		    	Tm = 1/((1/(delta_H
		   		    /
		   		    (delta_S + 1.9872 * Math.log(DNA_nM/4000000000.0)))) + correction) - T_KELVIN;
		    } 
		 }
		
		 return Tm;
	}
    



	/**
	 *  Return the melting temperature of a given sequence, 'seq', of any length.
	 * @param seq : The sequence.
	 * @param dna_conc : DNA concentration (nanomolar).
	 * @param salt_conc : Concentration of divalent cations (millimolar).
	 * @param divalent_conc : Concentration of divalent cations (millimolar)
	 * @param dntp_conc :  Concentration of dNTPs (millimolar)
	 * @param nn_max_len : The maximum sequence length for using the nearest neighbor model
	 * (as implemented in oligotm.  For sequences longer than this, seqtm uses the "GC%" formula implemented in long_seq_tm.
	 * @param tm_method :
	 * @param salt_corrections :
	 * @return
	 */
	static double sequenceTM(Sequence seq,  
          double dna_conc,  
          double salt_conc,  
          double divalent_conc, 
          double dntp_conc,    
          int  nn_max_len,  
          MeltingTemperatureMethod  tm_method,       
          SaltCorrectionMethod salt_corrections 
          ) {
		
		int len = seq.length();
		
		// No need for the check here
//		if (tm_method != MeltingTemperatureMethod.breslauer_auto
//			      && tm_method != MeltingTemperatureMethod.santalucia_auto)
//			return OLIGOTM_ERROR;
//		if (salt_corrections != SaltCorrectionMethod.schildkraut
//			      && salt_corrections != SaltCorrectionMethod.santalucia
//			      && salt_corrections != SaltCorrectionMethod.owczarzy)
//			return OLIGOTM_ERROR;
		
		
		 if (len > nn_max_len) {
			  return longSeqTM(seq, 0, len, salt_conc, divalent_conc, dntp_conc);
		  } else {
			  return oligoTM(seq, dna_conc, salt_conc, 
				      divalent_conc, dntp_conc, tm_method, salt_corrections);
		  }
	}
	
	public static double sequenceTM(char[] seq,  
	          double dna_conc,  
	          double salt_conc,  
	          double divalent_conc, 
	          double dntp_conc,    
	          int  nn_max_len,  
	          MeltingTemperatureMethod  tm_method,       
	          SaltCorrectionMethod salt_corrections 
	          ) {
		return sequenceTM(new Sequence(seq),  
		           dna_conc,  
		           salt_conc,  
		           divalent_conc, 
		           dntp_conc,    
		            nn_max_len,  
		            tm_method,       
		           salt_corrections 
		          ) ;
	}
	
	/** Return the delta G of disruption of oligo using the nearest neighbor model.
	   The length of seq should be relatively short, 
	   given the characteristics of the nearest
	   neighbor model.
	*/
	static double oligoDeltaG(Sequence seq,  MeltingTemperatureMethod tm_method  ) {
		
		// unwanted Check
//		if (tm_method != MeltingTemperatureMethod.breslauer_auto
//				&& tm_method != MeltingTemperatureMethod.santalucia_auto)
//			return OLIGOTM_ERROR;
				
		double dg = 0;
		SequenceIterator s = (SequenceIterator) seq.iterator();
		
		/* Use a finite-state machine (DFA) to calculate dg s. */
		//   c = *s; s++;
		char c = s.current(); s.next(); // or just c = s.next();
		HashMap<String, Integer> dg_Table = G;
		if(tm_method != MeltingTemperatureMethod.breslauer_auto)
		{	
			dg_Table = DG;
			dg=-1960; /* Initial dG */
			if(c == 'A' || c == 'T')  {
				 dg += -50; /* terminal AT penalty */
			}
		}
		
		char last = c;
		// TODO :: add try catch here to return OLIGOTM_ERROR in case if error
		while( s.hasNext() )
		{
			last = c;
			c = s.current(); s.next();
			dg += dg_Table.get(STATE(last, c));
			
		}
		
		 if(tm_method != MeltingTemperatureMethod.breslauer_auto) {
//		    int sym;
		    // --s; --s; c = *s;
		    c = s.previous();
		    if(c == 'A' || c == 'T')  {
		    	dg += -50; /* terminal AT penalty */
		    }
//		    sym = symmetry(s);
		    if(seq.symmetry())   {
		    	dg +=-430; /* symmetry correction for dG */
		     }
		   }
		
		
		 return dg / 1000.0;
	}

	/** Returns 1 if the sequence is self-complementary or symmetrical; 0
	   otherwise
	*/
//	int symmetry(char[] seq) {
//		return 0 ;
//	}
	

	/** Converts divalent salt concentration to monovalent salt concentration */
	static double divalentToMonovalent(double divalent, double dntp) {
		if(divalent==0) dntp=0;
		if(divalent<0 || dntp<0) return OLIGOTM_ERROR;
		/* According to theory, melting temperature does not depend on
		divalent cations */
		if(divalent<dntp) 
			divalent=dntp;  
		return 120*(Math.sqrt(divalent-dntp));
	}
	
	
	
	
	
	
	/*
	 * Two tables of nearest-neighbor parameters for di-nucleotide
	 * base pairs.
	 *
	 */
	
	/** 
	 * Table 1 (old parameters):
	 * See table 2 in the paper [Breslauer KJ, Frank R, Blöcker H and
	 * Marky LA (1986) "Predicting DNA duplex stability from the base
	 * sequence" Proc Natl Acad Sci 83:4746-50
	 * http://dx.doi.org/10.1073/pnas.83.11.3746]
	 */
	static HashMap<String, Integer> S= new HashMap<String, Integer>(); // Not sure what is the best way to represent this here
	static HashMap<String,Integer> H= new HashMap<String, Integer>();
	/**
	 *  Delta G's of disruption * 1000. 
	 */
	static HashMap<String,Integer> G= new HashMap<String, Integer>();

	/** Table 2, new parameters:
	 * Tables of nearest-neighbor thermodynamics for DNA bases, from the
	 * paper [SantaLucia JR (1998) "A unified view of polymer, dumbbell
	 * and oligonucleotide DNA nearest-neighbor thermodynamics", Proc Natl
	 * Acad Sci 95:1460-65 http://dx.doi.org/10.1073/pnas.95.4.1460]
	 */
	static HashMap<String,Integer> DS = new HashMap<String, Integer>();
	static HashMap<String,Integer> DH = new HashMap<String, Integer>();
	/**
	 *  Delta G's of disruption * 1000. 
	 */
	static HashMap<String,Integer> DG = new HashMap<String, Integer>();

	static public String STATE(char from , char to)
	{
		return from + "_" + to;
	}
	
	
	static {
		S.put("A_A", 240);
		S.put("A_C", 173);
		S.put("A_G", 208);
		S.put("A_T", 239);
		S.put("A_N", 215);

		S.put("C_A", 129);
		S.put("C_C", 266);
		S.put("C_G", 278);
		S.put("C_T", 208);
		S.put("C_N", 220);

		S.put("G_A", 135);
		S.put("G_C", 267);
		S.put("G_G", 266);
		S.put("G_T", 173);
		S.put("G_N", 210);

		S.put("T_A", 169);
		S.put("T_C", 135);
		S.put("T_G", 129);
		S.put("T_T", 240);
		S.put("T_N", 168);

		S.put("N_A", 168);
		S.put("N_C", 210);
		S.put("N_G", 220);
		S.put("N_T", 215);
		S.put("N_N", 203);
		
		
		H.put("A_A",91);
		H.put("A_C",65);
		H.put("A_G",78);
		H.put("A_T",86);
		H.put("A_N",80);

		H.put("C_A",58);
		H.put("C_C",110);
		H.put("C_G",119);
		H.put("C_T",78);
		H.put("C_N",91);

		H.put("G_A",56);
		H.put("G_C",111);
		H.put("G_G",110);
		H.put("G_T",65);
		H.put("G_N",85);

		H.put("T_A",60);
		H.put("T_C",56);
		H.put("T_G",58);
		H.put("T_T",91);
		H.put("T_N",66);

		H.put("N_A",66);
		H.put("N_C",85);
		H.put("N_G",91);
		H.put("N_T",80);
		H.put("N_N",80);
		
		
		/* Delta G's of disruption * 1000. */
		G.put("A_A",1900);
		G.put("A_C",1300);
		G.put("A_G",1600);
		G.put("A_T",1500);
		G.put("A_N",1575);

		G.put("C_A",1900);
		G.put("C_C",3100);
		G.put("C_G",3600);
		G.put("C_T",1600);
		G.put("C_N",2550);

		G.put("G_A",1600);
		G.put("G_C",3100);
		G.put("G_G",3100);
		G.put("G_T",1300);
		G.put("G_N",2275);

		G.put("T_A", 900);
		G.put("T_C",1600);
		G.put("T_G",1900);
		G.put("T_T",1900);
		G.put("T_N",1575);

		G.put("N_A",1575);
		G.put("N_C",2275);
		G.put("N_G",2550);
		G.put("N_T",1575);
		G.put("N_N",1994);
		
		DS.put("A_A",222);
		DS.put("A_C",224);
		DS.put("A_G",210);
		DS.put("A_T",204);
		DS.put("A_N",224);

		DS.put("C_A",227);
		DS.put("C_C",199);
		DS.put("C_G",272);
		DS.put("C_T",210);
		DS.put("C_N",272);

		DS.put("G_A",222);
		DS.put("G_C",244);
		DS.put("G_G",199);
		DS.put("G_T",224);
		DS.put("G_N",244);

		DS.put("T_A",213);
		DS.put("T_C",222);
		DS.put("T_G",227);
		DS.put("T_T",222);
		DS.put("T_N",227);

		DS.put("N_A",168);
		DS.put("N_C",210);
		DS.put("N_G",220);
		DS.put("N_T",215);
		DS.put("N_N",220);


		DH.put("A_A", 79);
		DH.put("A_C", 84);
		DH.put("A_G", 78);
		DH.put("A_T", 72);
		DH.put("A_N", 72);

		DH.put("C_A", 85);
		DH.put("C_C", 80);
		DH.put("C_G",106);
		DH.put("C_T", 78);
		DH.put("C_N", 78);

		DH.put("G_A",82);
		DH.put("G_C",98);
		DH.put("G_G",80);
		DH.put("G_T",84);
		DH.put("G_N",80);

		DH.put("T_A",72);
		DH.put("T_C",82);
		DH.put("T_G",85);
		DH.put("T_T",79);
		DH.put("T_N",72);

		DH.put("N_A",72);
		DH.put("N_C",80);
		DH.put("N_G",78);
		DH.put("N_T",72);
		DH.put("N_N",72);

		/* Delta G's of disruption * 1000. */
		DG.put("A_A",1000);
		DG.put("A_C",1440);
		DG.put("A_G",1280);
		DG.put("A_T",880);
		DG.put("A_N",880);

		DG.put("C_A",1450);
		DG.put("C_C",1840);
		DG.put("C_G",2170);
		DG.put("C_T",1280);
		DG.put("C_N",1450);

		DG.put("G_A",1300);
		DG.put("G_C",2240);
		DG.put("G_G",1840);
		DG.put("G_T",1440);
		DG.put("G_N",1300);

		DG.put("T_A", 580);
		DG.put("T_C",1300);
		DG.put("T_G",1450);
		DG.put("T_T",1000);
		DG.put("T_N", 580);

		DG.put("N_A", 580);
		DG.put("N_C",1300);
		DG.put("N_G",1280);
		DG.put("N_T", 880);
		DG.put("N_N", 580);
		
	}
	
	
	public static void main(String[] argv)
	{
		String s = "tattggtgaagcctcaggtagtgcagaatatgaaacttcaggatccagtgggcatgctactggtagtgctgccggccttacaggcattatggtggcaaagtcgacagagtttc".toUpperCase();
		Sequence seq = new Sequence(s.toCharArray());
		
		double mv = 100, d = 50;
		double dv = 0, n = 0;
		
		double tm = oligoTM(seq, d, mv, dv, n, MeltingTemperatureMethod.santalucia_auto, SaltCorrectionMethod.schildkraut);
		System.out.println("Tm = " + tm);
	}
	
	
}
