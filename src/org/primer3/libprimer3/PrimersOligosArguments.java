package org.primer3.libprimer3;

import org.primer3.p3_seq_lib.seq_lib;

/**
 * Arguments for individual oligos and/or primers 
 *
 *
 */
public class PrimersOligosArguments {

	public seq_lib       repeat_lib = null;
	public oligo_weights weights = new oligo_weights();
	
	private double optTm;
	private double minTm;
	private double maxTm;
	
	private double optGC;
	private double maxGC;
	private double minGC;

	
	private int    optSize;
	private int    minSize;
	private int    maxSize;
	
	
	/* 
	 * SALT_MONOVALENT
	 * The millimolar (mM) concentration of monovalent salt cations (usually KCl)
	 * Warning: also used for product Tm (TO DO: factor this out) */
	private double saltConcentration;

	/**
	 * SALT_DIVALENT 
	 * The millimolar concentration of divalent salt cations (usually MgCl^(2+))
	 * Primer3 converts concentration of divalent cations to concentration of monovalent cations using formula suggested in the paper 
	 * [Ahsen von N, Wittwer CT, Schutz E (2001) "Oligonucleotide Melting Temperatures under PCR Conditions: Nearest-Neighbor Corrections 
	 * for Mg^(2+), Deoxynucleotide Triphosphate, and Dimethyl Sulfoxide Concentrations with Comparison to Alternative Empirical Formulas", 
	 * Clinical Chemistry 47:1956-61 http://www.clinchem.org/cgi/content/full/47/11/1956].
	 * [Monovalent cations] = [Monovalent cations] + 120*(([divalent cations] - [dNTP])^0.5)

	 */
	private double divalentConcentration;
	/*
    DIVALENT_CONC and DNTP_CONC are both needed for enabling use of
    divalent cations for calculation of melting temperature of short
    and long oligos.  The formula for converting the divalent cations
    to monovalent cations is in the paper [Ahsen von N, Wittwer CT,
    Schutz E (2001) "Oligonucleotide Melting Temperatures under PCR
    Conditions: Nearest-Neighbor Corrections for Mg^2+,
    Deoxynucleotide Triphosphate, and Dimethyl Sulfoxide
    Concentrations with Comparision to Alternative Empirical
    Formulas", Clinical Chemistry 47:1956-61
    http://www.clinchem.org/cgi/content/full/47/11/1956] The default
    is 0.0.  (New in v. 1.1.0, added by Maido Remm and Triinu
    Koressaar.)
	 */

	/**
	 * The millimolar concentration of the sum of all deoxyribonucleotide triphosphates. 
	 * A reaction mix containing 0.2 mM ATP, 0.2 mM CTP, 0.2 mM GTP and 0.2 mM TTP would have a PRIMER_DNTP_CONC=0.8.
	 */
	private double dntpConcentration;

	/**
	 * A value to use as nanomolar (nM) concentration of each annealing oligo over the course the PCR. 
	 * used with oligotm
	 */
	private double dnaConcentration;
	
	
	
	private int    maxNumOfNsAccepted;
	
	

	
	/** 
	 * Maximum length of mononucleotide sequence in an oligo.
	 */
	private int    maxPolyX;      

	/**
	 * The minimum sequence quality (as specified by SEQUENCE_QUALITY) allowed within the 5' pentamer of a primer. 
	 */
	private int    minEndQuality;
	/**
	 *  Minimum quality permitted for oligo sequence.
	 */
	private int    minQuality;       

	
	/**
	 * only used when thermodynamic_oligo_alignment is 0
	 */
	private double maxSelfAny;  
	private double maxSelfEnd;
	
	
	/**
	 * used when thermodynamic_oligo_alignment is 1
	 * tendency of a primer to bind to itself
	 */
	private double maxSelfAnyTH;
	
	/**
	 * tries to bind the 3'-END to a identical primer.
	 */
	private double maxSelfEndTH;
	private double maxHairPinTH;
	
	/** 
	 * Acceptable complementarity with repeat
	 * sequences.
	 */
	private double maxRepeatCompl;   

	private double maxTemplateMispriming;
	private double maxTemplateMisprimingTH;

	// orginal char *
	char[] must_match_five_prime;
	// orginal char *
	char[] must_match_three_prime;
	/* Primers and Oligos must match this 5 prime and 3 prime sequences.
	     This allows to select a set of primer pair with an identical 3' end
	     to avoid primer dimers. On the 5 prime end bases quenching flourochromes can be avoided.
	     The sequence must be 5 nucletides long and can contain the following letters:

	     N Any nucleotide

	     A Adenine
	     G Guanine
	     C Cytosine
	     T Thymine

	     R Purine (A or G)
	     Y Pyrimidine (C or T)
	     W Weak (A or T)
	     S Strong (G or C)
	     M Amino (A or C)
	     K Keto (G or T)
	     B Not A (G or C or T)
	     H Not G (A or C or T)
	     D Not C (A or G or T)
	     V Not T (A or G or C)
	 */

	/**
	 * PRIMER_MIN_QUALITY
	 * @param datum
	 */
	public void setMinQuality(String datum) {
		this.setMinQuality(Integer.parseInt(datum));		
	}

	

	/**
	 * @return the minQuality
	 */
	public int getMinQuality() {
		return minQuality;
	}



	/**
	 * @param minQuality the minQuality to set
	 */
	public void setMinQuality(int minQuality) {
		this.minQuality = minQuality;
	}



	/**
	 * PRIMER_MAX_NS_ACCEPTED : Maximum number of unknown bases (N) allowable in any primer.
	 * @param datum
	 */
	public void setMaxNumOfNsAccepted(String datum) {
		this.setMaxNumOfNsAccepted(Integer.parseInt(datum));
	}

	/**
	 * @return the maxNumOfNsAccepted
	 */
	public int getMaxNumOfNsAccepted() {
		return maxNumOfNsAccepted;
	}

	/**
	 * @param maxNumOfNsAccepted the maxNumOfNsAccepted to set
	 */
	public void setMaxNumOfNsAccepted(int maxNumOfNsAccepted) {
		this.maxNumOfNsAccepted = maxNumOfNsAccepted;
	}

	/**
	 * PRIMER_DNA_CONC
	 * @param datum
	 */
	public void setDnaConcentration(String datum) {
		this.setDnaConcentration(Double.parseDouble(datum));		

	}

	/**
	 * @return the dnaConcentration
	 */
	public double getDnaConcentration() {
		return dnaConcentration;
	}

	/**
	 * @param dnaConcentration the dnaConcentration to set
	 */
	public void setDnaConcentration(double dnaConcentration) {
		this.dnaConcentration = dnaConcentration;
	}

	/**
	 * PRIMER_DNTP_CONC
	 * @param datum
	 */
	public void setDntpConcentration(String datum) {
		this.setDntpConcentration(Double.parseDouble(datum));		
	}

	/**
	 * @return the dntpConcentration : The millimolar concentration of the sum of all deoxyribonucleotide triphosphates
	 */
	public double getDntpConcentration() {
		return dntpConcentration;
	}

	/**
	 * @param dntpConcentration The millimolar concentration of the sum of all deoxyribonucleotide triphosphates to set
	 */
	public void setDntpConcentration(double dntpConcentration) {
		this.dntpConcentration = dntpConcentration;
	}

	/**
	 * PRIMER_SALT_DIVALENT
	 * @param datum
	 */
	public void setDivalentConcentration(String datum) {
		this.setDivalentConcentration(Double.parseDouble(datum));
	}

	/**
	 * @return the divalentConcentration : the millimolar concentration of divalent salt cations
	 */
	public double getDivalentConcentration() {
		return divalentConcentration;
	}

	/**
	 * @param divalentConcentration the millimolar concentration of divalent salt cations to set
	 */
	public void setDivalentConcentration(double divalentConcentration) {
		this.divalentConcentration = divalentConcentration;
	}

	/**
	 * SALT MONOVALENT Concentration
	 * @param datum
	 */
	public void setSaltConcentration(String datum) {
		this.setSaltConcentration(Double.parseDouble(datum));		

	}

	/**
	 * @return the SALT MONOVALENT Concentration
	 */
	public double getSaltConcentration() {
		return saltConcentration;
	}

	/**
	 * @param saltConcentration the SALT MONOVALENT Concentration to set
	 */
	public void setSaltConcentration(double saltConcentration) {
		this.saltConcentration = saltConcentration;
	}

	/**
	 * PRIMER_MAX_GC
	 * @param datum
	 */
	public void setMaxGC(String datum) {
		this.setMaxGC(Double.parseDouble(datum));		

	}

	/**
	 * @return the maxGC
	 */
	public  double getMaxGC() {
		return maxGC;
	}

	/**
	 * @param maxGC the maxGC to set
	 */
	public  void setMaxGC(double maxGC) {
		this.maxGC = maxGC;
	}

	/**
	 * PRIMER_MIN_GC
	 * @param datum
	 */
	public void setMinGC(String datum) {
		this.setMinGC(Double.parseDouble(datum));		
	}

	/**
	 * @return the minGC
	 */
	public double getMinGC() {
		return minGC;
	}

	/**
	 * @param minGC the minGC to set
	 */
	public void setMinGC(double minGC) {
		this.minGC = minGC;
	}

	/**
	 * PRIMER_MAX_TM
	 * @param datum
	 */
	public void setMaxTm(String datum) {
		this.maxTm = Double.parseDouble(datum);		
	}

	public void setMaxTm(double value) {
		this.maxTm = value;		
	}
	
	public double getMaxTm()
	{
		return this.maxTm;
	}
	
	/**
	 * PRIMER_MIN_TM
	 * @param datum
	 */
	public void setMinTm(String datum) {
		this.setMinTm(Double.parseDouble(datum));

	}

	/**
	 * @return the minTm
	 */
	public double getMinTm() {
		return minTm;
	}

	/**
	 * @param minTm the minTm to set
	 */
	public void setMinTm(double minTm) {
		this.minTm = minTm;
	}

	/**
	 * PRIMER_OPT_GC_PERCENT
	 * @param datum
	 */
	public void setOptGCContent(String datum) {
		this.setOptGC(Double.parseDouble(datum));		
	}

	/**
	 * @return the optGC
	 */
	public double getOptGC() {
		return optGC;
	}

	/**
	 * @param optGC the optGC to set
	 */
	public void setOptGC(double optGC) {
		this.optGC = optGC;
	}

	/**
	 * PRIMER_OPT_TM
	 * @param datum
	 */
	public void setOptTm(String datum) {
		this.optTm = Double.parseDouble(datum);		
	}

	public void setOptTm(double value) {
		this.optTm = value;		
	}
	
	public double getOptTm()
	{
		return this.optTm;
	}
	
	/**
	 * PRIMER_MAX_POLY_X
	 * @param datum
	 */
	public void setMaxPolyX(String datum) {
		this.setMaxPolyX(Integer.parseInt(datum));		

	}

	/**
	 * @return the maxPolyX
	 */
	public int getMaxPolyX() {
		return maxPolyX;
	}



	/**
	 * @param maxPolyX the maxPolyX to set
	 */
	public void setMaxPolyX(int maxPolyX) {
		this.maxPolyX = maxPolyX;
	}



	/**
	 * PRIMER_MAX_SIZE
	 * @param datum
	 */
	public void setMaxSize(String datum) {
		this.setMaxSize(Integer.parseInt(datum));		

	}

	/**
	 * @return the maxSize
	 */
	public int getMaxSize() {
		return maxSize;
	}

	/**
	 * @param maxSize the maxSize to set
	 */
	public void setMaxSize(int maxSize) {
		this.maxSize = maxSize;
	}

	/**
	 * PRIMER_MIN_SIZE
	 * @param datum
	 */
	public void setMinSize(String datum) {
		this.setMinSize(Integer.parseInt(datum));		

	}

	/**
	 * @return the minSize
	 */
	public int getMinSize() {
		return minSize;
	}

	/**
	 * @param minSize the minSize to set
	 */
	public void setMinSize(int minSize) {
		this.minSize = minSize;
	}

	/**
	 * PRIMER_OPT_SIZE
	 * @param datum
	 */
	public void setOptSize(String datum) {
		this.setOptSize(Integer.parseInt(datum));		
	}

	/**
	 * @return the optSize
	 */
	public int getOptSize() {
		return optSize;
	}

	/**
	 * @param optSize the optSize to set
	 */
	public void setOptSize(int optSize) {
		this.optSize = optSize;
	}

	/**
	 * PRIMER_MIN_END_QUALITY
	 * @param datum
	 */
	public void setMinEndQuality(String datum) {
		this.setMinEndQuality(Integer.parseInt(datum));		
	}

	/**
	 * @return the minEndQuality
	 */
	public int getMinEndQuality() {
		return minEndQuality;
	}



	/**
	 * @param minEndQuality the minEndQuality to set
	 */
	public void setMinEndQuality(int minEndQuality) {
		this.minEndQuality = minEndQuality;
	}



	/**
	 * PRIMER_MAX_SELF_ANY_TH
	 * @param datum
	 */
	public void setMaxSelfAnyTH(String datum) {
		this.setMaxSelfAnyTH(Double.parseDouble(datum));
	}

	/**
	 * @return the maxSelfAnyTH
	 */
	public double getMaxSelfAnyTH() {
		return maxSelfAnyTH;
	}



	/**
	 * @param maxSelfAnyTH the maxSelfAnyTH to set
	 */
	public void setMaxSelfAnyTH(double maxSelfAnyTH) {
		this.maxSelfAnyTH = maxSelfAnyTH;
	}



	/**
	 * PRIMER_MAX_SELF_END_TH
	 * @param datum
	 */
	public void setMaxSelfEndTH(String datum) {
		this.setMaxSelfEndTH(Double.parseDouble(datum));
	}

	/**
	 * @return the maxSelfEndTH
	 */
	public double getMaxSelfEndTH() {
		return maxSelfEndTH;
	}



	/**
	 * @param maxSelfEndTH the maxSelfEndTH to set
	 */
	public void setMaxSelfEndTH(double maxSelfEndTH) {
		this.maxSelfEndTH = maxSelfEndTH;
	}



	/**
	 * PRIMER_MAX_HAIRPIN_TH
	 * @param datum
	 */
	public void setMaxHairPinTH(String datum) {
		this.setMaxHairPinTH(Double.parseDouble(datum));
	}

	/**
	 * @return the maxHairPinTH
	 */
	public double getMaxHairPinTH() {
		return maxHairPinTH;
	}
	



	/**
	 * @param maxHairPinTH the maxHairPinTH to set
	 */
	public void setMaxHairPinTH(double maxHairPinTH) {
		this.maxHairPinTH = maxHairPinTH;
	}
	



	/**
	 * 
	 * @param datum
	 */
	public void setMaxRepeatCompl(String datum) {
		this.setMaxRepeatCompl(Double.parseDouble(datum));
	}

	/**
	 * @return the maxRepeatCompl
	 */
	public double getMaxRepeatCompl() {
		return maxRepeatCompl;
	}



	/**
	 * @param maxRepeatCompl the maxRepeatCompl to set
	 */
	public void setMaxRepeatCompl(double maxRepeatCompl) {
		this.maxRepeatCompl = maxRepeatCompl;
	}



	public void setMaxTemplateMispriming(String datum) {
		this.setMaxTemplateMispriming(Double.parseDouble(datum));
	}

	/**
	 * @return the maxTemplateMispriming
	 */
	public double getMaxTemplateMispriming() {
		return maxTemplateMispriming;
	}



	/**
	 * @param maxTemplateMispriming the maxTemplateMispriming to set
	 */
	public void setMaxTemplateMispriming(double maxTemplateMispriming) {
		this.maxTemplateMispriming = maxTemplateMispriming;
	}



	public void setMaxTemplateMisprimingTH(String datum) {
		this.setMaxTemplateMisprimingTH(Double.parseDouble(datum));
	}

	/**
	 * @return the maxTemplateMisprimingTH
	 */
	public double getMaxTemplateMisprimingTH() {
		return maxTemplateMisprimingTH;
	}



	/**
	 * @param maxTemplateMisprimingTH the maxTemplateMisprimingTH to set
	 */
	public void setMaxTemplateMisprimingTH(double maxTemplateMisprimingTH) {
		this.maxTemplateMisprimingTH = maxTemplateMisprimingTH;
	}



	public void set_must_match_five_prime(String datum) {
		this.must_match_five_prime = datum.toCharArray();
	}

	public void set_must_match_three_prime(String datum) {
		this.must_match_three_prime = datum.toCharArray();		
	}

	/**
	 * PRIMER_MAX_SELF_ANY
	 * @param datum
	 */
	public void setMaxSelfAny(String datum) {
		this.setMaxSelfAny(Double.parseDouble(datum));		
	}
	
	
	/**
	 * @return the maxSelfAny
	 */
	public double getMaxSelfAny() {
		return maxSelfAny;
	}



	/**
	 * @param maxSelfAny the maxSelfAny to set
	 */
	public void setMaxSelfAny(double maxSelfAny) {
		this.maxSelfAny = maxSelfAny;
	}



	/**
	 * PRIMER_MAX_SELF_END/PRIMER_INTERNAL_MAX_SELF_END
	 * @param datum
	 */
	public void setMaxSelfEnd(String datum) {
		this.setMaxSelfEnd(Double.parseDouble(datum));		
	}



	/**
	 * @return the maxSelfEnd
	 */
	public double getMaxSelfEnd() {
		return maxSelfEnd;
	}



	/**
	 * @param maxSelfEnd the maxSelfEnd to set
	 */
	public void setMaxSelfEnd(double maxSelfEnd) {
		this.maxSelfEnd = maxSelfEnd;
	}
}