package org.primer3.libprimer3;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import org.apache.commons.lang3.tuple.MutablePair;
import org.primer3.masker.MaskerParameters;
import org.primer3.masker.masker;
import org.primer3.masker.masking_direction;
import org.primer3.oligotm.MeltingTemperatureMethod;
import org.primer3.oligotm.SaltCorrectionMethod;

public class P3GlobalSettings {

	/* ================================================== */
	/* Arguments that control behavior of choose_primers() */

	public static P3GlobalSettings p3_create_global_settings(){
		P3GlobalSettings r = new P3GlobalSettings();
		r.pr_set_default_global_args_2();
		return r;
	} 

	public static P3GlobalSettings p3_create_global_settings(
			int default_version) {
		if (default_version == 1)
			return p3_create_global_settings_default_version_1();
		else if (default_version == 2) 
			return p3_create_global_settings();
		return null;
	}
	public static P3GlobalSettings p3_create_global_settings_default_version_1(){
		P3GlobalSettings r = new P3GlobalSettings();
		r.pr_set_default_global_args_1();
		return r;
	}
	/** 
	 * 2 if left primer only, 3 if right primer
	 * only, 4 if internal oligo only.  
	 */
	private P3Task   primerTask;

	private boolean    pickLeftPrimer;   

	private boolean    pickRightPrimer;  

	private boolean    pickInternalOligo;   

	/* TODO, See if this can be factored out. */
	int    file_flag; 

	/** 
	 * The index of the first base in the input
	 * sequence.  This parameter is ignored within
	 * pr_choice; pr_choice's caller must assure that
	 * all indexes are 0-based.  However, this
	 * parameter should used by output routines to
	 * adjust base indexes.
	 */
	private int    firstBaseIndex;    

	/** 
	 * If non-0 then turn characters other than
	 * [ATGCNatgcn] into N.
	 */
	private boolean    liberalBase;


	/** The number of best primer pairs to return. */
	private int    numReturn;
	/** Pick even if input primer or oligos 
	 * violate constraints. 
	 */
	private boolean    pickAnyway;

	/* ================================================== */
	/* Arguments for individual oligos and/or primers */

	/** If non-0, treat ambiguity codes in a mispriming/mishyb
	 * library as representing a consensus.  So, for example,
	 * S would match C or G.  N would match any nucleotide.
	 * It turns out that this _not_ what one normally wants,
	 * since many libraries contain strings of N, which then
	 * match every oligo (very bad).
	 */
	private boolean    libAmbiguityCodesConsensus;

	private int    qualityRangeMin;

	private int    qualityRangeMax;

	public PrimersOligosArguments primersArgs = new PrimersOligosArguments();  

	public PrimersOligosArguments oligosArgs = new PrimersOligosArguments(); 


	/* ================================================== */
	/* Start of arguments applicable to primers but not
     oligos */


	/** 
	 * tm_santalucia added by T.Koressaar for updated table
	 * thermodynamics.  Specifies details of melting temperature
	 * calculation.  (New in v. 1.1.0, added by Maido Remm and Triinu
	 * Koressaar.) 
	 * A value of 1 (recommended) directs primer3 to use the table of
	 * thermodynamic values and the method for melting temperature
	 * calculation suggested in the paper [SantaLucia JR (1998) "A
	 * unified view of polymer, dumbbell and oligonucleotide DNA
	 * nearest-neighbor thermodynamics", Proc Natl Acad Sci 95:1460-65
	 * http://dx.doi.org/10.1073/pnas.95.4.1460].
	 *  
	 * A value of 0 directs primer3 to a backward compatible calculation
	 * (in other words, the only calculation availble in previous
	 * version of primer3).
	 *  
	 * This backward compatible calculation uses the table of
	 * thermodynamic parameters in the paper [Breslauer KJ, Frank R,
	 * Bloecker H and Marky LA (1986) "Predicting DNA duplex stability
	 * from the base sequence" Proc Natl Acad Sci 83:4746-50
	 * http://dx.doi.org/10.1073/pnas.83.11.3746],
	 * and the method in the paper [Rychlik W, Spencer WJ and Rhoads
	 * RE (1990) "Optimization of the annealing temperature for DNA
	 * amplification in vitro", Nucleic Acids Res 18:6409-12
	 * http://www.pubmedcentral.nih.gov/articlerender.fcgi?tool=pubmed&pubmedid=2243783].
	 *   
	 * The default value is 0 only for backward compatibility.
	 */
	public MeltingTemperatureMethod mtMethod;


	/** 
     salt_corrections added by T.Koressaar for salt correction for Tm
     calculation.  A value of 1 (recommended) directs primer3 to use
     the salt correction formula in the paper [SantaLucia JR (1998) "A
     unified view of polymer, dumbbell and oligonucleotide DNA
     nearest-neighbor thermodynamics", Proc Natl Acad Sci 95:1460-65
     http://dx.doi.org/10.1073/pnas.95.4.1460]

     A value of 0 directs primer3 to use the the salt correction
     formula in the paper [Schildkraut, C, and Lifson, S (1965)
     "Dependence of the melting temperature of DNA on salt
     concentration", Biopolymers 3:195-208 (not available on-line)].
     This was the formula used in previous version of primer3.

     A value of 2 directs primer3 to use the salt correction formula
     in the paper [Owczarzy R, You Y, Moreira BG, Manthey JA, Huang L,
     Behlke MA and Walder JA (2004) "Effects of sodium ions on DNA
     duplex oligomers: Improved predictions of melting temperatures",
     Biochemistry 43:3537-54 http://dx.doi.org/10.1021/bi034621r].

     The default is 0 only for backward compatibility.
	 */
	private SaltCorrectionMethod saltCorrectionMethod;


	/** The maximum value allowed for the delta
	 * G of disruption for the 5 3' bases of
	 * a primer.
	 */
	private double maxEndStability;  

	/* ================================================== */
	/* Start of arguments related to primer and/or oligo
     location in the template. */

	/** The maximum number of allowed G or C
	 * for the 5 3' bases of a primer.
	 */
	private int    maxEndGC;

	/**
	 *  Required number of GCs at 3' end. 
	 */
	private int    gcClamp;

	private boolean lowercaseMasking; 
	/* 
     lowercase_masking added by T.Koressaar. Enables design of primers
     from lowercase masked template.  A value of 1 directs primer3 to
     reject primers overlapping lowercase a base exactly at the 3'
     end.

     This property relies on the assumption that masked features
     (e.g. repeats) can partly overlap primer, but they cannot overlap the
     3'-end of the primer.  In other words, lowercase bases at other
     positions in the primer are accepted, assuming that the masked
     features do not influence the primer performance if they do not
     overlap the 3'-end of primer.
	 */ 

	private SequencingParameters sequencingParameters = new SequencingParameters();  /* Used to calculate the position
					of sequencing primers */  

	/* ================================================== */
	/* Arguments for primer pairs and products. */

	/* Warning: Use p3_empty_gs_product_size_range and
     p3_add_to_gs_product_size_range to set these next
     three slots. */
	/** Minimum product sizes. */

	/** Multiply this value times the number of NTs
	 * from the 3' end to the the (unique) target to
	 * get the 'position penalty'.
	 * Meaningless if there are multiple targets
	 * or if the primer cannot be part of a pair
	 * that spans the target.
	 */
	private double outsidePenalty;
	/** Multiply this value times the number of NT
	 * positions by which the primer overlaps
	 * the (unique) target to the 'position penalty'.
	 * Meaningless if there are multiple targets
	 * or if the primer cannot be part of a pair
	 * that spans the target.
	 */
	private double insidePenalty;
	/**
	 * <Left : Minimum product sizes, Right Maximum product sizes>
	 */
	private List<MutablePair<Integer,Integer>> productSizeRanges = new ArrayList<MutablePair<Integer,Integer>>();
	//	int[]    pr_min = new int[LibPrimer3.PR_MAX_INTERVAL_ARRAY];
	/** Maximum product sizes. */
	//	int[]    pr_max = new int[LibPrimer3.PR_MAX_INTERVAL_ARRAY]; 
	/** 
	 * Number of product size intervals
	 * (i.e. number of elements in pr_min and
	 * pr_max)
	 */
	 //	int    num_intervals;        

	private int    productOptSize;
	
	private double productMaxTM;
	private double productMinTM;
	private double productOptTM;
	
	private double pairMaxTemplateMispriming;
	private double pairMaxTemplateMisprimingTH;
	private double pairRepeatCompl;
	private double pairComplAny;
	private double pairComplAnyTH;




	// refactor to true / false
//	public static final int USE_THERMODYNAMICS_ALIGNMENT = 1;
//	public static final int DONOT_USE_THERMODYNAMICS_ALIGNMENT = 0;

	private double pairComplEnd;
	private double pairComplEndTH;


	/** 
	 * Enables to use approach of thermodynamical alignment for dimer and hairpin calculations.
	 * <br>
	 * <code>DONOT_USE_THERMODYNAMICS_ALIGNMENT</code> 0 = Use alignment *not* based on thermodynamics - the only dimer
	 * calculation approach until primer3 2.2.0.  No hairpins are calculated.
	 * <br>
	 * <code>USE_THERMODYNAMICS_ALIGNMENT</code> 1 = Use alignment based on thermodynamics. Hairpins are calculated
	 * TODO :: PR_ASSERT ( thermodynamic_oligo_alignment == 0 || == 1 ) 
	 */
	private boolean thermodynamicOligoAlignment; 


	/** 
	 * Enables to use approach of thermodynamical alignment for template 
	 * mispriming calculations.
	 * 0 = Use alignment *not* based on thermodynamics.
	 * 1 = Use alignment based on thermodynamics.
	 */
	private boolean thermodynamicTemplateAlignment;

	/** Maximum allowed difference between temperature of primer and
	 * temperature of product.  Cannot be calculated until product is
	 * known. */
	private double maxDiffTm;
	
	private PairWeights  prPairWeights = new PairWeights();

	/** Minimum number of base pairs between the 3' ends of any two left
	 * or any two right primers when returning num_return primer pairs.
	 * The objective is get 'truly different' primer pairs.
	 * Please see the user documentation (primer3_manual.htm) for
	 * PRIMER_{LEFT,RIGHT}_MIN_THREE_PRIME_DISTANCE.
	 */
	private int    minLeft3PrimeDistance;   
	private int    minRight3PrimeDistance;

	/** The number of basepairs
	 * the primer has to overlap an overlap
	 * junction. 
	 */
	private int    min5PrimeOverlapOfJunction;

	private int    min3PrimeOverlapOfJunction;

	/** Turn on masking of the trimmed_orig_seq (added by M. Lepamets)*/
	private boolean    maskTemplate;
	private boolean    maskingParametersChanged;




	/** 
	 * a struct containing all masking parameters
	 */
	private MaskerParameters maskingParameters = new MaskerParameters();

	/** dump fields for global settings and seq args if dump == 1 */
	private boolean dump;

	


	/**
	 * @return the firstBaseIndex
	 */
	public int getFirstBaseIndex() {
		return firstBaseIndex;
	}





	/**
	 * @return the gcClamp
	 */
	public int getGcClamp() {
		return gcClamp;
	}

	/**
	 * @return the insidePenalty
	 */
	public double getInsidePenalty() {
		return insidePenalty;
	}


	/**
	 * @return the maskingParameters
	 */
	public MaskerParameters getMaskingParameters() {
		return maskingParameters;
	}

	/**
	 * @return the maxDiffTm
	 */
	public double getMaxDiffTm() {
		return maxDiffTm;
	}

	/**
	 * @return the maxEndGC
	 */
	public int getMaxEndGC() {
		return maxEndGC;
	}


	/**
	 * @return the maxEndStability
	 */
	public double getMaxEndStability() {
		return maxEndStability;
	}  

	/**
	 * @return the min3PrimeOverlapOfJunction
	 */
	public int getMin3PrimeOverlapOfJunction() {
		return min3PrimeOverlapOfJunction;
	}


	/**
	 * @return the min5PrimeOverlapOfJunction
	 */
	public int getMin5PrimeOverlapOfJunction() {
		return min5PrimeOverlapOfJunction;
	}


	/**
	 * @return the minLeft3PrimeDistance
	 */
	public int getMinLeft3PrimeDistance() {
		return minLeft3PrimeDistance;
	}
	/**
	 * 
	 * @param value
	 */
	public void set_min_left_three_prime_distance(int value) {
		this.minLeft3PrimeDistance = value;
	}


	/**
	 * @return the minRight3PrimeDistance
	 */
	public int getMinRight3PrimeDistance() {
		return minRight3PrimeDistance;
	}

	
	public void set_min_right_three_prime_distance(int value) {
		this.minRight3PrimeDistance = value;
	}
	
	/**
	 * @return the numReturn
	 */
	public int getNumReturn() {
		return numReturn;
	}

	/**
	 * @return the outsidePenalty
	 */
	public double getOutsidePenalty() {
		return outsidePenalty;
	}

	/**
	 * @return the pairComplAny
	 */
	public double getPairComplAny() {
		return pairComplAny;
	}

	/**
	 * @return the pairComplAnyTH
	 */
	public double getPairComplAnyTH() {
		return pairComplAnyTH;
	}

	/**
	 * @return the pairComplEnd
	 */
	public double getPairComplEnd() {
		return pairComplEnd;
	}

	/**
	 * @return the pairComplEndTH
	 */
	public double getPairComplEndTH() {
		return pairComplEndTH;
	}

	/**
	 * @return the pairMaxTemplateMispriming
	 */
	public double getPairMaxTemplateMispriming() {
		return pairMaxTemplateMispriming;
	}

	/**
	 * @return the pairMaxTemplateMisprimingTH
	 */
	public double getPairMaxTemplateMisprimingTH() {
		return pairMaxTemplateMisprimingTH;
	}

	/**
	 * @return the pairRepeatCompl
	 */
	public double getPairRepeatCompl() {
		return pairRepeatCompl;
	}

	public PrimersOligosArguments getPrimersArgs()
	{
		return primersArgs;
	}
	
	public PrimersOligosArguments getOligosArgs()
	{
		return oligosArgs;
	}
	/**
	 * @return the primerTask
	 */
	public P3Task getPrimerTask() {
		return primerTask;
	}

	/**
	 * @return the productMaxTM
	 */
	public double getProductMaxTM() {
		return productMaxTM;
	}

	/**
	 * @return the productMinTM
	 */
	public double getProductMinTM() {
		return productMinTM;
	}

	/**
	 * @return the productOptSize
	 */
	public int getProductOptSize() {
		return productOptSize;
	}

	/**
	 * @return the productOptTM
	 */
	public double getProductOptTM() {
		return productOptTM;
	}

	public MutablePair<Integer,Integer> getProductSizeRange(int i)
	{
		return this.productSizeRanges.get(i);
	}
	
	public void addProductSizeRange(int pr_min,int pr_max) {
		this.productSizeRanges.add(MutablePair.of(pr_min, pr_max));
	}

	
	/**
	 * this will would set the first entry in the list assuming only one will work
	 * @param minSize
	 */
	public void setProdcutMinSize(int minSize)
	{
		this.productSizeRanges.get(0).setLeft(minSize);
	}
	
	/**
	 * this will would set the first entry in the list assuming only one will work
	 * @param minSize
	 */
	public void setProdcutMaxSize(int maxSize)
	{
		this.productSizeRanges.get(0).setRight(maxSize);
	}
	
	
	public int getProdcutMaxSize()
	{
		return productSizeRanges.get(0).getRight();
	}
	
	public int getProdcutMinSize()
	{
		return productSizeRanges.get(0).getLeft();
	}
	
	public int getProductSizeRangesNumber()
	{
		return this.productSizeRanges.size();
	}
	/**
	 * @return the prPairWeights
	 */
	public PairWeights getPrPairWeights() {
		return prPairWeights;
	}
	/**
	 * @return the qualityRangeMax
	 */
	public int getQualityRangeMax() {
		return qualityRangeMax;
	}
	/**
	 * @return the qualityRangeMin
	 */
	public int getQualityRangeMin() {
		return qualityRangeMin;
	}

	/**
	 * @return the saltCorrectionMethod
	 */
	public SaltCorrectionMethod getSaltCorrectionMethod() {
		return saltCorrectionMethod;
	}

	/**
	 * @return the sequencingParameters
	 */
	public SequencingParameters getSequencingParameters() {
		return sequencingParameters;
	}

	public boolean isDefaultPositionPenalties() {
		return (LibPrimer3.PR_DEFAULT_INSIDE_PENALTY == this.insidePenalty  && LibPrimer3.PR_DEFAULT_OUTSIDE_PENALTY == this.outsidePenalty);
	}

	/**
	 * @return the dump
	 */
	public boolean isDump() {
		return dump;
	}

	/**
	 * @return the libAmbiguityCodesConsensus
	 */
	public boolean isLibAmbiguityCodesConsensus() {
		return libAmbiguityCodesConsensus;
	}
	/**
	 * @return the liberalBase
	 */
	public boolean isLiberalBase() {
		return liberalBase;
	}

	/**
	 * @return the lowercaseMasking
	 */
	public boolean isLowercaseMasking() {
		return lowercaseMasking;
	}

	/**
	 * @return the maskingParametersChanged
	 */
	public boolean isMaskingParametersChanged() {
		return maskingParametersChanged;
	}

	/**
	 * @return the maskTemplate
	 */
	public boolean isMaskTemplate() {
		return maskTemplate;
	}





	/**
	 * @return the pickAnyway
	 */
	public boolean isPickAnyway() {
		return pickAnyway;
	}

	/**
	 * @return the pickInternalOligo
	 */
	public boolean isPickInternalOligo() {
		return pickInternalOligo;
	}

	/**
	 * @return the pickLeftPrimer
	 */
	public boolean isPickLeftPrimer() {
		return pickLeftPrimer;
	}

	/**
	 * @return the pickRightPrimer
	 */
	public boolean isPickRightPrimer() {
		return pickRightPrimer;
	}

	/**
	 * @return the thermodynamicOligoAlignment
	 */
	public boolean isThermodynamicOligoAlignment() {
		return thermodynamicOligoAlignment;
	}

	/**
	 * @return the thermodynamicTemplateAlignment
	 */
	public boolean isThermodynamicTemplateAlignment() {
		return thermodynamicTemplateAlignment;
	}

	public boolean needPairTemplateMispriming() {
		return this.pairMaxTemplateMispriming >= 0
				|| this.prPairWeights.template_mispriming > 0.0;
	}

	public boolean needPairTemplateMisprimingTH()
	{   
		return this.pairMaxTemplateMisprimingTH >= 0
				|| this.prPairWeights.template_mispriming_th > 0.0;
	}



	public boolean needTemplateMispriming() {
		return this.primersArgs.getMaxTemplateMispriming() >= 0
				|| this.primersArgs.weights.template_mispriming > 0.0
				|| this.needPairTemplateMispriming();
	}
	public boolean needTemplateMisprimingTH() {
		return
				this.primersArgs.getMaxTemplateMisprimingTH() >= 0
				|| this.primersArgs.weights.template_mispriming_th > 0.0
				|| this.needPairTemplateMisprimingTH();
	}
	public void p3_print_args() {
		P3GlobalSettings p = this;
		int i;
		System.out.format("=============\n");
		System.out.format("BEGIN GLOBAL ARGS\n") ;
		System.out.format("  primer_task %s\n", p.getPrimerTask());
		System.out.format("  pick_left_primer %s\n", p.isPickLeftPrimer());
		System.out.format("  pick_right_primer %s\n", p.isPickRightPrimer());
		System.out.format("  pick_internal_oligo %s\n", p.isPickInternalOligo());
		System.out.format("  file_flag %s\n", p.file_flag) ;
		System.out.format("  first_base_index %s\n", p.firstBaseIndex);
		System.out.format("  liberal_base %s\n", p.liberalBase );
		System.out.format("  num_return %s\n", p.numReturn) ;
		System.out.format("  pick_anyway %s\n", p.pickAnyway);
		System.out.format("  lib_ambiguity_codes_consensus %s\n", p.libAmbiguityCodesConsensus) ;
		System.out.format("  quality_range_min %s\n", p.qualityRangeMin) ;
		System.out.format("  quality_range_max %s\n", p.qualityRangeMax) ;

		System.out.format("  tm_santalucia %s\n", p.mtMethod) ;
		System.out.format("  salt_corrections %s\n", p.getSaltCorrectionMethod()) ;
		System.out.format("  max_end_stability %f\n", p.maxEndStability) ;
		System.out.format("  gc_clamp %s\n", p.getGcClamp()) ;
		System.out.format("  max_end_gc %s\n", p.getMaxEndGC());
		System.out.format("  lowercase_masking %s\n", p.lowercaseMasking) ;
		System.out.format("  thermodynamic_oligo_alignment %s\n", p.thermodynamicOligoAlignment);
		System.out.format("  thermodynamic_template_alignment %s\n", p.thermodynamicTemplateAlignment);
		System.out.format("  outside_penalty %f\n", p.outsidePenalty) ;
		System.out.format("  inside_penalty %f\n", p.insidePenalty) ;
		System.out.format("  number of product size ranges: %d\n", p.productSizeRanges.size());
		System.out.format("  product size ranges:\n");
		for (i = 0; i < p.productSizeRanges.size(); i++) {
			System.out.format("  %d - %d \n", p.productSizeRanges.get(i).getLeft() ,
					p.productSizeRanges.get(i).getRight());
		}
		System.out.format("  product_opt_size %s\n", p.productOptSize) ;
		System.out.format("  product_max_tm %f\n", p.productMaxTM) ;
		System.out.format("  product_min_tm %f\n", p.productMinTM) ;
		System.out.format("  product_opt_tm %f\n", p.productOptTM) ;
		System.out.format("  pair_max_template_mispriming %f\n", p.pairMaxTemplateMispriming) ;
		System.out.format("  pair_max_template_mispriming_th %f\n", p.pairMaxTemplateMisprimingTH) ;
		System.out.format("  pair_repeat_compl %f\n", p.pairRepeatCompl) ;
		System.out.format("  pair_compl_any %f\n", p.pairComplAny) ;
		System.out.format("  pair_compl_end %f\n", p.pairComplEnd) ;
		System.out.format("  pair_compl_any_th %f\n", p.pairComplAnyTH) ;
		System.out.format("  pair_compl_end_th %f\n", p.pairComplEndTH) ;

		System.out.format("  min_left_three_prime_distance %s\n", p.minLeft3PrimeDistance) ;
		System.out.format("  min_right_three_prime_distance %s\n", p.minRight3PrimeDistance) ;
		System.out.format("  min_5_prime_overlap_of_junction %s\n", p.min5PrimeOverlapOfJunction);
		System.out.format("  min_3_prime_overlap_of_junction %s\n", p.min3PrimeOverlapOfJunction);
		System.out.format("  mask_template %s\n", p.maskTemplate);
		System.out.format("  failure_rate %f\n", p.getMaskingParameters().failure_rate);
		System.out.format("  nucl_masked_in_5p_direction %s\n", p.getMaskingParameters().nucl_masked_in_5p_direction);
		System.out.format("  nucl_masked_in_3p_direction %s\n", p.getMaskingParameters().nucl_masked_in_3p_direction);
		System.out.format("  list_prefix %s\n", p.getMaskingParameters().list_prefix);
		System.out.format("  dump %s\n", p.dump);

		p.prPairWeights.p3_print_args();
		System.out.format("\n\n");
		System.out.format("=============\n");
		System.out.format("BEGIN primer_args\n");
		System.out.format("begin oligo_weights\n");
		System.out.format("temp_gt %f\n", p.primersArgs.weights.temp_gt ) ;
		System.out.format("temp_gt %f\n", p.primersArgs.weights.temp_gt) ;
		System.out.format("temp_lt %f\n", p.primersArgs.weights.temp_lt) ;
		System.out.format("gc_content_gt %f\n", p.primersArgs.weights.gc_content_gt) ;
		System.out.format("gc_content_lt %f\n", p.primersArgs.weights.gc_content_lt) ;
		System.out.format("compl_any %f\n", p.primersArgs.weights.compl_any) ;
		System.out.format("compl_end %f\n", p.primersArgs.weights.compl_end) ;
		System.out.format("compl_any_th %f\n", p.primersArgs.weights.compl_any_th) ;
		System.out.format("compl_end_th %f\n", p.primersArgs.weights.compl_end_th) ;
		System.out.format("hairpin %f\n", p.primersArgs.weights.hairpin_th) ;
		System.out.format("repeat_sim %f\n", p.primersArgs.weights.repeat_sim) ;
		System.out.format("length_lt %f\n", p.primersArgs.weights.length_lt) ;
		System.out.format("length_gt %f\n", p.primersArgs.weights.length_gt) ;
		System.out.format("seq_quality %f\n", p.primersArgs.weights.seq_quality) ;
		System.out.format("end_quality %f\n", p.primersArgs.weights.end_quality) ;
		System.out.format("pos_penalty %f\n", p.primersArgs.weights.pos_penalty) ;
		System.out.format("end_stability %f\n", p.primersArgs.weights.end_stability) ;
		System.out.format("num_ns %f\n", p.primersArgs.weights.num_ns) ;
		System.out.format("template_mispriming %f\n", p.primersArgs.weights.template_mispriming) ;
		System.out.format("template_mispriming_th %f\n", p.primersArgs.weights.template_mispriming_th) ;
		System.out.format("failure_rate %f\n", p.primersArgs.weights.failure_rate) ;
		System.out.format("end oligo_weights\n") ;

		System.out.format("opt_tm %f\n", p.primersArgs.getOptTm()) ;
		System.out.format("min_tm %f\n", p.primersArgs.getMinTm()) ;
		System.out.format("max_tm %f\n", p.primersArgs.getMaxTm()) ;
		System.out.format("opt_gc_content %f\n", p.primersArgs.getOptGC()) ;
		System.out.format("max_gc %f\n", p.primersArgs.getMaxGC()) ;
		System.out.format("min_gc %f\n", p.primersArgs.getMinGC()) ;
		System.out.format("divalent_conc %f\n", p.primersArgs.getDivalentConcentration()) ;
		System.out.format("dntp_conc %f\n", p.primersArgs.getDntpConcentration()) ;
		System.out.format("dna_conc %f\n", p.primersArgs.getDnaConcentration()) ;
		System.out.format("num_ns_accepted %s\n", p.primersArgs.getMaxNumOfNsAccepted()) ;
		System.out.format("opt_size %s\n", p.primersArgs.getOptSize()) ;
		System.out.format("min_size %s\n", p.primersArgs.getMinSize()) ;
		System.out.format("max_size %s\n", p.primersArgs.getMaxSize()) ;
		System.out.format("max_poly_x %s\n", p.primersArgs.getMaxPolyX()) ;
		System.out.format("min_end_quality %s\n", p.primersArgs.getMinEndQuality()) ;
		System.out.format("min_quality %s\n", p.primersArgs.getMinQuality()) ;
		System.out.format("max_self_any %f\n", p.primersArgs.getMaxSelfAny()) ;
		System.out.format("max_self_end %f\n", p.primersArgs.getMaxSelfEnd()) ;
		System.out.format("max_self_any_th %f\n", p.primersArgs.getMaxSelfAnyTH()) ;
		System.out.format("max_self_end_th %f\n", p.primersArgs.getMaxSelfEndTH()) ;
		System.out.format("max_hairpin %f\n", p.primersArgs.getMaxHairPinTH()) ;
		System.out.format("max_repeat_compl %f\n", p.primersArgs.getMaxRepeatCompl()) ;
		System.out.format("max_template_mispriming %f\n", p.primersArgs.getMaxTemplateMispriming()) ;
		System.out.format("max_template_mispriming_th %f\n", p.primersArgs.getMaxTemplateMisprimingTH()) ;
		System.out.format("end primer args\n") ;

		System.out.format("begin internal oligo args (p.o_args.)\n") ;

		System.out.format("  begin internal oligo_weights (p.o_args.weights.)\n") ;
		System.out.format("    temp_gt %f\n", p.oligosArgs.weights.temp_gt) ;
		System.out.format("    temp_lt %f\n", p.oligosArgs.weights.temp_lt) ;
		System.out.format("    gc_content_gt %f\n", p.oligosArgs.weights.gc_content_gt) ;
		System.out.format("    gc_content_lt %f\n", p.oligosArgs.weights.gc_content_lt) ;
		System.out.format("    compl_any %f\n", p.oligosArgs.weights.compl_any) ;
		System.out.format("    compl_end %f\n", p.oligosArgs.weights.compl_end) ;
		System.out.format("    compl_any_th %f\n", p.oligosArgs.weights.compl_any_th) ;
		System.out.format("    compl_end_th %f\n", p.oligosArgs.weights.compl_end_th) ;
		System.out.format("    hairpin %f\n", p.oligosArgs.weights.hairpin_th) ;
		System.out.format("    repeat_sim %f\n", p.oligosArgs.weights.repeat_sim) ;
		System.out.format("    length_lt %f\n", p.oligosArgs.weights.length_lt) ;
		System.out.format("    length_gt %f\n", p.oligosArgs.weights.length_gt) ;
		System.out.format("    seq_quality %f\n", p.oligosArgs.weights.seq_quality) ;
		System.out.format("    end_quality %f\n", p.oligosArgs.weights.end_quality) ;
		System.out.format("    pos_penalty %f\n", p.oligosArgs.weights.pos_penalty) ;
		System.out.format("    end_stability %f\n", p.oligosArgs.weights.end_stability) ;
		System.out.format("    num_ns %f\n", p.oligosArgs.weights.num_ns) ;
		System.out.format("  end internal oligo_weights\n") ;

		System.out.format("  opt_tm %f\n", p.oligosArgs.getOptTm()) ;
		System.out.format("  min_tm %f\n", p.oligosArgs.getMinTm()) ;
		System.out.format("  max_tm %f\n", p.oligosArgs.getMaxTm()) ;
		System.out.format("  opt_gc_content %f\n", p.oligosArgs.getOptGC()) ;
		System.out.format("  max_gc %f\n", p.oligosArgs.getMaxGC()) ;
		System.out.format("  min_gc %f\n", p.oligosArgs.getMinGC()) ;
		System.out.format("  divalent_conc %f\n", p.oligosArgs.getDivalentConcentration()) ;
		System.out.format("  dntp_conc %f\n", p.oligosArgs.getDntpConcentration()) ;
		System.out.format("  dna_conc %f\n", p.oligosArgs.getDnaConcentration()) ;
		System.out.format("  num_ns_accepted %s\n", p.oligosArgs.getMaxNumOfNsAccepted()) ;
		System.out.format("  opt_size %s\n", p.oligosArgs.getOptSize()) ;
		System.out.format("  min_size %s\n", p.oligosArgs.getMinSize()) ;
		System.out.format("  max_size %s\n", p.oligosArgs.getMaxSize()) ;
		System.out.format("  max_poly_x %s\n", p.oligosArgs.getMaxPolyX()) ;
		System.out.format("  min_end_quality %s\n", p.oligosArgs.getMinEndQuality()) ;
		System.out.format("  min_quality %s\n", p.oligosArgs.getMinQuality()) ;
		System.out.format("  max_self_any %f\n", p.oligosArgs.getMaxSelfAny()) ;
		System.out.format("  max_self_end %f\n", p.oligosArgs.getMaxSelfEnd()) ;
		System.out.format("  max_repeat_compl %f\n", p.oligosArgs.getMaxRepeatCompl()) ;
		System.out.format("  end internal oligo args\n");
		System.out.format("\n");
		System.out.format("END GLOBAL ARGS\n");
		System.out.format("=============\n");
		System.out.format("\n");
	}


	/**
	 * PRIMER_INTERNAL_MAX_SELF_END - o_args.set_max_self_end
	 * @param datum
	 */
	public void p3_set_gs_primer_internal_oligo_self_end(String datum) {
		this.oligosArgs.setMaxSelfEnd(datum);
	}

	void p3_set_gs_primer_task(String task_tmp)
	{
		if (task_tmp.equalsIgnoreCase( "pick_pcr_primers")) {
			this.setPrimerTask(P3Task.GENERIC);
			this.setPickLeftPrimer(true);
			this.setPickRightPrimer(true);
			this.setPickInternalOligo(false);
		} else if (task_tmp.equalsIgnoreCase( "pick_pcr_primers_and_hyb_probe")) {
			this.setPrimerTask(P3Task.GENERIC);
			this.setPickLeftPrimer(true);
			this.setPickRightPrimer(true);
			this.setPickInternalOligo(true);
		} else if (task_tmp.equalsIgnoreCase( "pick_left_only")) {
			this.setPrimerTask(P3Task.GENERIC);
			this.setPickLeftPrimer(true);
			this.setPickRightPrimer(false);
			this.setPickInternalOligo(false);
		} else if (task_tmp.equalsIgnoreCase( "pick_right_only")) {
			this.setPrimerTask(P3Task.GENERIC);
			this.setPickLeftPrimer(false);
			this.setPickRightPrimer(true);
			this.setPickInternalOligo(false);
		} else if (task_tmp.equalsIgnoreCase( "pick_hyb_probe_only")) {
			this.setPrimerTask(P3Task.GENERIC);
			this.setPickLeftPrimer(false);
			this.setPickRightPrimer(false);
			this.setPickInternalOligo(true);
		} else if (task_tmp.equalsIgnoreCase( "generic")) {
			this.setPrimerTask(P3Task.GENERIC);
		} else if (task_tmp.equalsIgnoreCase( "pick_detection_primers")) {
			this.setPrimerTask(P3Task.GENERIC); /* Deliberate duplication for
					    backward compatibility. */
		} else if (task_tmp.equalsIgnoreCase( "pick_cloning_primers")) {
			this.setPrimerTask(P3Task.PICK_CLONING_PRIMERS);
		} else if (task_tmp.equalsIgnoreCase( "pick_discriminative_primers")) {
			this.setPrimerTask(P3Task.PICK_DISCRIMINATIVE_PRIMERS);
		} else if (task_tmp.equalsIgnoreCase( "pick_sequencing_primers")) {
			this.setPrimerTask(P3Task.PICK_SEQUENCING_PRIMERS);
		} else if (task_tmp.equalsIgnoreCase( "pick_primer_list")) {
			this.setPrimerTask(P3Task.PICK_PRIMER_LIST);
		} else if (task_tmp.equalsIgnoreCase("check_primers")) {
			this.setPrimerTask(P3Task.CHECK_PRIMERS);
		}
	}

	/**
	 * PRIMER_FIRST_BASE_INDEX
	 * @param datum
	 */
	public void p3_set_pa_first_base_index(String datum) {
		this.firstBaseIndex = Integer.parseInt(datum);
	}



	/**
	 * PRIMER_GC_CLAMP
	 * @param datum
	 */
	public void p3_set_pa_gc_clamp(String datum) {
		this.setGcClamp(Integer.parseInt(datum));
	}

	/**
	 * PRIMER_INSIDE_PENALTY
	 * @param datum
	 */
	public void p3_set_pa_inside_penalty(String datum) {
		this.insidePenalty = Double.parseDouble(datum);
	}

	/**
	 * PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS
	 * @param datum
	 */

	public void p3_set_pa_lib_ambiguity_codes_consensus(String datum) {
		if(Integer.parseInt(datum) != 0)
			this.libAmbiguityCodesConsensus = true;
		else 
			this.libAmbiguityCodesConsensus = false;
	}

	/**
	 * PRIMER_LIBERAL_BASE
	 * @param datum
	 */
	public void p3_set_pa_liberal_base(String datum) {

		if(Integer.parseInt(datum) != 0)
			this.liberalBase = true;
		else
			this.liberalBase = false;
	}

	/**
	 * PRIMER_LOWERCASE_MASKING
	 */
	public void p3_set_pa_lowercase_masking(boolean value) {
		this.lowercaseMasking = value;
	}

	/**
	 * PRIMER_LOWERCASE_MASKING
	 */
	public void p3_set_pa_lowercase_masking(String datum) {
		if(Integer.parseInt(datum)==0)
			this.lowercaseMasking = false;
		else
			this.lowercaseMasking = true;
	}

	/**
	 * 
	 * @param datum
	 */
	public boolean p3_set_pa_mask_template(String datum) {
		if(Integer.parseInt(datum)==0)
			this.maskTemplate = false;
		else
			this.maskTemplate = true;

		return this.maskTemplate;
	}

	/**
	 * PRIMER_PAIR_MAX_DIFF_TM
	 * @param datum
	 */
	public void p3_set_pa_max_diff_tm(String datum) {
		this.maxDiffTm = Double.parseDouble(datum);		
	}

	/**
	 * PRIMER_MAX_END_GC
	 * @param datum
	 */
	public void p3_set_pa_max_end_gc(String datum) {
		this.setMaxEndGC(Integer.parseInt(datum));		
	}

	/**
	 * PRIMER_MAX_END_STABILITY
	 * @param datum
	 */
	public void p3_set_pa_max_end_stability(String datum) {
		this.maxEndStability = Double.parseDouble(datum);
	}

	public void p3_set_pa_min_3_prime_overlap_of_junction(String datum) {
		this.min3PrimeOverlapOfJunction = Integer.parseInt(datum);
	}

	/**
	 * PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION
	 * @param datum
	 */
	public void p3_set_pa_min_5_prime_overlap_of_junction(String datum) {
		this.min5PrimeOverlapOfJunction = Integer.parseInt(datum);
	}

	public void p3_set_pa_min_left_three_prime_distance(String  datum) {
		this.minLeft3PrimeDistance = Integer.parseInt(datum);
	}

	public void p3_set_pa_min_right_three_prime_distance(String datum) {
		this.minRight3PrimeDistance =Integer.parseInt(datum);
	}

	/**
	 * PRIMER_NUM_RETURN
	 * @param datum
	 */
	public void p3_set_pa_num_return(String datum) {
		this.numReturn = Integer.parseInt(datum);
	}

	/**
	 * PRIMER_OUTSIDE_PENALTY
	 * @param datum
	 */
	public void p3_set_pa_outside_penalty(String datum) {
		this.outsidePenalty = Double.parseDouble(datum);
	}

	/**
	 * PRIMER_PAIR_MAX_COMPL_ANY
	 * @param datum
	 */
	public void p3_set_pa_pair_compl_any(String datum) {
		this.pairComplAny = Double.parseDouble(datum);
	}

	/**
	 * PRIMER_PAIR_MAX_COMPL_ANY_TH
	 * @param datum
	 */
	public void p3_set_pa_pair_compl_any_th(String datum) {
		this.pairComplAnyTH = Double.parseDouble(datum);
	}

	/**
	 * PRIMER_PAIR_MAX_COMPL_END
	 * @param datum
	 */
	public void p3_set_pa_pair_compl_end(String datum) {
		this.pairComplEnd = Double.parseDouble(datum);
	}

	/**
	 * PRIMER_PAIR_MAX_COMPL_END_TH
	 * @param datum
	 */
	public void p3_set_pa_pair_compl_end_th(String datum) {
		this.pairComplEndTH = Double.parseDouble(datum);
	}

	/**
	 * PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING
	 * @param datum
	 */

	public void p3_set_pa_pair_max_template_mispriming(String datum) {
		this.pairMaxTemplateMispriming = Double.parseDouble(datum);		
	}

	/**
	 * PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH
	 * @param datum
	 */

	public void p3_set_pa_pair_max_template_mispriming_th(String datum) {
		this.pairMaxTemplateMisprimingTH = Double.parseDouble(datum);
	}

	/**
	 * PRIMER_PAIR_MAX_LIBRARY_MISPRIMING
	 * @param datum
	 */
	public void p3_set_pa_pair_repeat_compl(String datum) {
		this.pairRepeatCompl = Double.parseDouble(datum);		
	}

	/**
	 * PRIMER_PICK_ANYWAY
	 * @param datum
	 */
	public void p3_set_pa_pick_anyway(String datum) {
		if(Integer.parseInt(datum) == 0 )
			this.pickAnyway = false;
		else
			this.pickAnyway = true;
	}

	/**
	 * PRIMER_PICK_INTERNAL_OLIGO
	 * @param datum
	 */
	public void p3_set_pa_pick_internal_oligo(String datum) {
		if(Integer.parseInt(datum) == 0)
			this.setPickInternalOligo(false);
		else
			this.setPickInternalOligo(true);

	}

	/**
	 * PRIMER_PICK_LEFT_PRIMER
	 * @param datum
	 */
	public void p3_set_pa_pick_left_primer(String datum) {
		if(Integer.parseInt(datum) == 0)
			this.setPickLeftPrimer(false);
		else
			this.setPickLeftPrimer(true);		
	}

	/**
	 * PRIMER_PICK_RIGHT_PRIMER
	 * @param datum
	 */
	public void p3_set_pa_pick_right_primer(String datum) {
		if(Integer.parseInt(datum) == 0)
			this.setPickRightPrimer(false);
		else
			this.setPickRightPrimer(true);

	}

	/**
	 * PRIMER_PRODUCT_MAX_TM
	 * @param datum
	 */
	public void p3_set_pa_product_max_tm(String datum) {
		this.productMaxTM = Integer.parseInt(datum);		
	}

	/**
	 * PRIMER_PRODUCT_MIN_TM
	 * @param datum
	 */
	public void p3_set_pa_product_min_tm(String datum) {
		this.productMinTM = Integer.parseInt(datum);
	}

	/**
	 * PRIMER_PRODUCT_OPT_SIZE
	 * @param datum
	 */
	public void p3_set_pa_product_opt_size(String datum) {
		this.productOptSize = Integer.parseInt(datum);
	}

	/**
	 * PRIMER_PRODUCT_OPT_TM
	 * @param datum
	 */
	public void p3_set_pa_product_opt_tm(String datum) {
		this.productOptTM = Integer.parseInt(datum);
	}

	/**
	 * PRIMER_PRODUCT_SIZE_RANGE
	 * @param datum
	 */
	public void p3_set_pa_product_size(String datum)
	{
		// TODO :: catch ex and clear values then throws it again
		String intervalSep = " ";
		String numSep = "-";
		String[] intervalStrs = datum.split(Pattern.quote(intervalSep));
		//		this.num_intervals=0;
		this.productSizeRanges.clear();
		for(String intervalStr : intervalStrs){
			//			if(this.num_intervals >=  LibPrimer3.PR_MAX_INTERVAL_ARRAY)
			//				throw new IndexOutOfBoundsException("Too many elements for tag ");

			String[] intervalNums = intervalStr.split(Pattern.quote(numSep));
			//			this.pr_min[this.num_intervals] = Integer.parseInt(intervalNums[0]);
			int pr_min = Integer.parseInt(intervalNums[0]);
			//			this.pr_max[this.num_intervals] = Integer.parseInt(intervalNums[1]);
			int pr_max = Integer.parseInt(intervalNums[1]);

			this.addProductSizeRange(pr_min, pr_max);
			//			this.num_intervals++;
		}		
	}

	/**
	 * PRIMER_QUALITY_RANGE_MAX
	 * @param datum
	 */
	public void p3_set_pa_quality_range_max(String datum) {
		this.qualityRangeMax = Integer.parseInt(datum);
	}

	public void p3_set_pa_quality_range_min(String datum) {
		this.qualityRangeMin = Integer.parseInt(datum);
	}

	/**
	 * PRIMER_SALT_CORRECTIONS
	 * @param datum
	 */
	public void p3_set_pa_salt_corrections(String datum) {
		int typeInt = Integer.parseInt(datum);

		if(typeInt == SaltCorrectionMethod.owczarzy.getValue() )
		{
			this.setSaltCorrectionMethod(SaltCorrectionMethod.owczarzy);
		} 
		else if(typeInt == SaltCorrectionMethod.santalucia.getValue() ) 
		{
			this.setSaltCorrectionMethod(SaltCorrectionMethod.santalucia);
		}
		else if(typeInt == SaltCorrectionMethod.schildkraut.getValue() ) 
		{
			this.setSaltCorrectionMethod(SaltCorrectionMethod.schildkraut);
		}
		else {
			// FIXME :: what to do in this case
			// TODO :: append a warnning to .... and used default value
		}
	}

	/**
	 * PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT
	 */
	public void p3_set_pa_thermodynamic_oligo_alignment(String datum) {
		if( Integer.parseInt(datum) == 1 )
			this.thermodynamicOligoAlignment = true ;
		else
			this.thermodynamicOligoAlignment = false ;
	}

	/**
	 * PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT
	 * @param datum
	 */

	public void p3_set_pa_thermodynamic_template_alignment(String datum) {
		if(Integer.parseInt(datum) == 1 )
			this.thermodynamicTemplateAlignment = true ;
		else 
			this.thermodynamicTemplateAlignment = false ;
	}

	/**
	 * PRIMER_TM_FORMULA
	 */
	public void p3_set_pa_tm_method_type(String datum) {
		int typeInt = Integer.parseInt(datum);
		if(typeInt == MeltingTemperatureMethod.santalucia_auto.getValue() )
		{
			this.mtMethod = MeltingTemperatureMethod.santalucia_auto;
		}
		else if(typeInt == MeltingTemperatureMethod.breslauer_auto.getValue() )
		{
			this.mtMethod = MeltingTemperatureMethod.breslauer_auto;
		}

	}

	/**
	 * Write the default values for default_values=1
	 */
	private void pr_set_default_global_args_1() {
		/* Arguments for primers ================================= */
		this.primersArgs.setOptSize(20);
		this.primersArgs.setMinSize(18);
		this.primersArgs.setMaxSize(27);

		this.primersArgs.setOptTm(60);
		this.primersArgs.setMinTm(57);
		this.primersArgs.setMaxTm(63);

		this.primersArgs.setMinGC(20.0);
		this.primersArgs.setOptGC(LibPrimer3.DEFAULT_OPT_GC_PERCENT);
		this.primersArgs.setMaxGC(80.0);
		this.primersArgs.setSaltConcentration(50.0);
		this.primersArgs.setDivalentConcentration(0.0);
		this.primersArgs.setDntpConcentration(0.0);

		this.primersArgs.setDnaConcentration(50.0);
		this.primersArgs.setMaxNumOfNsAccepted(0);
		this.primersArgs.setMaxSelfAny(8.0);
		this.primersArgs.setMaxSelfEnd(3.0);
		this.primersArgs.setMaxSelfAnyTH(47.0);
		this.primersArgs.setMaxSelfEndTH(47.0);
		this.primersArgs.setMaxHairPinTH(47.0);
		this.primersArgs.setMaxPolyX(5);
		this.primersArgs.setMaxRepeatCompl(12.0);
		this.primersArgs.setMinQuality(0);
		this.primersArgs.setMinEndQuality(0);
		this.primersArgs.setMaxTemplateMispriming(LibPrimer3.PR_UNDEFINED_ALIGN_OPT);
		this.primersArgs.setMaxTemplateMisprimingTH(LibPrimer3.PR_UNDEFINED_ALIGN_OPT);
		/* The following apply only to primers (and not to internal
		     oligos). */
		this.setGcClamp(0);
		this.setMaxEndGC(5);

		/* Weights for objective functions for oligos and pairs. */
		this.primersArgs.weights.compl_any     = 0;
		this.primersArgs.weights.compl_any_th  = 0;
		this.primersArgs.weights.compl_end     = 0;
		this.primersArgs.weights.compl_end_th  = 0;
		this.primersArgs.weights.end_quality   = 0;
		this.primersArgs.weights.end_stability = 0;
		this.primersArgs.weights.gc_content_gt = 0;
		this.primersArgs.weights.gc_content_lt = 0;
		this.primersArgs.weights.hairpin_th    = 0;
		this.primersArgs.weights.length_gt     = 1;
		this.primersArgs.weights.length_lt     = 1;
		this.primersArgs.weights.num_ns        = 0;
		this.primersArgs.weights.pos_penalty   = 1;
		this.primersArgs.weights.repeat_sim    = 0;
		this.primersArgs.weights.seq_quality   = 0;
		this.primersArgs.weights.temp_cutoff   = 5;
		this.primersArgs.weights.temp_gt       = 1;
		this.primersArgs.weights.temp_lt       = 1;
		this.primersArgs.weights.template_mispriming = 0.0;
		this.primersArgs.weights.template_mispriming_th = 0.0;
		this.primersArgs.must_match_five_prime  = null;
		this.primersArgs.must_match_three_prime = null;
		/* End of weights for objective functions for oligos and pairs. */

		/* End of arguments for primers =========================== */

		this.maxDiffTm         = 100.0;
		this.mtMethod       = MeltingTemperatureMethod.breslauer_auto;
		this.setSaltCorrectionMethod(SaltCorrectionMethod.schildkraut);
		this.pairComplAny      = 8.0;
		this.pairComplEnd      = 3.0;
		this.pairComplAnyTH   = 47.0;
		this.pairComplEndTH   = 47.0;
		this.thermodynamicOligoAlignment = false;
		this.thermodynamicTemplateAlignment = false;
		this.liberalBase        = false;
		this.setPrimerTask(P3Task.GENERIC);
		this.setPickLeftPrimer(true);
		this.setPickRightPrimer(true);
		this.setPickInternalOligo(false);
		this.firstBaseIndex    = 0;
		this.numReturn          = 5;
		this.addProductSizeRange(100, 300);
		//		this.pr_min[0]           = 100;
		//		this.pr_max[0]           = 300;
		//		this.num_intervals       = 1;
		this.pairRepeatCompl   = 24.0;
		this.qualityRangeMin   = 0;
		this.qualityRangeMax   = 100;
		this.outsidePenalty     = LibPrimer3.PR_DEFAULT_OUTSIDE_PENALTY;
		this.insidePenalty      = LibPrimer3.PR_DEFAULT_INSIDE_PENALTY;
		this.maxEndStability   = 100.0;
		this.lowercaseMasking   = false;
		this.productMaxTM      = LibPrimer3.PR_DEFAULT_PRODUCT_MAX_TM;
		this.productMinTM      = LibPrimer3.PR_DEFAULT_PRODUCT_MIN_TM;
		this.productOptTM      = LibPrimer3.PR_UNDEFINED_DBL_OPT;
		this.productOptSize    = LibPrimer3.PR_UNDEFINED_INT_OPT;
		this.pairMaxTemplateMispriming = LibPrimer3.PR_UNDEFINED_ALIGN_OPT;
		this.pairMaxTemplateMisprimingTH = LibPrimer3.PR_UNDEFINED_ALIGN_OPT;
		this.oligosArgs.setOptSize(20);
		this.oligosArgs.setMinSize(18);
		this.oligosArgs.setMaxSize(27);
		this.oligosArgs.setOptTm(60.0);
		this.oligosArgs.setMinTm(57.0);
		this.oligosArgs.setMaxTm(63.0);
		this.oligosArgs.setMinGC(20.0);
		this.oligosArgs.setMaxGC(80.0);
		this.oligosArgs.setOptGC(LibPrimer3.DEFAULT_OPT_GC_PERCENT);
		this.oligosArgs.setMaxPolyX(5);
		this.oligosArgs.setSaltConcentration(50.0);
		this.oligosArgs.setDivalentConcentration(0.0);
		this.oligosArgs.setDntpConcentration(0.0);
		this.oligosArgs.setDnaConcentration(50.0);
		this.oligosArgs.setMaxNumOfNsAccepted(0);
		this.oligosArgs.setMaxSelfAny(12.0);
		this.oligosArgs.setMaxSelfEnd(12.0);
		this.oligosArgs.setMaxSelfAnyTH(47.0);
		this.oligosArgs.setMaxSelfEndTH(47.0);
		this.oligosArgs.setMaxHairPinTH(47.0);
		this.oligosArgs.setMaxRepeatCompl(12.0);

		this.oligosArgs.setMinQuality(0);
		this.oligosArgs.setMinEndQuality(0);
		this.oligosArgs.setMaxTemplateMispriming(LibPrimer3.PR_UNDEFINED_ALIGN_OPT);
		this.oligosArgs.setMaxTemplateMisprimingTH(LibPrimer3.PR_UNDEFINED_ALIGN_OPT);
		this.oligosArgs.weights.temp_gt       = 1;
		this.oligosArgs.weights.temp_lt       = 1;
		this.oligosArgs.weights.length_gt     = 1;
		this.oligosArgs.weights.length_lt     = 1;
		this.oligosArgs.weights.gc_content_gt = 0;
		this.oligosArgs.weights.gc_content_lt = 0;
		this.oligosArgs.weights.compl_any     = 0;
		this.oligosArgs.weights.compl_end     = 0;
		this.oligosArgs.weights.compl_any_th  = 0;
		this.oligosArgs.weights.compl_end_th  = 0;
		this.oligosArgs.weights.hairpin_th    = 0;
		this.oligosArgs.weights.num_ns        = 0;
		this.oligosArgs.weights.repeat_sim    = 0;
		this.oligosArgs.weights.seq_quality   = 0;
		this.oligosArgs.weights.end_quality   = 0;
		this.oligosArgs.must_match_five_prime  = null;
		this.oligosArgs.must_match_three_prime = null;

		this.prPairWeights.primer_quality  = 1;
		this.prPairWeights.io_quality      = 0;
		this.prPairWeights.diff_tm         = 0;
		this.prPairWeights.compl_any       = 0;
		this.prPairWeights.compl_end       = 0;
		this.prPairWeights.compl_any_th    = 0;
		this.prPairWeights.compl_end_th    = 0;
		this.prPairWeights.temp_cutoff     = 5;
		this.prPairWeights.repeat_sim      = 0;

		this.prPairWeights.product_tm_lt   = 0;
		this.prPairWeights.product_tm_gt   = 0;
		this.prPairWeights.product_size_lt = 0;
		this.prPairWeights.product_size_gt = 0;

		this.libAmbiguityCodesConsensus   = true;
		/*  Set to 1 for backward compatibility. This _NOT_ what
		      one normally wants, since many libraries contain
		      strings of N, which then match every oligo (very bad).
		 */

		this.minLeft3PrimeDistance   = -1;
		this.minRight3PrimeDistance  = -1;

		this.sequencingParameters.setLead(50);
		this.sequencingParameters.setSpacing(500);
		this.sequencingParameters.setInterval(250);
		this.sequencingParameters.setAccuracy(20);

		this.min5PrimeOverlapOfJunction = 7;
		this.min3PrimeOverlapOfJunction = 4;

		this.maskTemplate                   = false;
		this.maskingParametersChanged      = false;
		this.getMaskingParameters().mdir                         = masking_direction.both_separately;
		this.getMaskingParameters().failure_rate                 = 0.1;
		this.getMaskingParameters().nucl_masked_in_5p_direction  = 1;
		this.getMaskingParameters().nucl_masked_in_3p_direction  = 0;
		this.getMaskingParameters().print_sequence               = false;
		this.getMaskingParameters().do_soft_masking              = true;
		this.getMaskingParameters().nlists                       = masker.DEFAULT_NLISTS;
		this.getMaskingParameters().list_prefix                  = masker.DEFAULT_LIST_FILE_PREFIX;
		this.getMaskingParameters().fp                           = null;
		this.getMaskingParameters().formula_intercept            = masker.DEFAULT_INTERCEPT;  		
	}

	/*
	 * Write the default values for default_values=2
	 */
	private void pr_set_default_global_args_2() {
		this.pr_set_default_global_args_1();
		this.mtMethod                    = MeltingTemperatureMethod.santalucia_auto;
		this.setSaltCorrectionMethod(SaltCorrectionMethod.santalucia);
		this.thermodynamicOligoAlignment    = true;
		this.thermodynamicTemplateAlignment = false;
		this.primersArgs.setDivalentConcentration(1.5);
		this.primersArgs.setDntpConcentration(0.6);
		this.libAmbiguityCodesConsensus    = false;
	}

	/**
	 * PRIMER_MASK_KMERLIST_PREFIX
	 */
	public void set_masking_parameters_KmerLis_Prefix(String datum) {
		this.getMaskingParameters().set_list_prefix(datum);
		this.maskingParametersChanged = true;
	}

	

	/**
	 * @param dump the dump to set
	 */
	public void setDump(boolean dump) {
		this.dump = dump;
	}

	/**
	 * @param firstBaseIndex the firstBaseIndex to set
	 */
	public void setFirstBaseIndex(int firstBaseIndex) {
		this.firstBaseIndex = firstBaseIndex;
	}

	/**
	 * @param gcClamp the gcClamp to set
	 */
	public void setGcClamp(int gcClamp) {
		this.gcClamp = gcClamp;
	}

	/**
	 * @param insidePenalty the insidePenalty to set
	 */
	public void setInsidePenalty(double insidePenalty) {
		this.insidePenalty = insidePenalty;
	}

	/**
	 * @param libAmbiguityCodesConsensus the libAmbiguityCodesConsensus to set
	 */
	public void setLibAmbiguityCodesConsensus(boolean libAmbiguityCodesConsensus) {
		this.libAmbiguityCodesConsensus = libAmbiguityCodesConsensus;
	}

	/**
	 * @param liberalBase the liberalBase to set
	 */
	public void setLiberalBase(boolean liberalBase) {
		this.liberalBase = liberalBase;
	}

	/**
	 * @param lowercaseMasking the lowercaseMasking to set
	 */
	public void setLowercaseMasking(boolean lowercaseMasking) {
		this.lowercaseMasking = lowercaseMasking;
	}

	/**
	 * @param maskingParameters the maskingParameters to set
	 */
	public void setMaskingParameters(MaskerParameters maskingParameters) {
		this.maskingParameters = maskingParameters;
	}

	/**
	 * @param maskingParametersChanged the maskingParametersChanged to set
	 */
	public void setMaskingParametersChanged(boolean maskingParametersChanged) {
		this.maskingParametersChanged = maskingParametersChanged;
	}

	/**
	 * @param maskTemplate the maskTemplate to set
	 */
	public void setMaskTemplate(boolean maskTemplate) {
		this.maskTemplate = maskTemplate;
	}

	/**
	 * @param maxDiffTm the maxDiffTm to set
	 */
	public void setMaxDiffTm(double maxDiffTm) {
		this.maxDiffTm = maxDiffTm;
	}

	/**
	 * @param maxEndGC the maxEndGC to set
	 */
	public void setMaxEndGC(int maxEndGC) {
		this.maxEndGC = maxEndGC;
	}

	/**
	 * @param maxEndStability the maxEndStability to set
	 */
	public void setMaxEndStability(double maxEndStability) {
		this.maxEndStability = maxEndStability;
	}

	/**
	 * @param min3PrimeOverlapOfJunction the min3PrimeOverlapOfJunction to set
	 */
	public void setMin3PrimeOverlapOfJunction(int min3PrimeOverlapOfJunction) {
		this.min3PrimeOverlapOfJunction = min3PrimeOverlapOfJunction;
	}

	/**
	 * @param min5PrimeOverlapOfJunction the min5PrimeOverlapOfJunction to set
	 */
	public void setMin5PrimeOverlapOfJunction(int min5PrimeOverlapOfJunction) {
		this.min5PrimeOverlapOfJunction = min5PrimeOverlapOfJunction;
	}

	/**
	 * @param minLeft3PrimeDistance the minLeft3PrimeDistance to set
	 */
	public void setMinLeft3PrimeDistance(int minLeft3PrimeDistance) {
		this.minLeft3PrimeDistance = minLeft3PrimeDistance;
	}

	/**
	 * @param minRight3PrimeDistance the minRight3PrimeDistance to set
	 */
	public void setMinRight3PrimeDistance(int minRight3PrimeDistance) {
		this.minRight3PrimeDistance = minRight3PrimeDistance;
	}

	/**
	 * @param numReturn the numReturn to set
	 */
	public void setNumReturn(int numReturn) {
		this.numReturn = numReturn;
	}

	/**
	 * @param outsidePenalty the outsidePenalty to set
	 */
	public void setOutsidePenalty(double outsidePenalty) {
		this.outsidePenalty = outsidePenalty;
	}

	/**
	 * @param pairComplAny the pairComplAny to set
	 */
	public void setPairComplAny(double pairComplAny) {
		this.pairComplAny = pairComplAny;
	}

	/**
	 * @param pairComplAnyTH the pairComplAnyTH to set
	 */
	public void setPairComplAnyTH(double pairComplAnyTH) {
		this.pairComplAnyTH = pairComplAnyTH;
	}

	/**
	 * @param pairComplEnd the pairComplEnd to set
	 */
	public void setPairComplEnd(double pairComplEnd) {
		this.pairComplEnd = pairComplEnd;
	}

	/**
	 * @param pairComplEndTH the pairComplEndTH to set
	 */
	public void setPairComplEndTH(double pairComplEndTH) {
		this.pairComplEndTH = pairComplEndTH;
	}

	/**
	 * @param pairMaxTemplateMispriming the pairMaxTemplateMispriming to set
	 */
	public void setPairMaxTemplateMispriming(double pairMaxTemplateMispriming) {
		this.pairMaxTemplateMispriming = pairMaxTemplateMispriming;
	}

	/**
	 * @param pairMaxTemplateMisprimingTH the pairMaxTemplateMisprimingTH to set
	 */
	public void setPairMaxTemplateMisprimingTH(double pairMaxTemplateMisprimingTH) {
		this.pairMaxTemplateMisprimingTH = pairMaxTemplateMisprimingTH;
	}

	/**
	 * @param pairRepeatCompl the pairRepeatCompl to set
	 */
	public void setPairRepeatCompl(double pairRepeatCompl) {
		this.pairRepeatCompl = pairRepeatCompl;
	}

	/**
	 * @param pickAnyway the pickAnyway to set
	 */
	public void setPickAnyway(boolean pickAnyway) {
		this.pickAnyway = pickAnyway;
	}

	/**
	 * @param pickInternalOligo the pickInternalOligo to set
	 */
	public void setPickInternalOligo(boolean pickInternalOligo) {
		this.pickInternalOligo = pickInternalOligo;
	}

	/**
	 * @param pickLeftPrimer the pickLeftPrimer to set
	 */
	public void setPickLeftPrimer(boolean pickLeftPrimer) {
		this.pickLeftPrimer = pickLeftPrimer;
	}

	/**
	 * @param pickRightPrimer the pickRightPrimer to set
	 */
	public void setPickRightPrimer(boolean pickRightPrimer) {
		this.pickRightPrimer = pickRightPrimer;
	}

	/**
	 * @param primerTask the primerTask to set
	 */
	public void setPrimerTask(P3Task primerTask) {
		this.primerTask = primerTask;
	}

	/**
	 * @param productMaxTM the productMaxTM to set
	 */
	public void setProductMaxTM(double productMaxTM) {
		this.productMaxTM = productMaxTM;
	}

	/**
	 * @param productMinTM the productMinTM to set
	 */
	public void setProductMinTM(double productMinTM) {
		this.productMinTM = productMinTM;
	}

	/**
	 * @param productOptSize the productOptSize to set
	 */
	public void setProductOptSize(int productOptSize) {
		this.productOptSize = productOptSize;
	}

	/**
	 * @param productOptTM the productOptTM to set
	 */
	public void setProductOptTM(double productOptTM) {
		this.productOptTM = productOptTM;
	}

	/**
	 * @param prPairWeights the prPairWeights to set
	 */
	public void setPrPairWeights(PairWeights prPairWeights) {
		this.prPairWeights = prPairWeights;
	}

	/**
	 * @param qualityRangeMax the qualityRangeMax to set
	 */
	public void setQualityRangeMax(int qualityRangeMax) {
		this.qualityRangeMax = qualityRangeMax;
	}

	/**
	 * @param qualityRangeMin the qualityRangeMin to set
	 */
	public void setQualityRangeMin(int qualityRangeMin) {
		this.qualityRangeMin = qualityRangeMin;
	}

	/**
	 * @param saltCorrectionMethod the saltCorrectionMethod to set
	 */
	public void setSaltCorrectionMethod(SaltCorrectionMethod saltCorrectionMethod) {
		this.saltCorrectionMethod = saltCorrectionMethod;
	}

	/**
	 * @param sequencingParameters the sequencingParameters to set
	 */
	public void setSequencingParameters(SequencingParameters sequencingParameters) {
		this.sequencingParameters = sequencingParameters;
	}

	/**
	 * @param thermodynamicOligoAlignment the thermodynamicOligoAlignment to set
	 */
	public void setThermodynamicOligoAlignment(boolean thermodynamicOligoAlignment) {
		this.thermodynamicOligoAlignment = thermodynamicOligoAlignment;
	}

	/**
	 * @param thermodynamicTemplateAlignment the thermodynamicTemplateAlignment to set
	 */
	public void setThermodynamicTemplateAlignment(
			boolean thermodynamicTemplateAlignment) {
		this.thermodynamicTemplateAlignment = thermodynamicTemplateAlignment;
	}
}