package com.primer3.oligotm;

/**
 * Melting Temperature Method 
 * For olgigotm() and seqtm()
 * Both functions return the melting temperature of the given oligo
 * calculated as specified by user, but oligotm _should_ only be used on
 * DNA sequences of length <= MAX_PRIMER_LENGTH (which is defined
 * elsewhere).  seqtm uses oligotm for sequences of length <=
 * MAX_PRIMER_LENGTH, and a different, G+C% based formula for longer
 * sequences.  For oligotm(), no error is generated on sequences
 * longer than MAX_PRIMER_LENGTH, but the formula becomes less
 * accurate as the sequence grows longer.  Caveat emptor.
 * 
 * We use the folowing typedefs:
 * If tm_method==santalucia_auto, then the table of
 * nearest-neighbor thermodynamic parameters and method for Tm
 * calculation in the paper [SantaLucia JR (1998) "A unified view of
 * polymer, dumbbell and oligonucleotide DNA nearest-neighbor
 * thermodynamics", Proc Natl Acad Sci 95:1460-65
 * http://dx.doi.org/10.1073/pnas.95.4.1460] is used.
 * THIS IS THE RECOMMENDED VALUE*.
 * If tm_method==breslauer_auto, then method for Tm
 * calculations in the paper [Rychlik W, Spencer WJ and Rhoads RE
 * (1990) "Optimization of the annealing temperature for DNA
 * amplification in vitro", Nucleic Acids Res 18:6409-12
 * http://www.pubmedcentral.nih.gov/articlerender.fcgi?tool=pubmed&pubmedid=2243783].
 * and the thermodynamic parameters in the paper [Breslauer KJ, Frank
 * R, Blöcker H and Marky LA (1986) "Predicting DNA duplex stability
 * from the base sequence" Proc Natl Acad Sci 83:4746-50
 * http://dx.doi.org/10.1073/pnas.83.11.3746], are is used.  This is
 * the method and the table that primer3 used up to and including
 * version 1.0.1
 */
public enum MeltingTemperatureMethod {
	breslauer_auto(0),
	santalucia_auto (1);

	private int id ;
	MeltingTemperatureMethod(int id ){
		this.id = id;
	}
	public int getValue() { return id; }



}
