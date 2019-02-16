package org.primer3.libprimer3;

import java.util.Formatter;


public class OligoStats {
	public int considered;          /* Total number of tested oligos of given type   */
	public int ns;                  /* Number of oligos rejected because of Ns       */
	public int target;              /* Overlapping targets.                          */
	public int excluded;            /* Overlapping excluded regions.                 */
	public int gc;                  /* Unacceptable GC content.                      */
	public int gc_clamp;            /* Don't have required number of GCs at 3' end.  */
	public int gc_end_high;         /* Too many G+Cs at the 3' end.                  */
	public int temp_min;            /* Melting temperature below t_min.              */
	public int temp_max;            /* Melting temperature more than t_max.          */
	public int size_min;            /* Primer shorter than minimal size.             */
	public int size_max;            /* Primer longer than minimal size.              */
	public int compl_any;           /* Self-complementarity too high.                */
	public int compl_end;           /* Self-complementarity at 3' end too high.      */
	public int hairpin_th;          /* Hairpin structure too stable in
			      thermodynamical approach                      */
	public int repeat_score;        /* Complementarity with repeat sequence too high.*/
	public int poly_x;              /* Long mononucleotide sequence inside.          */
	public int seq_quality;         /* Low quality of bases included.                */
	public int stability;           /* Stability of 5 3' bases too high.             */
	public int no_orf;              /* Would not amplify any of the specified ORF
                             (valid for left primers only).                 */
	public int template_mispriming; /* Template mispriming score too high.           */
	public int ok;                  /* Number of acceptable oligos.                  */
	public int gmasked;             /* Added by T. Koressaar, number of gmasked
			      oligos */
	public int must_match_fail;     /* Added by A. Untergasser, number of oligos 
			      failing must match */
	public int not_in_any_left_ok_region; /* Oligo not included in any of the
                                    left regions given in
                                    PRIMER_PAIR_OK_REGION_LIST. */
	public int not_in_any_right_ok_region;/* Oligo not included in any of the
                                    right regions given in
                                    PRIMER_PAIR_OK_REGION_LIST. */
	public String p3_oligo_explain_string() {
		StringBuilder  sbuf = new StringBuilder();
		Formatter sprintf = new Formatter(sbuf);
		sprintf.format("considered %d", considered);
		
		
//		sprintf.format("considered %d", this.considered);
		
		
		if(this.no_orf > 0)
			sprintf.format(", would not amplify any of the ORF %d", this.no_orf);
		if(this.ns != 0)
			sprintf.format(", too many Ns %d", this.ns);
		if(this.target != 0)
			sprintf.format(", overlap target %d", this.target);
		if(this.excluded != 0)
			sprintf.format(", overlap excluded region %d", this.excluded);
		if(this.gc != 0)
			sprintf.format(", GC content failed %d", this.gc);
		if(this.gc_clamp != 0)
			sprintf.format(", GC clamp failed %d", this.gc_clamp);
		if(this.temp_min != 0)
			sprintf.format(", low tm %d", this.temp_min);
		if( this.temp_max != 0)
			sprintf.format(", high tm %d", this.temp_max);
		if(this.compl_any != 0)
			sprintf.format(", high any compl %d", this.compl_any);
		if(this.compl_end != 0)
			sprintf.format(", high end compl %d", this.compl_end);
		if(this.hairpin_th != 0)
			sprintf.format(", high hairpin stability %d", this.hairpin_th);
		if(this.repeat_score != 0)
			sprintf.format(", high repeat similarity %d", this.repeat_score);
		if(this.poly_x != 0)
			sprintf.format(", long poly-x seq %d", this.poly_x);
		if(this.seq_quality != 0)
			sprintf.format(", low sequence quality %d", this.seq_quality);
		if(this.stability != 0)
			sprintf.format(", high 3' stability %d", this.stability);
		if(this.template_mispriming != 0)
			sprintf.format(", high template mispriming score %d", this.template_mispriming);
		if(this.gmasked != 0)
			sprintf.format(", lowercase masking of 3' end %d",this.gmasked);
		if(this.must_match_fail != 0)
			sprintf.format(", failed must_match requirements %d",this.must_match_fail);
		if(this.not_in_any_left_ok_region != 0)
			sprintf.format(", not in any ok left region %d",  this.not_in_any_left_ok_region);
		if(this.not_in_any_right_ok_region  != 0)
			sprintf.format(", not in any ok right region %d",this.not_in_any_right_ok_region);
		
		
		sprintf.format(", ok %d", this.ok);
		
		
		
		sprintf.close();
		return sbuf.toString();
	}
}