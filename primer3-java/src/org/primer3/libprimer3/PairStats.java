package org.primer3.libprimer3;

import java.util.Formatter;

public class PairStats {
	public int considered;          /* Total number of pairs or triples tested.      */
	public int product;             /* Pairs providing incorrect product size.       */
	public int target;              /* Pairs without any target between primers.     */
	public int temp_diff;           /* Melting temperature difference too high.      */
	public int compl_any;           /* Pairwise complementarity larger than allowed. */
	public int compl_end;           /* The same for 3' end complementarity.          */
	public int internal;            /* Internal oligo was not found.                 */
	public int repeat_sim;          /* Complementarity with repeat sequence too high.*/
	public int high_tm;             /* Product Tm too high.                          */
	public int low_tm;              /* Product Tm too low.                           */
	public int template_mispriming; /* Sum of template mispriming scores too high.   */

  /* Neither oligo in the pairs overlaps one of the "required sites".       */
	public int does_not_overlap_a_required_point;

  /* One of the oligos in the pair overlaps an oligo in a better_pair:       */
  public int overlaps_oligo_in_better_pair;

  /* The left and right oligos are not in any of the pair of regions given in
     PRIMER_PAIR_OK_REGION_LIST. */
  public int not_in_any_ok_region;

  /* Left primer to the right of right right primer. This can occur when
     the primers are provided by the caller. */
  public int reversed;

  public int ok;                  /* Number that were ok.                          */

  
  
  
  
  
  
	public String p3_pair_explain_string() {
		StringBuilder  sbuf = new StringBuilder();
		Formatter sprintf = new Formatter(sbuf);
		sprintf.format("considered %d", considered);
		
		
		if(this.target != 0)
			sprintf.format(", no target %d", this.target);
		if(this.product != 0)
			sprintf.format(", unacceptable product size %d", this.product);
		if(this.low_tm != 0)
			sprintf.format(", low product Tm %d", this.low_tm);
		if(this.high_tm != 0)
			sprintf.format(", high product Tm %d", this.high_tm);
		if(this.temp_diff != 0)
			sprintf.format(", tm diff too large %d",this.temp_diff);
		if(this.compl_any != 0)
			sprintf.format(", high any compl %d", this.compl_any);
		if(this.compl_end != 0)
			sprintf.format(", high end compl %d", this.compl_end);
		if(this.internal != 0)
			sprintf.format(", no internal oligo %d", this.internal);
		if( this.repeat_sim != 0)
			sprintf.format(", high mispriming library similarity %d",
		                  this.repeat_sim);
		if( this.does_not_overlap_a_required_point != 0)
			sprintf.format(", no overlap of required point %d",
		                  this.does_not_overlap_a_required_point);
		if(this.overlaps_oligo_in_better_pair != 0)
			sprintf.format(", primer in pair overlaps a primer in a better pair %d",
		                  this.overlaps_oligo_in_better_pair);
		if( this.template_mispriming != 0)
			sprintf.format(", high template mispriming score %d",
		                  this.template_mispriming);
		if( this.not_in_any_ok_region != 0)
			sprintf.format(", not in any ok region %d", 
		                  this.not_in_any_ok_region);
		if(this.reversed != 0)
			sprintf.format(", left primer to right of right primer %d",
				  this.reversed);

		sprintf.format(", ok %d", this.ok);
		
		
		
		sprintf.close();
		return sbuf.toString();
	} 
  
  
  
  
}