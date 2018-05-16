package com.primer3.libprimer3;

public class PairWeights {
	
	double primer_quality;
	double io_quality;
	double diff_tm;
	double compl_any;
	double compl_any_th;
	double compl_end;
	double compl_end_th;
	double temp_cutoff;
	double product_tm_lt;
	double product_tm_gt;
	double product_size_lt;
	double product_size_gt;
	double repeat_sim;
	
	
//	TODO :: PR_ASSERT(pa.pr_pair_weights.template_mispriming >= 0.0);
	double template_mispriming;

//	TODO :: PR_ASSERT(pa.pr_pair_weights.template_mispriming_th >= 0.0);
	double template_mispriming_th;
	
	public void set_template_mispriming_th(String datum) {
		// TODO Auto-generated method stub
		
	}

	public void set_template_mispriming(String datum) {
		// TODO Auto-generated method stub
		
	}

	public void set_repeat_sim(String datum) {
		// TODO Auto-generated method stub
		
	}

	public void set_product_size_lt(String datum) {
		// TODO Auto-generated method stub
		
	}

	public void set_product_size_gt(String datum) {
		// TODO Auto-generated method stub
		
	}

	public void set_product_tm_gt(String datum) {
		// TODO Auto-generated method stub
		
	}

	public void set_product_tm_lt(String datum) {
		// TODO Auto-generated method stub
		
	}

	public void set_compl_end_th(String datum) {
		// TODO Auto-generated method stub
		
	}

	public void set_compl_any_th(String datum) {
		// TODO Auto-generated method stub
		
	}

	public void set_compl_end(String datum) {
		// TODO Auto-generated method stub
		
	}

	public void set_compl_any(String datum) {
		// TODO Auto-generated method stub
		
	}

	public void set_diff_tm(String datum) {
		// TODO Auto-generated method stub
		
	}

	public void set_io_quality(String datum) {
		// TODO Auto-generated method stub
		
	}

	public void set_primer_quality(String datum) {
		// TODO Auto-generated method stub
		
	}



	public void p3_print_args() {
		System.out.format("  begin pr_pair_weights\n") ;
		System.out.format("    primer_quality %f\n", this.primer_quality) ;
		System.out.format("    io_quality %f\n", this.io_quality) ;
		System.out.format("    diff_tm %f\n", this.diff_tm) ;
		System.out.format("    compl_any %f\n", this.compl_any) ;
		System.out.format("    compl_end %f\n", this.compl_end) ;
		System.out.format("    compl_any_th %f\n", this.compl_any_th) ;
		System.out.format("    compl_end_th %f\n", this.compl_end_th) ;
		System.out.format("    product_tm_lt %f\n", this.product_tm_lt) ;
		System.out.format("    product_tm_gt %f\n", this.product_tm_gt) ;
		System.out.format("    product_size_lt %f\n", this.product_size_lt) ;
		System.out.format("    product_size_gt %f\n", this.product_size_gt) ;
		System.out.format("    repeat_sim %f\n", this.repeat_sim) ;
		System.out.format("    template_mispriming %f\n", this.template_mispriming) ;
		System.out.format("    template_mispriming_th %f\n", this.template_mispriming_th) ;
		System.out.format("  end pair_weights\n") ;		
	}
 }