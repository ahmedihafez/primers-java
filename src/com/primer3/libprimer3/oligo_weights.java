package com.primer3.libprimer3;

/** 
 * Arguments to the primer program as a whole.  Values for these arguments are
 * retained _across_ different input records.  (These are the so-called
 * "Global" arguments in the documentation.)
 */
public class oligo_weights {

	double compl_any;
	double compl_any_th;
	double compl_end;
	double compl_end_th;
	double end_quality;
	double end_stability;
	double gc_content_gt;
	double gc_content_lt;
	double hairpin_th;
	double length_gt;
	double length_lt;
    double num_ns;
    double pos_penalty;
    double repeat_sim;
    double seq_quality;
    double temp_cutoff;
    double temp_gt;
    double temp_lt;
    double template_mispriming;
    double template_mispriming_th;
    double failure_rate;
	
    public void set_template_mispriming_th(String datum) {
    	this.template_mispriming_th = Double.parseDouble(datum);
	}
    
	
    public void set_end_quality(String datum) {
    	this.end_quality = Double.parseDouble(datum);
	}
    
	public void set_seq_quality(String datum) {
		this.seq_quality =  Double.parseDouble(datum);
	}
	
	public void set_repeat_sim(String datum) {
		this.repeat_sim = Double.parseDouble(datum);
	}
	
	public void set_num_ns(String datum) {
		this.num_ns = Double.parseDouble(datum);
	}
	
	public void set_hairpin_th(String datum) {
		this.hairpin_th = Double.parseDouble(datum);	
	}
	
	public void set_compl_end_th(String datum) {
		this.compl_end_th = Double.parseDouble(datum);
	}
	
	public void set_compl_any_th(String datum) {
		this.compl_any_th = Double.parseDouble(datum);
	}
	
	public void set_compl_end(String datum) {
		this.compl_end = Double.parseDouble(datum);
	}
	
	public void set_compl_any(String datum) {
		this.compl_any = Double.parseDouble(datum);
	}
	
	public void set_length_gt(String datum) {
		this.length_gt = Double.parseDouble(datum);
	}
	
	public void set_length_lt(String datum) {
		this.length_lt = Double.parseDouble(datum);
	}
	
	public void set_gc_content_lt(String datum) {
		this.gc_content_lt = Double.parseDouble(datum);
	}
	
	public void set_gc_content_gt(String datum) {
		this.gc_content_gt = Double.parseDouble(datum);
	}
	
	public void set_temp_lt(String datum) {
		this.temp_lt = Double.parseDouble(datum);
	}
	
	public void set_temp_gt(String datum) {
		this.temp_gt = Double.parseDouble(datum);
	}
	public void set_failure_rate(String datum) {
		this.failure_rate = Double.parseDouble(datum);
	}
	public void set_template_mispriming(String datum) {
		this.template_mispriming = Double.parseDouble(datum);
	}
	public void set_end_stability(String datum) {
		this.end_stability = Double.parseDouble(datum);
	}
	public void set_pos_penalty(String datum) {
		this.pos_penalty = Double.parseDouble(datum);
	}

   }