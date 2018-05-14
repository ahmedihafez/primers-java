package com.primer3.libprimer3;

import com.primer3.thal.thal_alignment_type;
import com.primer3.thal.thal_args;

public class thal_arg_holder {
	
	thal_args any;
	thal_args end1;
	thal_args end2;
	thal_args hairpin_th;
	
	
	
	
	static public thal_arg_holder   create_thal_arg_holder (PrimersOligosArguments po_args) {
		thal_arg_holder h = new thal_arg_holder();
		
		h.any = new  thal_args();
		h.any.set_thal_default_args();
		h.any.type = thal_alignment_type.thal_any;
		h.any.mv = po_args.getSaltConcentration();
		h.any.dv = po_args.getDivalentConcentration();
		h.any.dntp = po_args.getDntpConcentration();
		h.any.dna_conc = po_args.getDnaConcentration();
		
		
		h.end1 = new thal_args();
		h.end1.set_thal_default_args();
		h.end1.type = thal_alignment_type.thal_end1;
		h.end1.mv = po_args.getSaltConcentration();
		h.end1.dv = po_args.getDivalentConcentration();
		h.end1.dntp = po_args.getDntpConcentration();
		h.end1.dna_conc = po_args.getDnaConcentration();
		
		h.end2 = new thal_args();
		h.end2.set_thal_default_args();
		h.end2.type = thal_alignment_type.thal_end2;
		h.end2.mv = po_args.getSaltConcentration();
		h.end2.dv = po_args.getDivalentConcentration();
		h.end2.dntp = po_args.getDntpConcentration();
		h.end2.dna_conc = po_args.getDnaConcentration();
		
		h.hairpin_th = new thal_args();
		h.hairpin_th.set_thal_default_args();
		h.hairpin_th.type = thal_alignment_type.thal_hairpin;
		h.hairpin_th.mv = po_args.getSaltConcentration();
		h.hairpin_th.dv = po_args.getDivalentConcentration();
		h.hairpin_th.dntp = po_args.getDntpConcentration();
		h.hairpin_th.dna_conc = po_args.getDnaConcentration();
		h.hairpin_th.dimer = 0;
		return h;
	}
}
