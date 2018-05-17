package org.primer3.libprimer3;

import org.primer3.thal.ThermodynamicAlignmentType;
import org.primer3.thal.ThermodynamicAlignmentArguments;

public class thal_arg_holder {
	
	ThermodynamicAlignmentArguments any;
	ThermodynamicAlignmentArguments end1;
	ThermodynamicAlignmentArguments end2;
	ThermodynamicAlignmentArguments hairpin_th;
	
	
	
	
	static public thal_arg_holder   create_thal_arg_holder (PrimersOligosArguments po_args) {
		thal_arg_holder h = new thal_arg_holder();
		
		h.any = new  ThermodynamicAlignmentArguments();
		h.any.set_thal_default_args();
		h.any.type = ThermodynamicAlignmentType.thal_any;
		h.any.mv = po_args.getSaltConcentration();
		h.any.dv = po_args.getDivalentConcentration();
		h.any.dntp = po_args.getDntpConcentration();
		h.any.dna_conc = po_args.getDnaConcentration();
		
		
		h.end1 = new ThermodynamicAlignmentArguments();
		h.end1.set_thal_default_args();
		h.end1.type = ThermodynamicAlignmentType.thal_end1;
		h.end1.mv = po_args.getSaltConcentration();
		h.end1.dv = po_args.getDivalentConcentration();
		h.end1.dntp = po_args.getDntpConcentration();
		h.end1.dna_conc = po_args.getDnaConcentration();
		
		h.end2 = new ThermodynamicAlignmentArguments();
		h.end2.set_thal_default_args();
		h.end2.type = ThermodynamicAlignmentType.thal_end2;
		h.end2.mv = po_args.getSaltConcentration();
		h.end2.dv = po_args.getDivalentConcentration();
		h.end2.dntp = po_args.getDntpConcentration();
		h.end2.dna_conc = po_args.getDnaConcentration();
		
		h.hairpin_th = new ThermodynamicAlignmentArguments();
		h.hairpin_th.set_thal_default_args();
		h.hairpin_th.type = ThermodynamicAlignmentType.thal_hairpin;
		h.hairpin_th.mv = po_args.getSaltConcentration();
		h.hairpin_th.dv = po_args.getDivalentConcentration();
		h.hairpin_th.dntp = po_args.getDntpConcentration();
		h.hairpin_th.dna_conc = po_args.getDnaConcentration();
		h.hairpin_th.dimer = 0;
		return h;
	}
}
