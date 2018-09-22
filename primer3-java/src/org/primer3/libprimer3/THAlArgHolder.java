package org.primer3.libprimer3;

import org.primer3.thal.ThermodynamicAlignmentType;
import org.primer3.thal.ThermodynamicAlignmentArguments;

public class THAlArgHolder {
	
	ThermodynamicAlignmentArguments any;
	ThermodynamicAlignmentArguments end1;
	ThermodynamicAlignmentArguments end2;
	ThermodynamicAlignmentArguments hairpin_th;
	
	
	
	
	static public THAlArgHolder   create_thal_arg_holder (PrimersOligosArguments po_args) {
		THAlArgHolder h = new THAlArgHolder();
		
		h.any = new  ThermodynamicAlignmentArguments();
		h.any.setThAlDefaultArgs();
		h.any.setAlignmentType(ThermodynamicAlignmentType.thal_any);
		h.any.setMonovalentConc(po_args.getSaltConcentration());
		h.any.setDivalentConc(po_args.getDivalentConcentration());
		h.any.setDntpConc(po_args.getDntpConcentration());
		h.any.setDnaConc(po_args.getDnaConcentration());
		
		
		h.end1 = new ThermodynamicAlignmentArguments();
		h.end1.setThAlDefaultArgs();
		h.end1.setAlignmentType(ThermodynamicAlignmentType.thal_end1);
		h.end1.setMonovalentConc(po_args.getSaltConcentration());
		h.end1.setDivalentConc(po_args.getDivalentConcentration());
		h.end1.setDntpConc(po_args.getDntpConcentration());
		h.end1.setDnaConc(po_args.getDnaConcentration());
		
		h.end2 = new ThermodynamicAlignmentArguments();
		h.end2.setThAlDefaultArgs();
		h.end2.setAlignmentType(ThermodynamicAlignmentType.thal_end2);
		h.end2.setMonovalentConc(po_args.getSaltConcentration());
		h.end2.setDivalentConc(po_args.getDivalentConcentration());
		h.end2.setDntpConc(po_args.getDntpConcentration());
		h.end2.setDnaConc(po_args.getDnaConcentration());
		
		h.hairpin_th = new ThermodynamicAlignmentArguments();
		h.hairpin_th.setThAlDefaultArgs();
		h.hairpin_th.setAlignmentType(ThermodynamicAlignmentType.thal_hairpin);
		h.hairpin_th.setMonovalentConc(po_args.getSaltConcentration());
		h.hairpin_th.setDivalentConc(po_args.getDivalentConcentration());
		h.hairpin_th.setDntpConc(po_args.getDntpConcentration());
		h.hairpin_th.setDnaConc(po_args.getDnaConcentration());
		h.hairpin_th.setCalcDimer(0);
		return h;
	}
}
