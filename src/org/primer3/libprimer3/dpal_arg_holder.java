package org.primer3.libprimer3;

import org.primer3.dpal.DPAlignmentArgs;
import org.primer3.dpal.DPAlignment;

public class dpal_arg_holder {
	DPAlignmentArgs local;
	DPAlignmentArgs end;
	DPAlignmentArgs local_end;
	DPAlignmentArgs local_ambig;
	DPAlignmentArgs local_end_ambig;
	
//	PRIMER_LEFT_0_SEQUENCE= GTGAAGCCTCAGGTAGTGCA
//	PRIMER_RIGHT_0_SEQUENCE= CTCTGTCGACTTTGCCACCA
	
	/* Create the dpal arg holder */
	public static dpal_arg_holder create_dpal_arg_holder (){
		dpal_arg_holder h = new dpal_arg_holder();
		
		
		h.local = new DPAlignmentArgs();
		h.local.set_dpal_args();
		h.local.flag = DPAlignment.DPAL_LOCAL;
		
		
		h.end = new DPAlignmentArgs();
		h.end.set_dpal_args();
		h.end.flag = DPAlignment.DPAL_GLOBAL_END;
		
		h.local_end = new DPAlignmentArgs();
		h.local_end.set_dpal_args();
		h.local_end.flag = DPAlignment.DPAL_LOCAL_END;
		
		h.local_ambig = new DPAlignmentArgs();
//		*h->local_ambig = *h->local;
		h.local_ambig.set_dpal_args();
		h.local_ambig.flag = DPAlignment.DPAL_LOCAL;
		h.local_ambig.dpal_set_ambiguity_code_matrix();
		
		h.local_end_ambig = new DPAlignmentArgs();
//		*h->local_end_ambig = *h->local_end;
		h.local_end_ambig.set_dpal_args();
		h.local_end_ambig.flag = DPAlignment.DPAL_LOCAL_END;
		h.local_end_ambig.dpal_set_ambiguity_code_matrix();
		return h;
	}
}
