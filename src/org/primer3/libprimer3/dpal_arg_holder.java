package org.primer3.libprimer3;

import org.primer3.dpal.dpal_args;
import org.primer3.dpal.dpallib;

public class dpal_arg_holder {
	dpal_args local;
	dpal_args end;
	dpal_args local_end;
	dpal_args local_ambig;
	dpal_args local_end_ambig;
	
//	PRIMER_LEFT_0_SEQUENCE= GTGAAGCCTCAGGTAGTGCA
//	PRIMER_RIGHT_0_SEQUENCE= CTCTGTCGACTTTGCCACCA
	
	/* Create the dpal arg holder */
	public static dpal_arg_holder create_dpal_arg_holder (){
		dpal_arg_holder h = new dpal_arg_holder();
		
		
		h.local = new dpal_args();
		h.local.set_dpal_args();
		h.local.flag = dpallib.DPAL_LOCAL;
		
		
		h.end = new dpal_args();
		h.end.set_dpal_args();
		h.end.flag = dpallib.DPAL_GLOBAL_END;
		
		h.local_end = new dpal_args();
		h.local_end.set_dpal_args();
		h.local_end.flag = dpallib.DPAL_LOCAL_END;
		
		h.local_ambig = new dpal_args();
//		*h->local_ambig = *h->local;
		h.local_ambig.set_dpal_args();
		h.local_ambig.flag = dpallib.DPAL_LOCAL;
		h.local_ambig.dpal_set_ambiguity_code_matrix();
		
		h.local_end_ambig = new dpal_args();
//		*h->local_end_ambig = *h->local_end;
		h.local_end_ambig.set_dpal_args();
		h.local_end_ambig.flag = dpallib.DPAL_LOCAL_END;
		h.local_end_ambig.dpal_set_ambiguity_code_matrix();
		return h;
	}
}
