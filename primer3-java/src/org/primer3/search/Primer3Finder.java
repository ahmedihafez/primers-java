package org.primer3.search;

import org.primer3.libprimer3.DPAlArgHolder;
import org.primer3.libprimer3.P3RetVal;
import org.primer3.libprimer3.THAlArgHolder;

public abstract class Primer3Finder {
	P3RetVal retval;
	DPAlArgHolder dpal_arg_to_use;
	THAlArgHolder thal_arg_to_use;
	THAlArgHolder thal_oligo_arg_to_use;
	
	
	public Primer3Finder (	P3RetVal retval,
	DPAlArgHolder dpal_arg_to_use,
	THAlArgHolder thal_arg_to_use,
	THAlArgHolder thal_oligo_arg_to_use)
	{
		this.retval  = retval;
		this. dpal_arg_to_use = dpal_arg_to_use;
		this. thal_arg_to_use = thal_arg_to_use;
		this. thal_oligo_arg_to_use = thal_oligo_arg_to_use;
	}
	
	
	abstract public void getNextResult() throws Exception;
	
}
