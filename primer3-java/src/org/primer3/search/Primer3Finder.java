package org.primer3.search;

import java.util.List;

import org.primer3.libprimer3.DPAlArgHolder;
import org.primer3.libprimer3.P3RetVal;
import org.primer3.libprimer3.PrimerPair;
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
	
	
	abstract protected void getLocalNextResult() throws Exception;
	
	
	public List<PrimerPair> getNextResult() throws Exception
	{
		// cache current set
		retval.best_pairs.cacheCurrent();
		this.getLocalNextResult();
		List<PrimerPair> newPairs = retval.best_pairs.pairs;
		retval.best_pairs.mergeBests();
		return newPairs;
	}
	
}
