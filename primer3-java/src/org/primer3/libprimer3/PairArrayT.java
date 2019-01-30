package org.primer3.libprimer3;

import java.util.ArrayList;
import java.util.List;

public class  PairArrayT {
//	  int         storage_size = 0;
	  public int         num_pairs = 0;
	  // TODO :: make sure it is init. somewhere
	  public List<PrimerPair> pairs = new ArrayList<PrimerPair>();
	  public PairStats  expl = new PairStats();
	  
	  
	  public String p3_get_pair_array_explain_string()
	  {
		  return expl.p3_pair_explain_string();
	  }


	
	  public void add_pair(PrimerPair the_best_pair) {

		  pairs.add(the_best_pair);
		  this.num_pairs++;
		  
	  }
	  
	  public void clearSet()
	  {
		  pairs.clear();
		  num_pairs = 0;
	  }


	List<PrimerPair> currentSet = null;
	public void cacheCurrent() {
		currentSet = this.pairs;
		this.pairs =new ArrayList<PrimerPair>();
		this.num_pairs = 0;
	}



	public void mergeBests() {
		
		currentSet.addAll(pairs);
		num_pairs =  currentSet.size();
		pairs = currentSet;
		currentSet = null;
		
	}
}