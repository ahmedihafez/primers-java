package com.primer3.libprimer3;

import java.util.ArrayList;
import java.util.List;

public class  pair_array_t {
//	  int         storage_size = 0;
	  public int         num_pairs = 0;
	  // TODO :: make sure it is init. somewhere
	  public List<primer_pair> pairs = new ArrayList<primer_pair>();
	  pair_stats  expl = new pair_stats();
	  
	  
	  public String p3_get_pair_array_explain_string()
	  {
		  return expl.p3_pair_explain_string();
	  }


	
	  public void add_pair(primer_pair the_best_pair) {

		  pairs.add(the_best_pair);
		  this.num_pairs++;
		  
	  }
}