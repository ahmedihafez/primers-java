package org.primer3.search;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;

import org.primer3.libprimer3.DPAlArgHolder;
import org.primer3.libprimer3.P3GlobalSettings;
import org.primer3.libprimer3.P3RetVal;
import org.primer3.libprimer3.PrimerPair;
import org.primer3.libprimer3.PrimerRecord;
import org.primer3.libprimer3.THAlArgHolder;
import org.primer3.p3_seq_lib.seq_lib;
import org.primer3.search.multiplex.MultiplexResult;

public class P3MultiplexSearch {

	public DPAlArgHolder dpal_arg_to_use;
	public THAlArgHolder thal_arg_to_use;
	public THAlArgHolder thal_oligo_arg_to_use;
	
	// is it normal pcr or real time pcr
	
	HashMap<String, ArrayList<P3RetVal> > targetsGroups = new HashMap<String, ArrayList<P3RetVal>>();
	HashMap<String, MultiplexResult > targetsGroupsResult = new HashMap<String, MultiplexResult>();
	public void addTarget(P3RetVal retval) {
		// TODO check retval is not null
		if (!targetsGroups.containsKey(retval.sa.mutliplexGroup)) {
			targetsGroups.put(retval.sa.mutliplexGroup , new ArrayList<P3RetVal>());
		}
		
		targetsGroups.get(retval.sa.mutliplexGroup).add(retval);
	}
	
	public boolean search() throws Exception {
		
		if(targetsGroups.size() > 0 )
		{
		
			for( Entry< String, ArrayList<P3RetVal>> group : targetsGroups.entrySet()) {
				if(group.getValue().size() == 1 )
				{
					// nothing to do here it is just one sequence
					return false;
				}
				MultiplexResult result = search(group.getKey(),group.getValue());
				targetsGroupsResult.put(group.getKey(), result);
				
			}
			
		}
		return true;
	}

	private MultiplexResult search(String gname ,ArrayList<P3RetVal> group) throws Exception {
		
		
		
		
		// how many pair do want get ??
		int globalNumToReturn = 1 ;
		P3GlobalSettings pa;
		pa = group.get(0).pa;
		pa.setNumReturn(globalNumToReturn);
		MultiplexResult result = new MultiplexResult(this,gname,pa);
		dpal_arg_to_use = DPAlArgHolder.create_dpal_arg_holder();// create_dpal_arg_holder();
		thal_arg_to_use = THAlArgHolder.create_thal_arg_holder(pa.primersArgs);// create_thal_arg_holder(&pa.p_args);
		thal_oligo_arg_to_use = THAlArgHolder.create_thal_arg_holder(pa.oligosArgs);

		
		
		
//		seq_lib targets = getTargets(group);
		// test mispriming against other targets
		for(P3RetVal retval : group )
		{
			HashSet<String> exculdeSet = new HashSet<String>();
			exculdeSet.add(retval.sa.getSequenceName());
			exculdeSet.add("reverse " + retval.sa.getSequenceName());
			// this should be for list and right and also internal
			// This can not be done here
			// we might need to end up will a partial list, make the test on the level of a result set together 
			for(PrimerRecord h : retval.fwd.oligo) {
				
			}
//			for(PrimerRecord h : retval.rev.oligo) {
//				h.oligo_repeat_library_mispriming(pa, targets, h.rec_type, 
//						retval.rev.expl, dpal_arg_to_use, new StringBuilder(),exculdeSet);
//				h.repeat_sim.score.clear();
//			}
		}
		
		Set<String> searchRound = new HashSet<String>();
		HashMap<String,Boolean> hasMorePairs = new HashMap<String, Boolean>();
		for(P3RetVal reval : group )
		{
			String targetName = reval.sa.getSequenceName();
			searchRound.add(targetName);
			
			result.addP3RetVal(reval);
			hasMorePairs.put(targetName, true);
//			reval.getNextRusult();
		}
		// TODO : this should be the same in all set

		result.calcSpecific();
		result.inisSearch(dpal_arg_to_use, thal_arg_to_use, thal_oligo_arg_to_use);
		boolean canContinue = false;
		while (true)
		{
		// retval is already have some : no 
			System.err.println("Next Round ");
			System.err.println(searchRound);
			

			for(P3RetVal reval : group )
			{
				String targetName = reval.sa.getSequenceName();
				if(searchRound.contains(targetName))
				{
					if(hasMorePairs.get(targetName))
					{
						List< PrimerPair> newPairs = result.getNextRusult(targetName);
						if(newPairs.size() == 0 )
						{	
							hasMorePairs.put(targetName,false);
						}
						else
						{
							
							result.addPairs(reval.sa.getSequenceName(),newPairs);
						}
					}
				}// reval.choose_pairs(dpal_arg_to_use, thal_arg_to_use, thal_oligo_arg_to_use);
				// init a set
			}
			
			searchRound.clear();
			searchRound = result.getNextSearchRound(hasMorePairs);
			System.err.println(hasMorePairs);
			if ( result.hasGoodSets())
			{
				break;
			}
			canContinue = false;
			for(boolean anyCanContinue : hasMorePairs.values())
				if(anyCanContinue)
					canContinue = true;
			if(!canContinue)
				break;
		}
		result.sort();
		return result;
		
	}

	public void print_boulder(int io_version) {
		for( Entry<String, MultiplexResult> groupResult : targetsGroupsResult.entrySet()) {
			groupResult.getValue().print_boulder(io_version);
		}
		
	}

	

}
