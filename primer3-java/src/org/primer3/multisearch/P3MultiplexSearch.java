package org.primer3.multisearch;

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
import org.primer3.libprimer3.THAlArgHolder;
import org.primer3.p3_seq_lib.seq_lib;
import org.primer3.primer.PrimerPair;
import org.primer3.primer.PrimerRecord;
import org.primer3.search.multiplex.MultiplexResult;

/**
 * CLI implementation of the search , Why ??
 * @author Ahmed Hafez
 *
 */
public class P3MultiplexSearch implements IP3MultiSearch {

	public DPAlArgHolder dpal_arg_to_use;
	/**
	 * th align args for primers
	 */
	public THAlArgHolder thal_primers_arg_to_use;
	/**
	 * th align args for probes
	 */
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
		
		dpal_arg_to_use = DPAlArgHolder.create_dpal_arg_holder();// create_dpal_arg_holder();
		thal_primers_arg_to_use = THAlArgHolder.create_thal_arg_holder(pa.primersArgs);// create_thal_arg_holder(&pa.p_args);
		thal_oligo_arg_to_use = THAlArgHolder.create_thal_arg_holder(pa.oligosArgs);
		final MultiTargetScanner multiTargetScan = new MultiTargetScanner(this,gname,pa);
		
		
		
//		seq_lib targets = getTargets(group);
		// test mispriming against other targets
//		for(P3RetVal retval : group )
//		{
//			HashSet<String> exculdeSet = new HashSet<String>();
//			exculdeSet.add(retval.sa.getSequenceName());
//			exculdeSet.add("reverse " + retval.sa.getSequenceName());
			// this should be for list and right and also internal
			// This can not be done here
			// we might need to end up will a partial list, make the test on the level of a result set together 
//			for(PrimerRecord h : retval.fwd.oligo) {
//				
//			}
//			for(PrimerRecord h : retval.rev.oligo) {
//				h.oligo_repeat_library_mispriming(pa, targets, h.rec_type, 
//						retval.rev.expl, dpal_arg_to_use, new StringBuilder(),exculdeSet);
//				h.repeat_sim.score.clear();
//			}
//		}
		
		for(P3RetVal reval : group )
		{
			String targetName = reval.sa.getSequenceName();
			
			
			multiTargetScan.addP3RetVal(reval);
//			reval.getNextRusult();
		}
		// TODO : this should be the same in all set
		multiTargetScan.prePare();
		
		final P3OptimzedMultiTargetFinder sFinderSp = new P3OptimzedMultiTargetFinder(multiTargetScan, dpal_arg_to_use, thal_primers_arg_to_use, thal_oligo_arg_to_use);
//		final P3OptimzedMultiTargetFinder sFindermLeft = new P3OptimzedMultiTargetFinder(multiTargetScan, dpal_arg_to_use, thal_primers_arg_to_use, thal_oligo_arg_to_use);
//		final P3OptimzedMultiTargetFinder sFindermRight = new P3OptimzedMultiTargetFinder(multiTargetScan, dpal_arg_to_use, thal_primers_arg_to_use, thal_oligo_arg_to_use);

		sFinderSp.initSearch(
//			new P3OptimzedMultiTargetFinder.TargetPrimerFiller() {
//			
//			@Override
//			public void fillTargetPrimer(int targetIndex) {
//				String targetName = multiTargetScan.targets.get(targetIndex);
//				sFinderSp.advanceSearchForProductOnly(targetIndex , multiTargetScan.targetsToSpFwd.get(targetName),multiTargetScan.targetsToSpRev.get(targetName));
//				// TODO :: need to be revised
////				sFinderSp.advanceSearchForProductOnly(targetIndex , multiTargetScan.multiFwd,multiTargetScan.targetsToSpRev.get(targetName));
//
//			}
//		}
		);
		
//		sFindermLeft.initSearch(new P3OptimzedMultiTargetFinder.TargetPrimerFiller() {
//			
//			@Override
//			public void fillTargetPrimer(int targetIndex) {
//				String targetName = multiTargetScan.targets.get(targetIndex);
//				sFindermLeft.advanceSearchForProductOnly(targetIndex , multiTargetScan.multiFwd,multiTargetScan.targetsToSpRev.get(targetName));
//			}
//		});
//		sFindermRight.initSearch(new P3OptimzedMultiFindermLeft.TargetPrimerFiller() {
//			
//			@Override
//			public void fillTargetPrimer(int targetIndex) {
//				String targetName = multiTargetScan.targets.get(targetIndex);
//				sFindermRight.advanceSearchForProductOnly(targetIndex , multiTargetScan.targetsToSpFwd.get(targetName), multiTargetScan.multiRev );
//			}
//		});
		
		Object ret = null; 
		while ( true )
		{
			if( sFinderSp.hasMore() )
				ret = sFinderSp.getLocalNextResult();
//			if( sFindermLeft.hasMore() )
//				ret = sFindermLeft.getLocalNextResult(result);
//			if( sFindermRight.hasMore() )
//				ret = sFindermRight.getLocalNextResult(result);
			if (sFinderSp.hasGoodSets())
			{
				break;
			}
			if( 	
					!sFinderSp.hasMore() 
//					&&
//					!sFindermLeft.hasMore() 
//					&& 
//					!sFindermRight.hasMore()
					)
				break;
		}
		
		return sFinderSp.getFinalResult();

		
	}

	public void print_boulder(int io_version) {
		for( Entry<String, MultiplexResult> groupResult : targetsGroupsResult.entrySet()) {
			groupResult.getValue().print_boulder(io_version);
		}
		
	}

	@Override
	public DPAlArgHolder getDPAL_ArgToUse() {
		return dpal_arg_to_use;
	}

	@Override
	public THAlArgHolder get_thal_primers_arg_to_use() {
		return thal_primers_arg_to_use;
	}

	

}
