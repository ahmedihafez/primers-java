package org.primer3.search.multiplex;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;

import org.primer3.dpal.AlignmentException;
import org.primer3.libprimer3.DPAlArgHolder;
import org.primer3.libprimer3.LibPrimer3;
import org.primer3.libprimer3.OligoType;
import org.primer3.libprimer3.P3GlobalSettings;
import org.primer3.libprimer3.P3RetVal;
import org.primer3.libprimer3.PrimerPair;
import org.primer3.libprimer3.PrimerRecord;
import org.primer3.libprimer3.THAlArgHolder;
import org.primer3.p3_seq_lib.seq_lib;
import org.primer3.search.P3MultiplexSearch;
import org.primer3.sequence.Sequence;
import org.primer3.thal.ThermodynamicAlignmentException;

public class MultiplexResult {
	
	class TargetCoverSet {
		public TargetCoverSet(String targetName) {
			this.targetName = targetName;
		}
		String targetName ;
		// Number of sets has a pair for this one
		int numOfCovers = 0;
		int numOfCoversInCompleteSet=0;
	}
	ArrayList<TargetCoverSet> coversSet = new ArrayList<MultiplexResult.TargetCoverSet>();
	HashMap<String ,TargetCoverSet> coversSetHash = new HashMap<String ,TargetCoverSet>();

	
	seq_lib targets = new seq_lib();
	// add just once
	HashMap<String,P3RetVal> sourcePairsResult = new HashMap<String, P3RetVal>();
	HashMap<String,List<P3RetVal>> pairsResult = new HashMap<String, List<P3RetVal>>();

	
	String groupName;
	ArrayList<MultiplexSet> multiplexSets = new ArrayList<MultiplexSet>();
	P3MultiplexSearch p3MultiplexSearcher = null;
	P3GlobalSettings pa ;
	
	
	
	public MultiplexResult(P3MultiplexSearch p3MultiplexSearch,String groupName,  P3GlobalSettings pa) {
		p3MultiplexSearcher = p3MultiplexSearch;
		this.pa = pa;
		multiplexSets.add(new MultiplexSet(this));
		this.groupName = groupName;
	}
	static final String REVERSE_PREFIX_TARGETS = "reverse " ;
	public void addP3RetVal(P3RetVal retVal) {
		
		// just add once
		String targetName = retVal.sa.getSequenceName();
		if(!sourcePairsResult.containsKey(targetName))
		{
			sourcePairsResult.put(targetName, retVal);
			TargetCoverSet newCoverSet =  new TargetCoverSet(targetName);
			coversSetHash.put(targetName,newCoverSet);
			coversSet.add(newCoverSet);
			expectedMaxNTargetsInSet++;
			targets.add(retVal.sa.getSequenceName(), retVal.sa.getTrimmedSequence(),  Sequence.p3_reverse_complement( retVal.sa.getTrimmedSequence()));
			targets.add(REVERSE_PREFIX_TARGETS + retVal.sa.getSequenceName() ,  retVal.sa.getUpcasedSeqRev(), retVal.sa.getUpcasedSeq());
		}
		
	}
	
	
	
	ArrayList<MultiplexSet> newAltsSets =  new ArrayList<MultiplexSet>();
	public void addPairs(String targetName, List<PrimerPair> newPairs) throws Exception {
		
		newAltsSets.clear();
		for(PrimerPair pair : newPairs)
		{
			
			boolean ifSingltonAddedOnce = false;
			for(MultiplexSet set : multiplexSets)
			{
				// check set has pair for this target
//				if(set.hasPair(targetName))
				{
					// if set has target pair
//					MultiplexSet newAltsSet = set.getAlt(targetName,pair);
//					if(newAltsSet != null) {
//						newAltsSets.add(newAltsSet);
//						coversSetHash.get(targetName).numOfCovers++;
//						if(set.getNTargets() == expectedMaxNTargetsInSet)
//						{
//							coversSetHash.get(targetName).numOfCoversInCompleteSet++;
//							
//						}
//
//					}
					// check if new pair could replace a better score than current one
				}
//				else
				{
					// check if pair can be added to the set
					if(set.addPair(targetName,pair))
					{
						// pair add no problems
						coversSetHash.get(targetName).numOfCovers++;
						if(set.getNTargets() == expectedMaxNTargetsInSet)
						{
						
							coversSetHash.get(targetName).numOfCoversInCompleteSet++;
						}
					}
					else
					{
						// pair can not be added and set does not have it
						// start a new set with this pair and add pairs from the current set that can be added to it 
						MultiplexSet newAltSet = new MultiplexSet(this);
						newAltSet.addPair(targetName, pair);
						for(Entry<String, PrimerPair> oldPair : set.set.entrySet())
						{
							if(newAltSet.addPair(oldPair.getKey(), oldPair.getValue()))
								coversSetHash.get(oldPair.getKey()).numOfCovers++;
						}
						if(set.hasPair(targetName))
						{
							// add only if it has a better score
							if(set.compareTo(newAltSet) < 0 )
							{
								newAltsSets.add(newAltSet);
								ifSingltonAddedOnce = true;
							}
						}
						else if(newAltSet.getNTargets() > 1  || !ifSingltonAddedOnce)
						{ 
							newAltsSets.add(newAltSet);
							ifSingltonAddedOnce = true;
						}
//						if(newAltSet.getNTargets() == 1 )
//							ifSingltonAddedOnce = true;
						coversSetHash.get(targetName).numOfCovers++;
						
					}
				}

			}	
		}
		multiplexSets.addAll(newAltsSets);
		
	}

	int numOfSetToReturn = 10;
	int numOfGoodSets = 0;
	int expectedMaxNTargetsInSet=0;
	int iter = 0;
	public boolean hasGoodSets() {
		numOfGoodSets = 0;
		for(MultiplexSet set : multiplexSets)
		{
			if(set.getNTargets() >= expectedMaxNTargetsInSet)
			{
				numOfGoodSets++;
			}
		}
		
		iter++;
		if(iter >= 4 )
			return true;
		
		return numOfGoodSets > numOfSetToReturn;
	}

	HashMap<PrimerRecord,HashMap<String,Double> > misprimingScores = new HashMap<PrimerRecord, HashMap<String,Double>>();
	public boolean calcMispriming(String targetName, PrimerPair pair, Collection<String> keySet) throws AlignmentException {
		
		ArrayList<String> targetsNameInLib =   new ArrayList<String>(keySet);
		for(String nName : keySet)
		{
			targetsNameInLib.add(REVERSE_PREFIX_TARGETS + nName);
		}
		
		PrimerRecord left  = pair.left , right = pair.right;
		boolean isValid  = calcMispriming(targetName,left,targetsNameInLib);
		if (isValid)
			isValid = calcMispriming(targetName,right,targetsNameInLib);
		if (isValid && pair.intl != null)
			isValid = calcMispriming(targetName,right,targetsNameInLib);
		
		return isValid;
	}
	
	/**
	 * return false if the given primer bind/has an alignment to any of the given seq targets
	 * @param targetName
	 * @param primer
	 * @param keySet
	 * @return
	 * @throws AlignmentException
	 */
	private boolean calcMispriming(String targetName, PrimerRecord primer, Collection<String> keySet) throws AlignmentException {
		double max_lib_compl;
		boolean  max_lib_compl_is_percent = false;
		/* First, check the oligo against the repeat library. */
		if (primer.rec_type == OligoType.OT_INTL) {
			max_lib_compl =  pa.oligosArgs.getMaxRepeatCompl();
			max_lib_compl_is_percent = pa.oligosArgs.maxRepeatComplIsPercent;
		} else {
			max_lib_compl =  pa.primersArgs.getMaxRepeatCompl();
			max_lib_compl_is_percent = pa.primersArgs.maxRepeatComplIsPercent;

		}
		double w = 0;
		for(int i = 0; i < targets.seq_lib_num_seq(); i++)
		{
			String libTargetSeqName = targets.getName(i);
			
			if(libTargetSeqName .equals(targetName) ||
				libTargetSeqName .equals( REVERSE_PREFIX_TARGETS+ targetName)
					||  !keySet.contains(libTargetSeqName)  )
				continue;
			if(!misprimingScores.containsKey(primer))
			{
				misprimingScores.put(primer, new HashMap<String, Double>());
			}
			if(misprimingScores.get(primer).containsKey(libTargetSeqName)){
				w = misprimingScores.get(primer).get(libTargetSeqName);
			}
			else {
				if(primer.rec_type == OligoType.OT_INTL)
				{
					w = LibPrimer3.align(primer.getOligoSeq(),
							targets.getSeq(i),
							(pa.isLibAmbiguityCodesConsensus() ?  p3MultiplexSearcher.dpal_arg_to_use.local_ambig
													: p3MultiplexSearcher.dpal_arg_to_use.local),
							pa.getMispriming3EndScore() );
				}
				else if(primer.rec_type == OligoType.OT_LEFT)
				{
					w = LibPrimer3.align(primer.getOligoSeq(),
							targets.getSeq(i),
							(pa.isLibAmbiguityCodesConsensus() ?  p3MultiplexSearcher.dpal_arg_to_use.local_end_ambig
													: p3MultiplexSearcher.dpal_arg_to_use.local_end),
							pa.getMispriming3EndScore() );
				}
				else if(primer.rec_type == OligoType.OT_RIGHT)
				{
					w = LibPrimer3.align(primer.getOligoRevSeq(),
							targets.getSeqRevCompl(i),
							(pa.isLibAmbiguityCodesConsensus() ?  p3MultiplexSearcher.dpal_arg_to_use.local_end_ambig
													: p3MultiplexSearcher.dpal_arg_to_use.local),
							pa.getMispriming3EndScore() );
				}
				if (max_lib_compl_is_percent)
					w = 100 *( w / primer.length);
				misprimingScores.get(primer).put(libTargetSeqName, w);
				if(w > max_lib_compl)
					return false;
			}
			
			
		}
		return true;
	}




	public static void main(String[] args)
	{
		HashMap<String, String> t = new HashMap<String, String>();
		t.put("a", "a");
		t.put("b", "c");
		t.put("c", "b");
		
		ArrayList<String> a = new ArrayList<String>(t.keySet());
		
		System.out.println(a);
		a.clear();
		System.out.println(a);
		System.out.println(t);
		
	}



	public void calcSpecific() throws AlignmentException
	{
		for( Entry<String,  P3RetVal> kValue : sourcePairsResult.entrySet() )
		{
			P3RetVal retval = kValue.getValue();
			List<P3RetVal> newRetVals = new ArrayList<P3RetVal>();
			P3RetVal spOnly = new P3RetVal(retval.pa, retval.sa);
			P3RetVal leftToSpRight = new P3RetVal(retval.pa, retval.sa);
			P3RetVal spLeftToRight = new P3RetVal(retval.pa, retval.sa);
			newRetVals.add(spOnly);
			newRetVals.add(leftToSpRight);
			newRetVals.add(spLeftToRight);
			
			pairsResult.put(kValue.getKey(), newRetVals);

			for(PrimerRecord h : retval.fwd.oligo) {
				 h.calcSpecific(targets,kValue.getKey(),  REVERSE_PREFIX_TARGETS , pa.isLibAmbiguityCodesConsensus() ?  p3MultiplexSearcher.dpal_arg_to_use.local_end_ambig
							: p3MultiplexSearcher.dpal_arg_to_use.local_end);
				 if(h.targetSpecificIndex.size() > 0 && h.isTargetSpecific)
				 {
					 leftToSpRight.fwd.add_oligo_to_oligo_array(h); 
				 }
				 if(h.isTargetSpecific && h.targetSpecificIndex.size() == 0) {
					 spLeftToRight.fwd.add_oligo_to_oligo_array(h); 
					 spOnly.fwd.add_oligo_to_oligo_array(h); 
					 
				 }
			}
			for(PrimerRecord h : retval.rev.oligo) {
				h.calcSpecific(targets,kValue.getKey(),  REVERSE_PREFIX_TARGETS, pa.isLibAmbiguityCodesConsensus() ?  p3MultiplexSearcher.dpal_arg_to_use.local_end_ambig
						: p3MultiplexSearcher.dpal_arg_to_use.local_end);
				if(h.targetSpecificIndex.size() > 0 && h.isTargetSpecific)
				 {
					 spLeftToRight.rev.add_oligo_to_oligo_array(h); 
				 }
				if(h.isTargetSpecific && h.targetSpecificIndex.size() == 0)
				 {
					 leftToSpRight .rev.add_oligo_to_oligo_array(h); 
					 spOnly.rev.add_oligo_to_oligo_array(h); 
				 }
			}
			System.err.println("Finish");
		}
	}
	
	

	public void print_boulder(int io_version) {
		System.out.println("Result Sets for " + this.groupName);
		for(int i = 0 ; i < multiplexSets.size();i++) {
			System.out.println("Sets #" + i);
			multiplexSets.get(i).print_buolder(io_version);
				
		}
		
	}




	public void sort() {
		
		multiplexSets.sort(new Comparator<MultiplexSet>() {

			@Override
			public int compare(MultiplexSet o1, MultiplexSet o2) {
				
				if(o1.getNTargets() == o2.getNTargets())
					return Double.compare(o1.getScore(), o2.getScore());
				return -Integer.compare(o1.getNTargets(), o2.getNTargets());
			}
		});
		
	}
	
	
	public Set<String> getNextSearchRound(HashMap<String, Boolean> hasMorePairs){
		HashSet<String> newRound = new HashSet<String>();
		
		coversSet.sort(new Comparator<TargetCoverSet>() {

			@Override
			public int compare(TargetCoverSet o1, TargetCoverSet o2) {
				// TODO Auto-generated method stub
				if(o1.numOfCoversInCompleteSet == o2.numOfCoversInCompleteSet)
					return Integer.compare(o1.numOfCovers, o2.numOfCovers);
				return Integer.compare(o1.numOfCoversInCompleteSet, o2.numOfCoversInCompleteSet);
			}
		});
		
		boolean zeroRound = false;
		expectedMaxNTargetsInSet = sourcePairsResult.size();
		for(TargetCoverSet c : coversSet) {
			
			if(!hasMorePairs.get(c.targetName))
			{
				if( c.numOfCoversInCompleteSet == 0 )
					expectedMaxNTargetsInSet --;
				continue;
			}
			
			if(c.numOfCoversInCompleteSet == 0 )
			{
					newRound.add(c.targetName);
					// break;
			}
			else
			{
				newRound.add(c.targetName);
			}

				
		}
		
		return newRound;
		
	}

	public List<PrimerPair> getNextRusult(String targetName) throws Exception {
		
		List<PrimerPair> newPairs = new ArrayList<PrimerPair>();
		
		for(P3RetVal retval : pairsResult.get(targetName))
		{
			newPairs.addAll(retval.getNextRusult());
		}
		
		return newPairs;
	}

	public void inisSearch(DPAlArgHolder dpal_arg_to_use, THAlArgHolder thal_arg_to_use,
			THAlArgHolder thal_oligo_arg_to_use) {
		

		
		for(List<P3RetVal> retvals : pairsResult.values() )
		{
			for(P3RetVal retval :  retvals)
			{
				retval.inisSearch( dpal_arg_to_use,  thal_arg_to_use,
						 thal_oligo_arg_to_use);
			}
		}
		for(P3RetVal retval : sourcePairsResult.values()){
			retval.inisSearch( dpal_arg_to_use,  thal_arg_to_use,
					 thal_oligo_arg_to_use);
		}
	}

	public PrimerRecord getLeftFor(String targetName, int altLeftStart, int len) {
		
		for(PrimerRecord left : this.sourcePairsResult.get(targetName).fwd.oligo)
		{
			if(left.start == altLeftStart && left.length == len )
				return left;
		}
		return null;
	}

}
