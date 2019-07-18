package org.primer3.search.multiplex;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;

import javax.sound.midi.Sequence;

import org.primer3.dpal.AlignmentException;
import org.primer3.libprimer3.LibPrimer3;
import org.primer3.primer.MultiTargetPrimerRecord;
import org.primer3.primer.PrimerPair;
import org.primer3.primer.PrimerRecord;
import org.primer3.thal.ThermodynamicAlignmentException;

public class MultiplexSet implements Comparable<MultiplexSet>  {

	boolean checkTmDiff = false;
	boolean checkProductsLen = true;
	double minTm = Double.MAX_VALUE , maxTm = Double.MIN_VALUE;
	MultiplexResult parentResult = null;
	
	HashMap<String, PrimerPair> set = new HashMap<String, PrimerPair>();
	
	HashMap<String , Integer> productsLen = new HashMap<String , Integer>();
	
	public MultiplexSet(MultiplexResult multiplexResult) {
		parentResult = multiplexResult;
	}

	public boolean hasPair(String targetName) {
		return set.containsKey(targetName);
	}

	public boolean addPair(String targetName, PrimerPair pair) throws Exception {
		// already have this target in my set ?
		if(hasPair(targetName)) 
			return false; // then return false
		// add new pair 
		// we already checked non target and target seq
		// now see 
		boolean status = true;
		
		
		
		
		
		
		
		
		
		
		
		PrimerRecord left  = pair.left , right = pair.right;
		// 1 product len/tm
		int newProdcutLen =  right.start - left.start+1;
		
		
		
		
		
		
		
		
		
		if (checkProductsLen)
		{
			
			if ( !checkSetProductsLen(newProdcutLen))
			{
				// not valid product can not add it here 
				// if got rejected becouse of len then see if the new pair has a
//				pair = checkUnversalPrimer( targetName,  pair);
//				if(pair == null )
				return false;
			}
		}
		// 1 TODO :: diff tm between primers.
		
		// 1 primers dimers between primers
		HashMap<String,Double> newPairMaxAny = new HashMap<String, Double>();
		
		
		double minTm = Double.min(this.minTm, left.temp);
		minTm = Double.min(minTm, right.temp);
		double maxTm = Double.max(this.maxTm, left.temp);
		maxTm = Double.max(maxTm, right.temp);
		
		double tmDiff = maxTm - minTm;
		if(tmDiff > this.parentResult.pa.getMaxDiffTm())
		{
			return false;
		}
		
		for(String  pTarget : set.keySet())
		{
			PrimerPair p = set.get(pTarget);
			

			
			PrimerRecord pLeft  = p.left , pRight = p.right;
			
			boolean leftIsUnversal = false , rightIsUnversal = false;
			if(org.primer3.sequence.Sequence.equals(pLeft.getOligoSeq(),left.getOligoSeq()))
				leftIsUnversal = true;
			if(org.primer3.sequence.Sequence.equals(pRight.getOligoSeq(),right.getOligoSeq()))
				rightIsUnversal = true;
			double maxAny = 0;
			if(!leftIsUnversal)
			{
				maxAny = checkDimers(pLeft, left);
				if(Double.isFinite(maxAny ))
				{
					newPairMaxAny.put(pTarget, maxAny);
				}
				else
				{
					return false;
				}
				maxAny = checkDimers(pLeft, right);
				if(Double.isFinite(maxAny ))
				{
					newPairMaxAny.put(pTarget, maxAny);
				}
				else
				{
					return false;
				}
			}
			if(!rightIsUnversal)
			{
				maxAny = checkDimers(pRight, left );
				if(Double.isFinite(maxAny ))
				{
					newPairMaxAny.put(pTarget, maxAny);
				}
				else
				{
					return false;
				}
				maxAny = checkDimers(pRight, right);
				if(Double.isFinite(maxAny ))
				{
					newPairMaxAny.put(pTarget, maxAny);
				}
				else
				{
					return false;
				}
			}
		}
		
		
		this.minTm = minTm;
		this.maxTm = maxTm;
		
		
		compl_any.put(targetName, newPairMaxAny);
		// update the rest
		for(String pTarget : newPairMaxAny.keySet())
		{
			compl_any.get(pTarget).put(targetName, newPairMaxAny.get(pTarget));
		}

		productsLen.put(targetName, newProdcutLen);
		set.put(targetName, pair);
		return status;
	}



//	private PrimerPair checkUnversalPrimer(String targetName, PrimerPair pair) throws Exception {
//		
//		// see first if it contain any universal primer
//		
//		// left is universal and right is specific
//		if(pair.right.targetSpecificIndex.size() == 0)
//		{
//			for(String existingTarget : set.keySet())
//			{
//				PrimerPair existingTargetPair = set.get(existingTarget);
//				if(    existingTargetPair.left.targetSpecificIndex.containsKey(targetName)) {
//					// existingTarget and targetName can share same left primer
////					System.err.println("existingTarget and targetName can share same left primer");
//					// TODO :: this is not the correct way to do that
//					PrimerRecord orgLeft = pair.left;
//					
//					// first check if the product is different then proceed 
//					
//					
//					int altLeftStart = existingTargetPair.left.targetSpecificIndex.get(targetName) ;
//					int newProductLen = pair.right.start - altLeftStart +1;
//					if(checkSetProductsLen(newProductLen))
//					{
//						PrimerPair newPair = new PrimerPair();
//						// new left  should be in the 
//						PrimerRecord newLeft = parentResult.getLeftFor(targetName,altLeftStart, existingTargetPair.left.length );
//						if(newLeft != null && newLeft.OK_OR_MUST_USE())
//						{
//							int pairStatus = newPair.characterize_pair(parentResult.sourcePairsResult.get(targetName),  newLeft, pair.right,  
//									parentResult.p3MultiplexSearcher.dpal_arg_to_use, parentResult.p3MultiplexSearcher.thal_primers_arg_to_use, false);
//							if(pairStatus == PrimerPair.PAIR_OK)
//							{
//								// TODO :: check internal 
//								
//								// TODO :: if ok
//								
//								newPair.obj_fn(parentResult.sourcePairsResult.get(targetName).pa);
//								return newPair;
//							}
//						}
//						
//					}
//				}
//			}
//		}
//		// right is universal and left is specific
//		if(pair.right.targetSpecificIndex.size() >= 0 && pair.left.targetSpecificIndex.size() == 0)
//		{
//			for(String existingTarget : set.keySet())
//			{
//				if(pair.right.targetSpecificIndex.containsKey(existingTarget)) {
//					// existingTarget and targetName can share same left primer
//					// System.err.println("existingTarget and targetName can share same right primer");
//				}
//			}
//		}
//		// also consider a unversal primer for a subgroup if we can not find any ?? TODO ::
//		return null;
//	}

	// compl_any as the new pair are added
	HashMap<String,HashMap<String,Double>>  compl_any = new HashMap<String,HashMap<String,Double>>();
	
	
//	double compl_end = 0;
	private double checkDimers_old(PrimerRecord p1, PrimerRecord p2) throws ThermodynamicAlignmentException {
		// TODO Auto-generated method stub
		char[] s1,s2 ,s1_rev, s2_rev;
		s1 = p1.getOligoSeq();
		s2 = p2.getOligoSeq();
		s1_rev = p1.getOligoRevSeq();
		s2_rev = p2.getOligoRevSeq();
		
		
		
		 double compl_any = LibPrimer3.align_thermod(s1, s2_rev, this.parentResult.p3MultiplexSearcher.get_thal_primers_arg_to_use().any);
	      if ( compl_any > this.parentResult.pa.getPairComplAnyTH()) {
//	         if (update_stats) {
//	            pair_expl.compl_any++; 
//	         }
	         return Double.POSITIVE_INFINITY;
	      }
	      double compl_end = 0 ;
	      
	      double compl_endTemp = LibPrimer3.align_thermod(s1, s2_rev, this.parentResult.p3MultiplexSearcher.get_thal_primers_arg_to_use().end1);
	      if (compl_endTemp > this.parentResult.pa.getPairComplEndTH())
	    	   return Double.POSITIVE_INFINITY;
	      compl_end = Double.max(compl_end, compl_endTemp);
	      
	      compl_endTemp = LibPrimer3.align_thermod(s1, s2_rev, this.parentResult.p3MultiplexSearcher.get_thal_primers_arg_to_use().end2); /* Triinu Please check */
	      if (compl_endTemp > this.parentResult.pa.getPairComplEndTH())
	    	  return Double.POSITIVE_INFINITY;
	      compl_end = Double.max(compl_end, compl_endTemp);	     
	      
	      
	      compl_endTemp = LibPrimer3.align_thermod(s1, s2_rev, this.parentResult.p3MultiplexSearcher.get_thal_primers_arg_to_use().end1);
	      if (compl_endTemp > this.parentResult.pa.getPairComplEndTH())
	    	  return Double.POSITIVE_INFINITY;
	      compl_end = Double.max(compl_end, compl_endTemp);
	      
	      compl_endTemp = LibPrimer3.align_thermod(s1, s2_rev, this.parentResult.p3MultiplexSearcher.get_thal_primers_arg_to_use().end2); /* Triinu Please check */
	      if (compl_endTemp > this.parentResult.pa.getPairComplEndTH())
	    	  return Double.POSITIVE_INFINITY;
	      compl_end = Double.max(compl_end, compl_endTemp);	
	      
	      
//	      if (compl_end > this.parentResult.pa.getPairComplEndTH()) {
//	         if (update_stats) {
//	            pair_expl.compl_end++; 
//	         }
//	         if (!must_use) return PAIR_FAILED;
//	         return false;
//	      }
	      
	     return compl_end;
	}
	
	
//	double compl_end = 0;
	private double checkDimers(PrimerRecord p1, PrimerRecord p2) throws ThermodynamicAlignmentException {
		// TODO Auto-generated method stub
		char[] s1,s2 ,s1_rev, s2_rev;
		s1 = p1.getOligoSeq();
		s2 = p2.getOligoSeq();
		s1_rev = p1.getOligoRevSeq();
		s2_rev = p2.getOligoRevSeq();
		
		
		
		 double compl_any = LibPrimer3.align_thermod(s1, s2_rev, this.parentResult.p3MultiplexSearcher.get_thal_primers_arg_to_use().any);
	      if ( compl_any > this.parentResult.pa.getPairComplAnyTH()) {
//	         if (update_stats) {
//	            pair_expl.compl_any++; 
//	         }
	         return Double.POSITIVE_INFINITY;
	      }
	      double compl_end = 0 ;
	      
	      double compl_endTemp = LibPrimer3.align_thermod(s1, s2_rev, this.parentResult.p3MultiplexSearcher.get_thal_primers_arg_to_use().end1);
	      if (compl_endTemp > this.parentResult.pa.getPairComplEndTH())
	    	   return Double.POSITIVE_INFINITY;
	      compl_end = Double.max(compl_end, compl_endTemp);
	      
	      compl_endTemp = LibPrimer3.align_thermod(s1, s2_rev, this.parentResult.p3MultiplexSearcher.get_thal_primers_arg_to_use().end2); /* Triinu Please check */
	      if (compl_endTemp > this.parentResult.pa.getPairComplEndTH())
	    	  return Double.POSITIVE_INFINITY;
	      compl_end = Double.max(compl_end, compl_endTemp);	     
	      
	      
	      compl_endTemp = LibPrimer3.align_thermod(s1, s2_rev, this.parentResult.p3MultiplexSearcher.get_thal_primers_arg_to_use().end1);
	      if (compl_endTemp > this.parentResult.pa.getPairComplEndTH())
	    	  return Double.POSITIVE_INFINITY;
	      compl_end = Double.max(compl_end, compl_endTemp);
	      
	      compl_endTemp = LibPrimer3.align_thermod(s1, s2_rev, this.parentResult.p3MultiplexSearcher.get_thal_primers_arg_to_use().end2); /* Triinu Please check */
	      if (compl_endTemp > this.parentResult.pa.getPairComplEndTH())
	    	  return Double.POSITIVE_INFINITY;
	      compl_end = Double.max(compl_end, compl_endTemp);	
	      
	      
//	      if (compl_end > this.parentResult.pa.getPairComplEndTH()) {
//	         if (update_stats) {
//	            pair_expl.compl_end++; 
//	         }
//	         if (!must_use) return PAIR_FAILED;
//	         return false;
//	      }
	      
	     return compl_end;
	}
	

	private boolean checkSetProductsLen(int newProdcutLen) {
		boolean isValid = true;
		
		ArrayList<Integer> productsLen = new ArrayList<Integer>() ;
		productsLen.addAll(this.productsLen.values());
		productsLen.add(newProdcutLen);
		productsLen.sort(new Comparator<Integer>() {

			@Override
			public int compare(Integer o1, Integer o2) {
				// TODO Auto-generated method stub
				return Integer.compare(o1, o2);
			}
		});
		for(int i = 1 ; i < productsLen.size() ; i ++ )
		{
			int len_p = productsLen.get(i-1);
			int len_n = productsLen.get(i);
			int diff =  Math.abs(len_n-len_p);
			
			int minLen = Integer.min(len_n, len_p);
			
			int diffMin = 45 ;
			
			
			if(minLen < 500  )
				diffMin = 45;
			else if(minLen < 850)
				diffMin = 90;
			else diffMin = 180;
			if( diff <= diffMin ) // not exactly should be
			{
				isValid = false;
				break;
			}
		}
		
		return isValid;
		
	}

	public MultiplexSet getAlt(String targetName, PrimerPair pair) throws Exception {

		// make sure you have a pair for targetName
		MultiplexSet newAltset = new MultiplexSet(parentResult);
		newAltset.productsLen = (HashMap<String, Integer>) productsLen.clone();
		newAltset.productsLen.remove(targetName);
//		for(Entry<String, Integer> entry :  productsLen.entrySet()  )
//		{
//			if(!targetName.equals(entry.getKey()))
//				newAltset.productsLen.put(entry.getKey(), entry.getValue());
//		}
		for(String target : compl_any.keySet()   )
		{
			if(!targetName.equals(target))
			{
				HashMap<String, Double> targetComplAnyClone =  (HashMap<String, Double>) compl_any.get(targetName).clone();
				targetComplAnyClone.remove(targetName);
				newAltset.compl_any.put(target, targetComplAnyClone);
			}
		}
		newAltset.set = (HashMap<String, PrimerPair>) this.set.clone();
		newAltset.set.remove(targetName);
		// now that we have acopy without the pair/target
		if(newAltset.addPair(targetName, pair))
		{
			return newAltset;
		}
		
		
		
		return null;
	}

	public void print_buolder( int io_version) {
		
		System.out.println("********");
		System.out.println("PRIMER_MULTIPLEX_NUM_TARGETS="+this.set.size());
		System.out.println("***");
		for(Entry<String, HashMap<String, Double>> mAnySp : compl_any.entrySet())
		{
			System.out.print(mAnySp.getKey() + " : ");
			for(Entry<String,  Double> mAny : mAnySp.getValue().entrySet())
			{
				System.out.print("( " + mAny.getKey() + " > " +  mAny.getValue() + " ),");
			}
			System.out.println("");
		}
		System.out.println("****");
		for(Entry<String,  PrimerPair >  pair : set.entrySet())
		{
			
			pair.getValue().print_boulder( parentResult.sourcePairsResult.get(pair.getKey()), pair.getKey(), io_version);
			System.out.println("***");
		}
		System.out.println("********");
	}

	
	double score = Double.MAX_VALUE;
	boolean scoreNeedUpdate = true;
	public double getScore() {
		if(scoreNeedUpdate)
		{
			score = 0 ;
			for(PrimerPair pair : set.values())
			{
				score += pair.pair_quality;
			}
		}
		return score;
	}

	public int getNTargets() {
		return set.size();
	}

	
	public HashMap<String, PrimerPair> getPairs ()
	{
		return set;
	}

	@Override
	public int compareTo(MultiplexSet o) {
		if(this.getNTargets() == o.getNTargets())
			return -Double.compare(this.getScore(), o.getScore());
		return Integer.compare(this.getNTargets(), o.getNTargets());
	}

}
