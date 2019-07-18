package org.primer3.multisearch;

import java.util.ArrayList;
import java.util.HashMap;

import org.primer3.primer.MultiTargetPrimerRecord;

/**
 * Hold and organize a pairs with a shared  a multi target primer this collection is should be obly for on primer
 * @author Ahmed Hafez
 *
 */
class SharedPotentialPairCollection {
		
		MultiTargetPrimerRecord sharedPrimer;
		HashMap<Integer,PotentialPairSet>[] pairsCollection ;
//		HashMap<Integer,Potenti> pairsCollection;
		public SharedPotentialPairCollection(MultiTargetPrimerRecord sharedPrimer , ArrayList<String> targetNames ) {
			this.sharedPrimer = sharedPrimer;
			pairsCollection = new HashMap[targetNames.size()];
			for(int i = 0 ; i < targetNames.size() ; i ++)
			{
				if( sharedPrimer.targetsToPrimer.containsKey(targetNames.get(i)))
					pairsCollection[i] = new HashMap<Integer, PotentialPairSet>();
			}
			
		}
		
		
		public void addPotentialPair(PotentialPair pair) {
			
			if(pairsCollection[pair.targetIndex] != null) {
				if(!pairsCollection[pair.targetIndex].containsKey(pair.productCriterion))
					pairsCollection[pair.targetIndex].put(pair.productCriterion,new  PotentialPairSet(pair.productCriterion));
				pairsCollection[pair.targetIndex].get(pair.productCriterion).add(pair);
			}
		}


		/**
		 * return true if there is one list is empty 
		 * call this at the end of filling the list with all target.
		 * @return
		 */
		public boolean hasMissingValues() {
			
			for (int i = 0; i < pairsCollection.length; i++) {
				if(pairsCollection[i]!= null && pairsCollection[i].size() == 0)
					return true;
			}
			
			return false;
		}

		/**
		 * return true if the list of targetIndex is empty
		 * call this after filling for the targetIndex
		 * @param targetIndex
		 * @return
		 */
		public boolean hasMissingValues(int targetIndex) {
			if(pairsCollection[targetIndex]!= null && pairsCollection[targetIndex].size() == 0)
				return true;
			return false;
		}
	}