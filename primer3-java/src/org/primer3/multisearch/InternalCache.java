package org.primer3.multisearch;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map.Entry;

import org.primer3.primer.MultiTargetPrimerRecord;
import org.primer3.primer.PrimerRecord;

class InternalCache {
		int mask ;
		int k ;
		MultiTargetScanner scanner;
		// this only at the leave level
		
		HashMap<Integer,PotentialPairSet>[] intCache ;
		HashMap<Integer,InternalCache> subCaches ;
		HashMap<MultiTargetPrimerRecord, SharedPotentialPairCollection> sharedCollection = new HashMap<MultiTargetPrimerRecord, SharedPotentialPairCollection>();

		
		/**
		 * top level node of cache
		 * @param MultiTargetScanner scanner 
		 */
		public InternalCache(MultiTargetScanner scanner  )
		{
			this.scanner = scanner;
			
			this.k = scanner.nTargets;
			intCache = new HashMap[k];
			for (int i = 0; i < this.intCache.length; i++) {
				intCache[i] = new HashMap<Integer, PotentialPairSet>();
			}
			subCaches = new HashMap<Integer,InternalCache>();
//			doneList.add(0);
			this.mask = 0;
			
			
//			createIndex();
		}
		InternalCache parentCache = null;
		public InternalCache( int mask , InternalCache parentCache) {
			this.mask = mask;
			this.k = parentCache.k;
			this.parentCache = parentCache;
			intCache = new HashMap[k];
			for (int i = 0; i < this.intCache.length; i++) {
				if(missing(i))
					intCache[i] = parentCache.intCache[i];
				else
					intCache[i] = new HashMap<Integer, PotentialPairSet>();
			}
			subCaches = new HashMap<Integer, InternalCache>();
			scanner =  parentCache.scanner;
		}
		
		// clean after filling targetIndex primers
		public void clean(int targetIndex) {
			System.err.println("Cleaning ");
			Iterator<Entry<MultiTargetPrimerRecord, SharedPotentialPairCollection>> it = sharedCollection.entrySet().iterator();
			while(it.hasNext())
			{
				Entry<MultiTargetPrimerRecord, SharedPotentialPairCollection> item = it.next();
				
				if(item.getValue().hasMissingValues(targetIndex))
					it.remove();
			}
			for (int baseIndex  : subCaches.keySet())  {
				subCaches.get(baseIndex).clean(targetIndex);
			}
			
		}



		HashSet<Integer> doneList = new HashSet<Integer>();
		Iterator<MultiTargetPrimerRecord> sharedNext ;
		public PotentialPairSetsCollection[] getNext() {
			// only once
			if(sharedNext == null)
				sharedNext = sharedCollection.keySet().iterator();
				
			if(!doneList.contains(0))
			{
				resetIntCache();
				PotentialPairSetsCollection[] res = getLocal();
				if(!sharedNext.hasNext()){
					doneList.add(0);
					sharedNext = sharedCollection.keySet().iterator();
					resetIntCache();
				}
				return res;
			}
			
			for (int baseIndex  : subCaches.keySet()) 
			{
				
				if(!doneList.contains(baseIndex))
				{
					while(true) {
						PotentialPairSetsCollection[] newCollection = subCaches.get(baseIndex).getNext();
						if(newCollection == null){
							if(!sharedNext.hasNext()){
								doneList.add(baseIndex);
								sharedNext = sharedCollection.keySet().iterator();
								resetIntCache();
								break;
							}
							else {
								resetIntCache();
								subCaches.get(baseIndex).reInitFromParent();
							}
						}
						else {
							return newCollection;
						}
									
					}
				}
			}
			return null;
		}

		private void reInitFromParent() {
			sharedNext = sharedCollection.keySet().iterator();
			for (int i = 0; i < this.intCache.length; i++) {
				if(missing(i))
					intCache[i] = parentCache.intCache[i];
			}
			
		}

		/** 
		 * make sure that there is next in sharedNext
		 */
		private void resetIntCache() {
			if(sharedCollection.size() == 0) return;
			MultiTargetPrimerRecord currentTrail = sharedNext.next();
			
			
			for (int i = 0; i < sharedCollection.get(currentTrail).pairsCollection.length; i++) {
				if(sharedCollection.get(currentTrail).pairsCollection[i] != null)
				{
					intCache[i] = sharedCollection.get(currentTrail).pairsCollection[i];
				}
			}
			
		}

		private PotentialPairSetsCollection[] getLocal() {
			PotentialPairSetsCollection[] productCalcPairs = new PotentialPairSetsCollection[k];
			for(int targetIndex =0 ; targetIndex < k;targetIndex++  )
			{
				productCalcPairs[targetIndex] = new PotentialPairSetsCollection();
				for (Integer pLen : intCache[targetIndex].keySet()) {
//					Integer pLen = (Integer) iterator.next();
//					calcPairsBucket.get(pLen).calcPairs[targetIndex].productCriterion = pLen;
					
					productCalcPairs[targetIndex].add(intCache[targetIndex].get(pLen));
				}
			}
			return productCalcPairs;
		}



//		private void createIndex() {
//			
//			int validK = (1 << k) - 1;
//			for (int i = 1; i <= validK; i++) {
//				// any numeber that can that does not have what i have can be added
//				if ( Integer.bitCount(i) == 1 )
//					continue;
//				if ( Integer.bitCount(i & mask ) == 0 )
//				{
//					if(!subCaches.containsKey(i))
//						subCaches.put(i, new InternalCache(i,this) );
//				}
//			}
//			
//		}




		private boolean missing(int i) {
			int validK = (1 << k) - 1;
			return (((1 << i ) & mask)) == 0;  
//			return false;
		}


		/**
		 * main entry to add a new pair 
		 * @param newPairInfo
		 */
		public void broadcast(PotentialPair newPairInfo) {
			
			PrimerRecord left = newPairInfo.left , right  = newPairInfo.right;
			
			int primerIndex = 0;
			int baseIndex = 0, complement = 0, validK = (1 << k) - 1;
			int nTargets = 0;
			if (left instanceof MultiTargetPrimerRecord || right instanceof MultiTargetPrimerRecord)
			{
				MultiTargetPrimerRecord selectedPrimer = (MultiTargetPrimerRecord) left;
				
				for (int i = 0; i < k; i++) {
					if(selectedPrimer.targetsToPrimer.containsKey(scanner.targets.get(i))) {
						baseIndex +=  (1 << i );
						nTargets++;
					}
				}
				
				
			}
			else
			{
				baseIndex = 1 << newPairInfo.targetIndex;
				// lets say that 0 is is the base for this case
				nTargets++;
			}
			complement = ~baseIndex & validK;
			int m1 = Integer.min(baseIndex, complement);
			int m2 = Integer.max(baseIndex, complement);
			if(nTargets == 1)
			{
				if(!intCache[newPairInfo.targetIndex].containsKey(newPairInfo.productCriterion))
					intCache[newPairInfo.targetIndex].put(newPairInfo.productCriterion, new PotentialPairSet(newPairInfo.productCriterion)  );
				intCache[newPairInfo.targetIndex].get(newPairInfo.productCriterion).add(newPairInfo);
			}
			else
			{
				if(!subCaches.containsKey(baseIndex))
					subCaches.put(baseIndex, new InternalCache(baseIndex,this) );
				subCaches.get(baseIndex).add(baseIndex,complement,newPairInfo);
			}
//			if(baseIndex > complement )
//			{
//				// if less do not do any thing
//				for (int i  : subCaches.keySet()) {
//					if(Integer.bitCount(i & baseIndex ) == 0)
//					{
//						subCaches.get(i).add(baseIndex, complement, newPairInfo);
//					}
//				}
//			}
		}
		
		public void fillMissing() {
			clean();
			int validK = (1 << k) - 1;
			for (int baseIndex  : subCaches.keySet()) {
				subCaches.get(baseIndex).clean();
				int complement = ~baseIndex & validK;
				if(baseIndex > complement )
				{
					// if less do not do any thing
					for (int i  : subCaches.keySet()) {
						if(Integer.bitCount(i & baseIndex ) == 0)
						{
							subCaches.get(i).complmentWith(subCaches.get(baseIndex));
						}
					}
				}
			}
			sharedNext = sharedCollection.keySet().iterator();
		}
		
		private void clean() {
			Iterator<Entry<MultiTargetPrimerRecord, SharedPotentialPairCollection>> it = sharedCollection.entrySet().iterator();
			while(it.hasNext())
			{
				Entry<MultiTargetPrimerRecord, SharedPotentialPairCollection> item = it.next();
				
				if(item.getValue().hasMissingValues())
					it.remove();
			}
			
		}

		private void complmentWith(InternalCache complCache) {
			
			
		}

		protected void add(int baseIndex,int complement, PotentialPair newPairInfo) {
			if (baseIndex == this.mask)
			{
				// add here
//				if(!intCache[newPairInfo.targetIndex].containsKey(newPairInfo.productCriterion))
//					intCache[newPairInfo.targetIndex].put(newPairInfo.productCriterion, new PotentialPairSet(newPairInfo.productCriterion)  );
//				intCache[newPairInfo.targetIndex].get(newPairInfo.productCriterion).add(newPairInfo);
				MultiTargetPrimerRecord keyPrimer = null;
				if(newPairInfo.left instanceof MultiTargetPrimerRecord)
				{
					keyPrimer = (MultiTargetPrimerRecord) newPairInfo.left;
				}
				if(newPairInfo.right instanceof MultiTargetPrimerRecord)
				{
					keyPrimer = (MultiTargetPrimerRecord) newPairInfo.right;
				}
				
				// do not add 
				if(keyPrimer == null)
					// FIXME :: need to debug this
					System.err.println("Keyprimer can not be null ?? ERR");
				
				// if this is not the first then we should have one in our collection other wise then ignore it since it is not going to be complemete collection
				int allPreMask = (1 << newPairInfo.targetIndex) - 1;
				if( (allPreMask & baseIndex) == 0 )
				{	
					sharedCollection.put(keyPrimer,new SharedPotentialPairCollection(keyPrimer, scanner.targets));
				}
				if(sharedCollection.containsKey(keyPrimer)) {
					sharedCollection.get(keyPrimer).addPotentialPair(newPairInfo);
				}
			}
			else
			{
				int validK = (1 << k) - 1;
				if(!subCaches.containsKey(baseIndex))
					subCaches.put(baseIndex, new InternalCache(baseIndex,this) );
				subCaches.get(baseIndex).add(baseIndex,complement,newPairInfo);
			}
			if(baseIndex > complement )
			{
				// if less do not do any thing
				for (int i  : subCaches.keySet()) {
					if(Integer.bitCount(i & baseIndex ) == 0)
					{
						subCaches.get(i).add(baseIndex, complement, newPairInfo);
					}
				}
			}
		}
		
	}