package org.primer3.multisearch;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.primer3.dpal.AlignmentException;
import org.primer3.libprimer3.DPAlArgHolder;
import org.primer3.libprimer3.DebugInfo;
import org.primer3.libprimer3.LibPrimer3;
import org.primer3.libprimer3.P3GlobalSettings;
import org.primer3.libprimer3.P3RetVal;
import org.primer3.libprimer3.P3Task;
import org.primer3.libprimer3.PairArrayT;
import org.primer3.libprimer3.PairStats;
import org.primer3.libprimer3.SeqArgs;
import org.primer3.libprimer3.THAlArgHolder;
import org.primer3.search.PotentialPairSetsCombination;
import org.primer3.search.multiplex.MultiplexResult;
import org.primer3.search.multiplex.MultiplexSet;
import org.primer3.multisearch.PotentialPair;
import org.primer3.primer.MultiTargetPrimerRecord;
import org.primer3.primer.PrimerPair;
import org.primer3.primer.PrimerRecord;
import org.primer3.thal.ThermodynamicAlignmentException;

import com.sun.xml.internal.bind.v2.runtime.unmarshaller.XsiNilLoader.Array;

// Stateful search of primers pairs 
public class P3OptimzedMultiTargetFinder {
	


	
	DPAlArgHolder dpal_arg_to_use;
	THAlArgHolder thal_arg_to_use;
	THAlArgHolder thal_oligo_arg_to_use;
	// Objects used in the Stateful search 
//	int[][] max_j_seen;
	
	MultiTargetScanner scanner = null;
	public P3OptimzedMultiTargetFinder(MultiTargetScanner scanner , DPAlArgHolder dpal_arg_to_use, THAlArgHolder thal_arg_to_use,
			THAlArgHolder thal_oligo_arg_to_use) {
		
		this. dpal_arg_to_use = dpal_arg_to_use;
		this. thal_arg_to_use = thal_arg_to_use;
		this. thal_oligo_arg_to_use = thal_oligo_arg_to_use;
		this.scanner = scanner;
//		max_j_seen = new int[scanner.nTargets][];
//		calcPairs = new TreeSet[scanner.nTargets];
//		calcPairsList = new ArrayList[scanner.nTargets];
//		pLens =  new List[scanner.nTargets];
//		productCalcPairs = new PotentialPairSetsCollection[scanner.nTargets];
//		intCache = new HashMap[scanner.nTargets];
		globalCache = new InternalCache(scanner);
//		calcPairsBucket = new TreeSet[scanner.nTargets];
//		for(int i = 0 ; i < scanner.nTargets;i++)
//		{
//			productCalcPairs[i] = new PotentialPairSetsCollection();
//		}
		
//		resetSearch();
		
	}
	
	
//	public void resetSearch()
//	{
//		for (int i = 0; i < scanner.nTargets; i++) {
//			for (int j = 0; j < max_j_seen[i].length; j++) {		
//				max_j_seen[i][j] = -1;
//			}
//		}
//	}
	

	static double ALPHA = 1;
	
	
	
	
//	class ijPairProductBucket {
//		PotentialPairSet[] calcPairs; 
//		int productSize;
//		public ijPairProductBucket(int k) {
//			calcPairs = new PotentialPairSet[k];
//			for (int i = 0; i < calcPairs.length; i++) {
//				calcPairs[i] = new PotentialPairSet();
//			}
//		}
//	}
//	HashMap<Integer, ijPairProductBucket> calcPairsBucket = new HashMap<Integer, P3OptimzedMultiTargetFinder.ijPairProductBucket>();
//	PotentialPairSetsCollection[] productCalcPairs ;
//	Comparator<ijPairProductBucket> ijPairProductBucketComparator = new Comparator<P3OptimzedMultiFinder.ijPairProductBucket>() {
//		
//		@Override
//		public int compare(ijPairProductBucket o1, ijPairProductBucket o2) {
//			return Integer.compare(o1.productSize, o2.productSize);
//		}
//	};
	public class mPairs {
		
		PotentialPair[] pairSet;
		
		int setStatus = PrimerPair.PAIR_UNCHARACTRIZED;
		
		public mPairs(PotentialPair[] pairSet)
		{
			if(pairSet != null )
			{
				this.pairSet = pairSet;
				estimatedQuilty = 0;
				for (int i = 0; i < pairSet.length; i++) {
					
					estimatedQuilty += pairSet[i].estimatedQuilty;
				}
			}
		}
		double estimatedQuilty = Double.MAX_VALUE;

		public int check() throws Exception {
			
			for (int p = 0; p < pairSet.length; p++) {
				PotentialPair pair = pairSet[p];
//				PrimerRecord left = scanner.multiFwd.get(pair.j);
//				PrimerRecord right = scanner.targetsToSpRev.get(scanner.targets.get(p)).get(pair.i);
				
				if(pair.pairStatus == PrimerPair.PAIR_OK)
				{
					
				}
				else if (pair.pairStatus == PrimerPair.PAIR_UNCHARACTRIZED) {
					pair.checkpair(dpal_arg_to_use, thal_arg_to_use);;
				}
				
				if(pair.pairStatus == PrimerPair.PAIR_FAILED) {	
					setStatus = PrimerPair.PAIR_FAILED;
					return setStatus;
				}
				
			}
//			for (int p = 0; p < pairSet.length; p++) {
//				PotentialPair pair = pairSet[p];
//				for (int i = 0; i < p; i++) {
//					if(!checkmTargetPrimers(pair, pairSet[i])) 
//					{
//						setStatus = PrimerPair.PAIR_FAILED;
//						return setStatus;
//					}
//				}
//				
//			}
			

			
			setStatus = PrimerPair.PAIR_OK;
			return setStatus;
		}
		
		private boolean checkmTargetPrimers(PotentialPair o1, PotentialPair o2) {
			
			if(o1.left instanceof MultiTargetPrimerRecord && o2.left instanceof MultiTargetPrimerRecord)
			{
//				// both are left
				MultiTargetPrimerRecord o1Left = (MultiTargetPrimerRecord) o1.left;
				MultiTargetPrimerRecord o2Left = (MultiTargetPrimerRecord) o2.left;
//				// if o1 left has the o2 taget then they should be the same left other wise no
				if(     o1Left.targetsToPrimer.containsKey(o2.getTargetName())  || 
						o2Left.targetsToPrimer.containsKey(o1.getTargetName()))
				{
					return o1Left == o2Left ;
				}
			}
			return true;
		}

	}
	
	
	
	
	
//	mPairs mPairs1 = new mPairs(null);
	

	// this
//	TreeSet<ijPair>[] calcPairs;
//	ArrayList<ijPair>[] calcPairsList;
	
//	ArrayList<HashMap<String,ijPair>> cachePiars = new ArrayList<HashMap<String,ijPair>>();
	PotentialPairSetsCombination mPairsIterator ;
	TreeSet<PotentialPairSetsCombination> allmPairIterators = new TreeSet<PotentialPairSetsCombination>(new Comparator<PotentialPairSetsCombination>() {

		@Override
		public int compare(PotentialPairSetsCombination o1, PotentialPairSetsCombination o2) {
			// TODO Auto-generated method stub
			int res = Double.compare(o1.getEstimatedScore(), o2.getEstimatedScore());
			if(res == 0 )
			{
				return o1 == o2 ? 0 : -1; 
			}
			return res ;
		}
	});
	public interface TargetPrimerFiller{
		void fillTargetPrimer(int targetIndex);
	}
	
	public void initSearch(TargetPrimerFiller filler) {
		for(int targetIndex =0 ; targetIndex < scanner.nTargets;targetIndex++  )
		{
//			String targetName = scanner.targets.get(targetIndex);
			
			filler.fillTargetPrimer(targetIndex);
			globalCache.clean(targetIndex);
//			advanceSearchForProductOnly(targetIndex , scanner.multiFwd,scanner.targetsToSpRev.get(targetName));
//			advanceSearchForProductOnly(targetIndex , scanner.targetsToSpFwd.get(targetName),scanner.targetsToSpRev.get(targetName));
//			advanceSearchForProductOnly(targetIndex , scanner.targetsToSpFwd.get(targetName), scanner.multiRev );
//			for (Integer pLen : intCache[targetIndex].keySet()) {
//				Integer pLen = (Integer) iterator.next();
//				calcPairsBucket.get(pLen).calcPairs[targetIndex].productCriterion = pLen;
//				productCalcPairs[targetIndex].add(intCache[targetIndex].get(pLen));
//			}
		}
		
		globalCache.fillMissing();
		
		// init values here
		while (true) {
			PotentialPairSetsCollection[] productCalcPairs = globalCache.getNext();
			if(productCalcPairs == null)
				break;
			PotentialPairSetsCombination newComp = new PotentialPairSetsCombination(productCalcPairs, PotentialPairSetsCombination.productLenConflictResolution );
			if(! (newComp.getEstimatedScore() >= Double.MAX_VALUE) )
				allmPairIterators.add(newComp);
		}
//		productCalcPairs = globalCache.getNext();
//		mPairsIterator = new PotentialPairSetsCombination(productCalcPairs);

//		mPairsIterator = allmPairIterators.pollFirst();
		
//		lastpLenComp = mPairsIterator.getComp();
		
//		allmPairIterators.add(mPairsIterator);
		
//		if(lastpLenComp == null) {
//			searchExausted = true;
//			return;
//		}
//		PotentialPairsLocalSearch newLocalSearch = new PotentialPairsLocalSearch(lastpLenComp);
//		localSearchmList.add(newLocalSearch);
//		lastpLenComp = null;
//		for(int i =0 ; i < scanner.nTargets;i++  )
//		{
//			advanceSearchForProductOnly(i , scanner.targetsToSpFwd.get(scanner.targets.get(i)),scanner.targetsToSpRev.get(scanner.targets.get(i)));
//		}
		
//		for(int targetIndex =0 ; targetIndex < scanner.nTargets;targetIndex++  ) {
//		
//		}
	}
	// use this to alternate search between different searcher
	int nMaxIterantion = 10000;
	
	
	InternalCache globalCache ;

//	HashMap<Integer,PotentialPairSet>[] intCache ;

	boolean searchExausted = false;
	PotentialPairSet[] lastpLenComp ;
	double topLastpLenCompScore = Double.MAX_VALUE;
	
	
	// one TreeSet for lastpLenComp
	TreeSet<PotentialPairsLocalSearch> localSearchmList = new TreeSet<PotentialPairsLocalSearch>( new Comparator<PotentialPairsLocalSearch>() {

		@Override
		public int compare(PotentialPairsLocalSearch o1, PotentialPairsLocalSearch o2) {
			return Double.compare(o1.getScore(), o2.getScore());
		}
	});
	public mPairs getLocalNextResult()  throws Exception {

		
//		advanceSearchForMultiplexSet();
		
//		int[] iS = new int[scanner.nTargets];
//		int[] jS = new int[scanner.nTargets];
//		int[] fwdSizes = new int[scanner.nTargets];
//		int[] revSizes = new int[scanner.nTargets];
//		for(int i =0 ; i < scanner.nTargets;i++  )
//		{
//			fwdSizes[i] = scanner.targetsToSpFwd.get(scanner.targets.get(i)).size();
//			revSizes[i] = scanner.targetsToSpRev.get(scanner.targets.get(i)).size();
//		}
//		
//		// who will clear the result ??
//		PairArrayT[] best_pairs =  new PairArrayT[scanner.nTargets]  ;//retval.best_pairs;
//		for (int i = 0; i < scanner.nTargets; i++) {
////			best_pairs = scanner.
//		}
//
////		int n_int; /* Index of the internal oligo */
//		boolean[] update_stats = new boolean[scanner.nTargets];  /* Flag to indicate whether pair_stats
//		                            should be updated. */
//		PrimerPair[] h = new PrimerPair[scanner.nTargets];             /* The current pair which is being evaluated. */
//		PrimerPair[] the_best_pair   = new PrimerPair[scanner.nTargets]  ;// new PrimerPair(); /* The best pair is being "remembered". */
//		PairStats[] pair_expl =   new PairStats[scanner.nTargets] ; // retval.best_pairs.expl; /* For statistics */
//		int[] product_sizes = new int[scanner.nTargets];
////		int the_best_i, the_best_j;
//
//
//		for(int i =0 ; i < scanner.nTargets;i++  )
//		{
//			pair_expl[i] = scanner.revals.get(i).best_pairs.expl;
//			best_pairs[i] = scanner.revals.get(i).best_pairs;
//		}

		

		
		// #############
		
//		ProductsLenComp3 iterator = new ProductsLenComp3(productCalcPairs);
		int nIterantion = nMaxIterantion;
		
		while(true) {
			
			// return a valid comp from available plens
//			if(nIterantion >= nMaxIterantion)
//				break;
			if(!searchExausted  && lastpLenComp == null ) {
				mPairsIterator = allmPairIterators.pollFirst();
				if(mPairsIterator!= null) {
					lastpLenComp = mPairsIterator.getComp();
					if(mPairsIterator.hasMore() )
						allmPairIterators.add(mPairsIterator);
				}
				if(mPairsIterator == null   ){
					searchExausted = true;
					topLastpLenCompScore = Double.MAX_VALUE;
				}
				if(lastpLenComp != null   ){
					topLastpLenCompScore = PotentialPairSet.getOrginalScore( lastpLenComp);
					if(topLastpLenCompScore >= Double.MAX_VALUE )
						System.err.println("Stop");
				}
				else
				{
					topLastpLenCompScore = Double.MAX_VALUE;
				}
			}

			if(lastpLenComp != null ) {

				boolean shouldAdd = false;
				if(localSearchmList.size() == 0 )
					shouldAdd = true;
				else {
//					topLastpLenCompScore = PotentialPairSet.calcScore( lastpLenComp);
					shouldAdd = topLastpLenCompScore < localSearchmList.first().getScore();
				}
				// FIXME :: How we could pull new combination faster without trying each time
				shouldAdd = true;
				if(shouldAdd)
				{
					PotentialPairsLocalSearch newLocalSearch = new PotentialPairsLocalSearch(lastpLenComp);
					if(newLocalSearch.hasMaxSet()) {
						localSearchmList.add(newLocalSearch);
					}
					lastpLenComp = null;
					
				}
			}
			if(searchExausted && localSearchmList.isEmpty())
			{
				break;
			}
			if(localSearchmList.isEmpty())
				continue;
//			while(true) {
				PotentialPairsLocalSearch lSearchs = localSearchmList.pollFirst();
//				if(lSearchs.calcPairs[0].productCriterion == 406)
//					System.err.println("Stop");
				PotentialPair[] newset = lSearchs.getNew();
				if(newset == null ) {
					// end if this lSearchs
					continue;
				}
				nIterantion--;
				mPairs newPairComp = new mPairs(newset);
				newPairComp.check();
				
				lSearchs.update();
				if(lSearchs.hasMaxSet()){
					// add it again after updating
					localSearchmList.add(lSearchs);
				}
				else
				{
					
					
				}
				if( newPairComp.setStatus == PrimerPair.PAIR_OK)
				{
					MultiplexSet res = new  MultiplexSet(scanner.result);
					for (int i = 0; i < newset.length; i++) {
						PotentialPair ijPair = newset[i];
						if(!res.addPair(scanner.targets.get(i), ijPair.pairRecord)) {
							System.err.println("Found One Fail !!");
							return null;
						}
					}
					scanner.result.addSet(res);
					System.err.println("Found One !!");
					return newPairComp;
				}
				
//			}
			

			//mPairsIterator.update(lastpLenComp);
//			lastpLenComp = null;
//				if(localSearchmList.size() > 10 )
				{
					while(localSearchmList.size() > 0 && ( localSearchmList.last().getScore() >= Double.MAX_VALUE ||  !localSearchmList.last().hasMaxSet()) )
						localSearchmList.pollLast();
				}
				if(nIterantion <= 0)
					break;
		}
		
		
//		System.err.println("Done");
		return null;
		
		//#############
		

	}

	
//	List<Integer>[] pLens = null;
	
	
	

	

	

	
	
	
	//	HashMap<Integer,ArrayList<ijPairs>> bins = new HashMap<Integer, ArrayList<ijPairs>>();
//	ArrayList<ijPairs>[][] bins = null;

	
	public void advanceSearchForProductOnly(final int targetIndex, ArrayList<PrimerRecord> leftList, ArrayList<PrimerRecord> rightList) {
		
		  
		
 //		HashMap<Integer, ijPairProductBucket> intCache = new HashMap<Integer, P3OptimzedMultiFinder.ijPairProductBucket>();
		for(int i = 0 ;i < rightList.size();i++)
		{
			PrimerRecord right = rightList.get(i);
			PrimerRecord selectedRight = right;
			if (right instanceof MultiTargetPrimerRecord) {
				selectedRight =  ((MultiTargetPrimerRecord)right).targetsToPrimer.get(scanner.targets.get(targetIndex));
				
			}
			
			if(selectedRight == null || !selectedRight.OK_OR_MUST_USE() )
				continue;
			
			
			for(int j = 0 ; j < leftList.size();j++) {
				PrimerRecord left = leftList.get(j);
				

				
				PrimerRecord selectedleft = left;
				if (left instanceof MultiTargetPrimerRecord) {
					selectedleft =  ((MultiTargetPrimerRecord)left).targetsToPrimer.get(scanner.targets.get(targetIndex));
				
//					// FIXME for debug only
//					if ( ((MultiTargetPrimerRecord)left).getNTarget() != 3 )
//						continue;
				}
				// does not contain this target;
				if(selectedleft == null || !selectedleft.OK_OR_MUST_USE() )
					continue;
				
				
				
				
				int product_size  = selectedRight.start - selectedleft.start+1;
				if (!scanner.pa.chackProductSizeRanges( product_size)) {
					scanner.revals.get(targetIndex).best_pairs.expl.product++;
					if(!selectedleft.must_use  && !selectedRight.must_use){
						scanner.revals.get(targetIndex).best_pairs.expl.product++;
						continue;
					}	
				}
				

				
				PotentialPair newPairInfo = new PotentialPair(scanner, left, right, scanner.pa.getPrPairWeights().primer_quality * 
						(left.quality + right.quality));
				

//				int binId = sizeToBins(product_size);
				
				newPairInfo.productCriterion = product_size;
				newPairInfo.targetIndex = targetIndex;
				
				globalCache.broadcast(newPairInfo);

//				if(!calcPairsBucket.containsKey(product_size))
//				{
//					ijPairProductBucket newBucket = new ijPairProductBucket(scanner.nTargets);
//					newBucket.productSize = product_size;
////					newBucket.calcPairs[targetIndex].productSize =  product_size;
//					calcPairsBucket.put(product_size, newBucket);
//					
//				}
//				if(!intCache[targetIndex].containsKey(product_size))
				{
//					calcPairsBucket.get(product_size).calcPairs[targetIndex].productSize = product_size;
//					productCalcPairs[targetIndex].add(calcPairsBucket.get(product_size).calcPairs[targetIndex]);
//					intCache[targetIndex].put(product_size, new PotentialPairSet(product_size) );
				}
//				calcPairs[targetIndex].add( newPairInfo);
//				intCache[targetIndex].get(product_size).add(newPairInfo);
//				calcPairsBucket.get(product_size).calcPairs[targetIndex].add(newPairInfo);

			}
			
		}
		
		
		

	}
	public boolean hasMore() {
		return ! ( searchExausted && localSearchmList.isEmpty());
	}


	public MultiplexResult getFinalResult() {
		scanner.result.sort();
		return scanner.result;
	}


	public boolean hasGoodSets() {
		return scanner.result.hasGoodSets();
	}



	
	
}
