/*
    This file is part of primer3 porting to java
    Primer3 and the libprimer3 library are free software;
    you can redistribute them and/or modify them under the terms
    of the GNU General Public License as published by the Free
    Software Foundation; either version 2 of the License, or (at
    your option) any later version.

	Original file are part of https://github.com/primer3-org/primer3
	Whitehead Institute for Biomedical Research, Steve Rozen
	(http://purl.com/STEVEROZEN/), Andreas Untergasser and Helen Skaletsky
	All rights reserved.



    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this software (file gpl-2.0.txt in the source
    distribution); if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
package org.primer3.multisearch;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.TreeSet;

import org.primer3.libprimer3.DPAlArgHolder;
import org.primer3.libprimer3.LibPrimer3;
import org.primer3.libprimer3.PCRType;
import org.primer3.libprimer3.THAlArgHolder;
import org.primer3.oligotm.OligoTMCalculator;
import org.primer3.primer.MultiTargetPrimerRecord;
import org.primer3.primer.PrimerPair;
import org.primer3.primer.PrimerRecord;
import org.primer3.search.ConflictResolution;
import org.primer3.search.PotentialPairSetsCombination;
import org.primer3.search.multiplex.MultiplexResult;
import org.primer3.search.multiplex.MultiplexSet;

// Stateful search of primers pairs 
public class P3OptimzedMultiTargetFinder {
		
	
// enum SearchStrategy {
	boolean 
		SpecificLeftSpecificRight = true,
		MultLeftSpecificRight = false,
		SpecificLeftMultiRight = false;
//	}
	boolean RealTimePCR = false;
	int realTimeReslotion = 10;
	DPAlArgHolder dpal_arg_to_use;
	THAlArgHolder thal_arg_to_use;
	THAlArgHolder thal_oligo_arg_to_use;
	// Objects used in the Stateful search 
//	int[][] max_j_seen;
	
	public MultiTargetScanner scanner = null;
	
	
	public P3OptimzedMultiTargetFinder(MultiTargetScanner scanner , DPAlArgHolder dpal_arg_to_use, THAlArgHolder thal_arg_to_use,
			THAlArgHolder thal_oligo_arg_to_use) {
		
		this. dpal_arg_to_use = dpal_arg_to_use;
		this. thal_arg_to_use = thal_arg_to_use;
		this. thal_oligo_arg_to_use = thal_oligo_arg_to_use;
		this.scanner = scanner;
		RealTimePCR = this.scanner.pa.getPcrType() == PCRType.Real_Time_PCR ;
		globalCache = new InternalCache(scanner);
		
	}	

	static double ALPHA = 1;
	
	
	
	

	
	
	
	
	

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
	
	// No need for this
//	public interface TargetPrimerFiller {
//		void fillTargetPrimer(int targetIndex);
//	}
	
	public void initSearch() {
		for(int targetIndex =0 ; targetIndex < scanner.nTargets;targetIndex++  )
		{
			String targetName = scanner.targets.get(targetIndex);
			
			
			this.advanceSearchForProductOnly(targetIndex , scanner.targetsToSpFwd.get(targetName),scanner.targetsToSpRev.get(targetName));
//			this.advanceSearchForProductOnly(targetIndex , scanner.multiFwd,scanner.targetsToSpRev.get(targetName));

			// filler.fillTargetPrimer(targetIndex);
			globalCache.clean(targetIndex);

		}
		
		globalCache.fillMissing();
		
		
		ConflictResolution<PotentialPairSet> conflictResolution = PotentialPairSetsCombination.productLenConflictResolution;
		if(RealTimePCR )
		{
			PotentialPairSetsCombination.productTmConflictResolution.setReslotion(realTimeReslotion);
			conflictResolution = PotentialPairSetsCombination.productTmConflictResolution;
		}
		
		// init values here
		while (true) {
			PotentialPairSetsCollection[] productCalcPairs = globalCache.getNext();
			if(productCalcPairs == null)
				break;
			PotentialPairSetsCombination newComp = new PotentialPairSetsCombination(productCalcPairs, conflictResolution );
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
				mPairs newPairComp = new mPairs(this, newset);
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
					MultiplexSet res = new  MultiplexSet(this);
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
					// scanner.revals.get(targetIndex).best_pairs.expl.product++;
					if(!selectedleft.must_use  && !selectedRight.must_use){
						scanner.revals.get(targetIndex).best_pairs.expl.product++;
						continue;
					}	
				}
				
				double product_tm =  OligoTMCalculator.longSeqTM(scanner.revals.get(targetIndex).sa.getTrimmedSequence(), selectedleft.start,
						selectedRight.start - selectedleft.start + 1,
						/* TO DO -- skewed, it would be better to not use p_args elements here */
						scanner.pa.primersArgs.getSaltConcentration(),
						scanner.pa.primersArgs.getDivalentConcentration(),
						scanner.pa.primersArgs.getDntpConcentration());
				
				if (scanner.pa.getProductMinTM() != LibPrimer3.PR_DEFAULT_PRODUCT_MIN_TM
						&& product_tm < scanner.pa.getProductMinTM()) {
					// scanner.revals.get(targetIndex).best_pairs.expl.low_tm++;
					if(!selectedleft.must_use  && !selectedRight.must_use){
						scanner.revals.get(targetIndex).best_pairs.expl.low_tm++;
						continue;
					}	
				}

				if (scanner.pa.getProductMaxTM() != LibPrimer3.PR_DEFAULT_PRODUCT_MAX_TM
						&& product_tm > scanner.pa.getProductMaxTM()) {
					if(!selectedleft.must_use  && !selectedRight.must_use){
						scanner.revals.get(targetIndex).best_pairs.expl.high_tm++;
						continue;
					}	
				}
				
				
				
				PotentialPair newPairInfo = new PotentialPair(scanner, left, right, scanner.pa.getPrPairWeights().primer_quality * 
						(left.quality + right.quality));
				

//				int binId = sizeToBins(product_size);
				
				if(RealTimePCR)
				{
					newPairInfo.productCriterion = (int)(product_tm*realTimeReslotion);
				}
				else
				{
					newPairInfo.productCriterion = product_size;
				}
				newPairInfo.targetIndex = targetIndex;
				
				globalCache.broadcast(newPairInfo);



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


	public boolean getIsRealTimePCR() {
		return RealTimePCR;
	}


	public int getRealTimePCRResoltion() {
		return realTimeReslotion;
	}



	
	
}
