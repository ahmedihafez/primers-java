/*
    This file is part of primer3 porting to java


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
package org.primer3.search;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.SortedSet;
import java.util.TreeSet;

import org.primer3.multisearch.PotentialPairSet;
import org.primer3.multisearch.PotentialPairSetsCollection;

/**
 * produce a combination between different potential set based on acceptable criteria
 * @author Ahmed Hafez
 *
 */
public class PotentialPairSetsCombination {

	
	static public final ConflictResolution<PotentialPairSet> productLenConflictResolution = new ConflictResolution<PotentialPairSet>() {
		
		@Override
		public boolean currentValid(PotentialPairSet o) {
			// TODO Auto-generated method stub
			return true;
		}
		
		@Override
		public boolean checkConflict(PotentialPairSet o1, PotentialPairSet o2) {
			 
			int productLen1 = o1.productCriterion;
			int productLen2 = o2.productCriterion;
			
				boolean isValid = true;
				
				
				int diff =  Math.abs(productLen1-productLen2);
				
				int minLen = Integer.min(productLen1,productLen2);
				
				int diffMin = 45 ;
				
				
				if(minLen < 500  )
					diffMin = 45;
				else if(minLen < 850)
					diffMin = 90;
				else diffMin = 180;
				if( diff <= diffMin ) // not exactly should be
				{
					isValid = false;
				}
			
			return isValid;
		}

		@Override
		public void setReslotion(int resoltionDiff) {
			// TODO Auto-generated method stub
		}
	} ;
	
	
	
	static public final ConflictResolution<PotentialPairSet> productTmConflictResolution = new ConflictResolution<PotentialPairSet>() {
		int tmResoltionDiff ;
		
		@Override
		public boolean currentValid(PotentialPairSet o) {
			// TODO Auto-generated method stub
			return true;
		}
		
		@Override
		public boolean checkConflict(PotentialPairSet o1, PotentialPairSet o2) {
			int productTm1 = o1.productCriterion;
			int productTm2 = o2.productCriterion;
			
			boolean isValid = true;
					
			int diff =  Math.abs(productTm1-productTm2);
				
				
			if( diff < tmResoltionDiff ) // not exactly should be
			{
				isValid = false;
			}
			
			return isValid;
		}

		@Override
		public void setReslotion(int resoltionDiff) {
			tmResoltionDiff = 	resoltionDiff;	
		}
	} ;
	
	
	ConflictResolution<PotentialPairSet> conflictResolution = null;
	PotentialPairSetsCollection[] productCalcPairs;
	int k;
	IterR rIt ;
	PotentialPairSet[] current = null;
	public PotentialPairSetsCombination(PotentialPairSetsCollection[] productCalcPairs, ConflictResolution<PotentialPairSet> conflictResolution) {
		this.productCalcPairs = productCalcPairs;
		this.k = productCalcPairs.length;
		this.conflictResolution  = conflictResolution;
		rIt = new IterR(productCalcPairs, k, 0);
//		Iterator<Integer>[] its = new Iterator[k];

	}

	class IterR {
		PotentialPairSetsCollection productCalcPairs;
		int kLevel ;
		Iterator<PotentialPairSet> it;
		IterR childNext ;
		int k;
		public IterR( PotentialPairSetsCollection[] productCalcPairs,int k, int level)
		{
			this.k = k;
			this.productCalcPairs = productCalcPairs[level];
			this.kLevel = level;
//			it = this.productCalcPairs.iterator();
			
			if( kLevel <  k - 1  )
			{
				childNext = new IterR( productCalcPairs,k,kLevel +1 );
			}
		
		}
		
		
		// each time a it is moving next it need to rest its child
		boolean reset(PotentialPairSet[] current) {
			if(kLevel == 0) {
				it = productCalcPairs.iterator();
			}
			else
			{
				// TODO :: offset here is assuming that the list is sorted based on the len :: otherwise it will fail
				// TODO :: ignore offset and use checkConflict to order based on score for greedy selection ???
				// 
//				SortedSet<CalcPairs> sub1 =  productCalcPairs.tailSet(offsetLarger(current[kLevel-1].productSize));
//				SortedSet<CalcPairs> sub1 =  productCalcPairs.headSet(offsetLess(current[kLevel-1].productSize));

				it = productCalcPairs.iterator(); // productCalcPairs.tailSet(offset(current[kLevel-1].productSize)).iterator();
			}
			if(!it.hasNext())
				return false;
//			current[kLevel] = it.next();
//			if(kLevel != 0){
				if(!advanceSaveMove(current))
				{
					return false;
				}
//			}
			if(childNext!= null){ 
				while(!childNext.reset(current))
				{
					if(!advanceSaveMove(current))
					{return false;}
				}
			}
			return true;
		}
		
		boolean advanceSaveMove(PotentialPairSet[] current)
		{
			do{
				if(it.hasNext())
				{
					current[this.kLevel] = it.next();
				}
				else
				{
					return false;
				}
			}while(current[kLevel].validPairs == 0|| current[kLevel].size() == 0 || checkConflict(current, kLevel)!=-1);
			return true;
		}
		
		private Integer offsetLarger(Integer len) {
			if(len < 500  )
				return len + 45;
			else if(len < 850)
				return len + 90;
			return len + 180;
		}
		private Integer offsetLess(Integer len) {
			if(len < 500  )
				return len - 45;
			else if(len < 850)
				return len - 90;
			return len - 180;
		}

		public boolean init(PotentialPairSet[] current)
		{
			if(this.reset(current))
			{
//				if(childNext != null)
//				{
//					return childNext.init(current);
//				}
				return true;
			}
			return false;	
			
//			if(it.hasNext()) {
//				current[kLevel] = it.next();
//					if(kLevel > 0 )
//					{
//						// check here and move fwd
//						while(checkConflict(current,kLevel) != -1)
//						{
//							if(it.hasNext())
//								current[kLevel] = it.next();
//							else
//								return false;
//						}
//					}
//				return true;
//			}
//			return false;
		}
	
	
		public boolean move(PotentialPairSet[] current)
		{
			if ( current[kLevel].validPairs == 0 || current[kLevel].size() == 0 ){
				it.remove();
				if(advanceSaveMove(current)) {
//					current[this.kLevel] = it.next();
					if(childNext!= null) {
						while(!childNext.reset(current))
						{
							if(!advanceSaveMove(current)){
								return false;
							}
						}
					}
					return true;
				}
				else
				{
					return false;
				}
			}
			
			
			if(childNext == null || !childNext.move(current))
			{
				if(advanceSaveMove(current)) {
//					current[this.kLevel] = it.next();
					if(childNext!= null)
						while(!childNext.reset(current))
						{
							if(!advanceSaveMove(current)){
								return false;
								}
						}					
					return true;
				}
				this.reset(current);
				return false;
				
			}
			
			return true;
		}
		
	
	}
	
	
	boolean hasMore = true;
//	public ArrayList<PotentialPairSet[]> getMoreComp(int n) {
//		ArrayList<PotentialPairSet[]> res = new ArrayList<PotentialPairSet[]>();
//		for (int i = 0; i < n; i++) {
//			if(hasMore)
//			{
//				res.add(current.clone());
//				if( !rIt.move(current))
//					hasMore = false;
//			}
//
//		}
//		return res;
//	}
	public PotentialPairSet[] getComp() {
	
			if(current == null)
			{
				// TODO :: search can not return a valid move FIXME in init
				current = new PotentialPairSet[k];
				if(!rIt.init(current)) {
					hasMore = false;
					return null;
				}
					
				return current.clone();
			}
			if(hasMore)
			{
				
				if(rIt.move(current))
				{
					PotentialPairSet[] res =  current.clone();
					return res;
				}
				hasMore = false;
				
			}
		return null;
	}
	
	private int checkConflict(PotentialPairSet[] current, int i) {
		
			for(int j = i-1 ; j >= 0;j-- )
			{
//				if(!checkSetProductsLen(current[i].productCriterion,current[j].productCriterion))
				if(! conflictResolution.checkConflict(current[i],current[j]))
				{
					return i;
				}
			}	
		return -1;
	}
//	private boolean checkSetProductsLen(int productLen1 , int productLen2) {
//		    boolean isValid = true;
//		
//		
//			int diff =  Math.abs(productLen1-productLen2);
//			
//			int minLen = Integer.min(productLen1,productLen2);
//			
//			int diffMin = 45 ;
//			
//			
//			if(minLen < 500  )
//				diffMin = 45;
//			else if(minLen < 850)
//				diffMin = 90;
//			else diffMin = 180;
//			if( diff <= diffMin ) // not exactly should be
//			{
//				isValid = false;
//			}
//		
//		return isValid;
//		
//	}
	
	
	
	
//	static  public void main (String[] args)
//	{
////		TreeSet<Integer> lens = new TreeSet<Integer>();
////		for (int i = 100; i <= 300; i+=1) {
////			lens.add(i);
////		}
////		
////		int k =4;
////		ProductsLenComp3 t = new ProductsLenComp3(lens,k);
////	
////		int p = 0;
////		for (Integer[] current : t.getMoreComp(10))
////		{
////			if(p == 0 )
////			{
////				for (int i = 0; i < k; i++) {
////					System.out.print(current[i] +" ");
////				}
////				System.out.print("\n");
////				p=0;
////			}
////			else
////				p--;
////
////		}
//		System.err.println("Done");
//		
//	}
	public void update(PotentialPairSet[] current) {
		
		for (int i = 0; i < current.length; i++) {
			 current[i].deleteUnwanted();
		}
		
	}
	public boolean hasMore() {
		
		return hasMore;
	}
	public double getEstimatedScore() {
		double score = 0;
		if(current != null)
			return PotentialPairSet.getOrginalScore( current);
		for (int i = 0; i < productCalcPairs.length; i++) {
			if(productCalcPairs[i].size() == 0 )
				return Double.MAX_VALUE;
			score += productCalcPairs[i].getEstimatedScore();
		}
		return score;
	}
}
