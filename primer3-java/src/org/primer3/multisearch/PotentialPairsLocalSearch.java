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

import java.util.HashMap;
import java.util.TreeSet;

import org.primer3.multisearch.PotentialPair;
import org.primer3.primer.MultiTargetPrimerRecord;
import org.primer3.libprimer3.P3GlobalSettings;
import org.primer3.search.ConflictResolution;
//import org.primer3.multisearch.P3OptimzedMultiTargetFinder.ijPairProductBucket;
//import org.primer3.multisearch.P3OptimzedMultiTargetFinder.mPairs;
//import org.primer3.primer.PrimerPair;
import org.primer3.search.MultiListIterator;

/**
 * search all primers combination between different potential pairSet of acceptable product criterion 
 * @author Ahmed Hafez
 *
 */
public class PotentialPairsLocalSearch {
	PotentialPairSet[] calcPairs; 

	public MultiListIterator<PotentialPair> pairComp ;
	
//	public MLocalSearch(Integer[] pLenComp, HashMap<Integer, ijPairProductBucket> calcPairsBucket) {
//		calcPairs = new CalcPairs[pLenComp.length];
//		for (int i = 0; i < calcPairs.length; i++) {
//			calcPairs[i] = calcPairsBucket.get(pLenComp[i]).calcPairs[i]; 
//		}
//		
//		pairComp = new ProductsLenComp2<ijPair>(calcPairs, new ConflictResolution<ijPair>() {
//			@Override
//			public boolean currentValid(ijPair o) {
//				
//				return o.pairStatus != 0;
//			}
//			@Override
//			public boolean checkConflict(ijPair o1, ijPair o2) {
//				if(o1.pairStatus == 0 || o2.pairStatus == 0)
//					return false;
//				{
//					
//				}
//				return true;
//			}
//		});
//	}

	public PotentialPairsLocalSearch(PotentialPairSet[] pLenComp) {
		this.calcPairs = pLenComp;
		pairComp = new MultiListIterator<PotentialPair>(calcPairs, new ConflictResolution<PotentialPair>() {

			@Override
			public boolean checkConflict(PotentialPair o1, PotentialPair o2) {
				//
				if(o1.pairStatus == 0 || o2.pairStatus == 0)
					return false;		
				if(o1.pairStatus == 1 && o2.pairStatus == 1) 
				{
					if( tmFail(o1,o2 , o1.scanner.pa)) {
						return false;
					}
				}
				
				
				
				// no need for this
//				if(o1.left instanceof MultiTargetPrimerRecord && o2.left instanceof MultiTargetPrimerRecord)
//				{
////					// both are left
////					// if o1 left has the o2 taget then they should be the same left other wise no
//					MultiTargetPrimerRecord o1Left = (MultiTargetPrimerRecord) o1.left;
//					MultiTargetPrimerRecord o2Left = (MultiTargetPrimerRecord) o2.left;
//					if(     o1Left.targetsToPrimer.containsKey(o2.getTargetName())  || 
//							o2Left.targetsToPrimer.containsKey(o1.getTargetName()))
//					{
//						return o1Left == o2Left ;
//					}
//				}
				
//				if(o1.pairStatus == 0 || o2.pairStatus == 0)
//					return false;
				{
					
				}
				return true;
			}



			@Override
			public boolean currentValid(PotentialPair o) {
				
				return o.pairStatus != 0;
			}



			@Override
			public void setReslotion(int resoltionDiff) {
				// TODO Auto-generated method stub
				
			}
		});
		
		
		lastScore = PotentialPairSet.calcScore(calcPairs);
	}

	public boolean hasMaxSet() {
		
//		boolean isValid = true;
		for (int i = 0; i < calcPairs.length; i++) {
			if (calcPairs[i].validPairs == 0 || calcPairs[i].size() == 0)
			{
				return false;
			}
		}
		
		return true;
	}

	
	// cache new selection
	double lastScore ;
	PotentialPair[] newSet = null;
	public PotentialPair[] getNew() {
		
		newSet = pairComp.getComp();
		return newSet;
	}

	public double getScore() {
		return lastScore;
	}
	@Override
	public String toString() {
		
		return String.format( "Score : %.2f" ,getScore() );
	}

	public void update() {
		pairComp.update();
		lastScore = 0;
		if(newSet != null)
		{
			for (int i = 0; i < newSet.length; i++) {
				if(newSet[i]== null)
				{
					lastScore = Double.MAX_VALUE;
					break;
				}
				lastScore += newSet[i].estimatedQuilty;
			}
		}
		lastScore = Double.max(lastScore,PotentialPairSet.calcScore(calcPairs));
		
	}
	
	
	public static boolean tmFail(PotentialPair pairRecord1, PotentialPair pairRecord2, P3GlobalSettings pa) {
		
		double minTm =  Double.min(pairRecord1.left.temp,pairRecord1.right.temp);
		double maxTm =  Double.max(pairRecord1.left.temp,pairRecord1.right.temp);
		
		minTm =  Double.min(minTm,pairRecord2.left.temp);
		minTm =  Double.min(minTm,pairRecord2.right.temp);

		
		maxTm =  Double.max(maxTm,pairRecord2.left.temp);
		maxTm =  Double.max(maxTm,pairRecord2.right.temp);

		if((maxTm-minTm) > pa.getMaxDiffTm() )
			return true;
		
		return false;
	}

	public boolean hasMore() {
		
		return this.pairComp.hasMore();
	}
	
}
