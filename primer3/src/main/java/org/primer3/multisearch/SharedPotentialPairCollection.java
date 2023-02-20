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