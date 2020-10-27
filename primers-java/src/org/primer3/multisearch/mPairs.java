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

import org.primer3.primer.MultiTargetPrimerRecord;
import org.primer3.primer.PrimerPair;

public class mPairs {

	/**
	 * 
	 */
	private final P3OptimzedMultiTargetFinder p3OptimzedMultiTargetFinder;

	PotentialPair[] pairSet;

	int setStatus = PrimerPair.PAIR_UNCHARACTRIZED;

	public mPairs(P3OptimzedMultiTargetFinder p3OptimzedMultiTargetFinder, PotentialPair[] pairSet)
	{
		this.p3OptimzedMultiTargetFinder = p3OptimzedMultiTargetFinder;
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
				pair.checkpair(this.p3OptimzedMultiTargetFinder.dpal_arg_to_use, this.p3OptimzedMultiTargetFinder.thal_arg_to_use);;
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