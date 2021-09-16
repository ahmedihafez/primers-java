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

import org.primer3.libprimer3.DPAlArgHolder;
import org.primer3.libprimer3.THAlArgHolder;
import org.primer3.primer.MultiTargetPrimerRecord;
import org.primer3.primer.PrimerPair;
import org.primer3.primer.PrimerRecord;

/**
 * Potential pair
 * @author Ahmed Hafez
 *
 */
public class PotentialPair {	
	PrimerRecord left,  right;
	MultiTargetScanner scanner;
	public int targetIndex;
	double estimatedQuilty = Double.MAX_VALUE;
	public PotentialPair(MultiTargetScanner scanner,PrimerRecord left,  PrimerRecord right , double q) {
//		this.i = i;
//		this.j = j;
		
		
		
		this.left = left;
		this.right =  right;
		this.estimatedQuilty = q;
		this.scanner = scanner;
	}
	
//	public double quality = Double.MAX_VALUE;
	public PrimerPair pairRecord =  null;
	public int pairStatus = PrimerPair.PAIR_UNCHARACTRIZED;
	
	
	/**
	 * pair critierial to select this could be the product len or the tm of the Product
	 */
	public int productCriterion;
	
	
	@Override
	public String toString() {
		return String.format("%.2f - %d " + 
	(pairStatus == PrimerPair.PAIR_FAILED ? "PAIR_FAILED" : pairStatus == PrimerPair.PAIR_OK ? "PAIR_OK" : "" )  , estimatedQuilty,productCriterion );
	}
	PotentialPairSet potentialPairParentSet;
	public void setParentSet(PotentialPairSet calcPairs) {
		this.potentialPairParentSet = calcPairs;
		
	}
	public void checkpair(	DPAlArgHolder dpal_arg_to_use, THAlArgHolder thal_arg_to_use) throws Exception {
		
		if(this.pairStatus != PrimerPair.PAIR_UNCHARACTRIZED)
			return;
		if (!left.OK_OR_MUST_USE() || !right.OK_OR_MUST_USE()) {
			this.pairStatus = PrimerPair.PAIR_FAILED;
			this.potentialPairParentSet.unValided(this);
			return;
		}

		
		PrimerRecord selectedLeft = left;
		PrimerRecord selectedRight = right;
		if (left instanceof MultiTargetPrimerRecord) {
			selectedLeft =  ((MultiTargetPrimerRecord)left).targetsToPrimer.get(scanner.targets.get(targetIndex));
			
		}
		if (right instanceof MultiTargetPrimerRecord) {
			selectedRight =  ((MultiTargetPrimerRecord)right).targetsToPrimer.get(scanner.targets.get(targetIndex));
			
		}

		pairRecord = new PrimerPair();
		this.pairStatus = pairRecord.characterize_pair(scanner.revals.get(targetIndex), selectedLeft  , selectedRight, dpal_arg_to_use, thal_arg_to_use, true);
		if (pairStatus == PrimerPair.PAIR_OK) {
			scanner.revals.get(targetIndex).best_pairs.expl.ok++;
			pairRecord.obj_fn(scanner.pa);
			scanner.revals.get(targetIndex).best_pairs.add_pair(pairRecord);
		}
		else
		{
//			System.err.println(selectedLeft.p3_get_ol_problem_string() +"/" +  selectedRight.p3_get_ol_problem_string());
			left.oligoProblems = selectedLeft.oligoProblems;
			right.oligoProblems = selectedRight.oligoProblems;
			this.potentialPairParentSet.unValided(this);
		}
			
	}
	public String getTargetName() {
		return scanner.targets.get(targetIndex);
	}
}
