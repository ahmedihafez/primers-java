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

import java.util.Comparator;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * Store and sort all different potential pairSet by the product criteria and average score of the primers in the list
 * for a specific target sequence
 * @author Ahmed Hafez
 *
 */
public class PotentialPairSetsCollection extends TreeSet<PotentialPairSet> {
	/**
	 * auto generated serial id 
	 */
	private static final long serialVersionUID = -5978203688435479112L;
		public PotentialPairSetsCollection() {
			super(new Comparator<PotentialPairSet>() {

				@Override
				public int compare(PotentialPairSet o1, PotentialPairSet o2) {
					int res =  Integer.compare(o1.productCriterion, o2.productCriterion);
					if (res == 0) 
						return res;
					return Double.compare(o1.getAverageScore(), o2.getAverageScore());
				}
			});
		}
		
		@Override
		public void clear() {
			// TODO Auto-generated method stub
			
			super.clear();
		}

		PotentialPairSet dummyObject = new PotentialPairSet(0);
		public SortedSet<PotentialPairSet> tailSet(Integer offset) {
			dummyObject.productCriterion = offset;
			return super.tailSet(dummyObject);
		}

		public SortedSet<PotentialPairSet> headSet(Integer offset) {
			dummyObject.productCriterion = offset;
			return super.headSet(dummyObject);
		}
		public PotentialPairSetsCollection execludeSet(Integer less, Integer larger)
		{
			dummyObject.productCriterion = larger;
			SortedSet<PotentialPairSet> tailSet = super.tailSet(dummyObject);
			dummyObject.productCriterion = less;
			SortedSet<PotentialPairSet> headSet = super.headSet(dummyObject);
			PotentialPairSetsCollection res = new PotentialPairSetsCollection();
			res.addAll(tailSet);
			res.addAll(headSet);
			return res;
		}

		public double getEstimatedScore() {
			if(size() == 0)
				return Double.MAX_VALUE;
			return first().getAverageScore();
		}
	}