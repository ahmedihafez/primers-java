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
package org.primer3.search.multiplex;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

import org.primer3.libprimer3.P3GlobalSettings;
import org.primer3.libprimer3.P3RetVal;
import org.primer3.multisearch.IP3MultiSearch;
import org.primer3.multisearch.P3MultiplexSearch;

public class MultiplexResult {
	String groupName;
	P3GlobalSettings pa;
	IP3MultiSearch p3MultiplexSearcher = null;
	// add just once
	HashMap<String,P3RetVal> sourcePairsResult = new HashMap<String, P3RetVal>();
	
	
	ArrayList<MultiplexSet> multiplexSets = new ArrayList<MultiplexSet>();
	
	
	
	
	public MultiplexResult(IP3MultiSearch p3MultiplexSearch,String groupName,  P3GlobalSettings pa) {
		p3MultiplexSearcher = p3MultiplexSearch;
		this.pa = pa;
		this.groupName = groupName;
	}
	public void addP3RetVal(P3RetVal retVal) {
		
		// just add once
		String targetName = retVal.sa.getSequenceName();
		if(!sourcePairsResult.containsKey(targetName))
		{
			sourcePairsResult.put(targetName, retVal);
			expectedMaxNTargetsInSet++;
		}
		
	}
	
	
	
	ArrayList<MultiplexSet> newAltsSets =  new ArrayList<MultiplexSet>();
	
	int numOfSetToReturn = 10;
	int numOfGoodSets = 0;
	int expectedMaxNTargetsInSet=0;
	public boolean hasGoodSets() {
		numOfGoodSets = 0;
		for(MultiplexSet set : multiplexSets)
		{
			if(set.getNTargets() >= expectedMaxNTargetsInSet)
			{
				numOfGoodSets++;
			}
		}
		return numOfGoodSets > numOfSetToReturn;
	}

	
	
	




	public static void main(String[] args)
	{
		HashMap<String, String> t = new HashMap<String, String>();
		t.put("a", "a");
		t.put("b", "c");
		t.put("c", "b");
		
		ArrayList<String> a = new ArrayList<String>(t.keySet());
		
		System.out.println(a);
		a.clear();
		System.out.println(a);
		System.out.println(t);
		
	}



	

	public void print_boulder(int io_version) {
		System.out.println("MULTIPLEX_GROUP=" + this.groupName);
		for(int i = 0 ; i < multiplexSets.size();i++) {
			System.out.println("MULTIPLEX_SET=" + i);
			multiplexSets.get(i).print_buolder(io_version);
				
		}
		System.out.println("MULTIPLEX_GROUP_END");

	}




	public void sort() {
		
		multiplexSets.sort(new Comparator<MultiplexSet>() {

			@Override
			public int compare(MultiplexSet o1, MultiplexSet o2) {
				
				if(o1.getNTargets() == o2.getNTargets())
					return Double.compare(o1.getScore(), o2.getScore());
				return -Integer.compare(o1.getNTargets(), o2.getNTargets());
			}
		});
		
	}
	
	



	public List<MultiplexSet> getMultiplexSet() {
		return multiplexSets;
	}
	

	public void addSet(MultiplexSet newSet) {
		multiplexSets.add(newSet);
	}

}
