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
import java.util.TreeSet;

@Deprecated 
public class ProductsLenComp {

	TreeSet<Integer> lens;
	int k;
	IterR rIt ;
	Integer[] current;
	public ProductsLenComp(TreeSet<Integer> lens, int k) {
		this.lens = lens;
		this.k = k;

		rIt = new IterR(lens, k, 0);
//		Iterator<Integer>[] its = new Iterator[k];
		current = new Integer[k];
		if(!rIt.init(current))
			hasMore = false;
	}

	class IterR {
		TreeSet<Integer> lens;
		int kLevel ;
		Iterator<Integer> it;
		IterR childNext ;
		int k;
		public IterR( TreeSet<Integer> lens,int k, int level)
		{
			this.k = k;
			this.lens = lens;
			this.kLevel = level;
//			it = lens.iterator();
			
			if( kLevel <  k - 1  )
			{
				childNext = new IterR( lens,k,kLevel +1 );
			}
		}
		
		
		// each time a it is moving next it need to rest its child
		boolean reset(Integer[] current) {
			if(kLevel == 0) {
				it = lens.iterator();
			}
			else
			{
				// TODO :: offset here is assuming that the list is sorted based on the len :: otherwise it will fail
				// TODO :: ignore offset and use checkConflict to order based on score for greedy selection ???
				// 
				it = lens.tailSet(offset(current[kLevel-1])).iterator();
			}
			if(!it.hasNext())
				return false;
			current[kLevel] = it.next();
//			if(kLevel != 0){
//				while(checkConflict(current, kLevel)!=-1)
//				{
//					if(it.hasNext())
//					{
//						current[this.kLevel] = it.next();
//					}
//					else
//					{
//						return false;
//					}
//				}
//			}
			if(childNext!= null)
				return childNext.reset(current);
			return true;
		}
		
		private Integer offset(Integer len) {
			if(len < 500  )
				return len + 45;
			else if(len < 850)
				return len + 90;
			return len + 180;
		}


		public boolean init(Integer[] current)
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
	
	
		public boolean move(Integer[] current)
		{
			
			if(childNext == null || !childNext.move(current))
			{
				if(it.hasNext()) {
					current[this.kLevel] = it.next();
					if(childNext!= null)
						return childNext.reset(current);
					return true;
				}
				this.reset(current);
				return false;
				
			}
			
			return true;
		}
		
	
	}
	
	
	boolean hasMore = true;
	public ArrayList<Integer[]> getMoreComp(int n) {
		ArrayList<Integer[]> res = new ArrayList<Integer[]>();
		for (int i = 0; i < n; i++) {
			if(hasMore)
			{
				res.add(current.clone());
				if( !rIt.move(current))
					hasMore = false;
			}

		}
		return res;
	}
	public Integer[] getComp() {
	
			if(hasMore)
			{
				Integer[] res =  current.clone();
				if( !rIt.move(current))
					hasMore = false;
				return res;
			}
		return null;
	}
	
	private int checkConflict(Integer[] current, int i) {
		
			for(int j = i-1 ; j >= 0;j-- )
			{
				if(!checkSetProductsLen(current[i],current[j]))
				{
					return i;
				}
			}	
		return -1;
	}
	private boolean checkSetProductsLen(int productLen1 , int productLen2) {
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
	
	
	
	
	static  public void main (String[] args)
	{
		TreeSet<Integer> lens = new TreeSet<Integer>();
		for (int i = 100; i <= 300; i+=1) {
			lens.add(i);
		}
		
		int k =4;
		ProductsLenComp t = new ProductsLenComp(lens,k);
	
		int p = 0;
		for (Integer[] current : t.getMoreComp(10))
		{
			if(p == 0 )
			{
				for (int i = 0; i < k; i++) {
					System.out.print(current[i] +" ");
				}
				System.out.print("\n");
				p=0;
			}
			else
				p--;

		}
		System.err.println("Done");
		
	}
}
