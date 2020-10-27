/*
    This file is part of primer3 porting to java


	Original file are part of https://github.com/primer3-org/primer3
	Whitehead Institute for Biomedical Research, Steve Rozen
	(http://purl.com/STEVEROZEN/), Andreas Untergasser and Helen Skaletsky
	All rights reserved to Primer3 authors.

    Primer3 and the libprimer3 library are free software;
    you can redistribute them and/or modify them under the terms
    of the GNU General Public License as published by the Free
    Software Foundation; either version 2 of the License, or (at
    your option) any later version.

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
package org.primer3.libprimer3;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

public class IntervalList {

	
	//	org.apache.commons.lang3.tuple.MutablePair<L, R>
	// refactor to a list array
//	private int[][] pairs = new int[LibPrimer3.PR_MAX_INTERVAL_ARRAY][2];

	private List<int[]> pairs = new ArrayList<int[]>();
	// 
	
	public List<int[]> getPairs() {
		return pairs;
	}
	
//	private int count = 0;
	public boolean oligoOverlapsInterval(int start, int len) {

		int i;
		int last = start + len - 1;
		for (i = 0; i <  pairs.size(); i++)
			if (!(last < pairs.get(i)[0]
					|| start > (pairs.get(i)[0] + pairs.get(i)[1] - 1)))
				return true;
		return false;
	}

	/* Add a pair of integers to an  array of intervals */
	public boolean addInterval( int i1, int i2)
	{
//		int c = this.count;
//		if (c >= LibPrimer3.PR_MAX_INTERVAL_ARRAY) return true;
		
		int[] pair = new int[]{i1,i2};
		
		this.pairs.add(pair);
		
//		this.pairs[c][0] = i1;
//		this.pairs[c][1] = i2;
//		this.count++;
		return  false;
	}


	static final String numSep = ",";
	static final String intervalSep = " ";
	
	
	public IntervalList append_interval(String datum) {
		
//		IntervalArrayT2 tar2 = new IntervalArrayT2();
		IntervalList tar2 = this;


		String[] intervalStrs = datum.split(Pattern.quote(intervalSep));
		for(String intervalStr : intervalStrs){
			String[] intervalNums = intervalStr.split(Pattern.quote(numSep));
			int i1,i2;
			i1 = i2  = -1;
			if(!intervalNums[0].isEmpty())
				i1 = Integer.parseInt(intervalNums[0]);
			if(!intervalNums[1].isEmpty())
				i2 = Integer.parseInt(intervalNums[1]);

			if(tar2.addInterval(i1, i2))
				throw new IndexOutOfBoundsException("Too many elements for tag ");
		}

		return tar2; 
	}


	/**
	 * Check intervals, and add any errors to err.
	 * Update the start of each interval to
	 * be relative to the start of the included region.
	 */
	public boolean checkAndAdjustInterval(String tag_name, int seq_len,
			int first_index, StringBuilder err, SeqArgs sa,
			StringBuilder warning, boolean empty_allowed) {
		return checkAndAdjustInterval(tag_name,this.pairs.size(),this.pairs,seq_len,first_index,err,sa,warning,empty_allowed);
	}

	/**
	 * Check intervals, and add any errors to err.
	 * Update the start of each interval to
	 * be relative to the start of the included region.
	 */
	public static boolean checkAndAdjustInterval(String tag_name,
			int count,
			int[][] pairs,
			int seq_len,
			int first_index, StringBuilder err, SeqArgs sa,
			StringBuilder warning, boolean empty_allowed) {
		int i;
		boolean outside_warning_issued = false;

		/* Subtract the first_index from the start positions in the interval
		     array */
		for (i = 0; i < count; i++) {
			if (empty_allowed && (pairs[i][0] == -1) && (pairs[i][1] == -1))
				continue;
			if (empty_allowed && (((pairs[i][0] == -1) && (pairs[i][1] != -1)) 
					|| ((pairs[i][0] != -1) && (pairs[i][1] == -1)))) {
				err.append(tag_name + " illegal interval");
				return true;
			}
			pairs[i][0] -= first_index;
		}

		for (i=0; i < count; i++) {
			if (empty_allowed && (pairs[i][0] == -1) && (pairs[i][1] == -1))
				continue;
			if (pairs[i][0] + pairs[i][1] > seq_len) {
				err.append(tag_name + " beyond end of sequence");
				return true;
			}
			/* Cause the interval start to be relative to the included region. */
			pairs[i][0] -= sa.getIncludedRegionStart();
			/* Check that intervals are within the included region. */
			if (pairs[i][0] < 0
					|| pairs[i][0] + pairs[i][1] > sa.getIncludedRegionLength()) {
				if (!outside_warning_issued) {
					warning.append(tag_name +  " outside of INCLUDED_REGION");
					outside_warning_issued = true;
				}
			}
			if (pairs[i][1] < 0) {
				err.append( "Negative " + tag_name + " length");
				return true;
			}
		}
		return false;
	}

	/**
	 * Check intervals, and add any errors to err.
	 * Update the start of each interval to
	 * be relative to the start of the included region.
	 */
	public static boolean checkAndAdjustInterval(String tag_name,
			int count,
			List<int[]> pairs,
			int seq_len,
			int first_index, StringBuilder err, SeqArgs sa,
			StringBuilder warning, boolean empty_allowed) {
		int i;
		boolean outside_warning_issued = false;

		/* Subtract the first_index from the start positions in the interval
		     array */
		for (i = 0; i < count; i++) {
			if (empty_allowed && (pairs.get(i)[0] == -1) && (pairs.get(i)[1] == -1))
				continue;
			if (empty_allowed && (((pairs.get(i)[0] == -1) && (pairs.get(i)[1] != -1)) 
					|| ((pairs.get(i)[0] != -1) && (pairs.get(i)[1] == -1)))) {
				err.append(tag_name + " illegal interval");
				return true;
			}
			pairs.get(i)[0] -= first_index;
		}

		for (i=0; i < count; i++) {
			if (empty_allowed && (pairs.get(i)[0] == -1) && (pairs.get(i)[1] == -1))
				continue;
			if (pairs.get(i)[0] + pairs.get(i)[1] > seq_len) {
				err.append(tag_name + " beyond end of sequence");
				return true;
			}
			/* Cause the interval start to be relative to the included region. */
			pairs.get(i)[0] -= sa.getIncludedRegionStart();
			/* Check that intervals are within the included region. */
			if (pairs.get(i)[0] < 0
					|| pairs.get(i)[0] + pairs.get(i)[1] > sa.getIncludedRegionLength()) {
				if (!outside_warning_issued) {
					warning.append(tag_name +  " outside of INCLUDED_REGION");
					outside_warning_issued = true;
				}
			}
			if (pairs.get(i)[1] < 0) {
				err.append( "Negative " + tag_name + " length");
				return true;
			}
		}
		return false;
	}
	
	public int getCount() {
		return pairs.size();
	}

	public int[] getInterval(int i) {
		return pairs.get(i);
	}

	public void addInterval(int[] interval) {
		this.pairs.add(interval);
	}


	public void remove(int i) {
		this.pairs.remove(i);
	}

}