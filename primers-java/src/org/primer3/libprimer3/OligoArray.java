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
import java.util.Comparator;
import java.util.List;

import org.primer3.primer.PrimerRecord;

public class OligoArray {

	 /* Array of oligo (primer) records. */
	 public List<PrimerRecord> oligo = new ArrayList<PrimerRecord>();

	 /* Number of initialized elements */
	 public int num_elem = 0;

	 /* Storage lengths of oligo */
//	 int storage_size = LibPrimer3.INITIAL_LIST_LEN;

	 /* Type of oligos in the array */
	 OligoType type;

	 /* Primers statistics. */
	 public OligoStats expl = new OligoStats();

	 
	// this is add here because it is been used a lot ?? : 
	public int extreme;


	public OligoArray(OligoType otType) {
		this.type = otType;
//		if(type == oligo_type.OT_LEFT)
//		{
//			extreme = Integer.MAX_VALUE;
//		} else if (type == oligo_type.OT_RIGHT) {
//			extreme = 0;
//		}
//		else
//		{
//			extreme = 0;
//		}
		if (type == OligoType.OT_RIGHT) {
			extreme = 0;
		}
		else
		{
			extreme = Integer.MAX_VALUE;
		}
		
		
	}

	/**
	 * add is updated to set and update extreme as new recorded are added to the array
	 * @param h
	 */
	public void add_oligo_to_oligo_array(PrimerRecord h) {
		
		oligo.add(h);
		num_elem++;
		
		// update extreme while adding
		/* Update the most extreme primer variable */
        if (( h.start < this.extreme) && (this.type != OligoType.OT_RIGHT))
        	this.extreme = h.start;
        /* Update the most extreme primer variable */
        if ((h.start > this.extreme) && (this.type == OligoType.OT_RIGHT))
        	this.extreme = h.start;
	}

	
	public String p3_get_oligo_array_explain_string() {
		return expl.p3_oligo_explain_string();
	}

	public void sort_primer_array() {
		oligo.sort(new Comparator<PrimerRecord>() {

			@Override
			public int compare(PrimerRecord a1, PrimerRecord a2) {
				if(a1.quality < a2.quality) return -1;
				if (a1.quality > a2.quality) return 1;

				/*
				 * We want primer_rec_comp to always return a non-0 result, because
				 * different implementations of qsort differ in how they treat "equal"
				 * elements, making it difficult to compare test output on different
				 * systems.
				 */
				if(a1.start > a2.start) return -1;
				if(a1.start < a2.start) return 1;
				if(a1.length < a2.length) return -1;
				return 1;
			}
		});
		
		
		
	}

	 
	 

}