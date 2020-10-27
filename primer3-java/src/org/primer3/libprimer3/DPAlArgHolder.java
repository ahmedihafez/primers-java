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

import org.primer3.dpal.DPAlignmentArgs;
import org.primer3.dpal.DPAlignment;

public class DPAlArgHolder {
	public DPAlignmentArgs local;
	public DPAlignmentArgs end;
	public DPAlignmentArgs local_end;
	public DPAlignmentArgs local_ambig;
	public DPAlignmentArgs local_end_ambig;
	
//	PRIMER_LEFT_0_SEQUENCE= GTGAAGCCTCAGGTAGTGCA
//	PRIMER_RIGHT_0_SEQUENCE= CTCTGTCGACTTTGCCACCA
	
	/* Create the dpal arg holder */
	public static DPAlArgHolder create_dpal_arg_holder (){
		DPAlArgHolder h = new DPAlArgHolder();
		
		
		h.local = new DPAlignmentArgs();
		h.local.set_dpal_args();
		h.local.flag = DPAlignment.DPAL_LOCAL;
		
		
		h.end = new DPAlignmentArgs();
		h.end.set_dpal_args();
		h.end.flag = DPAlignment.DPAL_GLOBAL_END;
		
		h.local_end = new DPAlignmentArgs();
		h.local_end.set_dpal_args();
		h.local_end.flag = DPAlignment.DPAL_LOCAL_END;
		
		h.local_ambig = new DPAlignmentArgs();
//		*h->local_ambig = *h->local;
		h.local_ambig.set_dpal_args();
		h.local_ambig.flag = DPAlignment.DPAL_LOCAL;
		h.local_ambig.dpal_set_ambiguity_code_matrix();
		
		h.local_end_ambig = new DPAlignmentArgs();
//		*h->local_end_ambig = *h->local_end;
		h.local_end_ambig.set_dpal_args();
		h.local_end_ambig.flag = DPAlignment.DPAL_LOCAL_END;
		h.local_end_ambig.dpal_set_ambiguity_code_matrix();
		return h;
	}
}
