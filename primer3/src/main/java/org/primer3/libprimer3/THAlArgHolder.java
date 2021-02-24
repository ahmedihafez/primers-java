/*
    This file is part of Primer3 porting to java (https://github.com/primer3-org/primer3)


	Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008
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

import org.primer3.thal.ThermodynamicAlignmentType;
import org.primer3.thal.ThermodynamicAlignmentArguments;

public class THAlArgHolder {
	
	public ThermodynamicAlignmentArguments any;
	public ThermodynamicAlignmentArguments end1;
	public ThermodynamicAlignmentArguments end2;
	public ThermodynamicAlignmentArguments hairpin_th;
	
	
	
	
	static public THAlArgHolder   create_thal_arg_holder (PrimersOligosArguments po_args) {
		THAlArgHolder h = new THAlArgHolder();
		
		h.any = new  ThermodynamicAlignmentArguments();
		h.any.setThAlDefaultArgs();
		h.any.setAlignmentType(ThermodynamicAlignmentType.thal_any);
		h.any.setMonovalentConc(po_args.getSaltConcentration());
		h.any.setDivalentConc(po_args.getDivalentConcentration());
		h.any.setDntpConc(po_args.getDntpConcentration());
		h.any.setDnaConc(po_args.getDnaConcentration());
		
		
		h.end1 = new ThermodynamicAlignmentArguments();
		h.end1.setThAlDefaultArgs();
		h.end1.setAlignmentType(ThermodynamicAlignmentType.thal_end1);
		h.end1.setMonovalentConc(po_args.getSaltConcentration());
		h.end1.setDivalentConc(po_args.getDivalentConcentration());
		h.end1.setDntpConc(po_args.getDntpConcentration());
		h.end1.setDnaConc(po_args.getDnaConcentration());
		
		h.end2 = new ThermodynamicAlignmentArguments();
		h.end2.setThAlDefaultArgs();
		h.end2.setAlignmentType(ThermodynamicAlignmentType.thal_end2);
		h.end2.setMonovalentConc(po_args.getSaltConcentration());
		h.end2.setDivalentConc(po_args.getDivalentConcentration());
		h.end2.setDntpConc(po_args.getDntpConcentration());
		h.end2.setDnaConc(po_args.getDnaConcentration());
		
		h.hairpin_th = new ThermodynamicAlignmentArguments();
		h.hairpin_th.setThAlDefaultArgs();
		h.hairpin_th.setAlignmentType(ThermodynamicAlignmentType.thal_hairpin);
		h.hairpin_th.setMonovalentConc(po_args.getSaltConcentration());
		h.hairpin_th.setDivalentConc(po_args.getDivalentConcentration());
		h.hairpin_th.setDntpConc(po_args.getDntpConcentration());
		h.hairpin_th.setDnaConc(po_args.getDnaConcentration());
		h.hairpin_th.setCalcDimer(0);
		return h;
	}
}
