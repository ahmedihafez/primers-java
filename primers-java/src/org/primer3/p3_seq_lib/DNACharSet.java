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
package org.primer3.p3_seq_lib;

import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.template.AbstractNucleotideCompoundSet;

public class DNACharSet extends AbstractNucleotideCompoundSet<NucleotideCompound> {

	

	public DNACharSet() {
		addNucleotideCompound("A", "T");
		addNucleotideCompound("C", "G");
		addNucleotideCompound("G", "C");
		addNucleotideCompound("T", "A");
		addNucleotideCompound("U", "A");

		addNucleotideCompound("B", "V");
		addNucleotideCompound("D", "H");
		addNucleotideCompound("H", "D");
		addNucleotideCompound("V", "B");
		addNucleotideCompound("R", "Y");
		addNucleotideCompound("Y", "R");
		addNucleotideCompound("K", "M");
		addNucleotideCompound("M", "K");
		addNucleotideCompound("S", "S");
		addNucleotideCompound("W", "W");
		
		addNucleotideCompound("N", "N");
		
		addNucleotideCompound("-", "-");
	}

	@Override
public NucleotideCompound newNucleotideCompound(String base, String complement, String... equivalents) {
		if(equivalents.length == 0) {
			return new NucleotideCompound(base, this, complement);
		}
		else {
			NucleotideCompound[] compounds = new NucleotideCompound[equivalents.length];
			for(int i=0; i<compounds.length; i++) {
				compounds[i] = getCompoundForString(equivalents[i]);
			}
			return new NucleotideCompound(base, this, complement, compounds);
		}
	}
}