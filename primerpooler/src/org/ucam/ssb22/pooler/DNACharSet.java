/*
  This file is part of Java porting of Primer Pooler (https://github.com/ssb22/PrimerPooler)
  Primer Pooler (c) Silas S. Brown.  For Wen.
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
    http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
package org.ucam.ssb22.pooler;

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