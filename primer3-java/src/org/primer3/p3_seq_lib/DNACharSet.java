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