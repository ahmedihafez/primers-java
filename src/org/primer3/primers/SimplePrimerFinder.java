package org.primer3.primers;

import java.util.ArrayList;
import java.util.List;

import javax.swing.text.AbstractDocument.LeafElement;

import org.primer3.libprimer3.LibPrimer3;
import org.primer3.libprimer3.P3GlobalSettings;
import org.primer3.libprimer3.PrimerRecord;
import org.primer3.libprimer3.oligo_type;
import org.primer3.libprimer3.p3_output_type;
import org.primer3.libprimer3.p3retval;
import org.primer3.libprimer3.seq_args;
import org.primer3.oligotm.OligoTMCalculator;
import org.primer3.thal.ThAlParameters;

public class SimplePrimerFinder {
	protected int default_version = 2;
	protected seq_args sarg  = new seq_args();
	protected P3GlobalSettings global_pa = P3GlobalSettings.p3_create_global_settings(default_version);
	
	
	public P3GlobalSettings getGlobal_pa()
	{
		return global_pa;
	}
	
	void init()
	{
		/* Check if any thermodynamical alignment flag was given */
		if ((global_pa.isThermodynamicOligoAlignment() ) || 
				(global_pa.isThermodynamicTemplateAlignment()))
			ThAlParameters.get_thermodynamic_values();
	}
	
	// TODO :: consider this http://www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html
	public void setSequence(char[] seq)
	{
		
	}
	public void setSequence(String seq,String seqName)
	{
		sarg.p3_set_sa_sequence(seq);
		sarg.p3_set_sa_sequence_name(seqName);
	}
	
	public String getSequence()
	{
		return new String( sarg.sequence);
	}
	
	
	


	

	
	

	public static void main(String[] args)
	{
		SimplePrimerFinder primerFinder = new SimplePrimerFinder();
		primerFinder.setSequence("tattggtgaagcctcaggtagtgcagaatatgaaacttcaggatccagtgggcatgctactggtagtgctgccggccttacaggcattatggtggcaaagtcgacagagttta","template 1");
		
	}


	

}
