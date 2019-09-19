package org.pooler;

import java.util.ArrayList;
import java.util.List;

import org.pooler.Fasta2Bit.SeqData;
import org.pooler.Fasta2Bit.SeqInfo;

public class Genome {

	public static List<SeqInfo> go_throgh_genome(Amplicons amplicons) {
		// TODO :: add read from fasta file
//		if(is_fasta(f)) return fasta_genome(f);
		try {
			// return fasta_genome(genomeFile);
		}
		catch(Exception ex) {
			// it is not fasta file
		}
		
		// try 2bit files
		Fasta2Bit genome2bit;
		ArrayList<SeqInfo> processedSeq = new ArrayList<>();
		try {
			genome2bit = new Fasta2Bit(amplicons.genomeFile);
			amplicons.allocateSeqs(genome2bit.nSeq);
			for(int seqNo=0;seqNo<genome2bit.nSeq;seqNo++) {
				SeqData seqData = genome2bit.parseNextSeq();
				// TODO :: please revise
				if (seqData.seq.seqname.contains("-") || seqData.seq.seqname.contains("_")) {
					continue;
				}
				else
				{
//					System.out.println(seqData.seq.seqname);
					
				}
				processedSeq.add(seqData.seq);
				processSequence(seqData,amplicons);
				// break;
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		

		return processedSeq;
	}

	static void processSequence(SeqData seqData, Amplicons amplicons) {
//		if( no problems ) {
//		if variant then skip
//		{
//			
//		}
//	else {
//			add this seq name to the retunr result
//		}
//	}
	
//	 if(isVariant || !seqNames) continue;
		
		while(seqData.next()) {
			

			amplicons.lookup(seqData.currentBuf,seqData.currentValid,seqData);
//			perform lookup here
			
//			 System.out.println( Long.toBinaryString(seqData.currentBuf));
			 //System.out.println( Long.toBinaryString(seqData.currentValid));
		}
	}
	
}
