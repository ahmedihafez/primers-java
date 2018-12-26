package org.primer3.p3_seq_lib;

import java.io.File;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.regex.Pattern;

import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.io.DNASequenceCreator;
import org.biojava.nbio.core.sequence.io.FastaReader;
import org.biojava.nbio.core.sequence.io.GenericFastaHeaderParser;

// TODO :: missing impl...
public class seq_lib {


	
	static double  PR_MAX_LIBRARY_WT = 100.0;

	
	/** An array of sequence names. */
	private List<String> names;
	/** An array of sequences. */
	private List<char[]> seqs;
	
	/** An array of reversed-complemented sequences.
	   x->rev_compl_seqs[i] is the reverse complement
	   of x->seqs[i], which lets us keep track of pairwise
	   mispriming.  See reverse_complement_seq_lib(). */
	private List<char[]> rev_compl_seqs;
	
	public List<Double> weight;
	
	
	StringBuilder warning = new StringBuilder();
	StringBuilder error = new StringBuilder();

	
	
	public seq_lib()
	{
		names = new ArrayList<String>();
		seqs = new ArrayList<char[]>();
		rev_compl_seqs = new ArrayList<char[]>();
		weight = new ArrayList<Double>();
	}
	
	
	
	/*  
	 * Reads any file in fasta format and returns a newly allocated
	 * seq_lib, lib.  Sets lib.error to a non-empty string on any error
	 * other than ENOMEM.  Returns NULL on ENOMEM.
	 */
	public static seq_lib read_and_create_seq_lib(String filename, String errfrag) throws Exception
	{
		seq_lib slib = new seq_lib();
		
		DNACharSet dnaSet = new DNACharSet();
		
		
		
		
		
		FastaReader<DNASequence, NucleotideCompound> fastaReader = new FastaReader<DNASequence, NucleotideCompound>(
				new FileInputStream(new File(filename)),
				new GenericFastaHeaderParser<DNASequence, NucleotideCompound>(),
				new DNASequenceCreator(dnaSet));
		
		while(true) {
			LinkedHashMap<String, DNASequence> dbFile = fastaReader.process(1);
			if(dbFile == null || dbFile.size() == 0 )
				break;
			for(Entry<String,DNASequence> seq : dbFile.entrySet())
			{
				//System.out.println(">" + seq.getKey());
				//System.out.println(seq.getValue().getSequenceAsString());
				slib.add(seq.getKey(),seq.getValue().getSequenceAsString().toCharArray(),
						seq.getValue().getReverseComplement().getSequenceAsString().toCharArray());
				
			}
		}
		
		int n = slib.names.size();
		for(int i = 0;i < n; i++)
		{
			slib.add("reverse " + slib.names.get(i), slib.rev_compl_seqs.get(i), slib.seqs.get(i));
		}
		
		return slib;
//		return null;
	}
	
	private void add(String key, char[] s, char[] srev) {
		
		// TODO :: still testing
		this.names.add(key);
		// this.names.add("reverse " + key);
//		char[] s = value.getSequenceAsString().toCharArray();
		this.seqs.add(s);
		//this.seqs.add(srev);
		this.rev_compl_seqs.add(srev);
		//this.rev_compl_seqs.add(s);
		this.weight.add(1.0);
		//this.weight.add(1.0);

	}



	/**
	 *  number of sequences in a seq_lib* */
	public int seq_lib_num_seq()
	{
		return this.seqs.size();
	}
	
	public  String seq_lib_warning_data()
	{
		return "";
	}
	
	private int add_seq_and_rev_comp_to_seq_lib(
		    char[] seq, 
		    String seq_id_plus, 
		    String errfrag) {
		
		return 0;
	}
	
	
	private int add_seq_to_seq_lib(char[] seq, 
		    String seq_id_plus, 
		    String errfrag) {
		
		names.add(seq_id_plus);
		double w = parse_seq_name(seq_id_plus);
		weight.add(w);
		if(w < 0 )
		{
			error.append("Illegal weight");
			return 1;
		}
		seqs.add(seq);
		if(seq.length == 0 )
		{
			error.append("Empty sequence in ");
			return 1;
		}
		
		char offender = upcase_and_check_char(seq);
		if(offender != 0){
			warning.append("Unrecognized character (" + offender + ") in " +errfrag + ", entry " + seq_id_plus);
		}
		
		return 0;
	}



	private char upcase_and_check_char(char[] seq) {
		// TODO :: 
		return 0;
	}



	private double parse_seq_name(String seq_id_plus) {
		
		String[] tokens = seq_id_plus.split(Pattern.quote("*"));
		if(tokens.length > 1 )
		{
			try
			{
				double val = Double.parseDouble(tokens[1]);
				if(val >  PR_MAX_LIBRARY_WT)
					return -1;
				return val;
			}
			catch(Exception ex)
			{
				return -1;
			}
		}
		return 1;
	}



	public String getErrors() {
		return error.toString();
	}



	public String getName(int i) {
		return this.names.get(i);
	}



	public char[] getSeq(int i) {
		return this.seqs.get(i);
	}



	public char[] getSeqRevCompl(int i) {
		return this.rev_compl_seqs.get(i);
	}
}
