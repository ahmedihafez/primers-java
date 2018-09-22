package org.primer3.p3_seq_lib;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import org.primer3.libprimer3.LibPrimer3;

// TODO :: missing impl...
public class seq_lib {

	
	static double  PR_MAX_LIBRARY_WT = 100.0;

	
	/** An array of sequence names. */
	public List<String> names;
	/** An array of sequences. */
	public List<char[]> seqs;
	/** An array of reversed-complemented sequences.
	   x->rev_compl_seqs[i] is the reverse complement
	   of x->seqs[i], which lets us keep track of pairwise
	   mispriming.  See reverse_complement_seq_lib(). */
	public List<char[]> rev_compl_seqs;
	
	public List<Double> weight;
	
	/* The number of names, sequences, and weights. */
	public int seq_num = 0;
	
	StringBuilder warning = new StringBuilder();
	StringBuilder error = new StringBuilder();

	
	
	public seq_lib()
	{
		names = new ArrayList<String>();
		seqs = new ArrayList<char[]>();
		rev_compl_seqs = new ArrayList<char[]>();
		weight = new ArrayList<Double>();
		seq_num = 0;
	}
	
	
	
	/*  
	 * Reads any file in fasta format and returns a newly allocated
	 * seq_lib, lib.  Sets lib.error to a non-empty string on any error
	 * other than ENOMEM.  Returns NULL on ENOMEM.
	 */
	public static seq_lib read_and_create_seq_lib(String filename, String errfrag) throws Exception
	{
		return null;
	}
	
	/**
	 *  number of sequences in a seq_lib* */
	static public int seq_lib_num_seq(seq_lib lib)
	{
		if (lib == null ) return 0;
		return lib.seq_num;
	}
	
	public  String seq_lib_warning_data()
	{
		return "";
	}
	
	public int add_seq_and_rev_comp_to_seq_lib(
		    char[] seq, 
		    String seq_id_plus, 
		    String errfrag) {
		
		return 0;
	}
	
	
	int add_seq_to_seq_lib(char[] seq, 
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
}
