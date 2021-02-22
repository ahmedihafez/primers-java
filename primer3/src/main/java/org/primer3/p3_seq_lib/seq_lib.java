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
package org.primer3.p3_seq_lib;

import java.io.File;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.HashMap;
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
	
	
	HashMap<String, Integer > seqNameToIndex = new HashMap<String, Integer>(); 
	
	StringBuilder warning = new StringBuilder();
	StringBuilder error = new StringBuilder();

	
	
	public seq_lib()
	{
		names = new ArrayList<String>();
		seqs = new ArrayList<char[]>();
		rev_compl_seqs = new ArrayList<char[]>();
		weight = new ArrayList<Double>();
	}
	
	public static seq_lib read_and_create_seq_lib(String filename, String libName) throws Exception
	{
		return read_and_create_seq_lib(filename, true);
	}
	
	/*  
	 * Reads any file in fasta format and returns a newly allocated
	 * seq_lib, lib.  Sets lib.error to a non-empty string on any error
	 * other than ENOMEM.  Returns NULL on ENOMEM.
	 */
	// Input args String errfrag
	public static seq_lib read_and_create_seq_lib(String filename, boolean addRev) throws Exception
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
		
		if(addRev)
		{
			int n = slib.names.size();
			for(int i = 0;i < n; i++)
			{
				slib.add("reverse " + slib.names.get(i), slib.rev_compl_seqs.get(i), slib.seqs.get(i));
			}
		}
		return slib;
//		return null;
	}
	
	public void add(String key, char[] s, char[] srev) {
		
		// TODO :: still testing
		seqNameToIndex.put(key, this.names.size());
		this.names.add(key);
		
		// this.names.add("reverse " + key);
//		char[] s = value.getSequenceAsString().toCharArray();
		this.seqs.add(s);
		//this.seqs.add(srev);
		this.rev_compl_seqs.add(srev);
		
		//this.rev_compl_seqs.add(s);
		
		// try to get the wieght
		double weight = 1;
		if(key.contains("*"))
		{
			String[] tokens = key.split("\\*");
			try {
				weight = Double.parseDouble(tokens[1]);
			}
			catch(Exception ex)
			{
				// ignore the wight
				System.err.println("Can not parse weight for "+ key +" Sequence in the sequence lib.");
			}
		}
		
		
		this.weight.add(weight);
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

	public char[] getSeq(String datum) {
		// TODO Auto-generated method stub
		Integer index= seqNameToIndex.get(datum);;
		if(index != null)
			return seqs.get(index);
		return  null;
	}
}
