package com.primer3.masker;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.List;

public class formula_parameters {
	
  	static int glistmaker_code_match = (int) ('G' << 24 | 'T' << 16 | '4' << 8 | 'C');

	
	
	/* If the list is created with GenomeTester4, 
	 * 210 char should be enough to contain the full list name */
	String list_file_name="";
	int oligo_length;
	
	/* binary mask is used for cutting the k-mer into the size of oligo_length */
	long binary_mask;
	
	/* number of unique k-mers in the given k-mer list */
	long words_in_list;
	
	/* pointer to k-mer list */
	ByteBuffer word_list;
	int header_size;
	byte[] pointer;
	int size;
	// size_t size;
	
	/* coefficients for all possible masking formula variables (and their squares) 
	 * concerning this k-mer list. 
	 * If certain variables are not used, their coefficiest are equal to 0 */
	double mm0;
	double mm1;
	double mm2;
	double mm0_2;
	double mm1_2;
	double mm2_2;
	
	
	// TODO :: replace with binary_mask =  ((long)Math.pow(2, this.oligo_length*2)-1)
	void create_binary_mask()
	{
		int i;
		long mask = 0L;

		for (i = 0; i < 2 * oligo_length; i++) {
			mask = (mask << 1) | 1;
		}
		binary_mask = mask;
	}
	
	
	
	
	void get_oligo_frequencies(oligo_counts oc,
			long word, int mm, int strand) {
		
		
		
		formula_parameters fp= this;
		int i, j;
		int count_0mm = 0;
		int count_1mm = 0;
		int count_2mm = 0;
		int mismatch;
		
		word &= fp.binary_mask;
		count_0mm = fp.get_frequency_of_canonical_oligo ( word);
		
		if (mm > 0) {
			for (i = 0; i < fp.oligo_length; i++) {
				for (mismatch = 1; mismatch < 4; mismatch++) {
					long mask = mismatch << (2 * i);
					count_1mm += fp.get_frequency_of_canonical_oligo ( word ^ mask);
					if (mm > 1) {
						for (j = i + 1; j < fp.oligo_length; j++) {
							long mask2 = mismatch << (2 * j);
							count_2mm += fp.get_frequency_of_canonical_oligo ( word ^ mask ^ mask2);
						}
					}
				}
			}
			
		}
		count_1mm += count_0mm;
		count_2mm += count_1mm;
		
		if (strand != masker.REV) {
			oc.count_mm0_fwd = count_0mm;
			oc.count_mm1_fwd = count_1mm;
			oc.count_mm2_fwd = count_2mm;
		}
		if (strand != masker.FWD) {
			oc.count_mm0_rev = count_0mm;
			oc.count_mm1_rev = count_1mm;
			oc.count_mm2_rev = count_2mm;		
		}		
	}
	
	
	private int get_frequency_of_canonical_oligo( long word) {
		
		formula_parameters fp =  this;
		int freq_fwd = 0, freq_rev = 0;
		freq_fwd = fp.binary_search (word);
		if (freq_fwd == 0 ) {
			freq_rev = fp.binary_search (formula_parameters.get_reverse_complement (word, fp.oligo_length));
			
			//fprintf (stderr, "rev %u\n", freq_rev);
			if(freq_rev==0){ /*heuristics is used, by default, for speed and memory issues, 
						kmer lists are generated with kmers freq >1*/
			   freq_rev=1;
			}
			return freq_rev;
		}
		
		//fprintf (stderr, "fwd %u\n", freq_fwd);
		if(freq_fwd==0){ /*heuristics is used, by default, for speed and memory issues,
				         kmer lists are generated with kmers freq >1*/
		       freq_fwd=1;
	        }
		return freq_fwd;
	}
	
	
	
	static public formula_parameters create_formula_parameters_from_list_file_name(String list_file_name) throws FormulaParametersException
	{
		formula_parameters fp = new formula_parameters();
		final int headerToRead = 40;
		
		byte[] data = readFile(list_file_name);
		fp.list_file_name = list_file_name;
//		long header_size;
		int magic;
		ByteBuffer wrapped = ByteBuffer.wrap(data,0,headerToRead);
		wrapped.order(ByteOrder.LITTLE_ENDIAN);
		magic  = wrapped.getInt();
		if(magic != glistmaker_code_match )
		{
			throw new FormulaParametersException("Given file is not a list file.");
		}
		fp.oligo_length = wrapped.getInt(12);
		fp.words_in_list = wrapped.getLong(16);
		fp.header_size = (int)wrapped.getLong(32);
		fp.size = data.length;
		fp.pointer = data;
		fp.word_list = ByteBuffer.wrap(data,fp.header_size,data.length-fp.header_size);
		fp.word_list.order(ByteOrder.LITTLE_ENDIAN);
		fp.create_binary_mask();
		
		System.out.println("words_in_list : " + fp.words_in_list);
		System.out.println("oligo_length : " + fp.oligo_length);
		System.out.println("binary_mask : " + fp.binary_mask + " == " + ((long)Math.pow(2, fp.oligo_length*2)-1));

		System.out.println("header_size : " + fp.header_size);
		
		
		fp.printtest();
		
		
		
		
		return fp;
	}
	
	
	private  void printtest() {
		System.out.println("Current " + word_list.position());
		System.out.println("capacity " + word_list.capacity());

		word_list.position(0);
		for(int i = 0 ; i < 50 ; i ++)
		{
			long currentWord  =  word_list.getLong(header_size + i*12);
			int freq  = word_list.getInt(header_size + (i*12) + 8);
			System.out.println("Current Word " + decode(currentWord)  + " : " + freq);
			
		}
		
	}




	private String decode(long currentWord) {
		// TODO Auto-generated method stub
		return null;
	}




	static public formula_parameters create_formula_parameters_from_list_file_prefix(String list_name_prefix,String kmer_lists_path, int word_length) throws FormulaParametersException{
		String list_file_name =  kmer_lists_path  + list_name_prefix + "_" + word_length + " .list"; 
		return create_formula_parameters_from_list_file_name(list_file_name);
		
	}
	
	static public List<formula_parameters> create_default_formula_parameters(String list_name_prefix, String kmer_lists_path){
		List<formula_parameters> fp = new ArrayList<formula_parameters>();
		try {
			formula_parameters newFP = create_formula_parameters_from_list_file_prefix(list_name_prefix, kmer_lists_path, masker.DEFAULT_WORD_LEN_1);
			fp.add(newFP);
			
			newFP = create_formula_parameters_from_list_file_prefix(list_name_prefix, kmer_lists_path, masker.DEFAULT_WORD_LEN_2);
			fp.add(newFP);
		} catch (FormulaParametersException e) {
//			e.printStackTrace();
			return null;
		}
		return fp;
	}
	
	
	
	static byte[] readFile(String list_file_name) throws FormulaParametersException
	{
//		log("Reading in binary file named : " + list_file_name);
	    File file = new File(list_file_name);
//	    log("File size: " + file.length());
	    byte[] result = new byte[(int)file.length()];
	    try {
	      InputStream input = null;
	      try {
	        int totalBytesRead = 0;
	        input = new BufferedInputStream(new FileInputStream(file));
	        while(totalBytesRead < result.length){
	          int bytesRemaining = result.length - totalBytesRead;
	          //input.read() returns -1, 0, or more :
	          int bytesRead = input.read(result, totalBytesRead, bytesRemaining); 
	          if (bytesRead > 0){
	            totalBytesRead = totalBytesRead + bytesRead;
	          }
	        }
	        /*
	         the above style is a bit tricky: it places bytes into the 'result' array; 
	         'result' is an output parameter;
	         the while loop usually has a single iteration only.
	        */
//	        log("Num bytes read: " + totalBytesRead);
	      }
	      finally {
//	        throw new FormulaParametersException("Closing input stream.");
	        input.close();
	      }
	    }
	    catch (FileNotFoundException ex) {
	    	throw new FormulaParametersException("File not found.");
	    }
	    catch (IOException ex) {
	    	throw new FormulaParametersException(ex.getMessage());
	    }
	    
	    return result;
	}
	
	





	public static void main(String[] args)
	{
		char[] c = new char[2];
		System.out.println("glistmaker_code_match : " + formula_parameters.glistmaker_code_match + " , " +c[0] );
		
		
		
		StringBuffer strBuffer = new StringBuffer(masker.MAX_BUFFER_SIZE);

		
		System.out.println(" asd :" +
		strBuffer.toString());
		
		
//		formula_parameters fp = new formula_parameters();
//		
//		try {
//			formula_parameters.create_formula_parameters_from_list_file_name("/data/softwares/primer3-primer3/kmer_lists/homo_sapiens_11.list");
//		} catch (FormulaParametersException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
	}

	void test(String list_file_name) throws FormulaParametersException
	{
		byte[] data = readFile(list_file_name);
		
		
		byte[] arr1 = { 'C', '4','T','G' };
		byte[] arr2 = { 'G','T', '4','C'};

		ByteBuffer wrapped = ByteBuffer.wrap(arr1); // big-endian by default
		int num = wrapped.getInt(); // 1
		System.out.println("magic1 : " + num);
		wrapped = ByteBuffer.wrap(arr2);
		num = wrapped.getInt(); // 1
		System.out.println("magic2 : " + num);
//		ByteBuffer dbuf = ByteBuffer.allocate(4);
//		dbuf.putShort(num);
//		byte[] bytes = dbuf.array(); // { 0, 1 }
//		
		
		
		
		wrapped = ByteBuffer.wrap(data, 0, 8);
		wrapped.order(ByteOrder.LITTLE_ENDIAN);
		int res = wrapped.getInt();
		
		System.out.println("magic : " + res);

		
		for(int i =0;i<20;i++)
			System.out.println((char)data[i]);
	}






	
	
	// TODO :: add this
	public static List<formula_parameters> read_formula_parameters_from_file(
			String lists_file_name, int nlist_parameters,
			parameters_builder pbuilder, double formula_intercept) {
		// TODO Auto-generated method stub
		return null;
	}






	public int binary_search(long word) {
		
		final int row_size = 8 + 4 ; // sizeof (long) + sizeof (int) 
		long current_word;
		int low, high, mid;
		int freq;
		low = 0;
		high = (int) (this.words_in_list - 1);
		mid = (low + high) / 2;
		while (low <= high) {
			
			current_word = this.word_list.getLong( header_size + (mid * row_size));
			if (current_word < word) {
				low = mid + 1;
			} else if (current_word > word) {
				if (mid == 0) break;
				high = mid - 1;
			} else {
				freq =  (this.word_list.getInt(header_size + (mid * row_size) + 8 ));
				return freq;
			}
			mid = (low + high) / 2;
		}
		return 0;
	}






	public static long get_reverse_complement(long word, int word_length) {
		int i;
		long mask, v, revcompl = 0L;

		word = ~word;
		mask = 3;
		for (i = 0; i < word_length; i++) {
			v = word & mask;
			revcompl <<= 2;
			revcompl |= v;
			word >>= 2;
		}
		return revcompl;
	}
}