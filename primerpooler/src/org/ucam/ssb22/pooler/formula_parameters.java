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

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;

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




	
	
	
	
	
	public static byte[] readFile(String list_file_name) throws IOException
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
	    	throw new IOException("File not found.");
	    }
	    catch (IOException ex) {
	    	throw ex;
	    }
	    
	    return result;
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