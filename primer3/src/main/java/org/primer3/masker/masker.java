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
package org.primer3.masker;


public class masker {
	/* ALPHABET defines which characters are considered as nucleotides */
	static public String ALPHABET  = "ACGTUacgtu";

	static public int FWD  = 1;
	static public int REV  = 2;
	static public int BOTH = 3;

	static public int MAX_BUFFER_SIZE = 5000;
	static public int MAX_SPLITS = 10;

	static public int NUCLEOTIDE = 0;
	static public int WHITESPACE = 1;
	static public int MASKED_CHAR = 2;

	static public boolean FLUSH_ALL = true;
	static public boolean KEEP_UNCERTAIN_REGION = false;

	static public double DEFAULT_COEF = 1.0;

	/* default formula parameter coefficients */
	static public int DEFAULT_NLISTS = 2;
	static public int DEFAULT_NLIST_PARAMETERS = 2;
	static public int DEFAULT_WORD_LEN_1 = 11;
	static public int DEFAULT_WORD_LEN_2 = 16;
	static public double DEFAULT_COEF_1 = 0.1772;
	static public double DEFAULT_COEF_2 = 0.239;
	static public double DEFAULT_INTERCEPT = -4.336;

	/* default masking parameter values */
	static public boolean PRINT_SEQUENCE = true; /* used only for commandline tool */
	static public double DEFAULT_FAILURE_RATE = 0.1;
	static public int DEFAULT_ABS_CUTOFF = 0 ;
	static public char DEFAULT_MASK_CHAR = 'N';
	static public masking_direction DEFAULT_MASKING_DIRECTION = masking_direction.both_on_same; /* used only for commandline tool */
	static public int DEFAULT_M5P = 1;
	static public int DEFAULT_M3P = 0;
	static public String DEFAULT_LIST_FILE_PREFIX =  "homo_sapiens";
	
	
	public static void delete_formula_parameters(formula_parameters[] fp,
			int nlists) {
		// TODO Auto-generated method stub
		
	}

	public static formula_parameters[] create_default_formula_parameters(
			char[] list_prefix, String kmer_lists_path,
			StringBuilder fatal_parse_err) {
		// TODO Auto-generated method stub
		return null;
	}

	public static long string_to_word(char[] oligo_seq, int length,
			int window_size) {
		// TODO Auto-generated method stub
		return 0;
	}

	

	
	
	
	
	
	
	public static void read_and_mask_sequence(input_sequence input_seq,
			output_sequence output_seq, MaskerParameters mp, boolean debug) {
		
		int is_header = 0;
		boolean init_round = true;
		long word_fwd = 0, word_rev = 0, nucl_value;
		int current_length = 0;
		int word_length = 0;
		long binary_mask = 0;
		masking_buffer mbuffer;
		long header_pos = 0, current_pos = 0;
		int i;
		
		for(formula_parameters fp : mp.fp)
		{
			if(fp.oligo_length > word_length )
			{
				word_length = fp.oligo_length;
				binary_mask = fp.binary_mask;
			}
		}
		mbuffer = new masking_buffer (word_length + mp.nucl_masked_in_3p_direction , mp);
		
		for(char c : input_seq.sequence_string) 
		{
			if (!init_round && mbuffer.wi == mbuffer.ri) {
				mbuffer.empty_buffer (output_seq, masker.KEEP_UNCERTAIN_REGION);
			}
			init_round = false;
			
			if (!ALPHABET.contains(c+"") && c > ' ') { 
				mbuffer.add_char_to_buffer (c,MASKED_CHAR);
				word_fwd = word_rev = 0;
				current_length = 0;
				continue;
			} else if (c <= ' ') {
				mbuffer.add_char_to_buffer (c, WHITESPACE);
				continue;
			} else {
				mbuffer.add_char_to_buffer (c, NUCLEOTIDE);
			}
			
			/* add next nucleotide to the word */
			nucl_value = get_nucl_value (c);
			if (mp.mdir != masking_direction.rev) {
				word_fwd <<= 2;
				word_fwd |= nucl_value;
			}
			if (mp.mdir != masking_direction.fwd) {
				word_rev >>= 2;
				word_rev |= ((~nucl_value & 3) << ((word_length - 1) * 2));
			}
			current_length += 1;
			
			if (current_length > word_length) {
				word_fwd &= binary_mask;
				word_rev &= binary_mask;
				current_length = word_length;
			}
			if (current_length == word_length) {
				oligo_pair h = new oligo_pair();
				h.fwd = word_fwd;
				h.rev = word_rev;
				
//				if (debug ) {
//					fprintf (stderr, "%llu %llu\n", h.fwd, h.rev);
//				}
				mask_oligo_region (h, mp, mbuffer, word_length, debug);
			}
			
		}
		
		mbuffer.empty_buffer (output_seq, FLUSH_ALL);

	}

	private static void mask_oligo_region(oligo_pair h, MaskerParameters mp,
			masking_buffer mbuffer, int word_length, boolean debug) {
		h.calculate_scores(mp,word_length);
		
		if (mp.mdir != masking_direction.rev && ((mp.failure_rate > 0 && h.score_fwd > mp.failure_rate) || 
				(mp.abs_cutoff > 0 && h.abs_score >= mp.abs_cutoff))) {
				int masked = 0, 
						i = mbuffer.TAKE_STEP_BACK(mbuffer.wi);
				while (masked < mp.nucl_masked_in_5p_direction) {
					if (!mbuffer.non_nucleotide_positions[i] && !mbuffer.mask_positions_fwd[i]) {
						mbuffer.mask_positions_fwd[i] = true;
						masked += 1;
					} else if (mbuffer.mask_positions_fwd[i]) {
						masked += 1;
					}
					i = mbuffer.TAKE_STEP_BACK(i);
				}
				mbuffer.mi = mp.nucl_masked_in_3p_direction;
			}
			
			if (mp.mdir != masking_direction.fwd && ((mp.failure_rate > 0 && h.score_rev > mp.failure_rate) || 
				(mp.abs_cutoff > 0 && h.abs_score >= mp.abs_cutoff))) {
				int masked = 0, i = mbuffer.TAKE_STEP_BACK(mbuffer.ei);		
				while (masked < mp.nucl_masked_in_5p_direction + mp.nucl_masked_in_3p_direction) {
					if (!mbuffer.non_nucleotide_positions[i] && !mbuffer.mask_positions_rev[i]) {
						mbuffer.mask_positions_rev[i] = true;
						masked += 1;
					} else if (mbuffer.mask_positions_rev[i]) {
						masked += 1;
					}
					i = mbuffer.TAKE_STEP_FORWARD(i);
				}
				
			}
		
	}

	static final long bit1 = 1 << 2;
	static final  long bit2 = 3 << 1;
	private static long get_nucl_value(char nucl) {

		if ((nucl & bit1) != 0) {
			return ((nucl >> 4) | 2) & 3;
		}
		return (nucl & bit2) >> 1;
	}

	
	
	
	
	
	
	
	
	
}
