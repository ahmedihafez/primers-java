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

/* 
 * masking_buffer: A round buffer for storing all the characters of the input sequence
 * as well as the masking and output information 
 */
class masking_buffer {
	/* the main array containing the characters that are read
	 * from the sequence input */
	char[] buffer = new char[masker.MAX_BUFFER_SIZE];
//	StringBuffer strBuffer = new StringBuffer(masker.MAX_BUFFER_SIZE);
	/* indicates the positions in the buffer that do not contain a
	 * nucleotide character */
	boolean[] non_nucleotide_positions = new boolean[masker.MAX_BUFFER_SIZE];
	
	/* indicates the nucleotide positions in the buffer that
	 * will be masked in he output */
	boolean[] mask_positions_fwd =  new  boolean[masker.MAX_BUFFER_SIZE];
	boolean[] mask_positions_rev =  new boolean[masker.MAX_BUFFER_SIZE];
	
	/* indices for reading from the buffer and 
	 writing to the buffer */
	int ri = 0;
	int wi = 0;
	
	/* index of the 5' end of the k-mer last entered plus the 
	 * number of nucleotides masked in 3' direction. This index is used
	 * to determine the first positions which can still be altered and therefore cannot be 
	 * flushed out of the buffer yet */
	int ei;
	
	/* mi indicates the number of nucleotides that are
	 * yet to be masked in 3' direction */
	int mi;
	
	
	public masking_buffer(int word_length , MaskerParameters mp)
	{
		this.ei = masker.MAX_BUFFER_SIZE - word_length +1;
		this.mp = mp;
	}
	
	
	
	MaskerParameters mp ;
	public void setMaskingParameters(MaskerParameters mp )
	{
		this.mp = mp ;
	}
	
	
	void add_char_to_buffer (char c, int char_type)
	{
		this.buffer[this.wi] = c;
		this.mask_positions_fwd[this.wi] = false;
		this.mask_positions_rev[this.wi] = false;
		this.non_nucleotide_positions[this.wi] = false;
		
		if (char_type != masker.WHITESPACE) {
			if (this.mi > 0) {
				this.mask_positions_fwd[this.wi] = true;
				this.mi -= 1;
			} else if (char_type == masker.MASKED_CHAR) {
				this.mask_positions_rev[this.wi] = true;
				this.mask_positions_fwd[this.wi] = true;
			}
			while (this.non_nucleotide_positions[this.ei] && !(this.mask_positions_fwd[this.ei])) {
				this.ei = TAKE_STEP_FORWARD(this.ei);
			}
			/* when it finds the first non-whitespace character it increases the ei one more time */
			this.ei = TAKE_STEP_FORWARD(this.ei);
		} 
		if (char_type == masker.WHITESPACE || char_type == masker.MASKED_CHAR) {
			this.non_nucleotide_positions[this.wi] = true;
		}
		
		this.wi = TAKE_STEP_FORWARD(this.wi);
	}


	public void empty_buffer (output_sequence output_seq, boolean flush_all)
	{
		int end = this.ei;
		if (flush_all) end = this.wi;
		while (this.ri != end) {
			if (COND0(this.ri)) {
				output_seq.write_char_to_output ( this.buffer[this.ri], this.buffer[this.ri], mp);
			} else {
				if (mp.do_soft_masking) {
					output_seq.write_char_to_output ( (char)(COND1(this.ri) ? this.buffer[this.ri] + 32 : this.buffer[this.ri]), 
							(char)   (COND2(this.ri) ? this.buffer[this.ri] + 32 : this.buffer[this.ri]), mp);
				} else {
					output_seq.write_char_to_output ( COND1(this.ri) ? mp.masking_char : this.buffer[this.ri], 
							      COND2(this.ri) ? mp.masking_char : this.buffer[this.ri], mp);				
				}
			}
			this.ri = TAKE_STEP_FORWARD(this.ri);
		}
		return;
	}
	
	/* next character is not nucleotide, or if nucleotide then already masked */
	boolean COND0(int i) 
	{
		return (this.non_nucleotide_positions[i] || (mp.do_soft_masking && this.buffer[i] >= 'a')) ;
	}
	/* output contains only one sequence and the next character is masked or
	 * output contains two sequences and the next character is masked on the forward strand */
	boolean COND1(int i) 
	{ 
	return 	 (
				(	
						mp.mdir != masking_direction.both_separately && 
						( this.mask_positions_fwd[i] ||  this.mask_positions_rev[i] ) 
				) || 
				(
						mp.mdir == masking_direction.both_separately && this.mask_positions_fwd[i])
				);
	}
		
	/* output contains two sequences and the next character is masked on the reverse strand */
	boolean COND2(int i){
		return (mp.mdir == masking_direction.both_separately && this.mask_positions_rev[i]);
	}
	
	
	public int TAKE_STEP_FORWARD(int i) {
		
		int BUFFER_SIZE = buffer.length;
		
		return (i == BUFFER_SIZE - 1) ? 0 : (i + 1);
	}
	
	public int TAKE_STEP_BACK(int i) {
		
		int BUFFER_SIZE = buffer.length;
		
		return (i == 0) ? (BUFFER_SIZE - 1) : (i - 1);
	}
}