/*
    This file is part of primer3 porting to java


	Original file are part of https://github.com/primer3-org/primer3
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
 * output_sequence: A wrapper for output sequence in case the output is given in a 
 * string variable
 */
public class output_sequence {
	
	
	// all can be converted into a StringBuffer or builder
	public char [] sequence;
	int pos;
	
	/* These will be used instead of char *sequence
	 * in case of both_separately masking direction */
	public char []sequence_fwd;
	public char []sequence_rev;
	
	
	

	public void write_char_to_output(char c, char c_other, MaskerParameters mp) {
		
		if (mp.print_sequence) {
			System.err.print(c);
		} 
		else {
			if (mp.mdir == masking_direction.both_separately) {
				this.sequence_fwd[this.pos] = c;
				this.sequence_rev[this.pos] = c_other;
			} else {
				this.sequence[this.pos] = c;
			}
			this.pos += 1;
		}
	}
	void  write_header_to_output (String header_name,MaskerParameters mp)
	{
		if (mp.print_sequence) 
			System.err.println(header_name);
		else {
			if (mp.mdir == masking_direction.both_separately) {
//				v = memcpy (output_seq.sequence_fwd + output_seq.pos, header_name, strlen(header_name));
//				if (v) v = memcpy (output_seq.sequence_rev + output_seq.pos, header_name, strlen(header_name));
			} else {
//				v = memcpy (output_seq.sequence + output_seq.pos, header_name, strlen(header_name));
			}	
//			this.pos += header_name.length();
		}
	}
	
	public static output_sequence create_output_sequence (int seq_len, masking_direction mdir)
	{
		output_sequence output_seq = new output_sequence();
		if (mdir == masking_direction.both_separately) {
			output_seq.sequence_fwd = new char [seq_len];
			output_seq.sequence_rev = new char [seq_len];
		} else {
			output_seq.sequence = new char [seq_len];
		}
		output_seq.pos = 0;
		return output_seq;
	}
	
}