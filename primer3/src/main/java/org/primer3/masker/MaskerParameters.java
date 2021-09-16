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

import java.util.List;


public class MaskerParameters {
	/* strand to mask */
	public masking_direction mdir;
	
	/* primer failure rate cutoff used in primer design, 
	 * potential locations in a sequence for primers with PCR 
	 * failure rate over the given cutoff are masked
	 *(see function calculate_scores() from masker.c) */
	public double failure_rate;
	
	/* absolute value cutoff, this can be used for masking all the k-mers in a sequence
	 * that have the frequency over abs_cutoff in a k-mer list */
	int abs_cutoff;
	
	/* number of nucleotides masked in 5' and 3' direction with respect
	 * to the 3' end of a primer */
	public int nucl_masked_in_5p_direction;
	public int nucl_masked_in_3p_direction;
	
	/* If masker is used as a separate application then always print_sequence=1, 
	 * i.e the output is sent to stdout.
	 * If print_sequence=0 the output is written in a string variable and can be forwarded
	 * to the next function */
	public boolean print_sequence;
	
	/* if do_soft_masking=1, masked nucleotides and converted to lower-case, else 
	 * masked nucleotide are converted to masking_char ('N' by default) */
	public boolean do_soft_masking;
	char masking_char;
	
	/* size of the masking window */
	public int window_size;
	
	/* number of k-mer lists used in the masking formula */
	public int nlists;
	/* k-mer lists and all their parameters which are used in the masking formula */
	public String list_prefix;
	public List<formula_parameters> fp;
	public double formula_intercept;
	
	public void set_nucl_masked_in_3p_direction(String datum) {
		this.nucl_masked_in_3p_direction = Integer.parseInt(datum);
	}
	
	public void set_list_prefix(String datum) {
		this.list_prefix = datum;
	}

	public void set_nucl_masked_in_5p_direction(String datum) {
		this.nucl_masked_in_5p_direction = Integer.parseInt(datum);
	}
	
	
	
	public void set_failure_rate(String datum) {
		this.failure_rate = Double.parseDouble(datum);		
	} 
	
	
	



}
