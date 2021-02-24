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

/**
 * oligo_pair: contains k-mer pairs and their calculated scores
 */
public class oligo_pair {
	/* fwd and rev k-mers */
	public long fwd = 0;
	public long rev = 0;
	
	/* a k-mer score is a failure rate of such a primer that 
	 *contains a given k-mer in its 3' end */
	public double score_fwd = 0;
	double score_rev = 0;
	
	/* abs_score is a number of given k-mers in a k-mer list */
	int abs_score = 0;

	public void calculate_scores(MaskerParameters mp, int word_length) {
		// TODO Auto-generated method stub
		for(formula_parameters fp : mp.fp)
		{
			oligo_counts oc = new oligo_counts();
			int mm = (fp.mm2 !=0 || fp.mm2_2 !=0) ? 2 : ((fp.mm1 !=0 || fp.mm1_2 !=0 ) ? 1 : 0);
			if ((mp.mdir == masking_direction.both_on_same || mp.mdir == masking_direction.both_separately || mp.abs_cutoff > 0 ) && word_length == fp.oligo_length) {
				double score = 0.0;
				int abs_score = 0;
				
				if (mp.mdir == masking_direction.rev) 
					fp.get_oligo_frequencies (oc,  this.rev, mm, masker.BOTH);
				else 
					fp.get_oligo_frequencies (oc,this.fwd, mm,  masker.BOTH);
				
				if (oc.count_mm0_fwd > 0 || oc.count_mm0_rev > 0) {
					int count = oc.count_mm0_fwd > 0 ? oc.count_mm0_fwd : oc.count_mm0_rev;
					score += fp.mm0 * Math.log(oc.count_mm0_fwd) + fp.mm0_2 *  Math.log(oc.count_mm0_fwd) *  Math.log(oc.count_mm0_fwd);
					abs_score = count;
				}
				if (oc.count_mm1_fwd > 0) {
					score += fp.mm1 *  Math.log(oc.count_mm1_fwd) + fp.mm1_2 *  Math.log(oc.count_mm1_fwd) *  Math.log(oc.count_mm1_fwd);
					abs_score = oc.count_mm1_fwd;
				}	
				if (oc.count_mm2_fwd > 0) {
					score += fp.mm2 *  Math.log(oc.count_mm2_fwd) + fp.mm2_2 *  Math.log(oc.count_mm2_fwd) *  Math.log(oc.count_mm2_fwd);
					abs_score = oc.count_mm2_fwd;
				}
				if (mp.abs_cutoff > 0) {
					this.abs_score = abs_score;
				} else {
					this.score_fwd += score;
					this.score_rev += score;
				}
			} else {
				if (mp.mdir != masking_direction.rev) {
					fp.get_oligo_frequencies (oc, this.fwd, mm, masker.FWD);
					
					//fprintf (stderr, "oc.count.fwd %u\n", oc.count_mm0_fwd);
					
					if (oc.count_mm0_fwd > 0) this.score_fwd += fp.mm0 * Math.log(oc.count_mm0_fwd) + fp.mm0_2 * Math.log(oc.count_mm0_fwd) * Math.log(oc.count_mm0_fwd);
					if (oc.count_mm1_fwd > 0) this.score_fwd += fp.mm1 * Math.log(oc.count_mm1_fwd) + fp.mm1_2 * Math.log(oc.count_mm1_fwd) * Math.log(oc.count_mm1_fwd);
					if (oc.count_mm2_fwd > 0) this.score_fwd += fp.mm2 * Math.log(oc.count_mm2_fwd) + fp.mm2_2 * Math.log(oc.count_mm2_fwd) * Math.log(oc.count_mm2_fwd);
					
					//fprintf (stderr, "1: sõnad: %s %s, skoorid %f %f, abs_skoor %u\n", word_to_string(h.fwd, word_length), word_to_string(h.rev, word_length), h.score_fwd, h.score_rev, h.abs_score);
				}
				if (mp.mdir != masking_direction.fwd) {
					fp.get_oligo_frequencies (oc, this.rev, mm, masker.REV);
					if (oc.count_mm0_rev > 0) this.score_rev += fp.mm0 * Math.log(oc.count_mm0_rev) + fp.mm0_2 * Math.log(oc.count_mm0_rev) * Math.log(oc.count_mm0_rev);
					if (oc.count_mm1_rev > 0) this.score_rev += fp.mm1 * Math.log(oc.count_mm1_rev) + fp.mm1_2 * Math.log(oc.count_mm1_rev) * Math.log(oc.count_mm1_rev);
					if (oc.count_mm2_rev > 0) this.score_rev += fp.mm2 * Math.log(oc.count_mm2_rev) + fp.mm2_2 * Math.log(oc.count_mm2_rev) * Math.log(oc.count_mm2_rev);
				}
			}
		}
		
		if (this.score_fwd > 0) this.score_fwd = Math.exp(this.score_fwd + mp.formula_intercept) / (1 + Math.exp(this.score_fwd + mp.formula_intercept));
		if (this.score_rev > 0) this.score_rev = Math.exp(this.score_rev + mp.formula_intercept) / (1 + Math.exp(this.score_rev + mp.formula_intercept));
		
	}

	

	
}
