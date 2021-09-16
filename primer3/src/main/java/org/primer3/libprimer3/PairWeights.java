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
package org.primer3.libprimer3;

public class PairWeights {
	
	public double primer_quality;
	public double io_quality;
	public double diff_tm;
	public double compl_any;
	public double compl_any_th;
	public double compl_end;
	public double compl_end_th;
	public double temp_cutoff;
	public double product_tm_lt;
	public double product_tm_gt;
	public double product_size_lt;
	public double product_size_gt;
	public double repeat_sim;
	
	
//	TODO :: PR_ASSERT(pa.pr_pair_weights.template_mispriming >= 0.0);
	public double template_mispriming;

//	TODO :: PR_ASSERT(pa.pr_pair_weights.template_mispriming_th >= 0.0);
	public double template_mispriming_th;
	
	public void set_template_mispriming_th(String datum) {
		this.template_mispriming_th = Double.parseDouble(datum);
	}

	public void set_template_mispriming(String datum) {
		this.template_mispriming = Double.parseDouble(datum);
		
	}

	public void set_repeat_sim(String datum) {
		this.repeat_sim = Double.parseDouble(datum);
		
	}

	public void set_product_size_lt(String datum) {
		this.product_size_lt = Double.parseDouble(datum);
		
	}

	public void set_product_size_gt(String datum) {
		this.product_size_gt = Double.parseDouble(datum);
		
	}

	public void set_product_tm_gt(String datum) {
		this.product_tm_gt = Double.parseDouble(datum);
		
	}

	public void set_product_tm_lt(String datum) {
		this.product_tm_lt = Double.parseDouble(datum);
		
	}

	public void set_compl_end_th(String datum) {
		this.compl_any_th = Double.parseDouble(datum);
		
	}

	public void set_compl_any_th(String datum) {
		this.compl_any_th = Double.parseDouble(datum);
		
	}

	public void set_compl_end(String datum) {
		this.compl_end = Double.parseDouble(datum);
		
	}

	public void set_compl_any(String datum) {
		this.compl_any = Double.parseDouble(datum);
		
	}

	public void set_diff_tm(String datum) {
		this.diff_tm = Double.parseDouble(datum);
		
	}

	public void set_io_quality(String datum) {
		this.io_quality = Double.parseDouble(datum);
		
	}

	public void set_primer_quality(String datum) {
		this.primer_quality = Double.parseDouble(datum);
		
	}



	public void p3_print_args() {
		System.out.format("  begin pr_pair_weights\n") ;
		System.out.format("    primer_quality %f\n", this.primer_quality) ;
		System.out.format("    io_quality %f\n", this.io_quality) ;
		System.out.format("    diff_tm %f\n", this.diff_tm) ;
		System.out.format("    compl_any %f\n", this.compl_any) ;
		System.out.format("    compl_end %f\n", this.compl_end) ;
		System.out.format("    compl_any_th %f\n", this.compl_any_th) ;
		System.out.format("    compl_end_th %f\n", this.compl_end_th) ;
		System.out.format("    product_tm_lt %f\n", this.product_tm_lt) ;
		System.out.format("    product_tm_gt %f\n", this.product_tm_gt) ;
		System.out.format("    product_size_lt %f\n", this.product_size_lt) ;
		System.out.format("    product_size_gt %f\n", this.product_size_gt) ;
		System.out.format("    repeat_sim %f\n", this.repeat_sim) ;
		System.out.format("    template_mispriming %f\n", this.template_mispriming) ;
		System.out.format("    template_mispriming_th %f\n", this.template_mispriming_th) ;
		System.out.format("  end pair_weights\n") ;		
	}
 }