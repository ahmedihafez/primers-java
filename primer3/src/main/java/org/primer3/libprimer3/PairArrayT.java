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

import java.util.ArrayList;
import java.util.List;

import org.primer3.primer.PrimerPair;

public class  PairArrayT {
//	  int         storage_size = 0;
	  public int         num_pairs = 0;
	  // TODO :: make sure it is init. somewhere
	  public List<PrimerPair> pairs = new ArrayList<PrimerPair>();
	  public PairStats  expl = new PairStats();
	  
	  
	  public String p3_get_pair_array_explain_string()
	  {
		  return expl.p3_pair_explain_string();
	  }


	
	  public void add_pair(PrimerPair the_best_pair) {

		  pairs.add(the_best_pair);
		  this.num_pairs++;
		  
	  }
	  
	  public void clearSet()
	  {
		  pairs.clear();
		  num_pairs = 0;
	  }


	List<PrimerPair> currentSet = null;
	public void cacheCurrent() {
		currentSet = this.pairs;
		this.pairs =new ArrayList<PrimerPair>();
		this.num_pairs = 0;
	}



	public void mergeBests() {
		
		currentSet.addAll(pairs);
		num_pairs =  currentSet.size();
		pairs = currentSet;
		currentSet = null;
		
	}
}