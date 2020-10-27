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
package org.primer3;

import java.util.HashMap;

public class Primer3Manual {

	public static final String PRIMER_PAIR_MAX_DIFF_TM = "PRIMER_PAIR_MAX_DIFF_TM";

	public static String PRIMER_PRODUCT_MAX_TM = "PRIMER_PRODUCT_MAX_TM";
	
	static HashMap<String,String> helpTxt;
	static HashMap<String,String> helpLongTxt;

	static {
		helpTxt = new HashMap<String, String>();
		helpLongTxt = new HashMap<String, String>();

		
		helpTxt.put(PRIMER_PRODUCT_MAX_TM, "The maximum allowed melting temperature of the amplicon.");
		helpLongTxt.put(PRIMER_PRODUCT_MAX_TM, "Primer3 calculates product Tm calculated using the formula from Bolton and McCarthy, PNAS 84:1390 (1962) as presented in Sambrook, Fritsch and Maniatis, Molecular Cloning, p 11.46 (1989, CSHL Press).");

		
		
		helpTxt.put(PRIMER_PAIR_MAX_DIFF_TM, "Maximum acceptable (unsigned) difference between the melting temperatures of the forward and reverse primers.");

	}
	
	
	public static String getHelp(String key) {
		return getHelp(key,false);
	}
	public static String getHelp(String key, boolean longVersion) {
		String shortHelp = helpTxt.get(key); 
		if(shortHelp == null)
			shortHelp = "";
		if(longVersion)
		{
			String longHelp = helpLongTxt.get(key);
			if(longHelp != null)
			{
				if(!shortHelp.endsWith("."))
					shortHelp += ".";
				shortHelp +=  " " + longHelp;
			}
		}
		return shortHelp;
	}
	
	
}
