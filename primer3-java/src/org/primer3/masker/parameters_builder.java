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

import java.util.ArrayList;
import java.util.List;

/*
 * parameters_builder: Helping structure for reading the formula parameters
 */
class parameters_builder {
	List<formula_parameters> fp_array = null;
	double intercept;
	List<String> used_lists = null;
//	int nslots;
//	int nfp;
	
	
	public int add_variable_to_formula_parameters(String[] list_values) throws FormulaParametersException {
		int i = 0;
		String list_name = null;
		formula_parameters fp = null;
		boolean add_parameters_to_existing_list = false;
		int add_position = -1 ;
		
		if(used_lists == null)
			used_lists = new ArrayList<String>();
		list_name = list_values[0];
		for( i = 0 ; i< used_lists.size();i++)
		{
			if(list_name.equals(used_lists.get(i)) ){
				add_parameters_to_existing_list = true;
				add_position = i;
				break;
			}
				
		}
		
			if(!add_parameters_to_existing_list)
			{
				try {
					fp = formula_parameters.create_formula_parameters_from_list_file_name(list_name);
					
					
					if(fp_array == null )
						fp_array = new ArrayList<formula_parameters>();
					
					used_lists.add(list_name);
					fp_array.add(fp);
					add_position = fp_array.size()-1;
				} catch (FormulaParametersException e) {
					return 1;
//					e.printStackTrace();
				}
			}
			else
			{
				fp = fp_array.get(add_position);
			}
			if (fp != null || add_parameters_to_existing_list)
			{
				boolean squared = false;
				int mm = 0;
				double coef = masker.DEFAULT_COEF;
				if(list_values.length > 1 )
				{
					try
					{
						coef = Double.parseDouble(list_values[1]);
					}
					catch(Exception coefEx)
					{
						throw new FormulaParametersException("Invalid coefficient value.");
					}
				}
				if(list_values.length > 2 )
				{
					try
					{
						mm = Integer.parseInt(list_values[2]);
					}
					catch(Exception mmEx)
					{
						throw new FormulaParametersException("Invalid mismatches value specified: " + list_values[2] +
								". Must be a positive integer less than 2.");
					}
				}
				if(list_values.length > 3 )
				{
					if(list_values[3].equals("sq"))
						squared = true;
				}
				switch(mm){
					case 0 :
						if(squared)
							fp.mm0_2 = coef;
						else
							fp.mm0   = coef;
						break;
					case 1 :
						if(squared)
							fp.mm1_2 = coef;
						else
							fp.mm1   = coef;
						break;
					case 2 :
						if(squared)
							fp.mm2_2 = coef;
						else
							fp.mm2   = coef;
						break;
					
				}
			}
		return 0;
	}
	
}