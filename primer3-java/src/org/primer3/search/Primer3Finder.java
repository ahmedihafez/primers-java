/*
    This file is part of primer3 porting to java


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
package org.primer3.search;

import java.util.List;

import org.primer3.libprimer3.DPAlArgHolder;
import org.primer3.libprimer3.P3RetVal;
import org.primer3.libprimer3.THAlArgHolder;
import org.primer3.primer.PrimerPair;

public abstract class Primer3Finder {
	P3RetVal retval;
	DPAlArgHolder dpal_arg_to_use;
	THAlArgHolder thal_arg_to_use;
	THAlArgHolder thal_oligo_arg_to_use;
	
	
	public Primer3Finder (	P3RetVal retval,
	DPAlArgHolder dpal_arg_to_use,
	THAlArgHolder thal_arg_to_use,
	THAlArgHolder thal_oligo_arg_to_use)
	{
		this.retval  = retval;
		this. dpal_arg_to_use = dpal_arg_to_use;
		this. thal_arg_to_use = thal_arg_to_use;
		this. thal_oligo_arg_to_use = thal_oligo_arg_to_use;
	}
	
	
	abstract protected void getLocalNextResult() throws Exception;
	
	
	public List<PrimerPair> getNextResult() throws Exception
	{
		// cache current set
		retval.best_pairs.cacheCurrent();
		this.getLocalNextResult();
		List<PrimerPair> newPairs = retval.best_pairs.pairs;
		retval.best_pairs.mergeBests();
		return newPairs;
	}
	
}
