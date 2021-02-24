/*
    This file is part of primer3 porting to java
    Primer3 and the libprimer3 library are free software;
    you can redistribute them and/or modify them under the terms
    of the GNU General Public License as published by the Free
    Software Foundation; either version 2 of the License, or (at
    your option) any later version.

	Original file are part of https://github.com/primer3-org/primer3
	Whitehead Institute for Biomedical Research, Steve Rozen
	(http://purl.com/STEVEROZEN/), Andreas Untergasser and Helen Skaletsky
	All rights reserved.



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
package org.primer3.multisearch;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.primer3.dpal.AlignmentException;
import org.primer3.libprimer3.P3GlobalSettings;
import org.primer3.libprimer3.P3RetVal;
import org.primer3.multisearch.PotentialPair;
import org.primer3.p3_seq_lib.seq_lib;
import org.primer3.primer.MultiTargetPrimerRecord;
import org.primer3.primer.PrimerRecord;
import org.primer3.search.multiplex.MultiplexResult;
import org.primer3.sequence.Sequence;

/**
 * helper class to scan multiple primer per sequence to group and sort then to multi target and one target specific primers
 * @author Ahmed Hafez
 *
 */
public class MultiTargetScanner {

//	HashMap<String, ArrayList<PrimerRecord>> targetsToLeft = new HashMap<String, ArrayList<PrimerRecord>>();
	static final String REVERSE_PREFIX_TARGETS = "reverse " ;
	private ArrayList<MultiTargetPrimerRecord> leftPrimers = new ArrayList<MultiTargetPrimerRecord>() ;
	private ArrayList<MultiTargetPrimerRecord> rightPrimers = new ArrayList<MultiTargetPrimerRecord>() ;
	int nTargets = 0;
	public ArrayList<String> targets =  new  ArrayList<String>(); 
	seq_lib targets_lib = new seq_lib();
	P3GlobalSettings pa;
	ArrayList<P3RetVal> revals = new ArrayList<P3RetVal>();
	
	IP3MultiSearch p3MultiplexSearcher;
	public MultiplexResult result;
	public MultiTargetScanner(IP3MultiSearch p3MultiplexSearcher, String gname ,  P3GlobalSettings pa)
	{
		this.p3MultiplexSearcher = p3MultiplexSearcher;
		this.pa = pa;
		this.result = new MultiplexResult(p3MultiplexSearcher,gname,pa);
	}
	
	public void addP3RetVal(P3RetVal reval) {
		
		
		this.revals.add(reval);
		String targetName = reval.sa.getSequenceName();
//		ArrayList<PrimerRecord> targetLeft =  new ArrayList<PrimerRecord>();
//		targetsToLeft.put(targetName, targetLeft);
		targetsToSpFwd.put(targetName, new ArrayList<PrimerRecord>());
		targetsToSpRev.put(targetName, new ArrayList<PrimerRecord>());
		targets.add(targetName);
		nTargets++;
		sortAndAdd(leftPrimers, targetName, reval.fwd.oligo);
		sortAndAdd(rightPrimers, targetName, reval.rev.oligo);
		targets_lib.add(reval.sa.getSequenceName(), reval.sa.getTrimmedSequence(),  Sequence.p3_reverse_complement( reval.sa.getTrimmedSequence()));
		targets_lib.add(REVERSE_PREFIX_TARGETS + reval.sa.getSequenceName() ,  reval.sa.getUpcasedSeqRev(), reval.sa.getUpcasedSeq());
	
		result.addP3RetVal(reval);
	}

	
	/**
	 * sort primer by sequence and group similar sequence together as one mutli-primer 
	 * TODO :: use bit pattern 
	 * @param multiTargetPrimers
	 * @param targetName
	 * @param oligo
	 */
	void sortAndAdd(ArrayList<MultiTargetPrimerRecord> multiTargetPrimers, String targetName,List<PrimerRecord> oligo) {
		for(int i = 0 ; i < oligo.size();i++)
		{
			PrimerRecord newLeft = oligo.get(i);
			int insertion = Collections.binarySearch(multiTargetPrimers, (Object)oligo.get(i).getOligoSeq(), new Comparator<Object>() {

				@Override
				public int compare(Object o1, Object o2) {
					if (o1 instanceof MultiTargetPrimerRecord && o2 instanceof char[]) {
						MultiTargetPrimerRecord r = (MultiTargetPrimerRecord) o1;
						char[] seq = (char[]) o2;
						return Sequence.compare(r.getOligoSeq(),seq);
					}
					return 0;
				}
			} );
			if(insertion < 0)
			{
				int insertionLocation = - (insertion + 1) ;
				multiTargetPrimers.add(insertionLocation,  new MultiTargetPrimerRecord(targetName,newLeft));
			}
			else
			{
				MultiTargetPrimerRecord multiTargetPrimerRecord = multiTargetPrimers.get(insertion);
				multiTargetPrimerRecord.addRecord(targetName, newLeft);
			}
		}
	}

	// I need one list for multitarget primers
	// one list for each target for specific primers
	
	/**
	 * multi target forward primer list
	 */
	public ArrayList<PrimerRecord> multiFwd = new ArrayList<PrimerRecord>();
	/**
	 * one target Specific forward primer list 
	 */
	ArrayList<PrimerRecord> spFwd = new ArrayList<PrimerRecord>();
	/**
	 * multi target reverse primer list
	 */
	ArrayList<PrimerRecord> multiRev = new ArrayList<PrimerRecord>();
	/**
	 * one target specific reverse primer
	 */
	ArrayList<PrimerRecord> spRev= new ArrayList<PrimerRecord>();
	/**
	 * specific forward primer for each target sequence
	 */
	public HashMap<String, ArrayList<PrimerRecord> > targetsToSpFwd = new HashMap<String, ArrayList<PrimerRecord>>();
	/**
	 * specific reverese primers for each target sequence
	 */
	public HashMap<String, ArrayList<PrimerRecord> > targetsToSpRev = new HashMap<String, ArrayList<PrimerRecord>>();
//	public Object sas;
	
	
	// three search stratagy 
	// 1 search specific left-to-right for one complex pair 
	
	
	
	


	public void prePare() throws AlignmentException {
		
		for(int i = 0 ;i < leftPrimers.size() ; i++)
		{
			if(leftPrimers.get(i).getNTarget() > 1 )
			{
				PrimerRecord  mPrimer = leftPrimers.get(i);
				
				mPrimer.calcSpecific(targets_lib,  REVERSE_PREFIX_TARGETS , pa.isLibAmbiguityCodesConsensus() ?  p3MultiplexSearcher.getDPAL_ArgToUse().local_end_ambig
						: p3MultiplexSearcher.getDPAL_ArgToUse().local_end);
				if(mPrimer.getIsTargetSpecific()) {
					multiFwd.add(mPrimer);
				}
			}
			else
			{
				PrimerRecord spPrimer = leftPrimers.get(i).getOrgPrimer();
				// not quite yet see if it target spific
				
				spPrimer.calcSpecific(targets_lib,spPrimer.getTargetName(),  REVERSE_PREFIX_TARGETS , pa.isLibAmbiguityCodesConsensus() ?  p3MultiplexSearcher.getDPAL_ArgToUse().local_end_ambig
							: p3MultiplexSearcher.getDPAL_ArgToUse().local_end);
				if(spPrimer.getIsTargetSpecific() && spPrimer.targetSpecificIndex.size() == 0) {
					spFwd.add(spPrimer);
					targetsToSpFwd.get(spPrimer.getTargetName()).add(spPrimer);
				}
			}
		}
		for(int i = 0 ;i < rightPrimers.size() ; i++)
		{
			
			if(rightPrimers.get(i).getNTarget() > 1 )
			{
				PrimerRecord  mPrimer = rightPrimers.get(i);
				
				mPrimer.calcSpecific(targets_lib,  REVERSE_PREFIX_TARGETS , pa.isLibAmbiguityCodesConsensus() ?  p3MultiplexSearcher.getDPAL_ArgToUse().local_end_ambig
						: p3MultiplexSearcher.getDPAL_ArgToUse().local_end);
				
				multiRev.add(rightPrimers.get(i));
			}
			else
			{
				PrimerRecord spPrimer = rightPrimers.get(i).getOrgPrimer();
				spPrimer.calcSpecific(targets_lib,spPrimer.getTargetName(),  REVERSE_PREFIX_TARGETS , pa.isLibAmbiguityCodesConsensus() ?  p3MultiplexSearcher.getDPAL_ArgToUse().local_end_ambig
						: p3MultiplexSearcher.getDPAL_ArgToUse().local_end);
				if(spPrimer.getIsTargetSpecific() && spPrimer.targetSpecificIndex.size() == 0) {
					spRev.add(spPrimer);
					targetsToSpRev.get(spPrimer.getTargetName()).add(spPrimer);
				}

			}
		}
		Comparator<PrimerRecord> primerComparator = new Comparator<PrimerRecord>() {

			@Override
			public int compare(PrimerRecord o1, PrimerRecord o2) {
				if (o1 instanceof MultiTargetPrimerRecord && o2 instanceof MultiTargetPrimerRecord ) {
					MultiTargetPrimerRecord  r1 = (MultiTargetPrimerRecord ) o1;
					MultiTargetPrimerRecord  r2 = (MultiTargetPrimerRecord ) o2;
					
					if(Math.abs(r1.quality - r2.quality ) <=  0.05 )
					{
						return -Integer.compare(r1.getNTarget(), r2.getNTarget());
					}
					return PrimerRecord.compare(r1, r2);
					
				}
				return PrimerRecord.compare(o1, o2);
			}
		};
		multiFwd.sort(primerComparator);
		multiRev.sort(primerComparator);
		
		for( ArrayList<PrimerRecord> spPrimerList : targetsToSpFwd.values())
		{
			spPrimerList.sort(primerComparator);
		}
		
		for( ArrayList<PrimerRecord> spPrimerList : targetsToSpRev.values())
		{
			spPrimerList.sort(primerComparator);
		}
		
		
		
		// TODEBUG ::
		// No need for them here
		leftPrimers.clear();;
		rightPrimers.clear();
		leftPrimers = null;
		rightPrimers = null;
		
	}

	
//	static private void incI(Iterator<Integer>[] its, Integer[] res, List<Integer>[] sources) {
//		for(int i = its.length-1; i > 0;i-- )
//		{
//			if(its[i].hasNext()  ){ 
//				res[i] = its[i].next();
//				return;
//			} 
//			else 
//			{
//				its[i] = sources[i].iterator();
//				res[i] = its[i].next();
//			}
//		}
//		res[0] = its[0].next();
//		
//	}
//	static private boolean ifAll(Iterator<Integer>[] its)
//	{
//		
//		for(int i = its.length-1; i >= 0;i-- )
//		{
//			if(its[i].hasNext())
//				return true;
//		}
//		return false;
//	}
//	static public void main(String[] argv)
//	{
//		int nTargets = 3;
//		int[] is = new int[nTargets];
//		int[] js = new int[nTargets];
//		int[] fwdSizes = new int[nTargets];
//		int[] revSizes = new int[nTargets];
//		List<Integer>[] sources = new List[nTargets];
//		Iterator<Integer>[] its = new Iterator[nTargets];
//		Integer[] res = new Integer[nTargets];
//		for(int i =0 ; i < nTargets;i++  )
//		{
//			sources[i]= new ArrayList<Integer>();
//			sources[i].add(1);
//			sources[i].add(2);
//			sources[i].add(3);
//			its[i] = sources[i].iterator();
//			res[i] = its[i].next();
//		}
//		for(; ifAll(its); incI(its,res,sources) )
//		{
//			for(int i =0 ; i < nTargets;i++)
//				System.err.print(res[i] + " ");
//			System.err.println();
//		}
//	}
	
}
