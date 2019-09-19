package org.pooler;

import java.io.PrintStream;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

import org.pooler.Fasta2Bit.SeqData;
import org.pooler.Fasta2Bit.SeqInfo;

public class Amplicons {

	public class PrimerToFind  implements java.lang.Comparable<Long> {
		long lhsBases; /* = p & minValid, for the bsearch (is actually the LAST len(minValid)/2 bases before the cursor, but we're reading backwards from the left see below) */
		public long p;
		public long valid; /* what we're looking for */
		int ampNo, onOrOff; /* what to do when we find it */
		int primerNo; 
		String name; /* for reports */
		@Override
		public int compareTo(Long k) {
			return Long.compareUnsigned( this.lhsBases,k);
			//					  if(k < this.lhsBases) return -1;
			//					  else return k > this.lhsBases;
			//			return 0;
		}
	}
	class AmpEvent {
		long baseEnd,baseStart; 
		int ampNo, onOrOff;
		String name; /* for reports */
	}
	HashMap<Integer, ArrayList<AmpEvent>> eventLists = new HashMap<>(); 
//	HashMap<String, Amplicon> amps = new HashMap<>(); 
	public List<Amplicon> amps = new ArrayList<>(); 


	AllPrimers ap;
	
	
	String genomeFile; 
	int maxAmpliconLen;
	String allAmpsFileName;
	// to report/output an input file for MultiPLX program
	// if 0 just write locations to the file
	int allAmpsIsMultiPLX;

	public int[] primerNoToAmpliconNo;
	public int nAmplicons;
	public Amplicons(AllPrimers ap, String genomeFile,
			int maxAmpliconLen,
			String multiplexFile, int allAmpsIsMultiPLX) {

		this.ap = ap;
		this.genomeFile = genomeFile;
		this.maxAmpliconLen = maxAmpliconLen;
		this.allAmpsIsMultiPLX = allAmpsIsMultiPLX;
		this.allAmpsFileName = multiplexFile;
	}
	PrintStream reportFileStream = null;
	public void setReportFile(PrintStream psReportFile) {
		this.reportFileStream= psReportFile;
		
	}


	PrimerToFind [] primerVariantsToFind;
	//	int nPrimerVariantsToFind;
	int[] didFind;
	public boolean[] getOverlappingAmplicons(PrintStream errStream) {
		
		PrintStream oldErrStream = System.err;
		if(errStream != null)
			System.setErr(errStream);
		
		
		Instant start = Instant.now();
		nAmplicons = amplicons_from_primer_names();
		boolean[] overlaps = null;
		if(nAmplicons > 0) {
			primerVariantsToFind = populatePF(nAmplicons);
			didFind = new int[ap.np];
			List<SeqInfo> names = Genome.go_throgh_genome(this);
			if(errStream != null)
				System.setErr(oldErrStream);

			checkPF(reportFileStream);
			overlaps = eventsToOverlaps(nAmplicons,
					names,
					reportFileStream,
					System.out,
					null);

		}
		return overlaps;
	}


	private PrimerToFind[] populatePF(int nAmp) {
		/* For how many primers do we want to search through the genome?
	     (I'm assuming degenerate bases are rare so we can write out all combinations here;
	     TODO: if this runs out of memory, try handling the degenerate primers with brute-force matching instead) */
		int ampNo,nToFind=0; 
		for(ampNo=0; ampNo<nAmp; ampNo++) 
			nToFind += ap.getNumPossibilities(ampliconNoToFwd[ampNo]) + ap.getNumPossibilities(ampliconNoToRev[ampNo]);
		nToFind *= 2;
		PrimerToFind[] pF=new PrimerToFind[nToFind];
		int outP = 0;
		boolean doneTruncate = false;
		for(ampNo=0; ampNo<nAmp; ampNo++) {
			int isFwd; 
			for(isFwd=0; isFwd<=1; isFwd++) {
				int primerNo = (isFwd == 1 ? ampliconNoToFwd:ampliconNoToRev )[ampNo];
				//			#ifdef Debug_AmpliconNo
				//			      if(Debug_AmpliconNo(ampNo)) {
				//			        fprintf(stderr,"p#=% 8d ",primerNo); printBasesMaybeD(ap,primerNo,stderr); fputs(":\n",stderr);
				//			      }
				//			#endif
				int nPo= ap.getNumPossibilities(primerNo);
				int i;
				boolean didTruncate=false;
				int reverseAndComplement;
				for(i=0; i<nPo; i++) {
					for(reverseAndComplement=0; reverseAndComplement<=1; reverseAndComplement++) {
						PrimerToFind newPrimerPoss =  new PrimerToFind();
						didTruncate = ap.Make2bit(primerNo,reverseAndComplement,reverseAndComplement,newPrimerPoss,i,nPo);
						//			#ifdef Debug_AmpliconNo
						//			          if(Debug_AmpliconNo(ampNo)) {
						//			            fprintf(stderr,"%d/%d (rc=%d): ",i,nPo,reverseAndComplement); debugPrnRTL(pF[outP].p,pF[outP].valid); fprintf(stderr," -> %d\n",outP); /* offset in primerVariantsToFind */
						//			          }
						//			#endif
						pF[outP]  = newPrimerPoss;
						pF[outP].ampNo = ampNo; 
						pF[outP].onOrOff = ((isFwd==reverseAndComplement)?-1:1)*(isFwd == 1 ?1:-1) * (1+ (isFwd == reverseAndComplement ? 1 :0)); // this should get +1,-1 for +ve strand 1st, or +2,-2 for -ve strand 1st
						pF[outP].primerNo = primerNo;
						pF[outP++].name= ap.names.get(primerNo);
					}
				} 
				if(didTruncate) {
					if(!doneTruncate) {
						System.err.println(""); 
						doneTruncate = true; 
					}
					System.err.format("%s\n",ap.names.get(primerNo));

				}
			}
		}
		//		fflush(stderr);
		minValid=~(long)0; int i;
		for(i=0; i<nToFind; i++) {
			if(pF[i].valid < minValid) {
				minValid = pF[i].valid;
			}
		}
		for(i=0; i<nToFind; i++)
			pF[i].lhsBases = pF[i].p & minValid;
		// qsort(pF,nToFind,sizeof(PrimerToFind),PF_cmp_func);
		Arrays.sort(pF, new Comparator<PrimerToFind>() {

			@Override
			public int compare(PrimerToFind o1, PrimerToFind o2) {
				return Long.compareUnsigned(o1.lhsBases, o2.lhsBases);
			}
		});
		// *nRet = nToFind; 
		//return pF;
		return pF;
	}



	int []ampliconNoToFwd = null;
	int []ampliconNoToRev = null;
	private int amplicons_from_primer_names() {
		ampliconNoToFwd=new int[(ap.np/2)];
		ampliconNoToRev=new int[(ap.np/2)];
		primerNoToAmpliconNo=new int[(ap.np)];
		for(int i=0;i < primerNoToAmpliconNo.length;i++) {
			primerNoToAmpliconNo[i]=-1;
		}
		int nAmp=0, foundForward=0, foundBackward=0;
		int i,jStart=0;
		for(i=0; i<ap.np; i++) {
			String n = ap.names.get(i); 
			int l = n.length();
			if(l > 0 && (n.charAt(l-1) =='F' || n.charAt(l-1)=='f')) {
				foundForward = 1;
				if(i==jStart) 
					jStart++; /* optimisation for if they're arranged as (f,r)(f,r)(f,r) in the input */

				for(int j=jStart; j<ap.np; j++) {
					String nAtJ = ap.names.get(j);
					if(j!=i && nAtJ.length() == l && 
							nAtJ.regionMatches(0, n, 0, l-1) &&
							"BbRr".indexOf(nAtJ.charAt(l-1)) >= 0 // && strchr("BbRr",ap.names[j][l-1]) &&
							//!ap.names[j][l]
							) {
						/* Do check for -B or -R, don't just take anything that's not -F, because
		             they might have TaqMan probes as -P or something (new in v1.36; previous
		             versions would be OK with this in the pooling stage but this overlap-check
		             stage could get confused) */
						foundBackward = 1;
						//		          #ifdef Debug_PrimerNamePrefix
						//		          if(!strncmp(ap.names[i],Debug_PrimerNamePrefix,sizeof(Debug_PrimerNamePrefix)-1))
						//		            fprintf(stderr,"Amplicon %d is P%d (%s) and %d (%s)\n",nAmp,i,ap.names[i],j,ap.names[j]);
						//		          #endif
						Amplicon newAmp = new Amplicon();
						(ampliconNoToFwd)[nAmp] = i;
						(ampliconNoToRev)[nAmp] = j;
						newAmp.fwdIndex = i ;
						newAmp.revIndex = j;
						newAmp.name = n.substring(0, l-2);
//						amps.put(newAmp.name, newAmp);
						amps.add( newAmp);
						(primerNoToAmpliconNo)[i] = (primerNoToAmpliconNo)[j] = nAmp++;
						if(j==jStart && jStart++ == i+1) 
							++i; /* NOT i = new jStart, because the loop counter will also incr it */
						break;
					}
				}
			}
		}
		//		if(!nAmp) { free(*ampliconNoToFwd); free(*ampliconNoToRev); }
		//		SetColour(Bright,Cyan,Black); 
		System.err.format("%d primer-sets found\n",nAmp); 
		//		ResetColour();
		boolean reported_missing = false;
		if(nAmp > 0 ) {
			for(i=0; i<ap.np; i++) {
				if((primerNoToAmpliconNo)[i] == -1) {
					if(!reported_missing) {
						System.err.format("The following primers were not matched into amplicon sets:\n(should you check the spelling and capitalisation?)\n");
						reported_missing = true;
					}
					System.err.format("%s\n",ap.names.get(i));
				}
			}
		} else if (foundForward == 0  || foundBackward == 0) {
			/* 2016-09: a user in India thought you had to end
		       primers with '-forward' and '-backward'.  Are my
		       instructions really that unclear?  :-(
		       (TODO: maybe we could parse the whole word after
		       the final hyphen.  but would need to ensure other
		       users understand if this happens.)  */
			System.err.format("Please end your forward primers with -F\nand your reverse primers with -R or -B as instructed\n");
		}
		return nAmp;
	}





	long minValid = 0;
	public void lookup(long currentBuf, long currentValid, SeqData seqData) {
		/* Called for every position in the genome:
	     find which (if any) of the primers we're looking for
	     ends at this position, and deal with it.
	     For ease of binary search, bases are shifted into
	     buf from the LEFT (not from the right as in 64.h).
	     We'll do a binary search on our SHORTEST primer len
	     and linear search the rest (saves having to launch a
	     separate binary search for each length: it's rare to
	     have a collection of primers where many share the
	     same tail that's as long as the shortest primer in
	     the set, so the linear search should be quite short)
		 */
		//		#ifdef Debug_BaseCheck
		//		  if(baseEnd==Debug_BaseCheck) { debugPrnRTL(buf,valid); fprintf(stderr,"\n"); }
		//		#endif
		long baseEnd = seqData.currentBaseNo;
		long lhsBases = seqData.currentBuf & minValid;
//		if (lhsBases == 6367499821904822272l)
//			System.out.println("6367499821904822272");
		int foundIndex = Arrays.binarySearch(primerVariantsToFind, lhsBases );

		//		PrimerToFind found = bsearch(lhsBases,primerVariantsToFind,nPrimerVariantsToFind,sizeof(*primerVariantsToFind),PF_srch_func);
		if (foundIndex < 0) 
			return; /* (as soon as possible) */
		/* (bsearch might not land at the START of a section with equal lhs bases) */
		//		while(foundIndex > primerVariantsToFind && lhsBases==found[-1].lhsBases) found--; 

		for(; foundIndex < primerVariantsToFind.length && primerVariantsToFind[foundIndex].lhsBases == lhsBases; 
				foundIndex++) 
		{
			if((currentBuf & primerVariantsToFind[foundIndex].valid) == primerVariantsToFind[foundIndex].p 
					&& (primerVariantsToFind[foundIndex].valid & currentValid) == primerVariantsToFind[foundIndex].valid) {
				/* OK, all up-to-32 bases match */
				didFind[primerVariantsToFind[foundIndex].primerNo]++;
				ArrayList<AmpEvent> l = eventLists.get(seqData.seq.index);
				if(l == null ) {
					l = new ArrayList<AmpEvent>();
					eventLists.put(seqData.seq.index,l);
				}
				//				assert(seqNo < numEventLists);
				//				if(l->ptr == l->size) {
				//					l->size = (l->size | 1) << 1;
				//					AmpEvent *events2 = (l->size > l->ptr) ? realloc(l->events,l->size*sizeof(AmpEvent)) : NULL; /* NULL if integer overflow (hopefully unlikely!) */
				//					if(memFail(events2,_memFail)) {
				//						fputs("Won't be able to record any more amplicon events\n",stderr); // TODO: and stop scanning the genome (but without making the error-handling code slow it down; use longjmp??)
				//						l->size = l->ptr;
				//					} else l->events = events2;
				//				} 
				//				if(l->ptr < l->size) {
				//					#ifdef Debug_AmpliconNo
				//					extern SeqName lastSequenceNameRead;
				//					if(Debug_AmpliconNo(found->ampNo))
				//						fprintf(stderr,"Event #%lu: amplicon %d state %d baseEnd=%s:%u \n",l->ptr,found->ampNo,found->onOrOff,lastSequenceNameRead,baseEnd);
				//					#endif
				AmpEvent newEvent = new AmpEvent();

				newEvent.baseEnd = baseEnd;
				newEvent.baseStart = baseEnd - (64-Long.numberOfTrailingZeros(primerVariantsToFind[foundIndex].valid))/2 + 1; /* +1 to bring into line with what UCSC browser does */
				newEvent.ampNo = primerVariantsToFind[foundIndex].ampNo;
				newEvent.onOrOff = primerVariantsToFind[foundIndex].onOrOff;
				newEvent.name = primerVariantsToFind[foundIndex].name;
				l.add(newEvent);
				//				}
			}
		}
	}


	/**
	 * Check if amplicons are in the genome and report to a file
	 * @param reportFile
	 */
	void checkPF(PrintStream reportFile) {
		int nPrimers = ap.np;
		//		  qsort(primerVariantsToFind,nPrimerVariantsToFind,sizeof(PrimerToFind),PF_by_didFind);
		Arrays.sort(primerVariantsToFind, new Comparator<PrimerToFind>() {

			@Override
			public int compare(PrimerToFind o1, PrimerToFind o2) {
				/* this is for the checkPF report at the end */
				int r = didFind[o1.primerNo] - didFind[o2.primerNo]; /* <0 if A comes first */
				if(r == 0) 
					r = o1.primerNo - o2.primerNo;
				return r;
			}
		});
		//		  FILE *reportFile = NULL;
		int i,nFound=0;
		for(i=0; i<nPrimers; i++) {
			if(didFind[i] != 0) nFound++;
		}
		//		  #ifndef Debug_ChromosomeCheck
		/* (no point printing all if reading 1 chromosome) */
		boolean reported=false;
		for(i=0; i< primerVariantsToFind.length; i++) {
			if( didFind[primerVariantsToFind[i].primerNo] == 0) {
				if(!reported) {
					reported = true; 
					System.err.println("The following primers were not found in the genome:");
					//		        reportFile = fopen(getReportFilename(),"w");
					if(reportFile!=null)
						reportFile.println("The following primers were not found in the genome:");
					// if(reportFile) fputs("The following primers were not found in the genome:\n",reportFile);
					/* Previously reported number found more than once
		           	as well, but this turned out not to be useful
		           	because some primers can occur frequently but
		           	not in their amplicons.  We could make didFind
		           	boolean rather than a counter. */
				}
				System.err.format("%s\n",primerVariantsToFind[i].name);
				if(reportFile!=null)
					reportFile.format("%s\n",primerVariantsToFind[i].name);
				//		      if(reportFile) fprintf(reportFile,"%s\n",primerVariantsToFind[i].name);
				didFind[primerVariantsToFind[i].primerNo] = 1; /* so we don't report this one a second time (as multiple entries in primerVariantsToFind map to the same primerNo) */
			}
		}
		//		  #endif
		if (nFound == nPrimers) {
			System.err.format("All %d primers were found in the genome\n",nPrimers);
			if(reportFile!=null)
				reportFile.format("All %d primers were found in the genome\n",nPrimers);
			//		  if(reportFile) fprintf(reportFile,"All %d primers were found in the genome\n",nPrimers);
		} 
		else {
			System.err.format("%d of %d primers were found in the genome\n",nFound,nPrimers);
			if(reportFile!=null)
				reportFile.format("%d of %d primers were found in the genome\n",nFound,nPrimers);
			//		  if(reportFile) fprintf(reportFile,"%d of %d primers were found in the genome\n",nFound,nPrimers);
		}
		//		  return reportFile;
	}


	public void allocateSeqs(int nSeq) {



	}
	public int nOverlaps = 0;
	boolean[] eventsToOverlaps( int nAmp,
			List<SeqInfo> names,
			PrintStream reportFileP,
			PrintStream allAmps,
			Object genome
			) {
		char[] ampsFound= new char[nAmp];
		boolean[] overlaps=  new boolean[nAmp*nAmp] ;// calloc(nAmp,nAmp); /* TODO: use triangle.h instead? */
		int[] inProgress= new int[nAmp];
		int[] inProgressI= new int [nAmp];
		//		  if(memFail(ampsFound,overlaps,inProgress,inProgressI,_memFail)) return NULL;
		nOverlaps = 0;
		if (allAmpsIsMultiPLX != 0) 
			ap.addTags(); // prior to v1.33 this was incorrectly placed below the start of the seqNo loop, resulting in additional copies of tags being added for each new chromosome (usually overflowing the selected bit size so you get only the last part of the 2nd tag)
		long maxLenFound = 0;
		int seqNo; 
		for(seqNo=0; seqNo <  eventLists.size()  ; seqNo++) {
			//		  memset(inProgress,0,nAmp*sizeof(int));
			//		  memset(inProgressI,0,nAmp*sizeof(int));
			for (int _i=0;_i<nAmp;_i++)
			{
				inProgress[_i] = inProgressI[_i] = 0;
			}
			ArrayList<AmpEvent> events = eventLists.get (  names.get(seqNo).index );
			
			//		  qsort(events,eventLists[seqNo].ptr,sizeof(AmpEvent),eventOrder);
			int i; 
			for(i=0;  events != null && i< events.size(); i++) {
				int onOrOff = events.get(i).onOrOff,
						ampNo = events.get(i).ampNo;
				if(onOrOff < 0) { /* it's an end event */
					if((inProgress[ampNo] & -onOrOff)  != 0)
						inProgress[ampNo] += onOrOff;
				} 
				else if((inProgress[ampNo] & onOrOff) != 0) {
					/* already in progress: ignore */
				} 
				else {
					int end = findEndEvent(events,i,maxAmpliconLen);
					if(end == 0) 
						continue;
					ampsFound[ampNo] = 1;
					inProgress[ampNo] += onOrOff;
					inProgressI[ampNo] = i;
					long ampLength = events.get(end).baseEnd + 1 - events.get(i).baseStart;
					if (ampLength>maxLenFound) maxLenFound= ampLength;
					if(allAmps != null ) {
						if (allAmpsIsMultiPLX != 0) {
//							/* TODO: what if duplicate ampsFound[ampNo] ? */
							AmpEvent sEvent = events.get(i);
							String sEventName = sEvent.name;
							int subIndex = sEventName.charAt(sEventName.length()-2) == '-' || sEventName.charAt(sEventName.length()-2) == '_'  ? 1 : 0 ;
							allAmps.format("%s",sEventName.substring(0, sEventName.length()-1-subIndex));
							// fprintf(allAmps,"%.*s",(int)strlen(events[i].name)-1-(strchr("-_",events[i].name[strlen(events[i].name)-2])!=0),events[i].name);
							allAmps.print("\t");
							//	fputc('\t',allAmps);
							ap.forward.get(ampliconNoToFwd[sEvent.ampNo]).printBases(allAmps);
							// printBasesMaybeD(ap,ampliconNoToFwd[events[i].ampNo],allAmps);
							allAmps.print("\t");
							// fputc('\t',allAmps);
							ap.forward.get(ampliconNoToRev[sEvent.ampNo]).printBases(allAmps);
							// printBasesMaybeD(ap,ampliconNoToRev[events[i].ampNo],allAmps);
							allAmps.print("\t");
							// fputc('\t',allAmps);
//							output_genome_segment(genome,seqNo,events[i].baseStart,ampLength,allAmps);
							allAmps.print("\n");
							// fputc('\n',allAmps);
						} 
						else {
//							fprintf(allAmps,"%s:%s (%s:%u%c%u)\n",events[i].name,events[end].name,&(names[seqNo][0]),events[i].baseStart,((onOrOff==1)?'+':'-'),events[end].baseEnd);
							allAmps.format("%s:%s (%s:%d%c%d)\n",events.get(i).name,events.get(end).name,names.get(seqNo).seqname , events.get(i).baseStart,((onOrOff==1)?'+':'-'),events.get(end).baseEnd);
						}
						
						Amplicon amp = amps.get(events.get(i).ampNo);
						amp.setLocation(names.get(seqNo).seqname, events.get(i).baseStart, events.get(end).baseEnd, onOrOff==1);
						
					}
					/* now check if there's any overlaps with other amplicons that are already running */
					int prnOver=0, j;
					for(j=0; j<nAmp; j++)
						if(j!=ampNo && inProgress[j] != 0 && overlaps[ampNo*nAmp+j] == false ) {
							overlaps[ampNo*nAmp+j] =
									overlaps[j*nAmp+ampNo] = true;
							if(prnOver != 0) 
								System.err.print(", "); 
							else {
								if(nOverlaps == 0) {
									System.err.print("Overlapping amplicons:\n");
//									if(!*reportFileP) *reportFileP = fopen(getReportFilename(),"w");
//									if(*reportFileP) fputs("Overlapping amplicons:\n",*reportFileP);
									if(reportFileP != null)
										reportFileP.print("Overlapping amplicons:\n");
								}
								System.err.format("%s:%s (%s:%d%c%d) / ",
										events.get(i).name,
										events.get(end).name,
										names.get(seqNo).seqname,
										events.get(i).baseStart,((onOrOff==1)?'+':'-'),
										events.get(end).baseEnd);
//								if(*reportFileP) fprintf(*reportFileP,"%s:%s (%s:%u%c%u) / ",events[i].name,events[end].name,&(names[seqNo][0]),events[i].baseStart,((onOrOff==1)?'+':'-'),events[end].baseEnd);
								if(reportFileP != null)
									reportFileP.format("%s:%s (%s:%d%c%d) / ",
											events.get(i).name,
											events.get(end).name,
											names.get(seqNo).seqname,
											events.get(i).baseStart,((onOrOff==1)?'+':'-'),
											events.get(end).baseEnd);
									
								prnOver=1;
							}
							nOverlaps++;
							AmpEvent overlapWithStart = events.get(inProgressI[j]),
							overlapWithEnd = events.get( findEndEvent(events,inProgressI[j],maxAmpliconLen));
							System.err.format("%s:%s (%d%c%d)",overlapWithStart.name,overlapWithEnd.name,overlapWithStart.baseStart,((overlapWithStart.onOrOff==1)?'+':'-'),overlapWithEnd.baseEnd);
							if(reportFileP != null)
								reportFileP.format("%s:%s (%d%c%d)",overlapWithStart.name,overlapWithEnd.name,overlapWithStart.baseStart,((overlapWithStart.onOrOff==1)?'+':'-'),overlapWithEnd.baseEnd);
							Amplicon amp = amps.get(events.get(i).ampNo);
							Amplicon overlappingAmp = amps.get(overlapWithStart.ampNo);
							amp.addOverlapingAmp(overlappingAmp);
//							if(*reportFileP) 
//								fprintf(*reportFileP,"%s:%s (%u%c%u)",overlapWithStart->name,overlapWithEnd->name,overlapWithStart->baseStart,((overlapWithStart->onOrOff==1)?'+':'-'),overlapWithEnd->baseEnd);
						}
					if(prnOver != 0) {
						System.err.print("\n");
						if(reportFileP != null)
							reportFileP.print("\n");
//						if(*reportFileP) fputs("\n",*reportFileP);
					}
				}
			}
		} /* finished - rest of this function is reporting */
		if (allAmpsIsMultiPLX !=  0) ap.removeTags();
		int numNotFound = 0, i;
		for(i=0; i<nAmp; i++)
			if( ampsFound[i] == 0) {
				if(numNotFound == 0) {
//					SetColour(Bright,Red,Black); 
					System.err.format("Amplicons not found in genome:\n"); 
//					ResetColour();
//					if(*reportFileP) fprintf(*reportFileP,"Amplicons not found in genome:\n");
					if(reportFileP != null)
						reportFileP.format("Amplicons not found in genome:\n"); 
						
				}
				numNotFound++;
//				#ifdef Debug_AmpliconNo
//				if(!Debug_AmpliconNo(i)) continue; /* don't let these get lost in the others */
//				#endif
				System.err.format("%s/%s\n",ap.names.get(ampliconNoToFwd[i]),
						ap.names.get(ampliconNoToRev[i]));
				if(reportFileP != null)
					reportFileP.format("%s/%s\n",ap.names.get(ampliconNoToFwd[i]),
							ap.names.get(ampliconNoToRev[i]));
//				if(*reportFileP) fprintf(*reportFileP,"%s/%s\n",ap.names[ampliconNoToFwd[i]],ap.names[ampliconNoToRev[i]]);
			}
		if(numNotFound != 0) {
//			SetColour(Bright,Red,Black); 
			System.err.format(maxAmpliconLen != 0 ? "%d of %d amplicons not found in genome (with length <=%d)\n":"%d of %d amplicons not found in genome\n",numNotFound,nAmp,maxAmpliconLen); 
//			ResetColour();
			if(reportFileP != null)
				reportFileP.format(maxAmpliconLen != 0 ? "%d of %d amplicons not found in genome (with length <=%d)\n":"%d of %d amplicons not found in genome\n",numNotFound,nAmp,maxAmpliconLen); 
//			if(*reportFileP) fprintf(*reportFileP,maxAmpliconLen?"%d of %d amplicons not found in genome (with length <=%d)\n":"%d of %d amplicons not found in genome\n",numNotFound,nAmp,maxAmpliconLen);
		} 
		else { 
//			SetColour(Dark,Green,Black); 
			System.err.print("All amplicons were found in the genome\n"); 
//			ResetColour(); 
			} 
//		free(ampsFound);
		System.err.format("Longest amplicon found: %d bp\n",maxLenFound);
		if(maxLenFound >= 100000 /* 50000 was longest on record and very difficult to achieve */) 
			/* relevant only if we're running on --amp-max=0 (v1.35+) */
			System.err.format("  - this might include coincidental matches on other chromosomes;\n    try setting a sensible --amp-max\n"); 
//		SetColour(Bright,Cyan,Black); 
		System.err.format("%d overlaps found\n",nOverlaps); 
//		ResetColour();
		if(reportFileP != null)
			reportFileP.format("%d overlaps found\n",nOverlaps);
//		if(*reportFileP) fprintf(*reportFileP,"%d overlaps found\n",nOverlaps);
//		free(inProgress); 
//		free(inProgressI); 
		return overlaps;
	}

	int findEndEvent(List<AmpEvent>  events,int startEvent,int maxAmpliconLen) {
//		  List<AmpEvent> events = eventLists.get(seq);
//		  #ifdef Debug_AmpliconNo
//		  if(Debug_AmpliconNo(events[startEvent].ampNo) && maxAmpliconLen)
//		    fprintf(stderr,"findEndEvent (%d/%d): baseStart=%u so max baseEnd=%u\n",events[startEvent].ampNo,events[startEvent].onOrOff,events[startEvent].baseStart,events[startEvent].baseStart+maxAmpliconLen);
//		  #endif
		  int i; 
		  for(i=startEvent+1; i< events.size(); i++) {
//		    #ifdef Debug_AmpliconNo
//		    if(Debug_AmpliconNo(events[startEvent].ampNo))
//		      fprintf(stderr," - checking #%d: %d/%d (ends %u)\n",i,events[i].ampNo,events[i].onOrOff,events[i].baseEnd);
//		    #endif
		    if(maxAmpliconLen != 0 && events.get(startEvent).baseStart + maxAmpliconLen < events.get(i).baseEnd) break; /* gone too far ahead */
		    if(events.get(i).ampNo ==  events.get(startEvent).ampNo) {
		      if (events.get(i).onOrOff == -events.get(startEvent).onOrOff) 
		    	  return i; /* found it */
		      else 
		    	  break; /* duplicate start??  abort now and take second (shorter) one (v1.35 added) */
		    }
		  }
//		  #ifdef Debug_AmpliconNo
//		  if(Debug_AmpliconNo(events[startEvent].ampNo))
//		    fprintf(stderr," - not found (lookahead=%d)\n",i-startEvent);
//		  #endif
		  return 0;
		}


	
}
