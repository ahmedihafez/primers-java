/*
  This file is part of Java porting of Primer Pooler (https://github.com/ssb22/PrimerPooler)
  Primer Pooler (c) Silas S. Brown.  For Wen.
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
    http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
package org.pooler;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map.Entry;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.io.DNASequenceCreator;
import org.biojava.nbio.core.sequence.io.FastaReader;
import org.biojava.nbio.core.sequence.io.GenericFastaHeaderParser;
import org.pooler.Amplicons.PrimerToFind;
import org.pooler.Triangle.TwoRanges;
import org.pooler.primers.CountResult;
import org.pooler.primers.IPrimer;
import org.pooler.primers.PrimerFactory;
import org.pooler.primers64.PrimerFactory64;
import org.primer3.p3_seq_lib.DNACharSet;

public class AllPrimers {

	static public  int  InvalidCombination=0x4000; 



	List<IPrimer> 	forward, // array of primers read "forward" (5'-3') 
	backward, // same primers read "backward" (3'-5')
	tags;  // array of tags (read 5'-3')
	List<Integer> whichTag; //  which tag no. applies to primer N (-1 for none) 
	List<String> names; // primer names

	int np; // total number of primers
	int maxLen = 0 , minLen = Integer.MAX_VALUE , maxTag= 0; // length of longest primer 

	PrimerFactory currentFactory = null;

	public void loadFasta(String fileName) throws Exception {

		whichTag = new ArrayList<Integer>();
		forward = new ArrayList<IPrimer>();
		backward = new ArrayList<IPrimer>();
		tags = new ArrayList<IPrimer>();
		names = new ArrayList<String>();

		// here decide which currentFactory = PrimerFactory.getCurrent();
		PrimerFactory.setCurrent( new PrimerFactory64());
		currentFactory = PrimerFactory.getCurrent();



		DNACharSet dnaSet = new DNACharSet();

		FastaReader<DNASequence, NucleotideCompound> fastaReader;
		try {
			fastaReader = new FastaReader<DNASequence, NucleotideCompound>(
					new FileInputStream(new File(fileName)),
					new GenericFastaHeaderParser<DNASequence, NucleotideCompound>(),
					new DNASequenceCreator(dnaSet));

			int p = 0 ; // number of seq so far
			HashMap<Character, Integer> lastByte_to_whichTag = new HashMap<Character, Integer>();
			HashMap<Character, Boolean> check_not_last = new HashMap<Character, Boolean>();
			int nextTag = 0;
			while(true) {
				LinkedHashMap<String, DNASequence> dbFile = fastaReader.process(1);
				if(dbFile == null || dbFile.size() == 0 )
					break;
				for(Entry<String,DNASequence> seq : dbFile.entrySet())
				{
					System.out.println(">" + seq.getKey());
					System.out.println(seq.getValue().getSequenceAsString());
					IPrimer mdp = currentFactory.createPrimer(seq.getValue().getSequenceAsString());
					String seqName =  seq.getKey();
					if (seqName.startsWith("tag") && seqName.length() == 4 ) {
						char tagType = Character.toUpperCase(seqName.charAt(3));
						if(!lastByte_to_whichTag.containsKey(tagType))
						{
							for(int i = 0 ; i < whichTag.size() ; i++)
							{
								if ( tagType ==  
										Character.toUpperCase( 
												names.get(i).charAt(names.get(i).length()-1) ) 
										) {
									whichTag.set(i,nextTag);
								}

							}
						} 
						else 
							check_not_last.put(tagType, true);
						lastByte_to_whichTag.put(tagType, nextTag);
						tags.add(mdp);
						nextTag++;
					}
					else {
						names.add(seqName);
						char tagType =  Character.toUpperCase(seqName.charAt(seqName.length()-1) );
						int whichTagValue = -1;
						if(lastByte_to_whichTag.containsKey(tagType))
							whichTagValue = lastByte_to_whichTag.get(tagType);
						whichTag.add(whichTagValue);
						check_not_last.put(tagType,false);
						forward.add(mdp);
					}


				}
			}

			p = 0;
			for(Entry<Character,Boolean> checkTag : check_not_last.entrySet() )
			{
				if(checkTag.getValue()) {
					if(p != 0) 
					{
						//fprintf(stderr,"WARNING: Same applies to >tag%c\n",nextTag);
					}
					else {
						p = 1;
						//fprintf(stderr,"\nWARNING: You have multiple >tag%c sequences\n         and the last one does not precede a >...%c primer.\n         This probably means you've made a mistake.\n         Apart from the first >tag%c, all >tag%c tags will apply to\n         >...%c primers AFTER the >tag%c (not before it).\n",nextTag,nextTag,nextTag,nextTag,nextTag,nextTag);
					}
				}
			}


			for(int i = 0 ;  i < forward.size();i++) {
				this.backward.add(forward.get(i).getReverse());
			}


			this.np = this.forward.size();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


	}

	
	public void loadSequences(ArrayList<String> seqNames,ArrayList<String> seqs) {

		whichTag = new ArrayList<Integer>();
		forward = new ArrayList<IPrimer>();
		backward = new ArrayList<IPrimer>();
		tags = new ArrayList<IPrimer>();
		names = new ArrayList<String>();

		// here decide which currentFactory = PrimerFactory.getCurrent();
		PrimerFactory.setCurrent( new PrimerFactory64());
		currentFactory = PrimerFactory.getCurrent();




		try {


			int p = 0 ; // number of seq so far
			HashMap<Character, Integer> lastByte_to_whichTag = new HashMap<Character, Integer>();
			HashMap<Character, Boolean> check_not_last = new HashMap<Character, Boolean>();
			int nextTag = 0;
			for(int seqIndex = 0 ; seqIndex < seqs.size();seqIndex++) {
				System.out.println(">" + seqNames.get(seqIndex));
				System.out.println(seqs.get(seqIndex));
				int seqLen = seqs.get(seqIndex).length();

				
				IPrimer mdp = currentFactory.createPrimer(seqs.get(seqIndex));
				String seqName = seqNames.get(seqIndex);
				if (seqName.startsWith("tag") && seqName.length() == 4 ) {
					
					
					if(maxTag < seqLen)
						maxTag = seqLen;
					char tagType = Character.toUpperCase(seqName.charAt(3));
					if(!lastByte_to_whichTag.containsKey(tagType))
					{
						for(int i = 0 ; i < whichTag.size() ; i++)
						{
							if ( tagType ==  
									Character.toUpperCase( 
											names.get(i).charAt(names.get(i).length()-1) ) 
									) {
								whichTag.set(i,nextTag);
							}

						}
					} 
					else 
						check_not_last.put(tagType, true);
					lastByte_to_whichTag.put(tagType, nextTag);
					tags.add(mdp);
					nextTag++;
				}
				else {
					
					if(minLen > seqLen)
						minLen = seqLen;
					if(maxLen < seqLen)
						maxLen = seqLen;
					names.add(seqName);
					char tagType =  Character.toUpperCase(seqName.charAt(seqName.length()-1) );
					int whichTagValue = -1;
					if(lastByte_to_whichTag.containsKey(tagType))
						whichTagValue = lastByte_to_whichTag.get(tagType);
					whichTag.add(whichTagValue);
					check_not_last.put(tagType,false);
					forward.add(mdp);
				}
			}
			

			p = 0;
			for(Entry<Character,Boolean> checkTag : check_not_last.entrySet() )
			{
				if(checkTag.getValue()) {
					if(p != 0) 
					{
						//fprintf(stderr,"WARNING: Same applies to >tag%c\n",nextTag);
					}
					else {
						p = 1;
						//fprintf(stderr,"\nWARNING: You have multiple >tag%c sequences\n         and the last one does not precede a >...%c primer.\n         This probably means you've made a mistake.\n         Apart from the first >tag%c, all >tag%c tags will apply to\n         >...%c primers AFTER the >tag%c (not before it).\n",nextTag,nextTag,nextTag,nextTag,nextTag,nextTag);
					}
				}
			}


			for(int i = 0 ;  i < forward.size();i++) {
				this.backward.add(forward.get(i).getReverse());
			}


			this.np = this.forward.size();

		}  catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


	}


	
	
	
	public int getNumberOfPrimers() {
		return np;
	}
	public int getNumberOfTags() {
		return tags.size();
	}
	public int getMaxPrimerLength() {
		return maxLen;
	}
	public int getMinPrimerLength() {
		return minLen;
	}
	public int getMaxTagLength() {
		return maxTag;
	}
	
	
	public static void main(String[] argv) {
		AllPrimers testAllPrimers = new AllPrimers();

		try {
			testAllPrimers.loadFasta("/data/softwares/pooler/example/example.txt");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}









	}




	public void addTags()  {
		for(int i=0; i<this.np; i++) {
			int whichTagIndex =  this.whichTag.get(i);
			if(whichTagIndex >=0) {
				IPrimer tag=  this.tags.get(whichTagIndex);
				MaybeDegeneratePrimerAddTag(i,whichTagIndex,true,true);

			}
		}
	}

	public void removeTags()  {
		for(int i=0; i<this.np; i++) {
			int whichTagIndex =  this.whichTag.get(i);
			if(whichTagIndex >=0) {
				IPrimer tag=  this.tags.get(whichTagIndex);
				MaybeDegeneratePrimerRmTag(i,whichTagIndex,true,true);

			}
		}
	}


	private void MaybeDegeneratePrimerAddTag(int primerIndex, int whichTagIndex, boolean addToForward, boolean addToBackward)  {

		IPrimer tag = this.tags.get(whichTagIndex);
		if(addToForward)
		{
			IPrimer primer = this.forward.get(primerIndex);
			// if needed upgradeToDegenerate64
			if(tag.isDegeneratePrimer() && ! primer.isDegeneratePrimer()) {
				primer = currentFactory.upgradeToDegeneratePrimer(primer);
				this.forward.set(primerIndex,primer);
			}

			if(primer.isDegeneratePrimer() && !tag.isDegeneratePrimer())
				primer.addTag(currentFactory.upgradeToDegeneratePrimer(tag));
			else	
				primer.addTag(tag);
		}
		if(addToBackward) {
			IPrimer primer = this.backward.get(primerIndex);
			// if needed upgradeToDegenerate64
			if(tag.isDegeneratePrimer() && ! primer.isDegeneratePrimer()) {
				primer = currentFactory.upgradeToDegeneratePrimer(primer);
				this.backward.set(primerIndex,primer);
			}

			if(primer.isDegeneratePrimer() && !tag.isDegeneratePrimer())
				primer.addTagBackward(currentFactory.upgradeToDegeneratePrimer(tag));
			else	
				primer.addTagBackward(tag);
		}

	}
	private void MaybeDegeneratePrimerRmTag(int primerIndex, int whichTagIndex, boolean rmFromForward, boolean rmFromBackward)  {

		IPrimer tag = this.tags.get(whichTagIndex);
		if(rmFromForward)
		{
			IPrimer primer = this.forward.get(primerIndex);
			// if needed upgradeToDegenerate64
			if(tag.isDegeneratePrimer() && ! primer.isDegeneratePrimer()) {
				primer = currentFactory.upgradeToDegeneratePrimer(primer);
				this.forward.set(primerIndex,primer);
			}

			if(primer.isDegeneratePrimer() && !tag.isDegeneratePrimer())
				primer.removeTag(currentFactory.upgradeToDegeneratePrimer(tag));
			else	
				primer.removeTag(tag);
		}
		if(rmFromBackward) {
			IPrimer primer = this.backward.get(primerIndex);
			// if needed upgradeToDegenerate64
			if(tag.isDegeneratePrimer() && ! primer.isDegeneratePrimer()) {
				primer = currentFactory.upgradeToDegeneratePrimer(primer);
				this.backward.set(primerIndex,primer);
			}

			if(primer.isDegeneratePrimer() && !tag.isDegeneratePrimer())
				primer.removeTagBackward(currentFactory.upgradeToDegeneratePrimer(tag));
			else	
				primer.removeTagBackward(tag);
		}

	}

	public  int[]  counts() {
		int countSize = currentFactory.getCurrentSize();

		int[] counts= new int[countSize];
		int maxS = 0;
		for(int i=0; i<np; i++) {
			IPrimer forwardPrimer = this.forward.get(i);
			for(int j=i; j<np; j++) {
				IPrimer backwordPrimer = this.backward.get(j);
				int score = currentFactory.getScore(forwardPrimer,backwordPrimer);
				counts[score]++; 
				if(score>maxS) 
					maxS=score;
			}
		}
		return counts;
	}




	public void dGandScoreCounts(AllPrimers ap, float[] table) {
		/* show how much the score can vary around a deltaG range
	     (in case anyone thinks it's more accurate than it is) */
		int[] counts=new int[0x4000];
		int[] maxScore= new int[0x4000];
		int[] minScore= new int[0x4000];
		for (int i =  0 ; i < minScore.length ;i++)
			minScore[i] = -1;
		//	  if(memFail(counts,minScore,maxScore,_memFail)) return;
		//	  memset(minScore,0xFF,0x4000*sizeof(int));
		Instant next = Triangle.t_ProgressStart("Precalculating dG... ");
		Instant start = next;		
		//	  #if defined(_OPENMP)
		//	  #pragma omp parallel
		//	  #endif
		{
			TwoRanges tr= Triangle.t_iBounds(np);
			int r,i,j,done=0;
			for(r=0; r<2; r++)
			{
				for( i= tr.r[r].start; i<tr.r[r].end; i++) {
					IPrimer forwardPrimer = this.forward.get(i);
					for(j=i; j<np; j++) {
						IPrimer backwordPrimer = this.backward.get(j);

						float dG = currentFactory.getDeltaG(forwardPrimer,backwordPrimer,table);
						int score = currentFactory.getScore(forwardPrimer,backwordPrimer);
						int bucket = HelperMethods.dGbucket(dG,0x4000-1);
						// #if defined(_OPENMP)
						// #pragma omp critical
						//#endif
						{
							counts[bucket]++;
							if(minScore[bucket]<0 || score<minScore[bucket])
								minScore[bucket] = score;
							if(score>maxScore[bucket]) 
								maxScore[bucket] = score;
						}
					}
					next = Triangle.t_Progress("Precalculating dG... ",tr,np,done,next);
					done+=np-i;
				}
			}
		}
		long numSecs = next.getEpochSecond() - start.getEpochSecond();
		System.err.format("\rPrecalculating dG... done (%d:%02d:%02d)",(int)(numSecs/3600L),(int)(numSecs%3600L)/60,(int)(numSecs%60L));


		for(int i=0; i<0x4000; i++) 
			if(counts[i] != 0) {
				System.out.format("%c%d.%d\t%d\t (score ",i != 0 ?'-':' ', i/2,(i%2) != 0 ? 5:0 ,counts[i]);
				if(minScore[i]==maxScore[i]) 
					System.out.format("%d)\n",minScore[i]);
				else System.out.format("%d-%d)\n",minScore[i],maxScore[i]);
			}

	}

	public class PS_cache {
		public int[] scores; 
		int[] primerMove_depends_on;
		int []fix_to_pool; 
		public int fix_min_pools;
		public void saturate_scores_of_overlapping_primers(boolean[] overlappingAmplicons, int[] primerNoToAmpliconNo,
				int nAmplicons, int np) {
			int i,j,p=0; 
			if(scores == null || nAmplicons == 0) return;
			//			  assert(overlappingAmplicons);
			for(i=0; i<np; i++) 
				for(j=i; j<np; j++) {
					//			      assert(*p<InvalidCombination);
					if(
							i!=j && 
							primerNoToAmpliconNo[i]!=-1 &&
							primerNoToAmpliconNo[j]!=-1 && 
							overlappingAmplicons[primerNoToAmpliconNo[i]*nAmplicons+primerNoToAmpliconNo[j]]
							)
						scores[p] = InvalidCombination;
					p++;
				}			
		}
	};

	PS_cache selfPreCalcCache = null;
	PS_cache PS_precalc(final float[] table)
	{
		if (selfPreCalcCache != null )
			return selfPreCalcCache;
		selfPreCalcCache =  PS_precalc(table,null,null,0);
		return selfPreCalcCache;
	}


	public PS_cache PS_precalc(final float[] table,
			final boolean[] overlappingAmplicons,
			final int[] primerNoToAmpliconNo,int nAmplicons) {

		if (overlappingAmplicons == null &&  nAmplicons == 0 && primerNoToAmpliconNo == null )
		{
			if (selfPreCalcCache != null )
				return selfPreCalcCache; 
		}

		PS_cache r = new PS_cache();
		this.addTags();
		r.scores = table != null ? this.dGtriangle(table) : this.triangle();
		this.removeTags();
		r.primerMove_depends_on = merge_scores_of_stuckTogether_primers(r.scores);
		r.fix_to_pool = pre_fix_primers_to_pools();
		//		  if(memFail(r.scores,r.primerMove_depends_on,r.fix_to_pool,_memFail)) r.scores = NULL;
		//		  else 
		{
			r.saturate_scores_of_overlapping_primers(overlappingAmplicons,primerNoToAmpliconNo,nAmplicons,np);
			r.fix_min_pools = 2;

			for(int i=0; i<np; i++) { 
				if(r.fix_to_pool[i]>=r.fix_min_pools) 
					r.fix_min_pools = r.fix_to_pool[i]+1;
			}
		}
		if (overlappingAmplicons == null &&  nAmplicons == 0 && primerNoToAmpliconNo == null )
		{
			selfPreCalcCache = r;
		}
		return r;
	}




	final int[] pre_fix_primers_to_pools() {
		int[] fix_to_pool=new int[this.np];
		for(int i=0; i< this.np; i++)
			fix_to_pool[i] = should_stick_to_pool(i);
		return fix_to_pool;
	}
	final int should_stick_to_pool(int i) {
		/* if the user wants some primers to be fixed to
		     specific pools (and we move the rest around) */
		String n = this.names.get(i);
		//		  if(*n == '@' && *(++n)>='0' && *n<='9') {
		//		    char *end; int pool=(int)strtol(n,&end,10);
		//		    if(*end==':') {
		//		      /* we have a valid @<pool number>: */
		//		      return pool-1; /* (internally start at 0) */
		//		    }
		//		  }
		return -1;
	}



	private int[] merge_scores_of_stuckTogether_primers(int[] scores) {
		int p=0; 
		if(scores == null ) 
			return null;
		int[] primerMove_depends_on =    new int[np];
		for(int i = 0 ; i < primerMove_depends_on.length ; i ++ )
			primerMove_depends_on[i]=-1;

		boolean[] pairedOK = new boolean[np];

		//		  memset(pairedOK,0,ap.np);
		boolean doneMerge = false;
		for(int i=0; i< np; i++) 
			for(int j=i; j<np; j++) {
				if( i!=j &&
						primerMove_depends_on[i] ==-1 &&
						primerMove_depends_on[j] ==-1 && 
						should_stick_together(i,j)) {
					/* For simplicity of pooling, we'll set it so:
		           - Interactions with i get maxed with those w.j
		           - Interactions with j itself "don't count"
		           - j is not allowed to be moved by itself
		           - j is always moved when i moves */
					scores[p] = 0; /* so S(i,j) = 0 */
					int k,kp=0; /* max S(k,i) with S(k,j): */
					int Sip=0; /* =0 to suppress compiler warning */
					for(k=0; k<j; k++) {
						if(k<i) {
							scores[kp + i-k] = Integer.max(scores[kp + i-k],scores[kp+j-k]);
							//		        	  updateMax( kp + i-k ,  scores[kp+j-k]); 
							scores[kp+j-k]=0;
						} else if(k==i) {
							Sip = kp+1; /* needed for S(i,k) */
						} else { /* max S(i,k) with S(k,j) */
							scores[Sip] = Integer.max(scores[Sip],scores[kp+j-k]);
							Sip++;
							//		        	  updateMax(Sip++,scores[kp+j-k]); 
							scores[kp+j-k]=0;
						}
						kp += (np-k);
					} k++; kp++; Sip++; /* ignore k==j */
					for(;k<np;k++) {
						/* max S(i,k) [=Sip] with S(j,k) [=kp] */
						//		          updateMax(Sip++,*kp); 
						scores[Sip] = Integer.max(scores[Sip],scores[kp]);
						Sip++;
						scores[kp++]=0;
					}
					primerMove_depends_on[j] = i;
					doneMerge = pairedOK[i] = pairedOK[j] = true;
				}
				p++;
			}
		if(doneMerge) {
			/* just check for lone primers, usually a bad sign */
			for(int i=0; i< np; i++) 
				if(!pairedOK[i]) 
					System.err.format("Warning: ungrouped primer %s\n",this.names.get(i));
		} else {
			/* same message as in amplicons.c (see comment there)
		       in case overlap-check was missed */
			System.err.format("WARNING: No primers are paired!\nPlease end your forward primers with -F\nand your reverse primers with -R or -B as instructed\n");
		}
		//		  free(pairedOK); 
		return primerMove_depends_on;
	}


	final boolean should_stick_together(int i,int j) {
		/* Names same except last letter = keep in same pool */
		String n1=this.names.get(i), n2=this.names.get(j);
		return   n1.substring(0, n1.length()-1) .equals(n2.substring(0, n2.length()-1)); /* TODO: case-insensitive? */
	}

	private int[] triangle() {
		int[] scores = new int[(np*(np+1)/2)];
		int i,j,p=0;
		for(i=0; i<np; i++) { 
			IPrimer forwardPrimer = this.forward.get(i);
			for(j=i; j<np; j++) {
				IPrimer backwordPrimer = this.backward.get(j);
				scores[p++] = (i==j ? 0 /* ignore self-interaction */ :
					currentFactory.getScore(forwardPrimer,backwordPrimer));
			}
		}
		return scores;
	}




	private int[] dGtriangle(float[] table) {
		/* To save doing a float version of the pool splitter, emulate 'score' by
	     using -dG*2, as bins of size 0.5 should be enough (*10 creates too many empty ones and split_into_pools would need changing) */
		Instant next = Triangle.t_ProgressStart("Precalculating dG... ");
		Instant start = next;	
		int[] scores = new int[(np*(np+1)/2)];
		//	  if(scores )
		//	    #if defined(_OPENMP)
		//	    #pragma omp parallel
		//	    #endif
		{
			TwoRanges tr= Triangle.t_iBounds(np);
			int r,i,j, p = 0 ,done=0;
			for(r=0; r<2; r++) {
				for(i= tr.r[r].start, p =  t_offset(np,i,i);
						i<tr.r[r].end; 
						i++) {
					IPrimer forwardPrimer = this.forward.get(i);
					for(j=i; j<np; j++) {
						IPrimer backwordPrimer = this.backward.get(j);
						scores[p++] = (i==j ? 0 : HelperMethods.dGbucket(currentFactory.getDeltaG(forwardPrimer,backwordPrimer,table),0x4000-1));
					}
					next = Triangle.t_Progress("Precalculating dG... ",tr,np,done,next);
					done+=np-i;
				}
			}
		} 
		long numSecs = next.getEpochSecond() - start.getEpochSecond();
		System.err.format("\rPrecalculating dG... done (%d:%02d:%02d)",(int)(numSecs/3600L),(int)(numSecs%3600L)/60,(int)(numSecs%60L));
		//	  prnSeconds((long)(time(NULL)-start)); fputs("\n",stderr); fflush(stderr);
		return scores;
	}
	static final int t_offset(int n,int i,int j) {
		if(j>=i) return _t_offset(n,i,j);
		else return _t_offset(n,j,i);
	}
	static final int _t_offset(int n,int i,int j) {
		/* Offset of the (i,j) pair assuming j >= i */
		/* line 0 starts at 0
		     line 1 starts at n
		     line 2 starts at n + n-1  = 2n - 1
		     line 3 starts at n + n-1 + n-2 = 3n - (1+2)
		     ...
		     line N starts at N*n - Tri(N-1)
		 */
		// return n*i - (i-1)*i/2 + j-i;
		// which is
		return (n-1)*i - (i-1)*i/2 + j;
	}


	public List<List<Pair<Double, Double>>> dGprintStats(int[] pools, int[] calcScores) {
		int[] counts= new int[0x4000+1];
		List<List<Pair<Double, Double>>> poolsCounts = new ArrayList<>();
		//		  if(memFail(counts,_memFail)) return;
		int i,j,nPools=0,pool;
		for(i=0; i<np; i++) 
			if(pools[i]>nPools) 
				nPools=pools[i];
		//		  const int *precalcScores2 = precalcScores; 
		int precalcScores2 = 0 , precalcScores = 0;
		nPools++;
		for(pool=0; pool<nPools; pool++) {
			List<Pair<Double, Double>> poolCounts = new ArrayList<>();
			poolsCounts.add(poolCounts);
			
			System.out.format("Pool %d:\n",pool+1);
			//		    memset(counts,0,(0x4000+1)*sizeof(int));
			precalcScores = precalcScores2;
			for(i=0; i<np; i++) 
				for(j=i; j<np; j++) 
					if(pools[i]==pools[j] && pools[i]==pool) {
						counts[calcScores[precalcScores++]]++;
					} else precalcScores++;
			int lines=0; 
			for(i=0x4000-1; i > 0; i--)
				if(counts[i] != 0 && ++lines==20) 
					break;
			for(; i<0x4000; i++) 
				if(counts[i] != 0) {
					System.out.format("%.3g\t%d\n",((float)(-i))/2.0,counts[i]);
//					poolCounts.add(new ImmutablePair<String, Integer>(
//							String.format("%.3g",((float)(-i))/2.0),
//							counts[i]
//							));
					poolCounts.add(
							new ImmutablePair<Double, Double>(((double)(-i))/2.0,(double) counts[i])
							);
				}
			int other=counts[0x4000];
			if(other != 0) 
				System.out.format("Overlaps\t%d\n",other);
		}
		
		return poolsCounts;

	}

	public List<List<Pair<Double, Double>>> printStats(int[] pools, int[] calcScores) {
		List<List<Pair<Double, Double>>> poolsCounts = new ArrayList<>();
		/* like pCounts64 but outputs per-pool */
		int i,j,nPools=0,pool;
		for(i=0; i<np; i++) 
			if(pools[i]>nPools) 
				nPools=pools[i];
		//		  const int *precalcScores2 = precalcScores; 
		int precalcScores = 0;
		nPools++; /* 1 higher than max */
		for(pool=0; pool<nPools; pool++) {
			List<Pair<Double, Double>> poolCounts = new ArrayList<>();
			poolsCounts.add(poolCounts);
			System.out.format("Pool %d:\n",pool+1);
			precalcScores = 0;
			int[] counts = new int[64];
			int maxS = 0, other=0;
			for(i=0; i<np; i++) 
				for(j=i; j<np; j++) 
					if(pools[i]==pools[j] && pools[i]==pool) {
						int score = calcScores != null  ? calcScores[precalcScores++] : currentFactory.getScore(this.forward.get(i), this.backward.get(j));
						if(score<64) {
							counts[score]++; if(score>maxS) maxS=score;
						} else other++;
					} else if(calcScores != null) 
						precalcScores++;
			for(i=0; i<=maxS; i++) { 
				System.out.format("%d\t%d\n",i,counts[i]);
				poolCounts.add(
						new ImmutablePair<Double, Double>((double)i,(double) counts[i])
						);
			}
			if(other != 0) 
				System.out.format("Overlaps\t%d\n",other);
		}
		//		  if(f!=stdout && f!=stderr) fclose(f);		
		return poolsCounts;
	}




	public int dGprintPooledCounts(final int[] bestPools, final int[] scores) {
		// TODO Auto-generated method stub
		return dGpCounts64(np,bestPools,scores);
	}
	static int dGpCounts64(int np,final int[] pools,
			final int[] scores) {
		  /* this is a combination of counts64 and pCounts64, for the delta-G variant.  pools can be NULL, but not precalcScores */
		  int i,j; 
		  int[] counts= new int[0x4000+1];
		  int precalcScores =0 ;
//		  if(!counts) return 0;
		  for(i=0; i<np; i++) 
			  for(j=i; j<np; j++) 
				  if(pools == null || pools[i]==pools[j]) {
		        counts[scores[precalcScores++]]++; /* this version has no precalcScores==NULL fallback; could put one in if don't mind passing around the extra parameters */
		      } 
				  else precalcScores++;
		  int lines=0; 
		  for(i=0x4000-1; i>0; i--)
		                 if(counts[i] != 0 && ++lines==20) 
		                	 break;
		  int first = 1;
		  for(; i<0x4000; i++) 
			  if(counts[i] != 0) {
		      System.err.format("%s%.3g\t%d",first != 0 ? "":"\n",((float)(-i))/2.0,counts[i]);
		      if ( first != 0) 
		    	  first =  0 ; 
		    }
		  int other=counts[0x4000]; 
//		  free(counts);
		  if(other != 0) { 
			  System.err.format("%sOverlaps\t%d",first != 0 ?"":"\n",other);
		      if ( first != 0) 
		    	  first =  0 ;
		  }
//		  if(f!=stdout && f!=stderr) 
//			  fclose(f);
		  return other;
		}




	public int printPooledCounts(int[] bestPools, int[] scores) {
		// TODO Auto-generated method stub
		return pCounts64(bestPools,scores);
	}
	
	int pCounts64(final int[] pools,final int[] scores) {
		  /* like counts64 but includes combinations only if they're in the same pool + count invalid/overlap scores */
		  int i,j;
		  int[] counts = new int[64];
		  int maxS = 0, other=0;
		  int precalcScores=0;
		  for(i=0; i<this.np; i++) 
			  for(j=i; j<this.np; j++) 
				  if(pools[i] == pools[j]) {
					  int score = scores != null ? scores[precalcScores++] : currentFactory.getScore( forward.get(i),backward.get(j));
		        if(score<64) {
		          counts[score]++; 
		          if(score>maxS) 
		        	  maxS=score;
		        } else 
		        	other++;
		      } else 
		    	  if(scores != null) 
		    		  precalcScores++;
		  boolean first = true;
		  for(i=0; i<=maxS; i++) {
		    System.err.format("%s%d\t%d",first  ?"":"\n",i,counts[i]);
		    if (first) 
		    	first=false;	
		  }
		  if(other != 0) 
		  {
			  System.err.format("%sOverlaps\t%d",first?"":"\n",other);
			    if (first) 
			    	first=false;
		  }
		  return other;
		}


	public void dGprintBonds(PrintStream out, double threshold, int[] pools, float[] table) {
		dGprintBonds(out,threshold,pools,table,-1);
		
	}

	public void dGprintBonds(PrintStream out, double threshold, int[] pools, float[] table, int onlyPoolI) {

		class DG_ScoreRecord {
			public DG_ScoreRecord(float dG, int i, int j) {
				this.i = i; this.j = j ; this.dG = dG;
			}
			int i , j;
			float dG;
		};

		Instant next = Triangle.t_ProgressStart("Sorting...");
		// DG_ScoreRecord[] sr = new DG_ScoreRecord[ Triangle.t_Nitems(np)];
		ArrayList<DG_ScoreRecord> sr = new ArrayList<DG_ScoreRecord>();
		// int srIndex = 0; 
		Instant start = Instant.now();
		// #if defined(_OPENMP)
		// #pragma omp parallel
		// #endif
		{
			TwoRanges tr= Triangle.t_iBounds(np);
			int r,i,j,done=0;
			for(r=0; r<2; r++)
				for(i=tr.r[r].start; i<tr.r[r].end; i++, next = Triangle.t_Progress("Sorting... ",tr,np,done,next),done+=np-i)
				{
					IPrimer forwardPrimer = forward.get(i);
					for(j=i; j<np; j++) {
						if( pools == null || pools[i]==pools[j]) {
							if(onlyPoolI != -1 && onlyPoolI != pools[i])
								// only then 
								continue;
							float dG =   currentFactory.getDeltaG(forwardPrimer, backward.get(j),table);
							if (dG <= threshold) {
								// #if defined(_OPENMP)
								// #pragma omp critical
								//#endif
								// if(sr) 
								{
									sr.add( new DG_ScoreRecord(dG, i, j)) ;
									// srIndex++;
								} 
								//	else 
								//	dGprint64MaybeD(forward[i],backward[j],names[i],names[j],dG,f,table);
							}
						}
					}
				}
		}
		//		  if(sr) 
		{
			sr.sort( new Comparator<DG_ScoreRecord>() {

				@Override
				public int compare(DG_ScoreRecord o1, DG_ScoreRecord o2) {
					return Float.compare(o1.dG,o2.dG) ;
				}
			}); 
			//		    qsort(sr,sr2-sr,sizeof(DG_ScoreRecord),dGhighestScore1st);
			System.err.format("\rSorting... done");
			HelperMethods.prnSeconds( Instant.now().getEpochSecond() - start .getEpochSecond()); 
			System.err.format("\n");
			//		    if(f!=stdout) 
			//		    { 
			//		    	fputs("Outputting... ",stderr); 
			//		    	start = time(NULL); 
			//		    	next = start + 2; 
			//		    }
			//		    fflush(stderr);
			int srIndex = 0;
			for(srIndex=0 ; srIndex<sr.size(); srIndex++) {
				//		      if(f!=stdout && time(NULL) > next) {
				//		        fprintf(stderr,"\rOutputting... (%d%%) ",100*(int)(s-sr)/(int)(sr2-sr)); fflush(stderr);
				//		        next = time(NULL) + 2;
				//		      }
				DG_ScoreRecord s =  sr.get(srIndex);
				dGprint64MaybeD(s.i,s.j,s.dG,out,table);
			}
			//		    free(sr);
			if(out != System.out) 
			{ 
				//		    	fputs("\rOutputting... done",stderr); 
				//		    	prnSeconds((long)(time(NULL)-start)); 
				//		    	fputs("\n",stderr); fflush(stderr); }
			}
		}
	}
	
	public void printBonds(PrintStream out, double threshold, int[] pools) {
		printBonds( out,  threshold, pools,-1);
	}
	public void printBonds(PrintStream out, double threshold, int[] pools , int onlyPoolI) {
		class ScoreRecord {
			public ScoreRecord(int score, int i, int j) {
				this.i = i; this.j = j ; this.score = score;
			}
			int i , j;
			int score;
		}
		ArrayList<ScoreRecord> sr = new ArrayList<ScoreRecord>();
		for(int i=0; i<np; i++) {
			IPrimer forwardPrimer = forward.get(i);
			for(int j=i; j<np; j++) {
				if( pools == null || pools[i]==pools[j])  {
					if(onlyPoolI != -1 && onlyPoolI != pools[i])
						// only then 
						continue;
					int score = currentFactory.getScore(forwardPrimer, backward.get(j));
					if (score >= threshold) {
						sr.add( new ScoreRecord(score, i, j)) ;
					}
				}
			}
		}
		sr.sort( new Comparator<ScoreRecord>() {

			@Override
			public int compare(ScoreRecord a, ScoreRecord b) {
				int res = Integer.compare(b.score,a.score) ;
				if(res != 0 )
					return res;
				res = Integer.compare(a.i,b.i) ;
				if(res != 0 )
					return res;
				return Integer.compare(a.j,b.j) ;
			}
		});
		for(ScoreRecord s :sr ) {
//			forward.get(s.i).print(backward.get( s.j), s.score, out);
			print64MaybeD(s.i,s.j,s.score, out);
		}
	}
	

	void dGprint64MaybeD(int fI ,int bI,float minDG,PrintStream out,final float[]table) {
		// if(!name1 || !*name1) name1="(no name)";
		// if(!name2 || !*name2) name2="(no name)";
		out.printf("%s versus %s\n",names.get(fI),names.get(bI));
		IPrimer forwardPrimer = this.forward.get(fI);
		IPrimer backwardPrimer = this.backward.get(bI);
		count64MaybeD(forwardPrimer,backwardPrimer,out);
		currentFactory.dGprint(forwardPrimer,backwardPrimer,minDG,out,table);
		out.println();
	}


	void print64MaybeD(int fI ,int bI,int maxScore,PrintStream out) {
		// if(!name1 || !*name1) name1="(no name)";
		// if(!name2 || !*name2) name2="(no name)";
		out.printf("%s versus %s\n",names.get(fI),names.get(bI));
		IPrimer forwardPrimer = this.forward.get(fI);
		IPrimer backwardPrimer = this.backward.get(bI);
		count64MaybeD(forwardPrimer,backwardPrimer,out);
		currentFactory.print(forwardPrimer,backwardPrimer,maxScore,out);
		out.println();
	}
	
	private void count64MaybeD(IPrimer forwardPrimer, IPrimer backwardPrimer, PrintStream out) {
		  /* Put clarification for any beginner users who haven't
	     been informed that we automatically try all positions
	     and print only the worst case */

		  CountResult result = currentFactory.getCount(forwardPrimer,backwardPrimer);
		  out.printf("Positions tried: %d\nBonding positions: %d%s\n",result.tried,result.count,(result.count>1)?"  (worst one shown here)":"");

	}








    // ampliconNoToFwd[ampNo],ampliconNoToRev[ampNo]
	public int getNumPossibilities(int i) {
		return forward.get(i).NumPossibilities_32bases();
	}




	public boolean Make2bit(int n, int useBackward, int doComplement, PrimerToFind newPrimerPoss,
			int possNo,int nPoss) {
		

	    IPrimer p = useBackward == 1 ? this.backward.get(n) : this.forward.get(n);
	    if(doComplement == 1) p = p.getComplement();
	    return currentFactory.upgradeToDegeneratePrimer(p).make2bit(newPrimerPoss,possNo,nPoss);
		

	}


	public String getPrimerName(int i) {
		return names.get(i);
	}






}
