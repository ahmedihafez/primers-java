/*
  This file is part of Java porting of Primer Pooler
  Original Primer Pooler (c) Silas S. Brown.  For Wen.
  Please refer to https://github.com/ssb22/PrimerPooler
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

import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;

import org.apache.commons.lang3.ArrayUtils;
import org.pooler.AllPrimers.PS_cache;

public class PoolSplit {
	static void randomise_pools(int np,
			final int[] primerMove_depends_on,
			final int[] fix_to_pool,
			final int[] scores,
			int nPools,
			int[] pools,
			long[] bContrib,
			int[] poolCounts,
			int maxCount) throws Exception {
		/* initialise to random distribution of pools, but note
		     primerMove_depends_on and maxCount when doing this.
		     Also initialise bContrib.  */
		int i; 
		for(i=0;i<poolCounts.length;i++)
			poolCounts[i]=0;
		//		  memset(poolCounts,0,nPools*sizeof(int));
		/* First set all fixed-pool primers in place,
		     before randomising the others around them */
		for(i=0; i<np; i++)
			if(primerMove_depends_on[i] == -1) {
				int pool = fix_to_pool[i];
				if(pool != -1) {
					if(maxCount != 0 && poolCounts[pool]==maxCount && !(maxCount==1 && nPools==np)) {
						/* (last part of that condition detects call by suggest_num_pools,
		             where it's OK if fixed-pool primers make us exceed 1 per pool) */
						System.err.format("randomise_pools ERROR: maxCount too small for fixed primer in pool %d\n",fix_to_pool[i]);
						throw new Exception("randomise_pools ERROR: maxCount too small for fixed primer in pool "+ fix_to_pool[i]);
					} 
					pools[i]=pool; 
					poolCounts[pool]++;
				}
			}
		for(i=0; i<np; i++)
			if(primerMove_depends_on[i] == -1 && fix_to_pool[i] == -1) {
				int pool = ThreadLocalRandom.current().nextInt(0, nPools*2) % nPools;
				int origPool = pool;
				while(maxCount != 0 && poolCounts[pool]>=maxCount) {
					pool++; /* not very random but it'll do for now */
					if(pool==nPools) pool=0;
					if(pool==origPool) {
						System.err.format( "randomise_pools ERROR: maxCount too small, can't fit\n");
						throw new Exception("randomise_pools ERROR: maxCount too small, can't fit");
					}
				} pools[i]=pool; poolCounts[pool]++;
			}
		for(i=0; i<np; i++)
			if(primerMove_depends_on[i]>-1)
				/* DON'T update poolCounts here (it's in pairs so moveTooLopsided doesn't have to account for this one) */
				pools[i]=pools[primerMove_depends_on[i]];
		//		  memset(bContrib,0,np*nPools*sizeof(ULL));
		for(i=0;i<bContrib.length;i++)
			bContrib[i]=0;
		for(i=0; i<np; i++) 
			badnessContrib(i,scores,np,nPools,pools,bContrib);
	}
	/* and here is code to set up & maintain that bContrib: */
	static void badnessContrib(int primer,
			final int[] scores,
			int np,
			int nPools,
			final int[] pools,
			long[] bContrib) {
		/* Assuming proposedPools[0:nPools] == 0 on entry, set proposedPools[0:nPools] to answer the Q: What contribution to the overall "badness" would primer make, assuming it were moved to (or left as-is in) proposedPools[n] and no other changes were made? */
		assert(primer>=0 && primer<np);
		//	  ULL *proposedPools = bContrib + primer*nPools;
		int proposedPools = primer*nPools;
		int scoresIndex = 0;

		int i;
		for(i=0; i<primer; i++) { /* ( <primer , primer ) */
			int pool = pools[i]; assert(pool>=0 && pool<nPools);
			bContrib[proposedPools+pool] = updateBadness( bContrib[proposedPools+pool],scores[scoresIndex+primer-i]); /* if we put 'primer' in the same pool as 'i' is, we'll get the badness of the interaction between i and primer */
			scoresIndex += (np-i);
		}
		++scoresIndex; /* i==primer, ignore interaction w. itself */
		for(++i; i<np; i++,scoresIndex++) { /* ( primer, >primer ) */
			int pool = pools[i]; assert(pool>=0 && pool<nPools);
			bContrib[proposedPools+pool] = updateBadness(bContrib[proposedPools+pool],scores[scoresIndex]);
		}
	}
	static final int maxScoreOfBadness(long badness) {
		return (int)(badness >>> 48);
	}
	final static int InvalidCombination=0x4000 ;
	static final long updateBadness(long badness,int score) {
		int max = maxScoreOfBadness(badness);
		assert(max <= InvalidCombination); 
		assert(score <= InvalidCombination);
		if (score > max) {
			if (score > max+2) {
				badness = ( ((long)score) << 48) | (((long)1) << 32);
				return badness;
			} else 
				while(score > max) 
				{ /* (TODO: could just write out the 2 cases of this loop, and change the below shifts into AND and add, IF profiling shows this needs it) */
					badness = ((badness >>> 16) & 0xFFFFFFFFL) | ((((badness)>>>48)+1)<<48); 
					max++;
				}
		} /* now score <= max */
		int lswScore = max - 2;
		if(score < lswScore) 
			return badness; /* this score is too insignificant to count at the current maximum level */
		int sL = (score-lswScore)*16;
		if(((badness>>sL)&0xFFFF)==0xFFFF) return badness; /* saturated */
		badness += (((long)1) << sL);
		assert(maxScoreOfBadness(badness) == max);
		return badness;
	}



	static final int primerAndPool_to_contribOffset(int primer,int pool,int nPools) {
		return primer * nPools + pool; 
	}
	static final int primerAndDest_to_moveNo(int primer,int newPool,int nPools,final int[] pools) {
		// works only if newPool != current pool
		return ((newPool+nPools-pools[primer]-1) % nPools) + primer*(nPools-1);
	}
	static public int suggest_num_pools(AllPrimers ap, PS_cache cache, float[] table) throws Exception {
		/* Apply a simple threshold-based allocation just for
	     suggesting a number of pools */
		int threshold = table != null ? 14 : 7; /* dG -7 or score 7.  TODO: customise?  but this function is for when the user is not sure, so perhaps we'd best hard-code the threshold */
		int nPools = ap.np; /* worst case is none of the primers are paired (unpaired primers could hang randomise_pools before v1.42 because this line said ap.np/2) */
		if(cache == null || cache.scores == null) 
			return 0;
		int[] scores = cache.scores; 
		int[] primerMove_depends_on = cache.primerMove_depends_on;
		int[] fix_to_pool = cache.fix_to_pool;

		long[] bContrib = new long[ap.np*nPools];
		int[] poolCounts=    new int[nPools];
		int[] pools = new int [ap.np];
		//	  if(memFail(bContrib,poolCounts,pools,_memFail))
		//	    return 0;
		PoolSplit.randomise_pools(ap.np,primerMove_depends_on,fix_to_pool,scores,nPools,pools,bContrib,poolCounts,1); /* puts 0 or 1 set in each pool (after the fixed ones) */
		int suggest_nPools = 1;
		int primer; for (primer=0; primer<ap.np; primer++) 
			if (primerMove_depends_on[primer]==-1) {
				if (fix_to_pool[primer]==-1) {
					int destPool; for (destPool=0; destPool < suggest_nPools; destPool++) 
						if(PoolSplit.maxScoreOfBadness(bContrib[primerAndPool_to_contribOffset(primer,destPool,nPools)]) <= threshold) 
							break; /* find first pool it will 'fit' in */
					if (destPool == suggest_nPools) 
						suggest_nPools++;
					if (pools[primer] != destPool) 
						make_a_move(primerAndDest_to_moveNo(primer,destPool,nPools,pools),ap.np,scores,primerMove_depends_on,nPools,pools,bContrib,poolCounts,ap.np);
				} else if (fix_to_pool[primer] >= suggest_nPools) {
					/* must have at least as many for the fixed-pool primers
	         (and fix_to_pool starts numbering at 0, so +1 of course) */
					suggest_nPools = fix_to_pool[primer] + 1;
				}
			}
		//	  free(bContrib); free(pools); free(poolCounts);
		return suggest_nPools;
	}
	static  void make_a_move(int m,
			int np,
			final int[] scores,
			final int[] primerMove_depends_on,
			int nPools,
			int[] pools,
			long[] bContrib,
			int[] poolCounts,
			int maxCount) {
		int primer = primerOfMove(m,nPools),
				oldPool = oldPoolOfMove(m,nPools,pools),
				newPool = poolOfMove(m,nPools,pools);
		assert(primer >= 0 && primer < np && oldPool>=0 && newPool>=0 && oldPool<nPools && newPool<nPools && oldPool != newPool);
		pools[primer] = newPool; /* 'm' changes meaning now */
		assert(poolCounts[oldPool] != 0);
		poolCounts[oldPool]--; poolCounts[newPool]++;
		assert(maxCount == 0 || poolCounts[newPool]<=maxCount);
		int i; for(i=0; i<np; i++) {
			if(primerMove_depends_on[i]==primer)
				pools[i] = newPool; /* see merge_scores_of_stuckTogether_primers (and DON'T need to update poolCounts here) */
			else badnessContribUpdate(i,scores,np,primer,oldPool,nPools,pools,bContrib,poolCounts);
		}
	}
	static final int primerOfMove(int m,int nPools) {
		return m/(nPools-1); }
	static final int oldPoolOfMove(int m,int nPools,final int[] pools) {
		return pools[primerOfMove(m,nPools)]; }
	static final int poolOfMove(int m,int nPools,final int[] pools) {
		return (oldPoolOfMove(m,nPools,pools)+1+(m % (nPools-1))) % nPools; 
	}

	static void badnessContribUpdate(int primer,
			final int[] scores,
			int np,
			int otherPrimer,
			int otherOldPool,
			int nPools,
			final int[] pools,
			long[] bContrib,
			final int[] poolCounts) {
		/* as above but just incrementally update primer's proposedPools in light of the fact that otherPrimer has just moved from otherOldPool to its current pool */
		assert(primer>=0 && primer<np && otherPrimer>=0 && otherPrimer<np && otherOldPool>=0 && otherOldPool<nPools);
		int proposedPools = primer*nPools ; // bContrib + primer*nPools;
		int  s, otherNewPool = pools[otherPrimer];
		int scoresIndex = 0 ;
		assert(otherNewPool>=0 && otherNewPool<nPools);
		if(otherPrimer == primer) return;
		s = scores[t_offset(np,primer,otherPrimer)]; /* the score-contribution of interaction between primer and otherPrimer */
		if(poolCounts[otherOldPool] == 0)
			bContrib[proposedPools+otherOldPool] = 0; /* like the loop below but also clears any saturation (since we know there won't be saturation if the pool was left empty) */
		else 
		{
			if(subtractBadness(bContrib , proposedPools+otherOldPool,s)) {
				/* oops, need to recalc the max of otherOldPool */
				bContrib[proposedPools+otherOldPool] = 0;
				int i;
				for(i=0; i<primer; i++) { /* ( <primer , primer ) */
					if(pools[i] == otherOldPool)
						bContrib[proposedPools+otherOldPool] = updateBadness(bContrib[proposedPools+otherOldPool],scores[scoresIndex+primer-i]);
					scoresIndex += (np-i);
				}
				++scoresIndex;
				for(++i; i<np; i++,scoresIndex++) {
					if(pools[i]==otherOldPool)
						bContrib[proposedPools+otherOldPool] = updateBadness(bContrib[proposedPools+otherOldPool],scores[scoresIndex]);
				}
			}
		}
		bContrib[proposedPools+otherNewPool] = updateBadness(bContrib[proposedPools+otherNewPool],s);
	}


	static final boolean subtractBadness(final long[] bContrib, int index ,int score) {
		/* for incremental updates.  Assume score has previously
		     been included in updateBadness, so we don't have to
		     worry about crossing 0 here.  */
		long  badness = bContrib[index];
		int max = maxScoreOfBadness(badness);
		assert(score <= max);
		int lswScore = max - 2;
		if(score < lswScore) return false; /* this score is too insignificant to affect the counters */
		int sL = (score-lswScore)*16;
		if(((badness>>sL)&0xFFFF)==0xFFFF) return false; /* if it was saturated, we'll have to leave it "stuck" there I'm afraid (unless return 1 to recalculate, but in many cases it would just saturate again) */
		badness -= (((long)1) << sL);
		bContrib[index] = badness;
		return ((badness>>32)&0xFFFF) == 0; /* recalc max if count(max)==0 */
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





	public static int[]  split_into_pools(AllPrimers ap,
			int nPools,
			int timeLimit,
			PS_cache cache,
			boolean seedless,
			final float[] table,
			int maxCount) throws Exception {
		if(cache == null ||  cache.scores == null) return null;  

		{
			if(nPools<cache.fix_min_pools) { 
				System.err.format("ERROR: @%d:primers need at least %d pools, but only got %d\n",cache.fix_min_pools,cache.fix_min_pools,nPools); 
				return null; 
			}
		}



		Splitter poolSplitter = new Splitter(ap,nPools,timeLimit,cache,seedless,table,maxCount);


		poolSplitter.runSplit();




		//		fflush(stderr);
		//		free(shared_moves); 
		return  poolSplitter.bestPools;
	}






	static int[] numInEachPool(final int[] pools,
			int np,
			int numPools,
			final int[] primerMove_depends_on) {
		/* for after everything has finished and the per-thread poolCounts has been freed */
		int[] counts= new int[numPools];
		//		  if(!counts) return NULL;
		int i; 
		for(i=0; i<np; i++)
			if(primerMove_depends_on[i]==-1)
				counts[pools[i]]++;
		Arrays.sort(counts );
		ArrayUtils.reverse(counts);
		return counts;
	}

	static void printNumInEachPool(final int[] poolCounts,int numPools) {
		System.err.format("\tPool sizes: "); 
		int i;
		for(i=0; i<numPools; i++) {
			if(i != 0 ) 
				System.err.format("|");
			System.err.format("%d",poolCounts[i] << 1); 
			/* TODO: this "<< 1" assumes countOf(primerMove_depends_on==-1) == np/2, 
			 * but that is almost certainly going to be the case, 
			 * unless somebody is doing something very strange, 
			 * and if the worst comes to the worst it's only an informational pool-size display going a bit wrong */
		} 
		System.err.format("\n");
	}

	static final long valueOfReduction(long from,long to) {
		/* for qsort etc.  We COULD return negative values as
		     LL, but that would be more computation and we might
		     as well just return 0, since we won't perform any
		     value-based moves that don't have positive reductions
		 */
		if(to > from) 
			return 0; 
		if( to == 0) 
			return from;
		if(maxScoreOfBadness(from) == maxScoreOfBadness(to)) {
			/* 'from-to' will be a 48-bit value, and if get here
		       then the high 16 bits of 'to' will be <= 'from',
		       but if the mid 16 bits are >= then we need to set
		       all of the low 32 bits to 0, and if the bottom
		       16 bits are >= then we need to set low 16 to 0.  */
			long hi = (from & ((long)0xFFFF<<32)) - (to & ((long)0xFFFF<<32)), /* do NOT factor out the & part! */
					mid1 = from & ((long)0xFFFF<<16),
					mid2 = to & ((long)0xFFFF<<16);
			if (mid2 > mid1) return hi;
			long lo1 = from & 0xFFFF, lo2 = to & 0xFFFF;
			if (lo2 > lo1) return hi | (mid1-mid2);
			return from-to;
		}
		/* maxScoreOfBadness(from) > maxScoreOfBadness(to) :
		     at the very least we need to 0-out the bottom 48 bits
		     (TODO: we might be able to add some back in if
		     maxScoreOfBadness(from) <= maxScoreOfBadness(to)+2,
		     but probably OK just to do this for now...)
		 */ return (from & ((long)0xFFFF<<48)) - (to & ((long)0xFFFF<<48));
	}
	static final int move_to_contribOffset(int m,int nPools,final int[] pools) {
		return primerAndPool_to_contribOffset(primerOfMove(m,nPools),poolOfMove(m,nPools,pools),nPools); 
	}
	
	static final long globalBadness(final int[] score,int np,final int[] pools) {
		  /* Measure across all primers in all pools.  We might
		     perhaps be able to optimise this by making use of
		     what we've already calculated in bContrib, but this
		     function is called only when we get to local maxima.
		  */
		  int i,j; 
		  long m=0;
		  int scoreIndex =0 ;
		  for(i=0; i<np; i++) 
			  for(j=i; j<np; j++) 
				  if(pools[i]==pools[j]) 
					  m = updateBadness(m,score[scoreIndex++]); 
				  else 
					  scoreIndex++;
		  return m;
		}
	
}
