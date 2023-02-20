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
package org.ucam.ssb22.pooler;

import java.io.PrintStream;
import java.time.Instant;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;

import org.ucam.ssb22.pooler.AllPrimers.PS_cache;

public class Splitter {

	boolean USE_QSORT = false;

	AllPrimers ap;
	int nPools; 
	int timeLimitInMins; 
	boolean seedless; 
	float[] table; 
	int maxCount ;

	int[] scores = null; 
	int[] primerMove_depends_on = null;
	int[] fix_to_pool = null;





	//  ##################
	List<Integer> shared_moves = null;
	int numMoves=0;
	public int[] bestPools = null;
	// ##################


	public Splitter(AllPrimers ap, int nPools, int timeLimitInMins,PS_cache cache, boolean seedless, float[] table, int maxCount) {
		this.ap = ap;
		this.nPools = nPools; 
		this.timeLimitInMins = timeLimitInMins; 
		this.seedless = seedless; 
		this.table = table; 
		this.maxCount = maxCount ;

		this.scores = cache.scores; 
		this.primerMove_depends_on = cache.primerMove_depends_on;
		this.fix_to_pool = cache.fix_to_pool;

	}

	public void runSplit(PrintStream redirErr) throws Exception {
		PrintStream orgErrStream = System.err;
		System.setErr(redirErr);

		try {
			runSplit();

		}
		catch (Exception e) {
			throw e;
		}
		finally {
			System.setErr(orgErrStream);
		}



	}


	public void runSplit() throws Exception {

		if(maxCount != 0 ) { 
			int denom=0; 
			for(int i=0; i<ap.np; i++) 
				if(primerMove_depends_on[i]!=-1) 
					denom++; 
			maxCount=maxCount*denom/ap.np; 
			if(maxCount == 0) 
				maxCount=1; 
		} /* pairs */

		//		int[] shared_moves;
		initMoves();

		/* =0 to stop warnings on old compilers */
		//		if(memFail(shared_moves,_memFail)) 
		//			return NULL;
		if(numMoves == 0) {
			System.err.format("Can't move anything!\n");
			//			free(shared_moves); 
			return;
		}

		bestPools = new int[ap.np];
		//		if(memFail(shared_moves,bestPools,_memFail)) return NULL;
		Instant start = Instant.now();
		//		ThreadLocalRandom.current().setSeed(seedless ? 1 : start.getEpochSecond());


		Instant limitTimeToRun = null;// Instant.ofEpochMilli(0); 
		if(timeLimitInMins != 0) { 
			limitTimeToRun = Instant.now().plusSeconds( timeLimitInMins*60 ); /* (timeLimit is in minutes) */
		}

		stop_state = 0 ; //  s_KeepGoing; 
		//		signal(SIGINT, intHandler);
		//		if(omp_get_max_threads() > 1) {
		//			if(seedless) {
		//				omp_set_num_threads(1); 
		//				fputs("NOT parallelising the pool trials, as you asked for predictability.\n",stderr); 
		//				}
		//			else 
		//				fprintf(stderr,"Parallelising pool trials: %d threads\n",omp_get_max_threads());
		//		}
		System.err.format("OK, here goes... (press Control-C to stop%s)\n",timeLimitInMins != 0 ?" early":""); 
		//		fflush(stderr);



		//		#if PARALLELIZE_POOLSPLIT && defined(_OPENMP)
		//		#pragma omp parallel
		//		#endif
		poolsplit_thread(
				//				shared_moves,
				//				ap,
				//				nPools,
				//				numMoves,
				//				primerMove_depends_on,
				//				fix_to_pool,scores,
				limitTimeToRun
				//				,
				//				bestPools,
				//				table, 
				//				&bestPools_init_yet,
				//				&gBadLast,
				//				&totalIterations,
				//				&lastOutputTime,
				//				&overlaps,
				//				&just_printed_counts,
				//				&threads_needing_to_reset_iter,
				//				maxCount
				);
		//		signal(SIGINT, SIG_DFL);
		if(just_printed_counts == 0) {
			System.err.format("... looks like this is the best I can do:\n");
			overlaps = table != null ? ap.dGprintPooledCounts(bestPools,scores) : ap.printPooledCounts(bestPools,scores);
			int[] counts = PoolSplit.numInEachPool(bestPools,ap.np,nPools,primerMove_depends_on);
			if(counts != null) {
				PoolSplit.printNumInEachPool(counts,nPools);
				//				free(counts);
			}
		}
		long numSecs = (Instant.now().getEpochSecond() - start.getEpochSecond());

		if(numSecs == 0) numSecs=1; /* so division doesn't crash */

		System.err.format("%d moves",totalIterations);

		HelperMethods.prnSeconds(numSecs); 

		System.err.format(" = %d/sec\n",totalIterations/numSecs);
		if(bestPools != null && overlaps != 0 ) {
			//			if(stop_state == s_tooManyIters) 
			//				System.err.format("WARNING: There are still overlaps in these pools,\neven after this number of moves.\nYou might need more pools.\n");
			//			else 
			//				System.err.format("WARNING: There are still overlaps in these pools.\nMaybe you should have let it run longer\nto see if these overlaps can be eliminated.\n");
		} 



	}

	// Thread related ??
	long threads_needing_to_reset_iter = 0;

	// shared between threads
	int bestPools_init_yet = 0; 
	long gBadLast=0;
	/* latter =0 to stop warnings on old compilers */
	long totalIterations = 0;
	Instant lastOutputTime = Instant.ofEpochMilli(0);
	/* so 1st maxima gets output no matter what (might be needed if break after) */
	int just_printed_counts = 0, overlaps = 0;


	// running stats to support interruption 
	int stop_state = 0;  // 0 s_KeepGoing , 1 s_ccPressed, 2 s_tooManyIters

	public void setCC() {
		stop_state = 1;
		
	}
	
	
	void poolsplit_thread(
			//			final List<Integer> shared_moves,
			//			AllPrimers ap,
			//			int nPools,
			//			int numMoves,
			//			final int[] primerMove_depends_on,
			//			final int[] fix_to_pool,
			//			final int[] scores,
			Instant limitTimeToRun
			//			,
			//			int *bestPools,
			//			final float[] table, 
			//			int* bestPools_init_yet,
			//			long[] gBadLast,
			//			long *totalIterations,
			//			Instant * lastOutputTime, 
			//			int *overlaps,
			//			int* just_printed_counts,
			//			ThreadMask* threads_needing_to_reset_iter,
			//			int maxCount
			) throws Exception {
		/* This is the inner part of split_into_pools.
		     Multiple instances may be called in parallel. */
		int iter = 0; 
		boolean willContinue= true;
		int shared_movesIndex = 0 ; // int *moves = (int*)shared_moves;
		long[] bContrib = new long[ap.np*nPools];
		int[] poolCounts = new int[(nPools)];
		int[] pools = null;
		//		  if(memFail(bContrib,poolCounts,_memFail))
		//		    willContinue = 0;
		//		  else 
		{
			pools = new int [ap.np];
			//		    if(memFail(pools,_memFail)) willContinue = 0;
			//		    else 
			if(USE_QSORT) { /* moves must be per-thread */
				//		      moves = malloc(numMoves*sizeof(int));
				//		      if(memFail(moves,_memFail)) 
				//		    	  willContinue = 0;
				//		      else 
				//		    	  wrapped_memcpy(moves,shared_moves,numMoves*sizeof(int));
			}
		}
		long myMask = ((long)1) <<  0 ;// omp_get_thread_num(); /* for threads_needing_to_reset_iter */
		if (myMask == 0) { // this is not guarantee 
			/* what, somebody's running us on >128 cores ?? (or >64 32-bit cores) */
			/* (versions below v1.16 would hit this after 32 cores and not detect it) */
			//		    #if defined(_OPENMP)
			//		    #pragma omp critical
			//		    #endif
			//		    fprintf(stderr,"Can't run thread number %d because ThreadMask type has only %d bits\n",omp_get_thread_num(),(int)sizeof(ThreadMask)*8);  /* If you hit this, I suggest you either find a wider ThreadMask type or else we'd better make it an array.  Haven't done it so far because I've tested only on a 4-core machine and I doubt the chances of being run on many more cores than that are particularly high in 2016 (future might be different) */
			//		    willContinue = 0;
		}
		int max_iterations = 10000000 ; /* TODO: customise? profile? (but low priority as we have an interrupt mechanism) */
		//		    / (omp_get_num_threads() > 10 ? 10 : omp_get_num_threads()); 
		/* TODO: customise this "10" as well? (it's maxMoves / minMoves) */
		while(willContinue ) {
			PoolSplit.randomise_pools(ap.np,primerMove_depends_on,
					fix_to_pool,
					scores,
					nPools,
					pools,
					bContrib,
					poolCounts,
					maxCount);
			for(; ; iter++) {
				//		      #if USE_QSORT
				//		      #if PARALLELIZE_POOLSPLIT && defined(_OPENMP)
				//		      #pragma omp critical
				//		      #endif
				{
					//		        qsort_pools = pools; qsort_nPools = nPools;
					//		        qsort_bContrib = bContrib;
					//		        qsort_poolCounts = poolCounts;
					//		        qsort_maxCount = maxCount;
					//		        qsort(moves,numMoves,sizeof(int),betterMoves1st);
				}
				//		      int bestMove = moves[0];
				//		      #else
				int bestMove = findBestMove(shared_movesIndex,
						//		    		  numMoves,
						nPools,pools,bContrib,poolCounts,maxCount);
				//		      #endif
				if( (threads_needing_to_reset_iter & myMask) != 0  ) {
					//		        #if PARALLELIZE_POOLSPLIT && defined(_OPENMP)
					//		        #pragma omp critical
					//		        #endif
					{
						threads_needing_to_reset_iter &= ~myMask;
						totalIterations += iter;
					} 
					iter = 0;
				}
				if(stop_state == 1)
					System.out.println("stop_state = " + stop_state);
				boolean timesUp = stop_state != 0 || (limitTimeToRun != null && Instant.now().isAfter(limitTimeToRun) );
				if(timesUp || valueOfMove(bestMove,nPools,pools,bContrib,poolCounts,maxCount) == 0) {
					/* looks like we're at a local maxima */
					willContinue = !timesUp;
					long gBad = PoolSplit.globalBadness(scores,ap.np,pools);
					boolean keep = bestPools_init_yet == 0 || gBad < gBadLast;
					if (keep)
					{
						//		          #if defined(_OPENMP)
						//		          #pragma omp critical
						//		          #endif
						keep = bestPools_init_yet == 0 || gBad < gBadLast ;
						if (keep) {
							bestPools_init_yet = 1;
							// copy from pools to bestPools
							//		        	  	wrapped_memcpy(bestPools,pools,ap.np*sizeof(int));
							bestPools = pools.clone();
							gBadLast = gBad; 
							totalIterations += iter;
						}
					}
					if(gBad < (long)1<<48) 
						willContinue = false; 
					// everything down to score 0 - can't very much improve on that (except for reducing # pools or size difference)
					if (keep) {
						iter = 0;
						if(Instant.now().isAfter( lastOutputTime.plusSeconds(2) ) ) {
							int should_print_counts = 0;
							//		            #if PARALLELIZE_POOLSPLIT && defined(_OPENMP)
							//		            #pragma omp critical
							//		            #endif
							if(Instant.now().isAfter( lastOutputTime.plusSeconds(2) )) {
								lastOutputTime = Instant.now();
								threads_needing_to_reset_iter = ~0;
								should_print_counts = 1;
							} 
							if(should_print_counts != 0) {
								overlaps = table != null ?  ap.dGprintPooledCounts(pools,scores) : ap.printPooledCounts(pools,scores);
								PoolSplit.printNumInEachPool(poolCounts,nPools);
								if(!willContinue) {
									just_printed_counts = 1; 
									break; 
								}
								System.err.format("Local maxima found after %d moves\nTrying to better it... (press Control-C to stop)\n",
										totalIterations+iter); /* TODO: what about the 'iter' values of other threads? (or just don't count them yet) */
								System.err.flush(); /* in case of broken Windows/WINE etc (see comments in user.c) */
							}
						}
					} else {
						/* this maxima doesn't beat the best we've seen */
						if(iter>max_iterations && stop_state == 0) {
							System.err.format("Too many moves without improvement: giving up\n");
							willContinue=false; /* and stop other threads: */
							stop_state = 2; // s_tooManyIters;
						}
					}
					if(!willContinue) break;
					if(keep) {
						/* already found a good local maxima, so just take a few random steps away from it... */
						int randomMoves = 5+ Math.abs(ThreadLocalRandom.current().nextInt()) %5;
						for(int i=0; i<randomMoves; i++) {
							int moveToMake = Math.abs(ThreadLocalRandom.current().nextInt()) % numMoves;
							if(maxCount != 0) 
								while(poolCounts[poolOfMove(shared_moves.get( shared_movesIndex+ moveToMake),nPools,pools)]==maxCount) 
									if(++moveToMake==numMoves) moveToMake=0;
							make_a_move(shared_moves.get( shared_movesIndex+ moveToMake),ap.np,scores,primerMove_depends_on,nPools,pools,bContrib,poolCounts,maxCount);
						} 
						continue; /* don't do the additional make_a_move below (we'd have to repeat the maxCount condition) */
					} else {
						/* local maximae getting worse...
		             get me out of here! */
						break;
					}
				}
				//		      #if USE_QSORT
				//		      int i = 0;
				//		      while(!(ThreadRand()%5) && i<numMoves-1 && valueOfMove(moves[i+1],nPools,pools,bContrib,poolCounts,maxCount)) 
				//		    	  ++i; 
				/* sometimes don't pick the best one, just in case (TODO: can we write code to "get the top N items" w/out a complete sort?) */
				//		      bestMove = moves[i];
				//		      #endif
				make_a_move(bestMove,ap.np,scores,primerMove_depends_on,nPools,pools,bContrib,poolCounts,maxCount);
			}
		}
		//		  if(bContrib) free(bContrib);
		//		  if(pools) free(pools);
		//		  free(poolCounts);
		//		  #if PARALLELIZE_POOLSPLIT && defined(_OPENMP)
		//		  #pragma omp critical
		//		  #endif
		totalIterations += iter;
	}



	private void initMoves() {
		if(nPools <= 1) return ;
		//		  int[] moves= new int[(np*(nPools-1))];
		shared_moves = new ArrayList<Integer>();
		//		  if(moves ) 
		{
			//		    int movesP =  0;
			for(int i=0; i<ap.np*(nPools-1); i++) {
				int primer = PoolSplit.primerOfMove(i,nPools);
				if(primerMove_depends_on[primer]==-1
						&& fix_to_pool[primer]==-1) 
					//		    	  moves[movesP++] = i;
					shared_moves.add(i);
			}
			//		    numMoves = movesP ;
			//		    moves = memTrim(moves,movesP);
		}
		numMoves = shared_moves.size();
	}

	static final void make_a_move(int m,int np,
			final int[] scores,
			final int[] primerMove_depends_on,
			int nPools,
			int[] pools,
			long[] bContrib,
			int[] poolCounts,
			int maxCount) {
		int primer = PoolSplit.primerOfMove(m,nPools),
				oldPool = PoolSplit.oldPoolOfMove(m,nPools,pools),
				newPool = poolOfMove(m,nPools,pools);
		assert(primer >= 0 && primer < np && oldPool>=0 && newPool>=0 && oldPool<nPools && newPool<nPools && oldPool != newPool);
		pools[primer] = newPool; /* 'm' changes meaning now */
		assert(poolCounts[oldPool] !=0 );
		poolCounts[oldPool]--; poolCounts[newPool]++;
		assert( maxCount == 0 || poolCounts[newPool]<=maxCount);
		int i; 
		for(i=0; i<np; i++) {
			if(primerMove_depends_on[i]==primer)
				pools[i] = newPool; /* see merge_scores_of_stuckTogether_primers (and DON'T need to update poolCounts here) */
			else 
				PoolSplit.badnessContribUpdate(i,scores,np,primer,oldPool,nPools,pools,bContrib,poolCounts);
		}
	}


	int findBestMove(
			int movesIndex,
			//			int numMoves,
			int nPools,
			final int[] pools,
			final long[] bContrib,
			final int[] poolCounts,
			int maxCount) {
		int bestMove = this.shared_moves.get(movesIndex);
		long bestVal = valueOfMove(bestMove,nPools,pools,bContrib,poolCounts,maxCount);
		//		  #if PARALLELIZE_BESTMOVE && defined(_OPENMP)
		//		  #pragma omp parallel
		//		  #endif
		{
			int priv_bestMove = bestMove;
			long priv_bestVal = bestVal;
			int i;
			//		    #if PARALLELIZE_BESTMOVE && defined(_OPENMP)
			//		    #pragma omp for schedule(static)
			//		    #endif
			for(i=1; i<numMoves; i++) {
				long thisVal = valueOfMove(shared_moves.get(movesIndex+i),nPools,pools,bContrib,poolCounts,maxCount);
				if(thisVal > priv_bestVal) {
					priv_bestVal = thisVal;
					priv_bestMove = shared_moves.get(movesIndex+i);
				}
			}
			if(priv_bestVal > bestVal) {
				//		      #if PARALLELIZE_BESTMOVE && defined(_OPENMP)
				//		      #pragma omp critical
				//		      #endif
				if (priv_bestVal > bestVal) {
					bestVal = priv_bestVal;
					bestMove = priv_bestMove;
				}
			}
		} 
		return bestMove;
	}
	static final long valueOfMove(int m,int nPools,
			final int[] pools,
			final long[] bContrib,
			final int[] poolCounts,
			int maxCount) {
		if(maxCount !=0  && poolCounts[PoolSplit.poolOfMove(m,nPools,pools)]==maxCount) 
			return 0;
		assert(maxCount == 0 || poolCounts[PoolSplit.poolOfMove(m,nPools,pools)] < maxCount);
		long from = bContrib[PoolSplit.primerAndPool_to_contribOffset(PoolSplit.primerOfMove(m,nPools),PoolSplit.oldPoolOfMove(m,nPools,pools),nPools)], /* what this primer was contributing to its old pool */
				to = bContrib[ PoolSplit.move_to_contribOffset(m,nPools,pools)]; /* what this primer will contribute to its new pool */
		assert( to == 0 || poolCounts[PoolSplit.poolOfMove(m,nPools,pools)] != 0  );
		return PoolSplit.valueOfReduction(from,to);
	}

	static final int poolOfMove(int m,int nPools,final int[] pools) 
	{ 
		return (PoolSplit.oldPoolOfMove(m,nPools,pools)+1+(m % (nPools-1))) % nPools;
	}

	




}
