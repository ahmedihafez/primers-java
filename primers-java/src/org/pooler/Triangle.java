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

import java.time.Instant;
import java.util.Date;

public class Triangle {

	public static class TwoRanges {
		Range[] r = {new Range(),new Range()};
	}

	public static class Range {
		int start, end;
	}
	public static  TwoRanges t_iBounds(int n) {
		/* distribute bounds to evenly schedule among threads.
		     Work units should look like:
		     0..1 and n..n
		     1..2 and n-1..n
		     2..3 and n-2..n-1
		     ...
		     n/2..n/2+1 and n-(n/2)..n-(n/2)+1 (if !=) */
		TwoRanges t = new TwoRanges();
		int nThreads =  1 ;// omp_get_num_threads(),
		int perThread = n/2/nThreads; /* NB it's rounded down */
		int tNum = 0;//omp_get_thread_num();
		if( perThread == 0) {
			/* unlikely but we should cover it */ 
			nThreads=n/2; 
			if(tNum>=nThreads)  { 
				t.r[0].start = t.r[0].end = t.r[1].start = t.r[1].end=0; 
				return t; 
			} 
			else 
				perThread=1; 
		}
		t.r[0].start = tNum*perThread;
		if(tNum==nThreads-1) {
			/* careful, due the above round-down */
			t.r[0].end = t.r[1].start = n/2; /* make SURE they meet in the middle and nothing gets left out */
		} else {
			t.r[0].end = (tNum+1)*perThread;
			t.r[1].start = n-(tNum+1)*perThread;
		}
		t.r[1].end = n-tNum*perThread;
		if(t.r[0].start==t.r[1].start) t.r[1].start=n;
		return t;
	}

	public static Instant t_ProgressStart(String p) {
		System.err.println(p);
		System.err.flush();
		return java.time.Instant.now();
	}
	
	static public Instant t_Progress(String p,TwoRanges tr,int n,int done,Instant next) {
//		  if(!omp_get_thread_num() /* only thread 0 prints */ && time(NULL) >= *next) {
		    Instant rnext =   Instant.now();
		    if(true && rnext.isAfter(next.plusSeconds(2) ) ) {
		    	System.err.format("\r%s%d%%",p,100*done/(n*(tr.r[0].end-tr.r[0].start))); 
		    	System.err.flush();
		    	return rnext;
		    }
		    return next;
		  }

	public static final int t_Nitems(int n) {
		  /* total number of items in an n*n triangle - this is
	     just normal triangle-numbers theory */
	  return n*(n+1)/2;
	}
	
}
