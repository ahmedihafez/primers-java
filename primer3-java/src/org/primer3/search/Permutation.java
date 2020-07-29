package org.primer3.search;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

@Deprecated
// TODO :: use loop instead of the recursion 
public class  Permutation<T> {

	// lets say it could be any thing
	T[] source;
	int k;
	IterR rIt ;
	int[] current;
	public Permutation(T[] source) {
		this.source = source;
		this.k = source.length;

		rIt = new IterR(source, k, 0);
//		Iterator<Integer>[] its = new Iterator[k];
		current = new int[k];
		if(!rIt.init(current))
			hasMore = false;
	}

	class IterR {
//		int[] source;
		int kLevel ;
		int it;
		IterR childNext ;
		int k;
		public IterR( T[] source,int k, int level)
		{
			this.k = k;
//			this.source = source;
			this.kLevel = level;
//			it = lens.iterator();
			
			if( kLevel <  k - 1  )
			{
				childNext = new IterR( source,k,kLevel +1 );
			}
		}
		
		
		// each time a it is moving next it need to rest its child
		boolean reset(int[] current) {
			if(kLevel == 0) {
				it = 0;
			}
			else
			{
				// TODO :: offset here is assuming that the list is sorted based on the len :: otherwise it will fail
				// TODO :: ignore offset and use checkConflict to order based on score for greedy selection ???
				// 
				it = 0  ;
			}
			if(it >= k)
				return false;
			current[kLevel] = it ;

			if(childNext!= null)
				return childNext.reset(current);
			return true;
		}
		



		public boolean init(int[] current)
		{
			if(this.reset(current))
			{
				return true;
			}
			return false;	
			

		}
	
	
		public boolean move(int[] current)
		{
			
			if(childNext == null || !childNext.move(current))
			{
				if(it < k -1 ) {
					current[this.kLevel] = ++it ;
					if(childNext!= null)
						return childNext.reset(current);
					return true;
				}
				this.reset(current);
				return false;
				
			}
			
			return true;
		}
		
	
	}
	
	
	boolean hasMore = true;
	public ArrayList<T[]> getMorePerm(int n) {
		ArrayList<T[]> resList = new ArrayList<T[]>();
		for (int i = 0; i < n; i++) {
//			if(hasMore)
//			{
//				boolean acceptedPerm = true;
//				Set<Integer> currentPerm = new HashSet<Integer>();
//				for (int j = 0; j < current.length; j++) {
//					if(!currentPerm.add(current[j]))
//					{
//						acceptedPerm = false;
//						break;
//					}
//				}
//				if(acceptedPerm)
//					res.add(current.clone());
//				if( !rIt.move(current))
//					hasMore = false;
//			}
			T[] res = getPerm();
			if(res == null )
				break;
			resList.add(res);

		}
		return resList;
	}
	public T[] getPerm() {
	
		while(true) {
			if(hasMore)
			{
				boolean acceptedPerm = true;
				Set<Integer> currentPerm = new HashSet<Integer>();
				for (int j = 0; j < current.length; j++) {
					if(!currentPerm.add(current[j]))
					{
						acceptedPerm = false;
						break;
					}
				}
				if(acceptedPerm) {
					T[] res =  (T[]) Array.newInstance(source[0].getClass(), k);; 
					for (int i = 0; i < res.length; i++) {
						res[i] = source[current[i]];
					}
					if( !rIt.move(current))
						hasMore = false;
					return res;
				}
				if( !rIt.move(current))
					hasMore = false;
			}
			else {
				break;
			}
		}
		return null;
	}
	
	
	
	
	 private static int permutate(int p, int r) {
	        int result = 1;
	        for (int i = p; r > 0; r--, p--) {
	            result *= p;
	        }

	        return result;
	    }

	    private static int permutate(int p) {
	        return permutate(p, p);
	    }
	
	static  public void main (String[] args)
	{
		
		
		int k =4;
		String[] lens = new String[] {"a","b","c","d"};
		
		
		Permutation<String> t = new Permutation<String>(lens);
	
		int p = 0;
		for (String[] current : t.getMorePerm(10000))
		{
			if(p == 0 )
			{
				for (int i = 0; i < k; i++) {
					System.out.print(current[i] +" ");
				}
				System.out.print("\n");
				p=0;
			}
			else
				p--;

		}
		for (int i = 0; i < 21; i++) {
			System.err.println(permutate(i));
		}
		
		
	}
}
