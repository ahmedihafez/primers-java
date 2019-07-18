package org.primer3.search;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Collection;
import java.util.ConcurrentModificationException;
import java.util.Iterator;
import java.util.List;
import java.util.NavigableSet;
import java.util.TreeSet;

/**
 * Iterator to iterate over all different possible combination in a given array of list 
 * list are sorted by the by a its own comparator so we could navigate easly  
 * Source are type of NavigableSet<T>
 * @author Ahmed Hafez
 *
 * @param <T>
 */
public  class  MultiListIterator<T> {


	static  public boolean flip = false;

	ConflictResolution<T> conflictResolution = null ;
	NavigableSet<T>[] sources;
	int[] lens;
	int k;
	//	IterR rIt ;
	Iterator<T>[] its;
	T[] current;


	boolean advanceSaveMove(int i)
	{
		do{
			if(its[i].hasNext())
			{
				current[i] = its[i].next();
			}
			else
			{
				return false;
			}
		}while( !this.conflictResolution.currentValid(current[i]) || 
				checkConflict(current,i) !=-1);
		return true;
	}

	boolean reset(int i) {
		// while reseting we could clear unwanted data
		its[i] = sources[i].iterator();

		if(!its[i].hasNext())
			return false;
		//		current[kLevel] = it.next();
		//		if(kLevel != 0){
		if(!advanceSaveMove(i))
		{
			return false;
		}
		//		}
		if(i < k - 1 ){ 
			while(!reset(i+1))
			{
				if(!advanceSaveMove(i)){
					return false;
				}
			}
		}
		return true;
	}

	public boolean init()
	{
		if(this.reset(0))
		{
			return true;
		}
		return false;
	}


	public boolean move(int i)
	{
//		boolean iDoNotMove = false;
//		if(conflictResolution.isFixing())
//		{
//			T currentIFixed = conflictResolution.isFixed(current,i);
//			if(currentIFixed != null)
//			{
//				current[i] = currentIFixed;
//				iDoNotMove= true;
//			}
//		}
		if ( !conflictResolution.currentValid(current[i]) ){
//			its[i].remove();
			if(advanceSaveMove(i)) {
				//				current[this.kLevel] = it.next();
				if(i < k - 1 ) {
					while(!reset(i+1))
					{
						if(!advanceSaveMove(i)){
							return false;
						}
					}
				}
				return true;
			}
			else
			{
				return false;
			}
		}


		if(i == (k-1) || !move(i+1))
		{
			if(advanceSaveMove(i)) {
				//				current[this.kLevel] = it.next();
				if(i < k - 1)
					while(!reset(i+1))
					{
						if(!advanceSaveMove(i)){
							return false;
						}
					}					
				return true;
			}
			this.reset(i);
			return false;

		}

		return true;
	}




	boolean isInit = false;
	public MultiListIterator(NavigableSet<T>[] sourceList , ConflictResolution<T> conflictResolution) {
		this.sources = sourceList;
		this.k = sourceList.length;
		this.its = new Iterator[this.k];
		this.lens = new int[this.sources.length];
		this.conflictResolution = conflictResolution;
	}

//	private boolean saveMove(int from , int to) {
//
//		for (int i = from; i <= to ; ) {
//			if(i == from && !its[i].hasNext() )
//				return false;
//			if(!its[i].hasNext()) {
//				its[i] = sources[i].iterator();
//				i--;
//				continue;
//			}
//
//			current[i] = its[i].next();
//
//			if(i > 0)
//			{
//				if(checkConflict(current, i)!=-1)
//				{
//					//					its[i] = lens[i].iterator();
//					//					i--;
//					continue;
//				}
//
//
//			}
//			i++;
//		}
//		return true;
//	}

//	 private boolean move() {
//
//		//		boolean reset = false;
//
//
//
//		for (int i = k-1; i >= 0; i--) {
//			// TODO :: check if i should be fixed based 
//			// if so do not move it and set it from outside
//			if(its[i].hasNext())
//			{
//				if(saveMove(i, k-1))
//				{
//					return true;
//				}
//			}
//			//			else
//			if(i>0)
//			{
//				its[i] = sources[i].iterator();
//			}
//
//		}
//		//		its[0] = lens[0].iterator();
//		if(!its[0].hasNext())
//			return false;
//		return true;
//	}


	boolean hasMore = true;
//	public ArrayList<T[]> getMoreComp(int n) {
//		ArrayList<T[]> res = new ArrayList<T[]>();
//		for (int i = 0; i < n; i++) {
//			if(hasMore)
//			{
//				res.add(current.clone());
//				if( !move())
//					hasMore = false;
//			}
//
//		}
//		return res;
//	}
	public T[] getComp() {

		if(!isInit)
		{
			isInit = true;
			if(sources[0].size() == 0 )
			{
				hasMore = false;
				return null;
			}
			for (int i = 0; i < sources.length; i++) {
				this.lens[i] = sources[i].size();
			}
			current = (T[]) Array.newInstance(sources[0].iterator().next().getClass(), k) ;
			if(!init())
			{	
				hasMore = false;
				return null;
			}
			
			
			return current.clone();


		}
		currentCheck();

//		for (int i = 0; i < sources.length; i++) {
//			if(sources[i].size() != lens[i] )
//			{
//				System.err.println("List has been changed");
//			}
//			
//			if ( !conflictResolution.currentValid(current[i]) ){
//				try  {
//					its[i].remove();
//					lens[i]--;
//				}
//				catch (IllegalStateException iex) {
//					System.err.println("IllegalStateException");
//				}
//				catch (ConcurrentModificationException ex )
//				{
//					System.err.println("ConcurrentModificationException");
//				}
//			}
//			// TODO :: CHECK if we need to break here and ??
////			if(lens[i].size() == 0)
////				hasMore = false;
//			
//		}

		if(hasMore)
		{

			if(move(0)) {
				T[] res =  current.clone();	
				return res;		
			}
			hasMore = false;

		}
		return null;
	}

	private int checkConflict(T[] current, int i) {		
		for(int j = i-1 ; j >= 0;j-- )
		{
			if(!conflictResolution.checkConflict(current[i],current[j]))
			{
				return i;
			}
		}	
		return -1;
	}



	public static final ConflictResolution<Integer> CheckSetProductsLen = new  ConflictResolution<Integer>() {
		@Override
		public boolean currentValid(Integer o) {
			if(o == 200 ) 
			{
				if(flip)
				{	
					//					flip = false;
					return false;
				}
				//				flip = true;
			}

			return true;
		}

		@Override
		public boolean checkConflict(Integer productLen1, Integer productLen2) {
			boolean isValid = true;


			int diff =  Math.abs(productLen1-productLen2);

			int minLen = Integer.min(productLen1,productLen2);

			int diffMin = 45 ;


			if(minLen < 500  )
				diffMin = 45;
			else if(minLen < 850)
				diffMin = 90;
			else diffMin = 180;
			if( diff <= diffMin ) // not exactly should be
			{
				isValid = false;
			}

			return isValid;
		}
	};



	static public int rt ()
	{
		int i ;
		
		return i = 5;
	}


	static  public void main (String[] args)
	{
		System.err.println(rt());
		
		NavigableSet<Integer> g1 = new TreeSet<Integer>();
		for (int i = 0; i < 10; i++) {
			g1.add(i);
		}
		
		Iterator<Integer> i1 = g1.tailSet(5).iterator();
		Iterator<Integer> i2 = g1.iterator();
		
		System.out.println(i2.next());
		System.out.println(i1.next());
		System.out.println(i1.next());
		i1.remove();
		
		System.out.println("after");
		i2 =g1.tailSet(0,false).iterator();
		System.out.println(i2.next());
		System.out.println(i2.next());
		System.out.println(i2.next());
		System.out.println(i2.next());
		System.out.println(i2.next());
		
		
		
		if(true )
			return;
		int k = 4;
		TreeSet<Integer> [] lens = new TreeSet[k];
		for (int i = 0; i < k ; i++) {
			lens[i] = new TreeSet<Integer>();
			for (int j = 100; j <= 500*(1+0); j+=100*(1+0)) {
				lens[i].add(j);
			}
		}



		//		int k =4;
		MultiListIterator<Integer> t = new MultiListIterator<Integer>(lens,MultiListIterator.CheckSetProductsLen );

		int p = 0;
		Integer[] current= null ;
		while( ( current = t.getComp() ) != null)
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
			MultiListIterator.flip = !MultiListIterator.flip;
		}
		System.err.println("Done");

	}

	
	
	public void currentCheck() {
		for (int i = 0; i < sources.length; i++) {
			if(sources[i].size() != lens[i] )
			{
//				System.err.println("List has been changed");
				its[i] = sources[i].tailSet(current[i],false).iterator();
				lens[i] = sources[i].size();
			}
			

			
		}
	}
	
	public void update() {
		for (int i = 0; i < sources.length; i++) {
			if(sources[i].size() != lens[i] )
			{
				System.err.println("List has been changed");
			}
			
			if ( !conflictResolution.currentValid(current[i]) ){
				try  {
					its[i].remove();
					lens[i]--;
				} 
				catch (IllegalStateException iex) {
					System.err.println("IllegalStateException");
				}
				catch (ConcurrentModificationException ex )
				{
					System.err.println("ConcurrentModificationException");
				}
			}
			// TODO :: CHECK if we need to break here and ??
//			if(lens[i].size() == 0)
//				hasMore = false;
			
		}
	}

	public boolean hasMore() {
		// TODO Auto-generated method stub
		return hasMore;
	}
	
	
	
	
}
