package org.primer3.search;

public interface ConflictResolution<T> {

	/**
	 * true of o is a valid object, this to avoid wasting time in trying different combination if the current one is not valid
	 * @param o
	 * @return
	 */
	boolean currentValid(T o );

	/**
	 * 
	 * @param o1
	 * @param o2
	 * @return true if it is a valid combination
	 */
	boolean checkConflict(T o1 , T o2);
	
	
	
	void setReslotion(int resoltionDiff);
}