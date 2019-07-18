package org.primer3.search.multiplex;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

import org.primer3.libprimer3.P3GlobalSettings;
import org.primer3.libprimer3.P3RetVal;
import org.primer3.multisearch.IP3MultiSearch;
import org.primer3.multisearch.P3MultiplexSearch;

public class MultiplexResult {
	String groupName;
	P3GlobalSettings pa;
	IP3MultiSearch p3MultiplexSearcher = null;
	// add just once
	HashMap<String,P3RetVal> sourcePairsResult = new HashMap<String, P3RetVal>();
	
	
	ArrayList<MultiplexSet> multiplexSets = new ArrayList<MultiplexSet>();
	
	
	
	
	public MultiplexResult(IP3MultiSearch p3MultiplexSearch,String groupName,  P3GlobalSettings pa) {
		p3MultiplexSearcher = p3MultiplexSearch;
		this.pa = pa;
		this.groupName = groupName;
	}
	public void addP3RetVal(P3RetVal retVal) {
		
		// just add once
		String targetName = retVal.sa.getSequenceName();
		if(!sourcePairsResult.containsKey(targetName))
		{
			sourcePairsResult.put(targetName, retVal);
			expectedMaxNTargetsInSet++;
		}
		
	}
	
	
	
	ArrayList<MultiplexSet> newAltsSets =  new ArrayList<MultiplexSet>();
	
	int numOfSetToReturn = 10;
	int numOfGoodSets = 0;
	int expectedMaxNTargetsInSet=0;
	public boolean hasGoodSets() {
		numOfGoodSets = 0;
		for(MultiplexSet set : multiplexSets)
		{
			if(set.getNTargets() >= expectedMaxNTargetsInSet)
			{
				numOfGoodSets++;
			}
		}
		return numOfGoodSets > numOfSetToReturn;
	}

	
	
	




	public static void main(String[] args)
	{
		HashMap<String, String> t = new HashMap<String, String>();
		t.put("a", "a");
		t.put("b", "c");
		t.put("c", "b");
		
		ArrayList<String> a = new ArrayList<String>(t.keySet());
		
		System.out.println(a);
		a.clear();
		System.out.println(a);
		System.out.println(t);
		
	}



	

	public void print_boulder(int io_version) {
		System.out.println("Result Sets for " + this.groupName);
		for(int i = 0 ; i < multiplexSets.size();i++) {
			System.out.println("Sets #" + i);
			multiplexSets.get(i).print_buolder(io_version);
				
		}
		
	}




	public void sort() {
		
		multiplexSets.sort(new Comparator<MultiplexSet>() {

			@Override
			public int compare(MultiplexSet o1, MultiplexSet o2) {
				
				if(o1.getNTargets() == o2.getNTargets())
					return Double.compare(o1.getScore(), o2.getScore());
				return -Integer.compare(o1.getNTargets(), o2.getNTargets());
			}
		});
		
	}
	
	



	public List<MultiplexSet> getMultiplexSet() {
		return multiplexSets;
	}
	

	public void addSet(MultiplexSet newSet) {
		multiplexSets.add(newSet);
	}

}
