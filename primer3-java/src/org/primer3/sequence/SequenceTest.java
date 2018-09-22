package org.primer3.sequence;

import org.primer3.sequence.Sequence.SequenceIterator;

public class SequenceTest {
	public static void main(String[] argv){
	
		
		

		{
			System.out.println(   1e-5 + " " + 0.000001);

		}
		
//		testIterator();
	}

	static void testIterator()
	{
		System.out.println("Sequence Test");

		Sequence seq = new Sequence("Test Sequence".toCharArray());
		System.out.println("Sequence  : normal Iterator :");
		for(char c : seq)
		{
			System.out.println(c);
		}
		
		
		
		
		System.out.println("using seqIterator current :");

		SequenceIterator seqIterator = (SequenceIterator) seq.iterator();
		char c = seqIterator.current();
		while(seqIterator.hasNext()  )
		{
			System.out.println(c);
			seqIterator.next();
			c = seqIterator.current();

		}
		System.out.println("After finish pre is " + seqIterator.previous() );

		System.out.println("using seqIterator [s,e] :");

		seqIterator = (SequenceIterator) seq.iterator(1,5);
		
		while(seqIterator.hasNext()  )
		{
			c = seqIterator.current();
			System.out.println(c);
			seqIterator.next();
		}

		System.out.println("End Test");
	}
	
}
