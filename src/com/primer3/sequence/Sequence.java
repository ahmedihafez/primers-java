package com.primer3.sequence;

import java.util.Iterator;
import java.util.ListIterator;
import java.util.NoSuchElementException;

import com.primer3.libprimer3.oligo_type;
import com.primer3.libprimer3.PrimerRecord;
import com.primer3.libprimer3.seq_args;

/**
 * Sequence Class not part of primer3 lib
 * It represent Java implementation to seqs as char* in c
 * @author ahafez
 *
 */
public class Sequence  implements Iterable<Character>{

	public static char NULL = '\0';
	
	
	char[] sequence;

	int len ;
	
	
	public Sequence(char[] seq)
	{
		this.sequence = seq;
		this.len = sequence.length;
		checkSymmetry();
		calcGCContent();
	}
	
	@Override
	public Iterator<Character> iterator() {
		// TODO Auto-generated method stub
		return new SequenceIterator();
	}
	/**
	 * return a sub view of the sequence in range [s,e]
	 * @param s : start index
	 * @param e : end index
	 * @return
	 */
	public Iterator<Character> iterator(int s, int e) {
		return new SequenceIterator(s,e);
	}
	
	public class SequenceIterator implements ListIterator<Character> {

		int startIndex,lastIndex ;
		int currentIndex = -1;
		
		public SequenceIterator()
		{
			// set Range Indices
			startIndex = 0;
			lastIndex = Sequence.this.sequence.length - 1;
			currentIndex = -1;
			
		}
		
		public SequenceIterator(int subStart, int subEnd)
		{
			// set Range Indices
			// TODO:: check range
			
			
			startIndex = subStart;
			lastIndex = subEnd;
			currentIndex = subStart-1;
			
		}
		
		@Override
		public boolean hasNext() {
			return currentIndex < lastIndex;
		}

		@Override
		public Character next() {
			if(hasNext())
			{
				return Sequence.this.sequence[++currentIndex];
			}
			throw new NoSuchElementException();
		}
		
		
		public Character current()
		{
			// if 
//			if(currentIndex == -1)
//				currentIndex = 0;
//			if(currentIndex<startIndex || currentIndex > lastIndex)
//				throw new NoSuchElementException();
//			return Sequence.this.sequence[currentIndex];
			if(hasNext())
			{
				return Sequence.this.sequence[currentIndex+1];
			}
			return NULL; // return '\0';
		}
		
		
		public void remove() {
            throw new UnsupportedOperationException();
        }

		@Override
		public boolean hasPrevious() {
			return currentIndex >= startIndex;
		}

		@Override
		public Character previous() {
			if(hasPrevious())
			{
				return Sequence.this.sequence[currentIndex--];
			}
			throw new NoSuchElementException();
		}

		@Override
		public int nextIndex() {
			return currentIndex+1;
		}

		@Override
		public int previousIndex() {
			return currentIndex;
		}

		@Override
		public void set(Character e) {
			throw new UnsupportedOperationException();			
		}

		@Override
		public void add(Character e) {
			throw new UnsupportedOperationException();			
		}
	}
	
	
	
	
	
	
	
	public static void main(String[] argv){
		testSymmetry();
//		testSubSeq();
	}

	static void testIterator()
	{
		System.out.println("Sequence Test");
		
		Sequence seq = new Sequence(new char[]{'a','b'});
		System.out.println("Sequence  1 :");
		for(char c : seq)
		{
			System.out.println(c);
		}
		seq = new Sequence("Test Sequence".toCharArray());
		System.out.println("Sequence  2 :");
		for(char c : seq)
		{
			System.out.println(c);
		}
		
		
		
		System.out.println("using seqIterator :");

		
		System.out.println("using seqIterator current :");

		SequenceIterator seqIterator = (SequenceIterator) seq.iterator();
		char c = seqIterator.current();
		while(c != '\0'  )
		{
			
			
			System.out.println(c);
//			System.out.println(  i + " : " + seqIterator.next());
			seqIterator.next();
			c = seqIterator.current();

		}
		System.out.println("After finish pre is " + seqIterator.previous() );

		
		seqIterator = (SequenceIterator) seq.iterator();
		while(seqIterator.hasNext())
		{
			int i = seqIterator.nextIndex();
			
			System.out.println(i + " : " + seqIterator.current());
//			System.out.println(  i + " : " + seqIterator.next());
			seqIterator.next();

			if(i == 5)
				break;
		}
//		System.out.println("Using seqIterator.Next one more time  :" + seqIterator.next());

		System.out.println("using seqIterator go back:");
		while(seqIterator.hasPrevious())
		{
			System.out.println(seqIterator.nextIndex() + " : " + seqIterator.current());
			System.out.println( seqIterator.previousIndex()  + " <- " + seqIterator.previous());

		}
		System.out.println(seqIterator.nextIndex() + " C " + seqIterator.current());

		System.out.println("End Test");
	}
	
	static void testSymmetry()
	{
		
		String s = "AATT";
		
		Sequence seq = new Sequence(s.toCharArray());
		System.out.println(s + " issymmetry =  " + seq.symmetry());
		
		s = "TTCA";
		seq = new Sequence(s.toCharArray());
		System.out.println(s + " issymmetry =  " + seq.symmetry());
	}
	static void testSubSeq()
	{
		
		String s = "AATTSSS[TTTTB]";
		
		Sequence seq = new Sequence(s.toCharArray());
		System.out.println(" seq =  " + seq);
		
		int start = 0;
		int end =  seq.len-1;
		System.out.println(" subSeq  =  ["+start+","+end+"] =  "+ seq.subSeqRange(start, end));
		
		int x = seq.len;
		start = x - 4;
		System.out.println("last 4 char subSeq  =  ["+start+","+end+"] =  "+ seq.subSeqRange(start, end));

	}
	/** Returns 1 if the sequence is self-complementary or symmetrical; 0
	   otherwise
	*/
	public boolean symmetry() {
		
		return isSymmetry;
	}

	boolean isSymmetry = false;
	private void checkSymmetry() {
		
		int mp = this.len/2;
		if(this.len%2 == 1)
		{
			isSymmetry = false;
			return;
		}
		int i =0;
		char s , e;
		while(i<mp)
		{
			s = this.sequence[i];
			e = this.sequence[this.len-i-1];
			if( (s=='A' && e!='T') ||
				(s=='T' && e!='A') || 
				(e=='A' && s!='T') || 
				(e=='T' && s!='A')) {
				isSymmetry = false;
				return;
			}
			if( (s=='C' && e!='G') || 
				(s=='G' && e!='C') || 
				(e=='C' && s!='G') || 
				(e=='G' && s!='C')) {
				isSymmetry = false;
				return;
			}
			i++;
		}
		
		isSymmetry = true;
	}
	
	public int length() {
		return this.len;
	}

	
	/**
	 * subSeq [start,end]
	 * @param start
	 * @param end
	 * @return
	 */
	public Sequence subSeqRange(int start, int end) {
		
		int subLen = end - start + 1;
		char[] subSeq = new char[subLen];
		for(int i = 0 ; i < subLen;i++)
		{
			subSeq[i] = this.sequence[start+i];
		}
		return new Sequence(subSeq);
	}
	
	
	public static char[] subSeqRange(char[] sequence, int start, int end) {
		
		int subLen = end - start + 1;
		char[] subSeq = new char[subLen];
		for(int i = 0 ; i < subLen;i++)
		{
			subSeq[i] = sequence[start+i];
		}
		return subSeq;
	}
	
	public static char[] subSeq(char[] sequence, int start, int subLen) {
		char[] subSeq = new char[subLen];
		for(int i = 0 ; i < subLen;i++)
		{
			subSeq[i] = sequence[start+i];
		}
		return subSeq;
	}
	
	public String toString()
	{
		return new String(this.sequence);
	}

	
	double gcPercent;
	int gcCount;
	public double getGCPrercent() {
		return gcPercent;
	}
	
	private void calcGCContent()
	{
		int gcCount = 0;
		
		for(char c : sequence)
		{
			if(c == 'C' || c == 'G')
			{
				gcCount++;
			}
		}
		this.gcCount =gcCount;
		this.gcPercent = (double)gcCount/(double)this.len;
	}

	public char get(int i) {
		
		return sequence[i];
	}

	public char[] getSequence() {
		return sequence;
	}

	public Sequence getReverse() {
		char[] revSeq = new char[this.sequence.length];
		for(int i = 0 ; i < len;i++)
		{
			revSeq[i] = this.sequence[len-i-1];
		}
		return new Sequence(revSeq);
	}
	
	
	public static char[] _pr_substr(char[] seq, int start, int len) {
		return subSeq(seq, start, len);
	}

	public static char[] p3_reverse_complement(char[] seq)
	{
		char[] revSeq = new char[seq.length];
		
		for(int i = 0 ; i < seq.length;i++)
		{
			 switch (seq[seq.length-i-1])
			 {
			 	case 'A': revSeq[i]='T'; break;
			 	case 'C': revSeq[i]='G'; break;
			 	case 'G': revSeq[i]='C'; break;
			 	case 'T': revSeq[i]='A'; break;
			 	case 'U': revSeq[i]='A'; break;

			 	case 'B': revSeq[i]='V'; break;
			 	case 'D': revSeq[i]='H'; break;
			 	case 'H': revSeq[i]='D'; break;
			 	case 'V': revSeq[i]='B'; break;
			 	case 'R': revSeq[i]='Y'; break;
			 	case 'Y': revSeq[i]='R'; break;
			 	case 'K': revSeq[i]='M'; break;
			 	case 'M': revSeq[i]='K'; break;
			 	case 'S': revSeq[i]='S'; break;
			 	case 'W': revSeq[i]='W'; break;

			 	case 'N': revSeq[i]='N'; break;

			 	case 'a': revSeq[i]='t'; break;
			 	case 'c': revSeq[i]='g'; break;
			 	case 'g': revSeq[i]='c'; break;
			 	case 't': revSeq[i]='a'; break;
			 	case 'u': revSeq[i]='a'; break;

			 	case 'b': revSeq[i]='v'; break;
			 	case 'd': revSeq[i]='h'; break;
			 	case 'h': revSeq[i]='d'; break;
			 	case 'v': revSeq[i]='b'; break;
			 	case 'r': revSeq[i]='y'; break;
			 	case 'y': revSeq[i]='r'; break;
			 	case 'k': revSeq[i]='m'; break;
			 	case 'm': revSeq[i]='k'; break;
			 	case 's': revSeq[i]='s'; break;
			 	case 'w': revSeq[i]='w'; break;

			 	case 'n': revSeq[i]='n'; break;
		      }
			
		}
		
		
		return revSeq;
	}
	
	
	/**
	 * 
	 * @param seq
	 * @param start
	 * @param len
	 * @return { gcPercent, nCount }
	 */
	public static double[]  calcGCandN(char[] seq , int start, int len)
	{
		
		int gcCount = 0;
		int nCount = 0 ;
		int gcatCount = 0;
		for( int i = start ; i < start + len ;i++)
		{
			char c = seq[i];
			if('N' == c)
			{
				nCount++;
			}
			else
			{
				gcatCount++;
				if(c == 'C' || c == 'G')
				{
					gcCount++;
				}
			}
		}
		double gcPercent = 0;
		if(gcatCount > 0 )
			gcPercent = 100.0 * ((double)gcCount)/(double)gcatCount;
		return new double[]{gcPercent,nCount};
	}

	
	
	/**
	 * 's' is the sequence in which to find the stop codon.
	 * 'start' is the position of a start codon.
	 *
	 * There are two modes depending on 'direction'.
	 *
	 * If direction is 1:
	 *
	 * Return the index of the first stop codon to the right of
	 * 'start' in 's' (in same frame).
	 * If there is no such stop codon return -1.
	 *
	 * If direction is -1:
	 *
	 * If 'start' is negative then return -1.
	 * Otherwise return the index the first stop codon to left of 'start'.
	 * If there is no such stop codon return -1.
	 *
	 * Note: we don't insist that the start codon be ATG, since in some
	 * cases the caller will not have the full sequence in 's', nor even
	 * know the postion of the start codon relative to s.
	 */
	public static int find_stop_codon(char[] s, int start,
			int direction) {

	

		  if(s == null)
			  throw new IllegalArgumentException("find_stop_codon : Sequence can not be null");
		  if( !(direction == 1 || direction == -1) ) 
			  throw new IllegalArgumentException("find_stop_codon : direction can be only 1 or -1");
		  
		  int len = (s.length);
		  if( !(len >= 3) )
			  throw new IllegalArgumentException("find_stop_codon : len can not be less than 3");
		  if( !(start <= (len - 3)))
			  throw new IllegalArgumentException("find_stop_codon : Start postion exceed len of the sequence");
		
		  int increment = 3 * direction;
		  
		  if (start < 0) {
			    if (direction == 1)
			      while (start < 0) start += increment;
			    else
			      return -1;
			  }
		  
		 int p , q;
		 for ( p = start; p>= 0 && p < (len-3) ;  p+=increment){
			 if(Character.toUpperCase(s[p]) != 'T') continue;
			 q = p + 1;
			 if(Character.toUpperCase(s[q]) == 'A') {
				 q++;
				 if(Character.toUpperCase(s[q]) == 'G' || Character.toUpperCase(s[q]) == 'A' ) {
					 return p;	 
				 }  
			 }
			 else if (Character.toUpperCase(s[q]) == 'G') {
				 q++;
				 if(Character.toUpperCase(s[q]) == 'A' ) {
					 return p;
				 }
			 }
		 }
		 
		  
		return -1;
	}

	

	
	
	
	
	
}
