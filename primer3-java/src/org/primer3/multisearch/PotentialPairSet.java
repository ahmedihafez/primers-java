package org.primer3.multisearch;

import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeSet;

/**
 * an Ordered set of potential pair sorted (Product or Tm according type of PCR used)
 * @author Ahmed Hafez
 *
 */
public class PotentialPairSet extends TreeSet<PotentialPair> {

	/**
	 * auto generated serialVerion ID 
	 */
	private static final long serialVersionUID = 6705010134823731730L;

	static Comparator<PotentialPair> potenialPairsComparator = new Comparator<PotentialPair>() {

		@Override
		public int compare(PotentialPair o1, PotentialPair o2) {
			if(o1 == o2)
				return 0;
//			if(o1.pairRecord == null || o2.pairRecord == null )
			{
				int res  = Double.compare(o1.productCriterion, o2.productCriterion);
				if (res == 0 )
				{
					res = Double.compare(o1.estimatedQuilty, o2.estimatedQuilty);
//					if(res == 0 )
//					{
//						res = Integer.compare(o1.i, o2.i);
//						if(res == 0 )
//							res = Integer.compare(o1.j, o2.j);
//					}
				}
				return res;
			}
//			int cRes2 =  LibPrimer3.compare_primer_pair(o1.pairRecord,o2.pairRecord);
//			return cRes2;
		}
	};
	
//		private final P3OptimzedMultiFinder p3OptimzedMultiFinder;
//		double minScore = Double.MAX_VALUE;
		public double validPairs = 0;
		double averageScore = 0;
		double orgMinScore = Double.MAX_VALUE;
//		double lastScore;
		
		// TODO :: this could be a range instead
		public int productCriterion;
		public PotentialPairSet(int productCriterion) {
			super(potenialPairsComparator);
			this.productCriterion = productCriterion;
//			this.p3OptimzedMultiFinder = p3OptimzedMultiFinder;
		}
		


		public void unValided(PotentialPair ijPair) {
//			validPairs--;
			toDelete.add(ijPair);

		}

//		@Override
//		public Iterator<ijPair> iterator() {
//			
//			this.deleteUnwanted();
//			return super.iterator();
//		}


		@Override
		public boolean add(PotentialPair e) {
			// TODO Auto-generated method stub
//			if(size() == 0)
			averageScore += e.estimatedQuilty;
			
//			else
//				averageScore = averageScore + ALPHA * (e.estimatedQuilty - averageScore);
			if(super.add(e))
			{
				e.setParentSet(this);
				validPairs++;
				if(orgMinScore >= e.estimatedQuilty)
					orgMinScore = e.estimatedQuilty;
				return true;
			}
			return false;
		}
		
		@Override
		public boolean remove(Object o) {
			
			if(o instanceof PotentialPair )
			{
				averageScore -=  ((PotentialPair) o).estimatedQuilty;
			}

			return super.remove(o);
		}
		
		HashSet<PotentialPair> toDelete = new HashSet<PotentialPair>();
		public void markDeletion(PotentialPair removePair)
		{
			toDelete.add(removePair);
		}
		
		public void deleteUnwanted(Iterator<PotentialPair> saveDelete) {
			
		}
		public void deleteUnwanted() {
			for(PotentialPair p : toDelete)
			{
				this.remove(p);
			}
			toDelete.clear();
		}
		@Override
		public String toString() {
			// TODO Auto-generated method stub
			return this.productCriterion + " " + super.toString();
		}

		public double getminScore() {
			if(size() == 0 )
				return Double.MAX_VALUE;
			return first().estimatedQuilty;
		}
		public double getOrgMinScore() {
			if(size() == 0 )
				return Double.MAX_VALUE;
			return orgMinScore;
		}
		public double getAverageScore() {
			// TODO Auto-generated method stub
			if(size() == 0 )
				return Double.MAX_VALUE;
			double scoreDiff =  ( this.last().estimatedQuilty -  this.first().estimatedQuilty ) /4;
//			return  this.first().estimatedQuilty + scoreDiff;
			return this.first().estimatedQuilty + scoreDiff ;// +( averageScore/size());
		}
		
		
		
		static public double calcScore(PotentialPairSet[] pairsList )
		{
			double score  = 0;
			for (int i = 0; i < pairsList.length; i++) {
				if(pairsList[i].size() == 0 )
					return Double.MAX_VALUE;
				score += pairsList[i].getminScore();
			}
			return score;
		}



		public static double getOrginalScore(PotentialPairSet[] pairsList) {
			double score  = 0;
			for (int i = 0; i < pairsList.length; i++) {
				if(pairsList[i].size() == 0 )
					return Double.MAX_VALUE;
				score += pairsList[i].getOrgMinScore();
			}
			return score;
		}
		
		@Override
		public Iterator<PotentialPair> iterator() {
			// TODO Auto-generated method stub
			return super.iterator();
		}
		
	}