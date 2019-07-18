package org.primer3.multisearch;

import java.util.Comparator;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * Store and sort all different potential pairSet by the product criteria and average score of the primers in the list
 * for a specific target sequence
 * @author Ahmed Hafez
 *
 */
public class PotentialPairSetsCollection extends TreeSet<PotentialPairSet> {
	/**
	 * auto generated serial id 
	 */
	private static final long serialVersionUID = -5978203688435479112L;
		public PotentialPairSetsCollection() {
			super(new Comparator<PotentialPairSet>() {

				@Override
				public int compare(PotentialPairSet o1, PotentialPairSet o2) {
					int res =  Integer.compare(o1.productCriterion, o2.productCriterion);
					if (res == 0) 
						return res;
					return Double.compare(o1.getAverageScore(), o2.getAverageScore());
				}
			});
		}
		
		@Override
		public void clear() {
			// TODO Auto-generated method stub
			
			super.clear();
		}

		PotentialPairSet dummyObject = new PotentialPairSet(0);
		public SortedSet<PotentialPairSet> tailSet(Integer offset) {
			dummyObject.productCriterion = offset;
			return super.tailSet(dummyObject);
		}

		public SortedSet<PotentialPairSet> headSet(Integer offset) {
			dummyObject.productCriterion = offset;
			return super.headSet(dummyObject);
		}
		public PotentialPairSetsCollection execludeSet(Integer less, Integer larger)
		{
			dummyObject.productCriterion = larger;
			SortedSet<PotentialPairSet> tailSet = super.tailSet(dummyObject);
			dummyObject.productCriterion = less;
			SortedSet<PotentialPairSet> headSet = super.headSet(dummyObject);
			PotentialPairSetsCollection res = new PotentialPairSetsCollection();
			res.addAll(tailSet);
			res.addAll(headSet);
			return res;
		}

		public double getEstimatedScore() {
			if(size() == 0)
				return Double.MAX_VALUE;
			return first().getAverageScore();
		}
	}