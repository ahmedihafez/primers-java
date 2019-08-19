package org.primer3.multisearch;

import org.primer3.primer.MultiTargetPrimerRecord;
import org.primer3.primer.PrimerPair;

public class mPairs {

	/**
	 * 
	 */
	private final P3OptimzedMultiTargetFinder p3OptimzedMultiTargetFinder;

	PotentialPair[] pairSet;

	int setStatus = PrimerPair.PAIR_UNCHARACTRIZED;

	public mPairs(P3OptimzedMultiTargetFinder p3OptimzedMultiTargetFinder, PotentialPair[] pairSet)
	{
		this.p3OptimzedMultiTargetFinder = p3OptimzedMultiTargetFinder;
		if(pairSet != null )
		{
			this.pairSet = pairSet;
			estimatedQuilty = 0;
			for (int i = 0; i < pairSet.length; i++) {

				estimatedQuilty += pairSet[i].estimatedQuilty;
			}
		}
	}
	double estimatedQuilty = Double.MAX_VALUE;

	public int check() throws Exception {

		for (int p = 0; p < pairSet.length; p++) {
			PotentialPair pair = pairSet[p];
			//				PrimerRecord left = scanner.multiFwd.get(pair.j);
			//				PrimerRecord right = scanner.targetsToSpRev.get(scanner.targets.get(p)).get(pair.i);

			if(pair.pairStatus == PrimerPair.PAIR_OK)
			{

			}
			else if (pair.pairStatus == PrimerPair.PAIR_UNCHARACTRIZED) {
				pair.checkpair(this.p3OptimzedMultiTargetFinder.dpal_arg_to_use, this.p3OptimzedMultiTargetFinder.thal_arg_to_use);;
			}

			if(pair.pairStatus == PrimerPair.PAIR_FAILED) {	
				setStatus = PrimerPair.PAIR_FAILED;
				return setStatus;
			}

		}
		//			for (int p = 0; p < pairSet.length; p++) {
		//				PotentialPair pair = pairSet[p];
		//				for (int i = 0; i < p; i++) {
		//					if(!checkmTargetPrimers(pair, pairSet[i])) 
		//					{
		//						setStatus = PrimerPair.PAIR_FAILED;
		//						return setStatus;
		//					}
		//				}
		//				
		//			}



		setStatus = PrimerPair.PAIR_OK;
		return setStatus;
	}

	private boolean checkmTargetPrimers(PotentialPair o1, PotentialPair o2) {

		if(o1.left instanceof MultiTargetPrimerRecord && o2.left instanceof MultiTargetPrimerRecord)
		{
			//				// both are left
			MultiTargetPrimerRecord o1Left = (MultiTargetPrimerRecord) o1.left;
			MultiTargetPrimerRecord o2Left = (MultiTargetPrimerRecord) o2.left;
			//				// if o1 left has the o2 taget then they should be the same left other wise no
			if(     o1Left.targetsToPrimer.containsKey(o2.getTargetName())  || 
					o2Left.targetsToPrimer.containsKey(o1.getTargetName()))
			{
				return o1Left == o2Left ;
			}
		}
		return true;
	}

}