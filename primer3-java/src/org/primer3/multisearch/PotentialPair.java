package org.primer3.multisearch;

import org.primer3.libprimer3.DPAlArgHolder;
import org.primer3.libprimer3.THAlArgHolder;
import org.primer3.primer.MultiTargetPrimerRecord;
import org.primer3.primer.PrimerPair;
import org.primer3.primer.PrimerRecord;

/**
 * Potential pair
 * @author Ahmed Hafez
 *
 */
public class PotentialPair {	
	PrimerRecord left,  right;
	MultiTargetScanner scanner;
	public int targetIndex;
	double estimatedQuilty = Double.MAX_VALUE;
	public PotentialPair(MultiTargetScanner scanner,PrimerRecord left,  PrimerRecord right , double q) {
//		this.i = i;
//		this.j = j;
		
		
		
		this.left = left;
		this.right =  right;
		this.estimatedQuilty = q;
		this.scanner = scanner;
	}
	
//	public double quality = Double.MAX_VALUE;
	public PrimerPair pairRecord =  null;
	public int pairStatus = PrimerPair.PAIR_UNCHARACTRIZED;
	
	
	/**
	 * pair critierial to select this could be the product len or the tm of the Product
	 */
	public int productCriterion;
	
	
	@Override
	public String toString() {
		return String.format("%.2f - %d " + 
	(pairStatus == PrimerPair.PAIR_FAILED ? "PAIR_FAILED" : pairStatus == PrimerPair.PAIR_OK ? "PAIR_OK" : "" )  , estimatedQuilty,productCriterion );
	}
	PotentialPairSet potentialPairParentSet;
	public void setParentSet(PotentialPairSet calcPairs) {
		this.potentialPairParentSet = calcPairs;
		
	}
	public void checkpair(	DPAlArgHolder dpal_arg_to_use, THAlArgHolder thal_arg_to_use) throws Exception {
		
		if(this.pairStatus != PrimerPair.PAIR_UNCHARACTRIZED)
			return;
		if (!left.OK_OR_MUST_USE() || !right.OK_OR_MUST_USE()) {
			this.pairStatus = PrimerPair.PAIR_FAILED;
			this.potentialPairParentSet.unValided(this);
			return;
		}

		
		PrimerRecord selectedLeft = left;
		PrimerRecord selectedRight = right;
		if (left instanceof MultiTargetPrimerRecord) {
			selectedLeft =  ((MultiTargetPrimerRecord)left).targetsToPrimer.get(scanner.targets.get(targetIndex));
			
		}
		if (right instanceof MultiTargetPrimerRecord) {
			selectedRight =  ((MultiTargetPrimerRecord)right).targetsToPrimer.get(scanner.targets.get(targetIndex));
			
		}

		pairRecord = new PrimerPair();
		this.pairStatus = pairRecord.characterize_pair(scanner.revals.get(targetIndex), selectedLeft  , selectedRight, dpal_arg_to_use, thal_arg_to_use, true);
		if (pairStatus == PrimerPair.PAIR_OK) {
			scanner.revals.get(targetIndex).best_pairs.expl.ok++;
			pairRecord.obj_fn(scanner.pa);
			scanner.revals.get(targetIndex).best_pairs.add_pair(pairRecord);
		}
		else
		{
//			System.err.println(selectedLeft.p3_get_ol_problem_string() +"/" +  selectedRight.p3_get_ol_problem_string());
			left.oligoProblems = selectedLeft.oligoProblems;
			right.oligoProblems = selectedRight.oligoProblems;
			this.potentialPairParentSet.unValided(this);
		}
			
	}
	public String getTargetName() {
		return scanner.targets.get(targetIndex);
	}
}
