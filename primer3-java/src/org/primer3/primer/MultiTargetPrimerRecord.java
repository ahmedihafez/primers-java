package org.primer3.primer;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import org.primer3.dpal.AlignmentException;
import org.primer3.dpal.DPAlignmentArgs;
import org.primer3.libprimer3.OligoType;
import org.primer3.libprimer3.SeqArgs;
import org.primer3.p3_seq_lib.seq_lib;
import org.primer3.sequence.Sequence;

public class MultiTargetPrimerRecord  extends PrimerRecord {

	public HashMap<String, PrimerRecord> targetsToPrimer = new HashMap<String, PrimerRecord>();
	HashMap<String, Integer> primerStarts = new HashMap<String, Integer>();

	PrimerRecord orgPrimer;



	public MultiTargetPrimerRecord(String target , PrimerRecord newRecord) {
		super(newRecord.rec_type);
		addRecord( target ,  newRecord);
		this.orgPrimer = newRecord;
		this.oligoSeq = newRecord.getOligoSeq();
		this.oligoRevSeq = newRecord.getOligoRevSeq();

		// pairwise copy
		// this is target specific value, we are keeping it in a list and storing the smallest
		this.start = newRecord.start;
		
		this.length = newRecord.length;
		this.quality= newRecord.quality; 

		this.temp = newRecord.temp; 

		this.gc_content = newRecord.gc_content;

		this.position_penalty= newRecord.position_penalty;
		


		this.end_stability= newRecord.end_stability;
		


		this.seq_quality= newRecord.seq_quality; 
		this.seq_end_quality= newRecord.seq_end_quality; 

		this.self_any= newRecord.self_any; 

		this.self_end= newRecord.self_end; 

		this.hairpin_th= newRecord.hairpin_th; 

		this. template_mispriming= newRecord.template_mispriming;
	
		this. template_mispriming_r= newRecord.template_mispriming_r;

		this.length= newRecord.length; 
		this. num_ns= newRecord.num_ns; 

		this. must_use= newRecord.must_use; 
		this. overlaps= newRecord.overlaps; 

		this. oligoProblems = newRecord.oligoProblems;
		this. overlaps_overlap_position= newRecord.overlaps_overlap_position;

		this. template_mispriming_ok= newRecord.template_mispriming_ok; 

		this.failure_rate= newRecord.failure_rate; 

	}
	boolean isRepeatsInTarget = false;
	public void addRecord(String target, PrimerRecord newRecord) {
		// ensure integrity
		// we could also detect repeats here ??
		// TODO :: if we have same target then this is a repeat unless the start is differnt
		if(!targetsToPrimer.containsKey(target) ) {
			targetsToPrimer.put(target, newRecord);
			primerStarts.put(target, newRecord.start);

			// take the smallest start
			if(this.start >= newRecord.start)
				this.start = newRecord.start;
		}
		else
		{
			// TODO :: FIXME :: check this case
			if(newRecord.start != primerStarts.get(target))
				isRepeatsInTarget = true;
		}
	}
	public int getNTarget() {

		return this.primerStarts.size();
	}
	public PrimerRecord getOrgPrimer() {
		return orgPrimer;
	}
	@Override
	public String toString() {
		return  String.copyValueOf(oligoSeq) + " Q: " + this.quality + " NT : " + this.getNTarget() ;
	}
	
	
	@Override
	public void calcSpecific(seq_lib targets_lib, String reversePrefixTargets, DPAlignmentArgs dpAlignmentArgs) throws AlignmentException {
		
		calcSpecific(targets_lib, targetsToPrimer.keySet()  , reversePrefixTargets, dpAlignmentArgs);
	}
	@Override
	public char[] pr_oligo_sequence(SeqArgs sa) {
		return getOligoSeq();
	}
	@Override
	public char[] pr_oligo_rev_c_sequence(SeqArgs sa) {
		return getOligoRevSeq();
	}
	
}
