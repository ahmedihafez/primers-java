package org.primer3.libprimer3;

/**
 *  Enum to define tasks primer3 can do
 */
public enum P3Task { 
	/*  For backward compatibility, equivalent to
	 * generic
	 * plus pick_left_primer =1
	 * plus pick_right_primer = 1
	 * plus pick_internal_oligo = 0 
	 */
	
	PICK_PCR_PRIMERS              (0),

	PICK_PCR_PRIMERS_AND_HYB_PROBE (1),
	PICK_LEFT_ONLY                 (2),
	PICK_RIGHT_ONLY                (3),
	PICK_HYB_PROBE_ONLY            ( 4),
	GENERIC                        ( 5),
	PICK_CLONING_PRIMERS           ( 6),
	PICK_DISCRIMINATIVE_PRIMERS    ( 7),    
	PICK_SEQUENCING_PRIMERS        ( 8),
	PICK_PRIMER_LIST               ( 9),
	CHECK_PRIMERS                  ( 10);

	private int id ;
	P3Task(int id ){
		this.id = id;
	}
	public int getValue() { return id; }
}