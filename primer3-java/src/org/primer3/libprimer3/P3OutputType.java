package org.primer3.libprimer3;

/**
 *  Enum explaining if output are pairs 
 */
public enum P3OutputType {
	primer_pairs    (0),
	primer_list     (1);

	private int id ;
	P3OutputType(int id ){
		this.id = id;
	}
	public int getValue() { return id; }
}