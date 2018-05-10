package com.primer3.libprimer3;

public enum oligo_type { 
	OT_LEFT (0),
	OT_RIGHT (1),
	OT_INTL (2) ;
	private int id ;
	oligo_type(int id ) {
    	 this.id = id;
    }
    public int getValue() { return id; }	
}