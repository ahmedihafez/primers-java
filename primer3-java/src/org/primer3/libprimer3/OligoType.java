package org.primer3.libprimer3;

public enum OligoType { 
	OT_LEFT (0),
	OT_RIGHT (1),
	OT_INTL (2) ;
	private int id ;
	OligoType(int id ) {
    	 this.id = id;
    }
    public int getValue() { return id; }	
}