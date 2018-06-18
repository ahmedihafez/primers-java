package org.primer3;

enum P3FileType {
    all_parameters(0),
    sequence        (1),
    settings        (2);
    
	private int id ;
	P3FileType(int id ) {
		this.id = id;
    }
    public int getValue() { return id; }
}