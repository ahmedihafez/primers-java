package com.primer3;

enum p3_file_type {
    all_parameters(0),
    sequence        (1),
    settings        (2);
    
	private int id ;
	p3_file_type(int id ) {
		this.id = id;
    }
    public int getValue() { return id; }
}