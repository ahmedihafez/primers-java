package org.primer3.masker;

public enum masking_direction {
	 both_on_same(0),
	 both_separately (1),
	 fwd(2),
	 rev(3);
	 private int id ;
	 masking_direction(int id ) {
		 	this.id = id;
	 }
	 public int getValue() { return id; }
 }