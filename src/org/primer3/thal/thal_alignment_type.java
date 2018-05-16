package org.primer3.thal;

public enum thal_alignment_type {
  thal_any (1),
  thal_end1 (2),
  thal_end2 (3),
  thal_hairpin (4);
  private int id ;
  thal_alignment_type(int id ){
    	 this.id = id;
  }
  public int getValue() { return id; }
}