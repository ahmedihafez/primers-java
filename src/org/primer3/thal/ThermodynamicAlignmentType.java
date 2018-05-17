package org.primer3.thal;

public enum ThermodynamicAlignmentType {
  thal_any (1),
  thal_end1 (2),
  thal_end2 (3),
  thal_hairpin (4);
  private int id ;
  ThermodynamicAlignmentType(int id ){
    	 this.id = id;
  }
  public int getValue() { return id; }
}