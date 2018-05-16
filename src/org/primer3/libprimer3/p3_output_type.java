package org.primer3.libprimer3;

/**
    *  Enum explaining if output are pairs 
    */
   public enum p3_output_type {
     primer_pairs    (0),
     primer_list     (1);
     
     private int id ;
     p3_output_type(int id ){
    	 this.id = id;
     }
     public int getValue() { return id; }
   }