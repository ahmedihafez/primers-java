package com.primer3.libprimer3;

/* Enum to define tasks primer3 can do */
   public enum task { 
     pick_pcr_primers              (0),
     /*  For backward compatibility, equivalent to
         generic
         plus pick_left_primer =1
         plus pick_right_primer = 1
         plus pick_internal_oligo = 0 */
     pick_pcr_primers_and_hyb_probe (1),
     pick_left_only                 (2),
     pick_right_only                (3),
     pick_hyb_probe_only            ( 4),
     generic                        ( 5),
     pick_cloning_primers           ( 6),
     pick_discriminative_primers    ( 7),    
     pick_sequencing_primers        ( 8),
     pick_primer_list               ( 9),
     check_primers                  ( 10);
     
     private int id ;
     task(int id ){
    	 this.id = id;
     }
     public int getValue() { return id; }
   }