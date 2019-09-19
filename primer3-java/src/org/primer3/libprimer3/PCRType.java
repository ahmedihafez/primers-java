package org.primer3.libprimer3;

/**
  
 */
public enum PCRType {
	 
	
	NORMAL_PCR(0,"Conventional PCR"),
	QPCR(1,"Real Time PCR");
	
	
	private int id ;
	private String name;
	PCRType(int id,String name ) {
		this.id = id;
		this.name = name;
    }
    public int getValue() { return id; }
     
    public String getName()
    {
    	return name;
    }
    
    
    public String toString()
    {
    	return name;
    }
}
