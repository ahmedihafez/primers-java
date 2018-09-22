package org.primer3.oligotm;

/**
   If salt_corrections==schildkraut, then formula for
   salt correction in the paper [Schildkraut, C, and Lifson, S (1965)
   "Dependence of the melting temperature of DNA on salt
   concentration", Biopolymers 3:195-208 (not available on-line)] is
   used.  This is the formula that primer3 used up to and including
   version 1.0.1.

   If salt_corrections==santalucia, then formula for
   salt correction suggested by the paper [SantaLucia JR (1998) "A
   unified view of polymer, dumbbell and oligonucleotide DNA
   nearest-neighbor thermodynamics", Proc Natl Acad Sci 95:1460-65
   http://dx.doi.org/10.1073/pnas.95.4.1460] is used.

   *THIS IS THE RECOMMENDED VALUE*. 
  
   If salt_corrections==owczarzy, then formula for
   salt correction in the paper [Owczarzy, R., Moreira, B.G., You, Y., 
   Behlke, M.A., and Walder, J.A. (2008) "Predicting stability of DNA 
   duplexes in solutions containing magnesium and monovalent cations", 
   Biochemistry 47:5336-53 http://dx.doi.org/10.1021/bi702363u] is used.
 */
public enum SaltCorrectionMethod {
	 
	
	schildkraut(0,"schildkraut"),
	santalucia(1,"Santalucia 1998"),
	owczarzy(2,"owczarzy");
	
	
	private int id ;
	private String name;
	SaltCorrectionMethod(int id,String name ) {
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
