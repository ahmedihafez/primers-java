package org.pooler.primers64;

import org.pooler.primers.DegeneratePrimer;
import org.pooler.primers.IPrimer;
import org.pooler.primers.PrimerFactory;

/**
 * Current Implementation is 64 bit 
 * @author Ahmed Hafez
 *
 */
public class PrimerFactory64 extends PrimerFactory {

	
	public IPrimer createPrimer(String primerSeq)  {
		
		boolean isDegn = false;
		if(!primerSeq.matches("[ACTGactg]*"))
		{
			isDegn = true;
		}
		if(isDegn)
		{
			return new DegeneratePrimer64(primerSeq);
		}
		else
		{
			return new Primer64(primerSeq);

		}
	}
	
	public IPrimer createPrimer(char[] primerSeq)  {
//		throw new Exception("Not Emplemented Yet");
		// TODO :: implement it
		return null;
	}
	
	
	
	public static void main(String[] argv) {
		
		System.out.println("ssssss".matches("[ATCGatcgBDHKMNRSVWYbdhkmnrsvwy]+"));
	}
	public  DegeneratePrimer upgradeToDegeneratePrimer(IPrimer primer) {
		if(primer.isDegeneratePrimer() )
		{
			return (DegeneratePrimer) primer;
		}
			
		Primer64 p = (Primer64) primer;
		return new DegeneratePrimer64(p.AorT,p.GorT,p.valid);
	}
	public int getCurrentSize() {
		return 64;
	}

}
