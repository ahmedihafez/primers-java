/*
  This file is part of Java porting of Primer Pooler (https://github.com/ssb22/PrimerPooler)
  Primer Pooler (c) Silas S. Brown.  For Wen.
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
    http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
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
