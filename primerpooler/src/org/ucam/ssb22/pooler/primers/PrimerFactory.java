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
package org.ucam.ssb22.pooler.primers;

import java.io.PrintStream;

import org.ucam.ssb22.pooler.primers64.PrimerFactory64;

public abstract class PrimerFactory {

	static PrimerFactory currentInstance = null;
	static public PrimerFactory getCurrent() {
		return currentInstance;
	}
	static public void setCurrent(PrimerFactory primerFactory) {
		currentInstance = primerFactory;
	}
	public abstract IPrimer createPrimer(String primerSeq) ;
	
	abstract public IPrimer createPrimer(char[] primerSeq) ;
	
	
	
	abstract public  DegeneratePrimer upgradeToDegeneratePrimer(IPrimer primer) ;
	abstract public int getCurrentSize() ;
	
	
	public int getScore(IPrimer forwardPrimer, IPrimer backwordPrimer) {
		if(forwardPrimer.isDegeneratePrimer() && !backwordPrimer.isDegeneratePrimer())
		{
			return forwardPrimer.calcScore(upgradeToDegeneratePrimer(backwordPrimer));
		}
		// TODO :: here make sure that the order of primer does not matter here
		if(!forwardPrimer.isDegeneratePrimer() && backwordPrimer.isDegeneratePrimer())
		{
			return backwordPrimer.calcScore(upgradeToDegeneratePrimer(forwardPrimer));
		}
		return forwardPrimer.calcScore(backwordPrimer);
	}

	public float getDeltaG(IPrimer forwardPrimer, IPrimer backwordPrimer, float[] table) {
		if(forwardPrimer.isDegeneratePrimer() && !backwordPrimer.isDegeneratePrimer())
		{
			return forwardPrimer.calcDeltaG(upgradeToDegeneratePrimer(backwordPrimer),table);
		}
		if(!forwardPrimer.isDegeneratePrimer() && backwordPrimer.isDegeneratePrimer())
		{
			return upgradeToDegeneratePrimer(forwardPrimer).calcDeltaG(backwordPrimer,table);
		}
		return forwardPrimer.calcDeltaG(backwordPrimer,table);
	}
	
	public  CountResult getCount(IPrimer forwardPrimer, IPrimer backwardPrimer) {
		if(forwardPrimer.isDegeneratePrimer() && !backwardPrimer.isDegeneratePrimer())
		{
			return forwardPrimer.calcCount(upgradeToDegeneratePrimer(backwardPrimer));
		}
		if(!forwardPrimer.isDegeneratePrimer() && backwardPrimer.isDegeneratePrimer())
		{
			return upgradeToDegeneratePrimer(forwardPrimer).calcCount(backwardPrimer);
		}
		return forwardPrimer.calcCount(backwardPrimer);
	}
	public void dGprint(IPrimer forwardPrimer, IPrimer backwardPrimer, float minDG, PrintStream out, float[] table) {
		if(forwardPrimer.isDegeneratePrimer() && !backwardPrimer.isDegeneratePrimer())
		{
			forwardPrimer.dGprint(upgradeToDegeneratePrimer(backwardPrimer) , minDG, out,table);
			return ;
		}
		if(!forwardPrimer.isDegeneratePrimer() && backwardPrimer.isDegeneratePrimer())
		{
			upgradeToDegeneratePrimer(forwardPrimer).dGprint(backwardPrimer, minDG, out,table);
			return ;
		}
		forwardPrimer.dGprint(backwardPrimer, minDG, out,table);		
	}
	public void print(IPrimer forwardPrimer, IPrimer backwardPrimer, int maxScore, PrintStream out) {
		if(forwardPrimer.isDegeneratePrimer() && !backwardPrimer.isDegeneratePrimer())
		{
			forwardPrimer.print(upgradeToDegeneratePrimer(backwardPrimer) , maxScore, out);
			return ;
		}
		if(!forwardPrimer.isDegeneratePrimer() && backwardPrimer.isDegeneratePrimer())
		{
			upgradeToDegeneratePrimer(forwardPrimer).print(backwardPrimer, maxScore, out);
			return ;
		}
		forwardPrimer.print(backwardPrimer, maxScore, out);
	}

}
