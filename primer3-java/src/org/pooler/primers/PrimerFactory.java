package org.pooler.primers;

import java.io.PrintStream;

import org.pooler.primers64.PrimerFactory64;

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

}
