package org.pooler.primers;

public abstract class NonDegeneratePrimer implements IPrimer {

	
	@Override
	public final boolean isDegeneratePrimer() {
		return false;
	}
	@Override
	public final int NumPossibilities_32bases() {
		return 1;
	}
}
