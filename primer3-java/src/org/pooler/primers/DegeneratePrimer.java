package org.pooler.primers;

import org.pooler.Amplicons.PrimerToFind;

public abstract class DegeneratePrimer implements IPrimer {

	@Override
	public final boolean isDegeneratePrimer() {
		return true;
	}

	public abstract boolean make2bit(PrimerToFind newPrimerPoss, int possNo, int nPoss);
}
