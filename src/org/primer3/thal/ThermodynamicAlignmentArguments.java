package org.primer3.thal;


public class ThermodynamicAlignmentArguments {
	int debug; /* if non zero, print debugging info to stderr */

	/* one of the
    1 THAL_ANY, (by default)
    2 THAL_END1,
    3 THAL_END2,
    4 THAL_HAIRPIN */
	private ThermodynamicAlignmentType alignmentType; 

	private int maxLoop;  /* maximum size of loop to consider; longer than 30 bp are not allowed */
	private double monovalentConc; /* concentration of monovalent cations */
	private double divalentConc; /* concentration of divalent cations */
	private double dntpConc; /* concentration of dNTP-s */
	private double dnaConc; /* concentration of oligonucleotides */
	private double temperature; /* temperature from which hairpin structures will be calculated */
	private boolean tempOnly; /* if non zero, print only temperature to stderr */
	private int calcDimer; /* if non zero, dimer structure is calculated */


	/**
	 * @return the alignmentType
	 */
	public ThermodynamicAlignmentType getAlignmentType() {
		return alignmentType;
	}

	/**
	 * @return the calcDimer
	 */
	public int getCalcDimer() {
		return calcDimer;
	}

	/**
	 * @return the divalentConc
	 */
	public double getDivalentConc() {
		return divalentConc;
	}

	/**
	 * @return the dnaConc
	 */
	public double getDnaConc() {
		return dnaConc;
	}

	/**
	 * @return the dntpConc
	 */
	public double getDntpConc() {
		return dntpConc;
	}

	/**
	 * @return the maxLoop
	 */
	public int getMaxLoop() {
		return maxLoop;
	}

	/**
	 * @return the monovalentConc
	 */
	public double getMonovalentConc() {
		return monovalentConc;
	}

	/**
	 * @return the temperature
	 */
	public double getTemperature() {
		return temperature;
	}

	/**
	 * @return the tempOnly
	 */
	public boolean isTempOnly() {
		return tempOnly;
	}

	/**
	 * @param alignmentType the alignmentType to set
	 */
	public void setAlignmentType(ThermodynamicAlignmentType alignmentType) {
		this.alignmentType = alignmentType;
	}

	/**
	 * @param calcDimer the calcDimer to set
	 */
	public void setCalcDimer(int calcDimer) {
		this.calcDimer = calcDimer;
	}

	/**
	 * @param divalentConc the divalentConc to set
	 */
	public void setDivalentConc(double divalentConc) {
		this.divalentConc = divalentConc;
	}

	/**
	 * @param dnaConc the dnaConc to set
	 */
	public void setDnaConc(double dnaConc) {
		this.dnaConc = dnaConc;
	}

	/**
	 * @param dntpConc the dntpConc to set
	 */
	public void setDntpConc(double dntpConc) {
		this.dntpConc = dntpConc;
	}

	/**
	 * @param maxLoop the maxLoop to set
	 */
	public void setMaxLoop(int maxLoop) {
		this.maxLoop = maxLoop;
	}

	/**
	 * @param monovalentConc the monovalentConc to set
	 */
	public void setMonovalentConc(double monovalentConc) {
		this.monovalentConc = monovalentConc;
	}

	/**
	 * @param temperature the temperature to set
	 */
	public void setTemperature(double temperature) {
		this.temperature = temperature;
	}

	/**
	 * @param tempOnly the tempOnly to set
	 */
	public void setTempOnly(boolean tempOnly) {
		this.tempOnly = tempOnly;
	}

	public void setThAlDefaultArgs() {
		this.debug = 0;
		this.alignmentType = ThermodynamicAlignmentType.thal_any; /* thal_alignment_type THAL_ANY */
		this.maxLoop = ThermodynamicAlignment.MAX_LOOP;
		this.monovalentConc = 50; /* mM */
		this.divalentConc = 0.0; /* mM */
		this.dntpConc = 0.8; /* mM */
		this.dnaConc = 50; /* nM */
		this.temperature = ThermodynamicAlignment.TEMP_KELVIN; /* Kelvin */
		this.tempOnly = true; /* return only melting temperature of predicted structure */
		this.calcDimer = 1; /* by default dimer structure is calculated */
	}

	/**
	 * Set default args for oligo
	 */
	public void setThAlOligoDefaultArgs(){
		this.debug = 0;
		this.alignmentType = ThermodynamicAlignmentType.thal_any; /* thal_alignment_type THAL_ANY */
		this.maxLoop = ThermodynamicAlignment.MAX_LOOP;
		this.monovalentConc = 50; /* mM */
		this.divalentConc = 0.0; /* mM */
		this.dntpConc = 0.0; /* mM */
		this.dnaConc = 50; /* nM */
		this.temperature = ThermodynamicAlignment.TEMP_KELVIN; /* Kelvin */
		this.tempOnly = true; /* return only melting temperature of predicted structure */
		this.calcDimer = 1; /* by default dimer structure is calculated */
	}


}
