package org.primer3;

import java.util.HashMap;

public class Primer3Manual {

	public static final String PRIMER_PAIR_MAX_DIFF_TM = "PRIMER_PAIR_MAX_DIFF_TM";

	public static String PRIMER_PRODUCT_MAX_TM = "PRIMER_PRODUCT_MAX_TM";
	
	static HashMap<String,String> helpTxt;
	static HashMap<String,String> helpLongTxt;

	static {
		helpTxt = new HashMap<String, String>();
		helpLongTxt = new HashMap<String, String>();

		
		helpTxt.put(PRIMER_PRODUCT_MAX_TM, "The maximum allowed melting temperature of the amplicon.");
		helpLongTxt.put(PRIMER_PRODUCT_MAX_TM, "Primer3 calculates product Tm calculated using the formula from Bolton and McCarthy, PNAS 84:1390 (1962) as presented in Sambrook, Fritsch and Maniatis, Molecular Cloning, p 11.46 (1989, CSHL Press).");

		
		
		helpTxt.put(PRIMER_PAIR_MAX_DIFF_TM, "Maximum acceptable (unsigned) difference between the melting temperatures of the forward and reverse primers.");

	}
	
	
	public static String getHelp(String key) {
		return getHelp(key,false);
	}
	public static String getHelp(String key, boolean longVersion) {
		String shortHelp = helpTxt.get(key); 
		if(shortHelp == null)
			shortHelp = "";
		if(longVersion)
		{
			String longHelp = helpLongTxt.get(key);
			if(longHelp != null)
			{
				if(!shortHelp.endsWith("."))
					shortHelp += ".";
				shortHelp +=  " " + longHelp;
			}
		}
		return shortHelp;
	}
	
	
}
