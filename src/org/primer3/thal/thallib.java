package org.primer3.thal;

import java.io.IOException;
import java.io.InputStream;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.io.IOUtils;
import org.primer3.sequence.Sequence;


public class thallib {
	
	
	private static final boolean debug = false;
	static double[][] atpS = new double[5][5]; /* AT penalty */
	static double[][] atpH = new double[5][5]; /* AT penalty */
	
	static double[][][] dangleEntropies3 = new double[5][5][5]; /* thermodynamic paramteres for 3' dangling ends */
	static double[][][] dangleEnthalpies3 = new double[5][5][5]; /* ther params for 3' dangling ends */
	static double[][][] dangleEntropies5 = new double[5][5][5];  /* ther params for 5' dangling ends */
	static double[][][] dangleEnthalpies5 = new double[5][5][5]; /* ther params for 5' dangling ends */
	static double[][][][] stackEntropies = new double[5][5][5][5]; /* ther params for perfect match pairs */
	static double[][][][] stackEnthalpies = new double[5][5][5][5]; /* ther params for perfect match pairs */
	static double[][][][] stackint2Entropies = new double[5][5][5][5]; /*ther params for perfect match and internal mm */
	static double[][][][] stackint2Enthalpies = new double[5][5][5][5]; /* ther params for perfect match and internal mm*/
	static double[] interiorLoopEntropies = new double[30]; /* interior loop params according to length of the loop */
	static double[] bulgeLoopEntropies = new double[30]; /* bulge loop params according to length of the loop */
	static double[] hairpinLoopEntropies = new double[30]; /* hairpin loop params accordint to length of the loop */
	static double[] interiorLoopEnthalpies = new double[30]; /* same as interiorLoopEntropies but values of entropy */
	static double[] bulgeLoopEnthalpies = new double[30]; /* same as bulgeLoopEntropies but values of entropy */
	static double[] hairpinLoopEnthalpies = new double[30]; /* same as hairpinLoopEntropies but values of entropy */
	static double[][][][] tstackEntropies = new double[5][5][5][5]; /* ther params for terminal mismatches */
	static double[][][][] tstackEnthalpies = new double[5][5][5][5]; /* ther params for terminal mismatches */
	static double[][][][] tstack2Entropies = new double[5][5][5][5]; /* ther params for internal terminal mismatches */
	static double[][][][] tstack2Enthalpies = new double[5][5][5][5]; /* ther params for internal terminal mismatches */

	
	static HashMap<Integer,Double> triloopEntropies = new HashMap<Integer, Double>(); /* ther penalties for given triloop seq-s */
	static HashMap<Integer,Double> triloopEnthalpies = new HashMap<Integer, Double>(); /* ther penalties for given triloop seq-s */
	static HashMap<Integer,Double> tetraloopEntropies = new HashMap<Integer, Double>(); /* ther penalties for given tetraloop seq-s */
	static HashMap<Integer,Double> tetraloopEnthalpies = new HashMap<Integer, Double>(); /* ther penalties for given tetraloop seq-s */
	
//	
//	static HashMap<Integer,triloop> triloopEntropies = null; /* ther penalties for given triloop seq-s */
//	static HashMap<Integer,triloop> triloopEnthalpies = null; /* ther penalties for given triloop seq-s */
//	static HashMap<Integer,tetraloop> tetraloopEntropies = null; /* ther penalties for given tetraloop seq-s */
//	static HashMap<Integer,tetraloop> tetraloopEnthalpies = null; /* ther penalties for given tetraloop seq-s */

	/* Read the thermodynamic values (parameters) from the parameter files
	   in the directory specified by 'path'.  Return 0 on success and -1
	   on error. The thermodynamic values are stored in multiple static
	   variables. */
	public static int  get_thermodynamic_values(String path, thal_results o)
	{
		if(path == null || path.isEmpty())
		{
			
		}
		
		getStack();
		getStackint2();
		getDangle();
		getLoop();
		getTstack();
		getTstack2();
		getTriloop();
		getTetraloop();
		tableStartATS();
		tableStartATH();
		/* getting the AT-penalties */

		return 0;
	}

	
	
	public static int  get_thermodynamic_values()
	{
		return get_thermodynamic_values(null,null);
	}
	static void tableStartATS()
	{
		for (int i = 0; i < 5; ++i)
			for (int j = 0; j < 5; ++j)
				atpS[i][j] = 0.00000000001;
		atpS[0][3] = atpS[3][0] = thal.AT_S;
	}


	static void 
	tableStartATH()
	{
		//
		for (int i = 0; i < 5; ++i)
			for (int j = 0; j < 5; ++j)
				atpH[i][j] = 0.0;
		atpH[0][3] = atpH[3][0] = thal.AT_H;
	}
	
	static void  getStack()
	{
		
		List<String> sFileList = loadResource("stack.ds");
		List<String> hFileList = loadResource("stack.dh");
		
		Iterator<String> sFile =  sFileList.iterator();
		Iterator<String> hFile =  hFileList.iterator();

		// stackEntropies, stackEnthalpies
		for (int i = 0; i < 5; ++i) {
			for (int ii = 0; ii < 5; ++ii) {
				for (int j = 0; j < 5; ++j) {
					for (int jj = 0; jj < 5; ++jj) {
						if (i == 4 || j == 4 || ii == 4 || jj == 4) {
							stackEntropies[i][ii][j][jj] = -1.0;
							stackEnthalpies[i][ii][j][jj] = Double.POSITIVE_INFINITY;
						} else {
							stackEntropies[i][ii][j][jj]  = readDouble(sFile);
							stackEnthalpies[i][ii][j][jj] = readDouble(hFile);
							if (!Double.isFinite(stackEntropies[i][ii][j][jj]) || !Double.isFinite(stackEnthalpies[i][ii][j][jj])) {
								stackEntropies[i][ii][j][jj] = -1.0;
								stackEnthalpies[i][ii][j][jj] = Double.POSITIVE_INFINITY;
							}
						}
					}
				}
			}
		}
	}
	
	static void getStackint2()
	{
		List<String> sFileList = loadResource("stackmm.ds");
		List<String> hFileList = loadResource("stackmm.dh");
		
		Iterator<String> sFile =  sFileList.iterator();
		Iterator<String> hFile =  hFileList.iterator();
		for (int i = 0; i < 5; ++i) {
			for (int ii = 0; ii < 5; ++ii) {
				for (int j = 0; j < 5; ++j) {
					for (int jj = 0; jj < 5; ++jj) {
						if (i == 4 || j == 4 || ii == 4 || jj == 4) {
							stackint2Entropies[i][ii][j][jj] = -1.0;
							stackint2Enthalpies[i][ii][j][jj] = Double.POSITIVE_INFINITY;
						} else {
							stackint2Entropies[i][ii][j][jj] = readDouble(sFile);
							stackint2Enthalpies[i][ii][j][jj] = readDouble(hFile);
							if (!Double.isFinite(stackint2Entropies[i][ii][j][jj]) || !Double.isFinite(stackint2Enthalpies[i][ii][j][jj])) {
								stackint2Entropies[i][ii][j][jj] = -1.0;
								stackint2Enthalpies[i][ii][j][jj] = Double.POSITIVE_INFINITY;
							}
						}
					}
				}
			}
		}
	}
	
	
	static void getDangle()
	{
		List<String> sFileList = loadResource("dangle.ds");
		List<String> hFileList = loadResource("dangle.dh");
		
		Iterator<String> sFile =  sFileList.iterator();
		Iterator<String> hFile =  hFileList.iterator();
		

		for (int i = 0; i < 5; ++i)
		{
			for (int j = 0; j < 5; ++j) 
			{
				for (int k = 0; k < 5; ++k) 
				{
					if (i == 4 || j == 4) {
						dangleEntropies3[i][k][j] = -1.0;
						dangleEnthalpies3[i][k][j] = Double.POSITIVE_INFINITY;
					} else if (k == 4) {
						dangleEntropies3[i][k][j] = -1.0;
						dangleEnthalpies3[i][k][j] = Double.POSITIVE_INFINITY;
					} else {
						dangleEntropies3[i][k][j] = readDouble(sFile);
						dangleEnthalpies3[i][k][j] = readDouble(hFile);
						if(!Double.isFinite(dangleEntropies3[i][k][j]) || !Double.isFinite(dangleEnthalpies3[i][k][j])) {
							dangleEntropies3[i][k][j] = -1.0;
							dangleEnthalpies3[i][k][j] = Double.POSITIVE_INFINITY;	     
						}
					}
				}
			}
		}
		
		
		for (int i = 0; i < 5; ++i)
		{
			for (int j = 0; j < 5; ++j) {
				for (int k = 0; k < 5; ++k) {
					if (i == 4 || j == 4) {
						dangleEntropies5[i][j][k] = -1.0;
						dangleEnthalpies5[i][j][k] = Double.POSITIVE_INFINITY;
					} else if (k == 4) {
						dangleEntropies5[i][j][k] = -1.0;
						dangleEnthalpies5[i][j][k] = Double.POSITIVE_INFINITY;
					} else {
						dangleEntropies5[i][j][k] = readDouble(sFile);
						dangleEnthalpies5[i][j][k] = readDouble(hFile);
						if(!Double.isFinite(dangleEntropies5[i][j][k]) || !Double.isFinite(dangleEnthalpies5[i][j][k])) {
							dangleEntropies5[i][j][k] = -1.0;
							dangleEnthalpies5[i][j][k] = Double.POSITIVE_INFINITY;
						}
					}
				}
			}
		}
	}
	static void getLoop() {
		List<String> sFileList = loadResource("loops.ds");
		List<String> hFileList = loadResource("loops.dh");
		
		Iterator<String> sFile =  sFileList.iterator();
		Iterator<String> hFile =  hFileList.iterator();
		for (int k = 0; k < 30; ++k) {
			double[] loopEntropies		= readLoop(sFile);
			interiorLoopEntropies[k]	= loopEntropies[0]; 
			bulgeLoopEntropies[k]    	= loopEntropies[1];
			hairpinLoopEntropies[k]		= loopEntropies[2];
			double[] loopEnthalpies 	= readLoop(hFile );
			interiorLoopEnthalpies[k] 	= loopEnthalpies[0]; 
			bulgeLoopEnthalpies[k]		= loopEnthalpies[1];
			hairpinLoopEnthalpies[k]	= loopEnthalpies[2];
		}
	}
	static void  getTstack()
	{
		List<String> sFileList = loadResource("tstack_tm_inf.ds");
		List<String> hFileList = loadResource("tstack.dh");
		
		Iterator<String> sFile =  sFileList.iterator();
		Iterator<String> hFile =  hFileList.iterator();
		for (int i1 = 0; i1 < 5; ++i1)
		{
			for (int i2 = 0; i2 < 5; ++i2)
			{
				for (int j1 = 0; j1 < 5; ++j1)
				{
					for (int j2 = 0; j2 < 5; ++j2)
					{
						if (i1 == 4 || j1 == 4) {
							tstackEnthalpies[i1][i2][j1][j2] = Double.POSITIVE_INFINITY;
							tstackEntropies[i1][i2][j1][j2] = -1.0;
						} else if (i2 == 4 || j2 == 4) {
							tstackEntropies[i1][i2][j1][j2] = 0.00000000001;
							tstackEnthalpies[i1][i2][j1][j2] = 0.0;
						} else {
							tstackEntropies[i1][i2][j1][j2] = readDouble(sFile);
							tstackEnthalpies[i1][i2][j1][j2] = readDouble(hFile);
							if (!Double.isFinite(tstackEntropies[i1][i2][j1][j2]) || !Double.isFinite(tstackEnthalpies[i1][i2][j1][j2])) {
								tstackEntropies[i1][i2][j1][j2] = -1.0;
								tstackEnthalpies[i1][i2][j1][j2] = Double.POSITIVE_INFINITY;
							}
						}
					}
				}
			}
		}
	}
	
	static void  getTstack2()
	{
		List<String> sFileList = loadResource("tstack2.ds");
		List<String> hFileList = loadResource("tstack2.dh");
		
		Iterator<String> sFile =  sFileList.iterator();
		Iterator<String> hFile =  hFileList.iterator();
		for (int i1 = 0; i1 < 5; ++i1)
		{
			for (int i2 = 0; i2 < 5; ++i2)
			{
				for (int j1 = 0; j1 < 5; ++j1)
				{
					for (int j2 = 0; j2 < 5; ++j2)
						if (i1 == 4 || j1 == 4)  {
							tstack2Enthalpies[i1][i2][j1][j2] = Double.POSITIVE_INFINITY;
							tstack2Entropies[i1][i2][j1][j2] = -1.0;
						} else if (i2 == 4 || j2 == 4) {
							tstack2Entropies[i1][i2][j1][j2] = 0.00000000001;
							tstack2Enthalpies[i1][i2][j1][j2] = 0.0;
						} else {
							tstack2Entropies[i1][i2][j1][j2] = readDouble(sFile);
							tstack2Enthalpies[i1][i2][j1][j2] = readDouble(hFile);
							if (!Double.isFinite(tstack2Entropies[i1][i2][j1][j2]) || !Double.isFinite(tstack2Enthalpies[i1][i2][j1][j2])) {
								tstack2Entropies[i1][i2][j1][j2] = -1.0;
								tstack2Enthalpies[i1][i2][j1][j2] = Double.POSITIVE_INFINITY;
							}
						}
				}
			}
		}
	}
	
	static void getTriloop()
	{
		List<String> sFileList = loadResource("triloop.ds");
		List<String> hFileList = loadResource("triloop.dh");
		
//		Iterator<String> sFile =  sFileList.iterator();
//		Iterator<String> hFile =  hFileList.iterator();
		
		for(String triLoopLine : sFileList )
		{
			String[] tokens = triLoopLine.split("\t");
			int hashKey = hashLoop(tokens[0].toCharArray());
			double value = Double.parseDouble(tokens[1]);
			triloopEntropies.put(hashKey, value);
		}
		for(String triLoopLine : hFileList )
		{
			String[] tokens = triLoopLine.split("\t");
			int hashKey = hashLoop(tokens[0].toCharArray());
			double value = Double.parseDouble(tokens[1]);
			triloopEnthalpies.put(hashKey, value);
		}
		
	}
	static void getTetraloop()
	{
		List<String> sFileList = loadResource("tetraloop.ds");
		List<String> hFileList = loadResource("tetraloop.dh");
		
//		Iterator<String> sFile =  sFileList.iterator();
//		Iterator<String> hFile =  hFileList.iterator();
		
		for(String triLoopLine : sFileList )
		{
			String[] tokens = triLoopLine.split("\t");
			int hashKey = hashLoop(tokens[0].toCharArray());
			double value = Double.parseDouble(tokens[1]);
			tetraloopEntropies.put(hashKey, value);
		}
		for(String triLoopLine : hFileList )
		{
			String[] tokens = triLoopLine.split("\t");
			int hashKey = hashLoop(tokens[0].toCharArray());
			double value = Double.parseDouble(tokens[1]);
			tetraloopEnthalpies.put(hashKey, value);
		}
		
	}
	public static int hashLoop(char[] loop)
	{
		int hash = 0;
		int len = loop.length;
		int base = 1;
		for(int i = len-1 ; i >= 0; i--)
		{
			hash += base * (str2int(loop[i]) +1);
			base *= 10;
		}
		
		return hash;
	}
	public static int hashLoop(int[] loop)
	{
		
		
		int hash = 0;
		int len = loop.length;
		int base = 1;
		for(int i = len-1 ; i >= 0; i--)
		{
			hash += base * (loop[i]+1);
			base *= 10;
		}
		
		return hash;
	}
	
	/**
	 *  converts DNA sequence to int; 0-A, 1-C, 2-G, 3-T, 4-whatever 
	 */
	public static int str2int(char c) {
		switch (c) {
			case 'A': case '0':
				return 0;
			case 'C': case '1':
				return 1;
			case 'G': case '2':
				return 2;
			case 'T': case '3':
				return 3;
		}
		return 4;
	}


	private static double[] readLoop(Iterator<String> sFile) {

		String line = sFile.next();
		String[] lineTokens =  line.trim().split("\t");
		double[] parsedList = new double[3];
		
		parsedList[0] = Double.parseDouble(lineTokens[1].replace("inf", "Infinity"));
		parsedList[1] = Double.parseDouble(lineTokens[2].replace("inf", "Infinity"));
		parsedList[2] = Double.parseDouble(lineTokens[3].replace("inf", "Infinity"));

		return parsedList;
	}


	private static double  readDouble(Iterator<String> sFile) {
		String line = sFile.next();
		line = line.replace("inf", "Infinity").trim();
		return Double.parseDouble(line);
	}


	static  List<String> loadResource(String res_name) {
		String resourceFolder = "resources/";
		ClassLoader classLoader = thallib.class .getClassLoader();
		String resourceName = resourceFolder + res_name;		
		InputStream file = classLoader.getResourceAsStream(resourceName);
		List<String> fileText = null;
		try {
			fileText = IOUtils.readLines(file); 
		} catch (IOException e) {
			e.printStackTrace();
		}
		return fileText;
	}
	
	
	
	public static void main(String[] args)
	{
		
		String[] testCases = new String[]{
				"-s2 AACCCGCTTGCTAGCAAAAAAAAAACCCGCTTGCTAGCAAAAAAAAAACCCGCTTGCTAGCAAAAAAAAAACCCGCTTGCTAGCAAAAAAAAAACCCGCTTGCTAGCAA -a HAIRPIN",
				"-s1 ccgcagtaagctgcgg -a HAIRPIN",
				"-s1 AAAACCCGCTTTGCTAGCTACG -s2 AACCCGCTTGCTAGCAAAAAAAA -a ANY -t 50",
				"-s1 ccgcagtaagctgcgg -s2 ccgcagtaagctgcgg -a ANY -maxloop 1",
				"-dv 1.5 -n 0.6 -s1 TGGTGGCAAAGTCGACAGAG -s2 TGCACTACCTGAGGCTTCAC -a END2",
				"-dv 1.5 -n 0.6 -s1 AACTTCAGGATCCAGTGGGC  -a HAIRPIN"
		}; 
		
		
		args = testCases[5].split(" ");
		
		String s1 = null, s2 = null;
		
		CommandLine commandLine = setArgs(args);
		if(commandLine != null)
		{
			thal_args a = new thal_args();
			a.set_thal_default_args();
			a.temponly=0;
			if(commandLine.hasOption("mv"))
			{
				a.mv = Double.parseDouble(commandLine.getOptionValue("mv"));
			}
			if(commandLine.hasOption("dv"))
			{
				a.dv = Double.parseDouble(commandLine.getOptionValue("dv"));
				if(a.dv < 0 )
				{
					System.err.println("divalent_conc can not be less than 0");
					System.exit(-1);
				}
			}
			if(commandLine.hasOption("n"))  /* concentration of dNTPs */
			{
				a.dntp = Double.parseDouble(commandLine.getOptionValue("n"));
				if(a.dntp < 0 )
				{
					System.err.println("(n) concentration of dNTPs can not be less than 0");
					System.exit(-1);
				}
			}
			if(commandLine.hasOption("maxloop"))  /* maximum size of loop calculated; 
						      this value can not be larger than 30 */
			{
				a.maxLoop = Integer.parseInt(commandLine.getOptionValue("maxloop"));
				 if(a.maxLoop > thal.MAX_LOOP ) {
					 a.maxLoop = thal.MAX_LOOP;
				 }  else if(a.maxLoop < thal.MIN_LOOP) {	 
					 a.maxLoop = thal.MIN_LOOP;
				 }
			}
			if(commandLine.hasOption("a"))
			{
				String type = commandLine.getOptionValue("a");
				if(type != null )
				{
					if(type.equals("ANY"))
					{
						a.type =thal_alignment_type.thal_any;
					} else if(type.equals("END1")) { 
						a.type =thal_alignment_type.thal_end1;
					} else if(type.equals("END2")) {
						a.type =thal_alignment_type.thal_end2;
					} else if(type.equals("HAIRPIN")) {
						a.type =thal_alignment_type.thal_hairpin;
						 a.dimer = 0;
					}
				}	
			}
			if(commandLine.hasOption("d")) // dna conc
			{
				a.dna_conc = Double.parseDouble(commandLine.getOptionValue("d"));
				if(a.dna_conc <= 0 )
				{
					System.err.println("dna conc can not be less than or eqaul 0");
					System.exit(-1);
				}
					
			}
			if(commandLine.hasOption("t")) // temperature at which sec str are calculated in C
			{
				a.temp = Double.parseDouble(commandLine.getOptionValue("t")) + thal.ABSOLUTE_ZERO;	
			}
			if(commandLine.hasOption("s1")) //s1
			{
				s1 = commandLine.getOptionValue("s1");	
			}
			if(commandLine.hasOption("s2")) //s1
			{
				s2 = commandLine.getOptionValue("s2");	
			}
			
			
			 /* check the input correctness */
			if(a.dimer == 1 && (s2==null || s1==null)) { /* if user wants to calculate structure 
									       of dimer then two sequences must be defined*/
				System.err.println("two sequences must be defined");
				System.exit(-1);
			   }
			   if(a.dimer==0 && (s2==null && s1==null)) { /* if user wants to calculate structure
									       of monomer then only one sequence must be defined */
					System.err.println("two sequences must be defined");
					System.exit(-1);
			   }			
			
			
			thal_results o = null;
			int thal_trace = 1;
			// run 
			get_thermodynamic_values(null,null);
			try {
				if(a.dimer == 0 && s1 != null)
				{
					Sequence oligo1 = new Sequence(s1.toCharArray());
					thal thal_1 = new thal(oligo1, oligo1, a);
					o = thal_1.calc_thal();
					if (thal_trace != 0) {
						System.out.format(  "thal, thal_args, type=%s maxLoop=%d mv=%f dv=%f dntp=%f dna_conc=%f, temp=%f, temponly=%d dimer=%d\n",
								a.type, a.maxLoop, a.mv, a.dv, a.dntp, a.dna_conc, 
								a.temp, a.temponly, a.dimer);
						System.out.format(  "thal: s1=%s s2=%s temp=%f msg=%s end1=%d end2=%d\n", 
							s1, s2, o.temp, o.msg, o.align_end_1, o.align_end_2);
					}
					 
				} else if (a.dimer == 0 && s1 == null && s2 != null)
				{
					Sequence oligo2 = new Sequence(s2.toCharArray());
					thal thal_1 = new thal(oligo2, oligo2, a);
					o = thal_1.calc_thal();
				} else {
					Sequence oligo1 = new Sequence(s1.toCharArray());
					Sequence oligo2 = new Sequence(s2.toCharArray());
					thal thal_1 = new thal(oligo1, oligo2, a);
					o = thal_1.calc_thal();
					
					if (thal_trace != 0) {
						System.out.format(  "thal, thal_args, type=%s maxLoop=%d mv=%f dv=%f dntp=%f dna_conc=%f, temp=%f, temponly=%d dimer=%d\n",
								a.type, a.maxLoop, a.mv, a.dv, a.dntp, a.dna_conc, 
								a.temp, a.temponly, a.dimer);
						System.out.format(  "thal: s1=%s s2=%s temp=%f msg=%s end1=%d end2=%d\n", 
							s1, s2, o.temp, o.msg, o.align_end_1, o.align_end_2);
					}
					
					
				}
			} catch (ThermodynamicAlignmentException e) {
				
				
				System.err.println("Error: " + e.getMessage());
				
				if(debug)
				{
					e.printStackTrace();
				}
			}
		}
	}
	
	static String string(char[] s2) {
		if(s2 == null)
			return "";
		return new String(s2);
	}
	static CommandLine setArgs(String[] args)
	{
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		Option  argOption   = Option.builder("mv").argName( "monovalent_conc" )
                .hasArg()
                .desc("concentration of monovalent cations in mM, by default 50 mM" )
                
                .build();
		options.addOption(argOption);
		
		argOption   = Option.builder("dv").argName( "divalent_conc" )
                .hasArg()
                .desc("concentration of divalent cations in mM, by default 0 mM" )
                .build();
		options.addOption(argOption);
		
		argOption   = Option.builder("n").argName( "dNTP_conc" )
                .hasArg()
                .desc("concentration of deoxynycleotide triphosphate in mM, by default 0 mM" )
                .build();
		options.addOption(argOption);
		
		argOption   = Option.builder("d").argName( "dna_conc" )
                .hasArg()
                .desc("concentration of DNA strands in nM, by default 50 nM")
                .build();
		options.addOption(argOption);
		
		argOption   = Option.builder("a").argName( "mode" )
                .hasArg()
                .desc("alignment type, END1, END2, ANY and HAIRPIN, by default ANY (when duplex)")
                .build();
		options.addOption(argOption);
		
		argOption   = Option.builder("t").argName( "temp" )
                .hasArg()
                .desc("temperature at which duplex is calculated, by default 37C")
                .build();
		options.addOption(argOption);
		
		argOption   = Option.builder("r").argName( "TempOnly" )
                .desc("causes the alignment NOT to be displayed on stderr, _only_ Tm is printed" )
                .build();
		options.addOption(argOption);
		
		argOption   = Option.builder("maxloop").argName( "maxloop" )
				.hasArg()
                .desc("the maximum size of secondary structures loops.\n" + 
                		"                       Default is 30 (this is maximum allowed length, currently)." )
                .build();
		options.addOption(argOption);
		
		argOption   = Option.builder("s1").argName( "DNA_oligomer1" )
                .hasArg()
                .desc("DNA_oligomer 1")
                .build();
		options.addOption(argOption);
		argOption   = Option.builder("s2").argName( "DNA_oligomer2" )
                .hasArg()
                .desc("DNA_oligomer2" )
                .build();
		
		options.addOption(argOption);
		
		Option helpOption = new Option( "help", "Print Usage" );

		options.addOption(helpOption);

		// automatically generate the help statement
		HelpFormatter formatter = new HelpFormatter();
		formatter.setOptionComparator(null);
		try {
			// parse the command line arguments
		    CommandLine line = parser.parse( options, args );
		    
		    if(line.getOptions().length == 0 )
		    {
				formatter.printHelp( "ntthal", options );
		    }
		    
		    return line;		    
		}
		catch( ParseException exp ) {
			// oops, something went wrong
		    System.err.println( "Parsing failed.  Reason: " + exp.getMessage() );
			formatter.printHelp( "ntthal", options );
		}
		
		return null;
	}
	
	
	// 5
	class triloop {
		  char[] loop;
		  double value; 
	};

	// 6
	class tetraloop {
		  char[] loop;
		  double value; 
	};
}
