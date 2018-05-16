package com.primer3.masker;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class maskerMain {
	static String pr_programme_name = "primer3_masker";
	/* maximum number of variables allowed in the probability formula */
	static int MAX_VARIABLES = 100;

	/* flags for choosing between hard and soft masking */
	static boolean HARD_MASKING = false;
	static boolean SOFT_MASKING = true;
	
	public static void main(String[] args){
		
		
//		int idx;
		// char *end;
		boolean debug = false;
		
		/* sequence file name specified by the user, otherwise STDIN will be used */
		String sequence_file_name = null;
		String lists_file_name = null;
		String kmer_lists_path = null;
		
		/* data structure for all k-mer lists used for masking */
		int nlists = 0;
		int nlist_parameters = 0;
		int npos = 0;
		int[] list_pos = new int[MAX_VARIABLES], list_components = new int[MAX_VARIABLES];
		
		MaskerParameters mp = new MaskerParameters();
		parameters_builder pbuilder = new parameters_builder();
		input_sequence input_seq = null;
		output_sequence output_seq = new output_sequence();
		StringBuilder parse_err = new StringBuilder();
		StringBuilder warnings = new StringBuilder();
		
		

		/* fill mp with default parameters */
		mp.mdir = masker.DEFAULT_MASKING_DIRECTION;
		mp.failure_rate = masker.DEFAULT_FAILURE_RATE;
		mp.abs_cutoff = masker.DEFAULT_ABS_CUTOFF;
		mp.nucl_masked_in_5p_direction = masker.DEFAULT_M5P;
		mp.nucl_masked_in_3p_direction = masker.DEFAULT_M3P;
		mp.print_sequence = masker.PRINT_SEQUENCE;
		mp.do_soft_masking = HARD_MASKING;
		mp.masking_char = masker.DEFAULT_MASK_CHAR;
		mp.list_prefix = masker.DEFAULT_LIST_FILE_PREFIX;
		kmer_lists_path = "../kmer_lists/";
		
		
		
		
		CommandLine commandLine = setArgs(args);
		
		if(commandLine != null)
		{
		
			
			// parse parameters here
			
			List<String[]> listOfParameters = new ArrayList<String[]>();
			listOfParameters.add("/data/softwares/primer3-primer3/kmer_lists/homo_sapiens_16.list".split(" "));
			try
			{
				
				
				
				input_seq = input_sequence.create_input_sequence_from_file_name (sequence_file_name);
				if (lists_file_name != null) {
					/* if lists are given in a text file */
					mp.fp = formula_parameters.read_formula_parameters_from_file (lists_file_name, nlist_parameters, pbuilder, mp.formula_intercept);
//					nlists = pbuilder.nfp;
				}
				
				if (listOfParameters.size() > 0) {
					/* if lists are given by commandline arguments (can be added to the ones given in a file) */
					pbuilder.fp_array = mp.fp;

					for (String[] values : listOfParameters) {
						
						
						if (values.length == 1) {
							try
							{
								mp.formula_intercept =  Double.parseDouble(values[0]);
								continue;
							}
							catch(Exception ex )
							{
								// do nothing it is filename
							}
						}
						pbuilder.add_variable_to_formula_parameters (values );
					}		
//					nlists = pbuilder.nfp;
					mp.fp = pbuilder.fp_array;
					if(mp.fp != null)
						nlists = mp.fp.size();

				} else if (mp.fp == null && lists_file_name == null) {
					/* if there are no lists specified use the default formula */
					mp.fp = formula_parameters.create_default_formula_parameters (mp.list_prefix, kmer_lists_path);
					mp.formula_intercept = masker.DEFAULT_INTERCEPT;
					nlists =  masker.DEFAULT_NLISTS;
					nlist_parameters =  masker.DEFAULT_NLIST_PARAMETERS;
				} 

				mp.nlists = nlists;
				
				print_parameters(mp);
				
				masker.read_and_mask_sequence(input_seq,output_seq,mp,debug);
				
				
				
			}
			catch(Exception ex)
			{
				ex.printStackTrace();
			}

			
		}
	}
	
	static void print_parameters (MaskerParameters mp)
	{
		int i;
		System.err.format("Current masker parameters:\n");
		System.err.format( "    masking_direction = %s\n", mp.mdir == masking_direction.both_on_same ? "both" : (mp.mdir == masking_direction.fwd ? "fwd" : "rev"));
		System.err.format( "    failure_rate = %f\n", mp.failure_rate);
		System.err.format( "    abs_cutoff = %d\n", mp.abs_cutoff);
		System.err.format( "    m5p = %d\n", mp.nucl_masked_in_5p_direction);
		System.err.format( "    m3p = %d\n", mp.nucl_masked_in_3p_direction);
		System.err.format( "    print_sequence = %b\n", mp.print_sequence);
		System.err.format( "    do_soft_masking = %b\n", mp.do_soft_masking);
		System.err.format( "    masking_char = %c\n", mp.masking_char);
		System.err.format( "    nlists = %d\n", mp.nlists);
		System.err.format( "    intercept = %f\n", mp.formula_intercept);
		for (i = 0; i < mp.nlists; i++) {
			System.err.format( "    LIST: %d\n", i);
			formula_parameters fp = mp.fp.get(i);
			System.err.format( "        name = %s\n", fp.list_file_name);
			System.err.format( "        oligo_length = %d\n", fp.oligo_length);
			System.err.format( "        words_in_list = %d\n", fp.words_in_list);
			System.err.format( "        mm0 %f mm1 %f mm2 %f mm0_2 %f mm1_2 %f mm2_2 %f\n", fp.mm0, fp.mm1, fp.mm2, fp.mm0_2, fp.mm1_2, fp.mm2_2);
		}
		return;
		
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
				formatter.printHelp( pr_programme_name, options );
		    }
		    
		    return line;		    
		}
		catch( ParseException exp ) {
			// oops, something went wrong
		    System.err.println( "Parsing failed.  Reason: " + exp.getMessage() );
			formatter.printHelp( pr_programme_name, options );
		}
		
		return null;
	}
}
