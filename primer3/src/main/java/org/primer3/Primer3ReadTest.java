package org.primer3;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.primer3.libprimer3.LibPrimer3;
import org.primer3.libprimer3.OligoArray;
import org.primer3.libprimer3.OligoStats;
import org.primer3.libprimer3.OligoType;
import org.primer3.libprimer3.P3GlobalSettings;
import org.primer3.libprimer3.P3OutputType;
import org.primer3.libprimer3.P3RetVal;
import org.primer3.libprimer3.P3Task;
import org.primer3.libprimer3.PairArrayT;
import org.primer3.libprimer3.PairStats;
import org.primer3.libprimer3.SeqArgs;
import org.primer3.masker.formula_parameters;
import org.primer3.masker.masker;
import org.primer3.masker.masking_direction;
import org.primer3.p3_seq_lib.seq_lib;
import org.primer3.primer.PrimerPair;
import org.primer3.primer.PrimerRecord;

public class Primer3ReadTest {
	

	static String pr_release = LibPrimer3.libprimer3_release();
	static String pr_program_name = "" ;
	static HelpFormatter formatter;
	static Options options;

	static public void main(String[] args)
	{

		/* Setup the input data structures handlers */
		boolean format_output = false;
		boolean strict_tags = false;
		boolean echo_settings = false;
		int io_version = 4;
		int default_version = 2;
		boolean dump_args = false;	
		int about = 0;
		int compat = 0;
		int invalid_flag = 0;
		String fastaInputFile = null;
		
		if(args.length > 3) {
			fastaInputFile = args[3];
		}

		String p3_settings_path = null;
		String output_path = args[2];
		String error_path = null;
		StringBuilder fatal_parse_err = new StringBuilder();
		StringBuilder nonfatal_parse_err = new StringBuilder();
		StringBuilder warnings = new StringBuilder();



		P3GlobalSettings global_pa;
		SeqArgs sarg;




		P3RetVal retval = null;
		int input_found=0;

		/* Get the program name for correct error messages */
		// pr_release = libprimer3_release();
		pr_program_name = args[0];
		LibPrimer3.p3_set_program_name(pr_program_name);

		if(output_path != null && !output_path.isEmpty())
		{
			try
			{
				PrintStream outStream = outputFile(output_path);
				System.setOut(outStream);
			}
			catch(FileNotFoundException fnEx)
			{
				System.err.println( "Error creating file : " + output_path );
				System.exit(-1);
			}
		}

		if (about == 1) {
			System.out.println(pr_release + "\n");
			System.exit(0);
		}
		if ((io_version == -1) || (invalid_flag == 1) || (default_version == -1)) {
			print_usage();
			System.exit(-1);
		}




		if(args.length > 0)
		{
			String inputFileName = args[1];
			try {
				System.setIn(new FileInputStream(new File(inputFileName)));

			} catch (FileNotFoundException e) {
				System.err.println("Error opening file " + inputFileName);
				System.exit(-1);
			}
		}

		Scanner scan  = new Scanner(System.in);

		/* Allocate the space for global settings and fill in default parameters */
		global_pa = P3GlobalSettings.p3_create_global_settings(default_version);

		if(global_pa == null)
		{
			print_usage();
			System.exit(-1);
		}

		global_pa.setDump(dump_args) ;

		// will assume that the input will come from a fasta file
		if (fastaInputFile != null)
		{
			try {
				seq_lib inputFasta = seq_lib.read_and_create_seq_lib(fastaInputFile, false);
				global_pa.setInputFasta(inputFasta);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				System.exit(-1);
			}
		}


		/* Settings files have to be read in just below, and
	    	the functions need a temporary sarg */
		sarg = new SeqArgs();
//		read_boulder_record_res = new read_boulder_record_results();
		/* Read data from the settings file until a "=" line occurs.  Assign
	     parameter values for primer picking to pa and sa. */
		if(p3_settings_path != null && !p3_settings_path.isEmpty())
		{
			boulder.read_p3_file(p3_settings_path, P3FileType.settings,
					echo_settings && !format_output, 
					strict_tags, 
					global_pa, sarg, 
					fatal_parse_err,
					nonfatal_parse_err, warnings);

			/* Check if any thermodynamical alignment flag was given */
			if ((global_pa.isThermodynamicOligoAlignment() ) || 
					(global_pa.isThermodynamicTemplateAlignment() ))
				read_thermodynamic_parameters();
			/* Check if masking template flag was given */
			if (global_pa.isMaskTemplate())
				validate_kmer_lists_path();


		}

		/* We also need to print out errors here because the loop erases all
	     errors at start. If there are fatal errors, write the proper
	     message and exit */
		if (!fatal_parse_err.toString().isEmpty()) {
			if (format_output) {
				boulder.format_error(sarg.getSequenceName(), fatal_parse_err.toString());
			} else {
				boulder.print_boulder_error(fatal_parse_err.toString());
			}
			System.err.format("%s: %s\n",pr_program_name, fatal_parse_err.toString());
			System.exit(-4);
		}

		/* If there are nonfatal errors, write the proper message
		 * and skip to the end of the loop */
		if (!nonfatal_parse_err.toString().isEmpty()) {
			if (format_output) {
				boulder.format_error( sarg.getSequenceName(), nonfatal_parse_err.toString());
			} else {
				boulder.print_boulder_error(nonfatal_parse_err.toString());
			}
		}
		/* The temporary sarg is not needed any more */
		//	    destroy_seq_args(sarg);
		sarg = null ;

//		read_boulder_record_res = new read_boulder_record_results();

		// for new part in multiplex 
		ArrayList<P3RetVal> lateResult = new ArrayList<P3RetVal>();

		List<SeqArgs> SargsList = new ArrayList<SeqArgs>(); //A LIST FOR ALL THE SEQUENCE ARGUMENTS
		List<P3GlobalSettings> PaList = new ArrayList<P3GlobalSettings>(); //A LIST FOR ALL THE GLOBAL SETTINGS
		try {
			//THIS LOOP READ AND CREATE ALL THE SEQUENCE ARGUMENTS AND GLOBAL SETTINGS THAT WILL BE USED IN RETVAL READ AND PRINT
			while(true)
			{
			
				global_pa = P3GlobalSettings.p3_create_global_settings(default_version);
				
				if(global_pa == null)
				{
					print_usage();
					System.exit(-1);
				}

				global_pa.setDump(dump_args) ;

				// will assume that the input will come from a fasta file
				if (fastaInputFile != null)
				{
					try {
						seq_lib inputFasta = seq_lib.read_and_create_seq_lib(fastaInputFile, false);
						global_pa.setInputFasta(inputFasta);
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
						System.exit(-1);
					}
				}
				
				
				/* Create and initialize a seq_args data structure. sa (seq_args *) is
				 * initialized here because Values are _not_ retained across different
				 * input records. */
				sarg = new SeqArgs();

				/* Reset all errors handlers and the return structure */
				fatal_parse_err = new StringBuilder();
				nonfatal_parse_err= new StringBuilder();
				warnings = new StringBuilder();
				System.gc();
				retval = null;

				/* There were no more boulder records */
				//!format_output on false
				if(! boulder.read_boulder_record(scan,
						strict_tags,
						io_version,
						false,
						P3FileType.all_parameters,
						global_pa,
						sarg,
						fatal_parse_err,
						nonfatal_parse_err,
						warnings))
				{
					break; /* There were no more boulder records */
				}
				//Adding to the list
				SargsList.add(sarg);
				
				//PA other modifications
				if(global_pa.isMaskTemplate()){
					global_pa.setLowercaseMasking(global_pa.isMaskTemplate());
				}
				/* Check if any thermodynamical alignment flag was given and the
			    path to the parameter files changed - we need to reread them */
				if (((global_pa.isThermodynamicOligoAlignment() ) ||
						(global_pa.isThermodynamicTemplateAlignment() ))
						&& (thermodynamic_path_changed ))
					read_thermodynamic_parameters();

				/* Check if template masking flag was given */
				if (global_pa.isMaskTemplate() )
					validate_kmer_lists_path();


				/* Check that we found the thermodynamic parameters in case any thermodynamic flag was set to 1. */
				if (((global_pa.isThermodynamicOligoAlignment() ) ||
						(global_pa.isThermodynamicTemplateAlignment() ))
						&& (thermodynamic_params_path == null)) {
					/* no parameter directory found, error */
					System.out.println("PRIMER_ERROR=thermodynamic approach chosen, but path to thermodynamic parameters not specified\n=\n");
					System.exit(-1);
				}

				/* Check that we found the kmer lists in case masking flag was set to 1. */
				if (global_pa.isMaskTemplate()  && kmer_lists_path == null){
					System.out.println("PRIMER_ERROR=masking template chosen, but path to kmer lists not specified\n=\n");
					System.exit(-1);
				}

				/* Set up some masking parameters */
				/* edited by M. Lepamets */
				if (global_pa.isMaskTemplate() ) {
					global_pa.getMaskingParameters().window_size = masker.DEFAULT_WORD_LEN_2;

					if (global_pa.isPickRightPrimer() ) global_pa.getMaskingParameters().mdir = masking_direction.fwd;
					else if (global_pa.isPickLeftPrimer() ) global_pa.getMaskingParameters().mdir = masking_direction.rev;
					/* Check if masking parameters (k-mer list usage) have changed */
					if (global_pa.isMaskingParametersChanged()) {
						//	            masker.delete_formula_parameters (global_pa.mp.fp, global_pa.mp.nlists);
						global_pa.getMaskingParameters().fp = formula_parameters.create_default_formula_parameters (global_pa.getMaskingParameters().list_prefix, kmer_lists_path);
						global_pa.setMaskingParametersChanged(false);
					}
				}

				input_found = 1;
				if ((global_pa.getPrimerTask() == P3Task.GENERIC)
						&& (global_pa.isPickInternalOligo() )){
					PR_ASSERT(global_pa.isPickInternalOligo());
				}
				
				//Adding to the list
				PaList.add(global_pa);
			} //END SARG_PA WHILE
			
			
			//WE READ ALL THE RETVALS FROM THE FILE 				
			List<P3RetVal> retvalList = read_from_file(args[0], PaList,SargsList);
			
			//WE PRINT AGAIN THE RETVALS
			for(int g = 0; g < SargsList.size();g++){
				
				if(SargsList.get(g).isMultiplex() )
					lateResult.add(retvalList.get(g));

				if (global_pa.isPickAnyway() && format_output) {
					if (SargsList.get(g).getLeftInput() != null) {
						retvalList.get(g).add_must_use_warnings( "Left primer", retvalList.get(g).fwd);
					}
					if (SargsList.get(g).getRightInput() != null) {
						retvalList.get(g).add_must_use_warnings("Right primer", retvalList.get(g).rev);
					}
					if (SargsList.get(g).getInternalInput() != null) {
						retvalList.get(g).add_must_use_warnings("Hybridization probe", retvalList.get(g).intl);
					}
				}
	
	
				if (retvalList.get(g).glob_err.length() == 0 && retvalList.get(g).per_sequence_err.length() == 0  ) {

					if (PaList.get(g).getFileFlag() != 0) {

						LibPrimer3.p3_print_oligo_lists(retvalList.get(g), SargsList.get(g), PaList.get(g),
								retvalList.get(g).per_sequence_err,
								SargsList.get(g).getSequenceName());
					}
				}
				
				//PRINTING AGAIN
				if(! SargsList.get(g).isMultiplex() )
				if (format_output) {
					boulder.print_format_output(io_version, PaList.get(g),
							SargsList.get(g), retvalList.get(g), pr_release,
							PaList.get(g).getExplainFlag());
				} else {
					/* Use boulder output */
					
					retvalList.get(g).print_boulder(io_version, 
							PaList.get(g).getExplainFlag() != 0);
				}


			}// end foreach


		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}



	private static void PR_ASSERT(boolean pick_internal_oligo) {
		// TODO Auto-generated method stub
		assert pick_internal_oligo;
	}



	static boolean   thermodynamic_path_changed = true;/* if this is set to 1, we need to re-read the thermodynamic parameters from new path */
	static String thermodynamic_params_path; /* path to thermodynamic parameter files */
	static String kmer_lists_path = null; /* path to kmer lists files */
	private static void validate_kmer_lists_path() {
		// TODO Auto-generated method stub

	}
	private static void read_thermodynamic_parameters() {

		//		ThAl.get_thermodynamic_values();

		// FIXME :: remove this
		thermodynamic_params_path = "";
		thermodynamic_path_changed = false;
	}
	private static void print_usage() {
		// TODO Auto-generated method stub

	}


	protected static PrintStream outputFile(String name) throws FileNotFoundException {
		return new PrintStream(new BufferedOutputStream(new FileOutputStream(name)), true);
	}
	
	//Custom object for storing the data of each retval in a structurized way
	//normal_data: data found before the start of the primer atributes list
	//primer_data: a list of atributes of each primer found in a retval
	static public class P3RetValData{
		List<String> normal_data;
		List<List<String>> primer_data;
		
		public P3RetValData() {
			normal_data = new ArrayList<String>();
			primer_data = new ArrayList<List<String>>();
		}
	}
	
	public static List<P3RetVal> read_from_file(String path, List<P3GlobalSettings> pa,List<SeqArgs> sa){
		List<P3RetVal> result = new ArrayList<P3RetVal>();
		Pattern p = Pattern.compile("[0-9]+");
		List<P3RetValData> listOResults = new ArrayList<P3RetValData>(); //list of the generated retvals
		P3RetValData dat_obj = new P3RetValData();
		List<String> primer_data = null;//initialization later
		int index = -1;
		int tmp_index;
		
		try {
			File file = new File(path);
			Scanner rfile = new Scanner(file);
			
			// This loop fills the list listOResults with the data of the differents retvals
			while (rfile.hasNextLine()) {
				String data = rfile.nextLine();
				if (!data.equals("=")) {
					//left right pair internal pair?
					// get the primers lines
					if (data.matches("PRIMER_[A-Z]+_[0-9]+.*")) {
						//get the primer number
						//System.out.println(data);
						Matcher matcher = p.matcher(data);
						matcher.find();
						tmp_index = Integer.parseInt(matcher.group());
						// if the index has changed we store the array and start a new one
						if (index < tmp_index) {
							index = tmp_index;
							if(primer_data != null) {dat_obj.primer_data.add(primer_data);}
							primer_data = new ArrayList<String>();
						}
						primer_data.add(data);
					} else {
						dat_obj.normal_data.add(data);
					}
				} else {
					//a new retval starts
					dat_obj.primer_data.add(primer_data);
					listOResults.add(dat_obj);
					index= -1;
					dat_obj = new P3RetValData();
					primer_data = null;
				}
			}
			rfile.close();
			
		} catch (FileNotFoundException e){
			System.out.println("An error occurred.");
		    e.printStackTrace();
		    return null;
		}

		int indi = 0;
		for ( P3RetValData obj : listOResults) {
			P3RetVal new_retval = new P3RetVal(pa.get(indi), sa.get(indi));
			
			if (pa.get(indi).isPickLeftPrimer() && pa.get(indi).isPickRightPrimer()) {
				new_retval.output_type = P3OutputType.primer_pairs;
			} else {
				new_retval.output_type = P3OutputType.primer_list;
			}
			if (	pa.get(indi).getPrimerTask() == P3Task.PICK_PRIMER_LIST ||
					pa.get(indi).getPrimerTask() == P3Task.PICK_SEQUENCING_PRIMERS) {
				new_retval.output_type = P3OutputType.primer_list;
			}
			
			
			int print_fwd = 0;
			int print_rev = 0;
			int print_int = 0;
			new_retval.fwd = new OligoArray(OligoType.OT_LEFT);
			new_retval.rev = new OligoArray(OligoType.OT_RIGHT);
			new_retval.intl = new OligoArray(OligoType.OT_INTL);
			new_retval.best_pairs = new PairArrayT();
			for ( String line : obj.normal_data) {
			
				//Printing 3
				if(line.contains("PRIMER_STOP_CODON_POSITION")) { 
					new_retval.stop_codon_pos = Integer.parseInt(line.replace("PRIMER_STOP_CODON_POSITION=", ""));
					new_retval.set_retval_both_stop_codons();
				}
				if(line.contains("PRIMER_LEFT") && line.contains("_NUM_RETURNED=")) {
					print_fwd = Integer.parseInt(line.substring(line.indexOf("=")+1));
				}
				if(line.contains("PRIMER_RIGHT") && line.contains("_NUM_RETURNED=")) {
					print_rev = Integer.parseInt(line.substring(line.indexOf("=")+1));
				}
				if(line.contains("PRIMER_INTERNAL") && line.contains("_NUM_RETURNED=")) {
					print_int = Integer.parseInt(line.substring(line.indexOf("=")+1));
				}
				//STATS
				if(line.contains("PRIMER_LEFT_EXPLAIN")) {
					new_retval.fwd.expl = oligostats_from_explain_string(line.substring(line.indexOf("=")+1));
				}
				if(line.contains("PRIMER_RIGHT_EXPLAIN")) {
					new_retval.rev.expl = oligostats_from_explain_string(line.substring(line.indexOf("=")+1));
				}
				if(line.contains("PRIMER_INTERNAL_EXPLAIN")) {
					new_retval.intl.expl = oligostats_from_explain_string(line.substring(line.indexOf("=")+1));
				}
				if(line.contains("PRIMER_PAIR_EXPLAIN")) {
					new_retval.best_pairs.expl = pairstats_from_explain_string(line.substring(line.indexOf("=")+1));
				}
				
				
			}
			
			//print the primer data
			for (List<String> datalist : obj.primer_data) {
				PrimerPair pp = new PrimerPair();
				
				PrimerRecord left = new PrimerRecord(OligoType.OT_LEFT);
				PrimerRecord right = new PrimerRecord(OligoType.OT_RIGHT);
				PrimerRecord intl = new PrimerRecord(OligoType.OT_INTL);
				
				if(datalist != null) {
					for (String line : datalist) {
						//QUALITY
						if(line.contains("PRIMER_PAIR") && line.contains("_PENALTY=")) {
							pp.pair_quality = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						
						if(line.contains("PRIMER_LEFT") && line.contains("_PENALTY=")) {
							left.quality = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_RIGHT") && line.contains("_PENALTY=")) {
							right.quality = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_INTERNAL") && line.contains("_PENALTY=")) {
							intl.quality = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						
						//PROBLEMS
						
						//SEQUENCES
						//useless
						
						//START AND LENGHT
						if(line.matches("PRIMER_LEFT_[0-9]+=[0-9]+,[0-9]+")) {
							List<String> start_length = Arrays.asList(line.substring(line.indexOf("=")+1).split(","));
							left.start = Integer.parseInt(start_length.get(0)) - sa.get(indi).getIncludedRegionStart() - pa.get(indi).getFirstBaseIndex();
							left.length = Integer.parseInt(start_length.get(1));
						}
						if(line.matches("PRIMER_RIGHT_[0-9]+=[0-9]+,[0-9]+")) {
							List<String> start_length = Arrays.asList(line.substring(line.indexOf("=")+1).split(","));
							right.start = Integer.parseInt(start_length.get(0)) - sa.get(indi).getIncludedRegionStart() - pa.get(indi).getFirstBaseIndex();
							right.length = Integer.parseInt(start_length.get(1));
						}
						if(line.matches("PRIMER_INTERNAL_[0-9]+=[0-9]+,[0-9]+")) {
							List<String> start_length = Arrays.asList(line.substring(line.indexOf("=")+1).split(","));
							intl.start = Integer.parseInt(start_length.get(0)) - sa.get(indi).getIncludedRegionStart() - pa.get(indi).getFirstBaseIndex();
							intl.length = Integer.parseInt(start_length.get(1));
						}
						
						//TM
						if(line.contains("PRIMER_LEFT") && line.contains("_TM=")) {
							left.temp = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_RIGHT") && line.contains("_TM=")) {
							right.temp = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_INTERNAL") && line.contains("_TM=")) {
							intl.temp = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						
						//GC CONTENT
						if(line.contains("PRIMER_LEFT") && line.contains("_GC_PERCENT=")) {
							left.gc_content = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_RIGHT") && line.contains("_GC_PERCENT=")) {
							right.gc_content = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_INTERNAL") && line.contains("_GC_PERCENT=")) {
							intl.gc_content = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						
						//SELF ANY
						if(line.contains("PRIMER_LEFT") && line.contains("_SELF_ANY=")) {
							left.self_any = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_RIGHT") && line.contains("_SELF_ANY=")) {
							right.self_any = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_INTERNAL") && line.contains("_SELF_ANY=")) {
							intl.self_any = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_INTERNAL") && line.contains("_SELF_ANY_TH=")) {
							intl.self_any = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_LEFT") && line.contains("_SELF_ANY_TH=")) {
							left.self_any = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_RIGHT") && line.contains("_SELF_ANY_TH=")) {
							right.self_any = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						
						//SELF END
						if(line.contains("PRIMER_LEFT") && line.contains("_SELF_END=")) {
							left.self_end = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_RIGHT") && line.contains("_SELF_END=")) {
							right.self_end = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_INTERNAL") && line.contains("_SELF_END=")) {
							intl.self_end = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_INTERNAL") && line.contains("_SELF_END_TH=")) {
							intl.self_end = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_LEFT") && line.contains("_SELF_END_TH=")) {
							left.self_end = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_RIGHT") && line.contains("_SELF_END_TH=")) {
							right.self_end = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						
						//HAIRPIN TH
						if(line.contains("PRIMER_LEFT") && line.contains("_HAIRPIN_TH=")) {
							left.hairpin_th = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_RIGHT") && line.contains("_HAIRPIN_TH=")) {
							right.hairpin_th = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_INTERNAL") && line.contains("_HAIRPIN_TH=")) {
							intl.hairpin_th = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						
						//mispriming
						if(line.contains("PRIMER_LEFT") && line.contains("_LIBRARY_MISPRIMING=")) {
							List<String> max_name = Arrays.asList(line.substring(line.indexOf("=")+1).split(","));
							left.repeat_sim.score.add(Double.parseDouble(max_name.get(0)));
							left.repeat_sim.name = max_name.get(1);
						}
						if(line.contains("PRIMER_RIGHT") && line.contains("_LIBRARY_MISPRIMING=")) {
							List<String> max_name = Arrays.asList(line.substring(line.indexOf("=")+1).split(","));
							right.repeat_sim.score.add(Double.parseDouble(max_name.get(0)));
							right.repeat_sim.name = max_name.get(1);
						}
						if(line.contains("PRIMER_PAIR") && line.contains("_LIBRARY_MISPRIMING=")) {
							List<String> max_name = Arrays.asList(line.substring(line.indexOf("=")+1).split(","));
							pp.repeat_sim = Double.parseDouble(max_name.get(0));
							pp.rep_name = max_name.get(1);
						}
						if(line.contains("PRIMER_INTERNAL") && line.contains("_LIBRARY_MISHYB=")) {
							List<String> max_name = Arrays.asList(line.substring(line.indexOf("=")+1).split(","));
							intl.repeat_sim.score.add(Double.parseDouble(max_name.get(0)));
							intl.repeat_sim.name = max_name.get(1);
						}
						
						//_MIN_SEQ_QUALITY
						if(line.contains("PRIMER_LEFT") && line.contains("_MIN_SEQ_QUALITY=")) {
							left.seq_quality = Integer.parseInt(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_RIGHT") && line.contains("_MIN_SEQ_QUALITY=")) {
							right.seq_quality = Integer.parseInt(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_INTERNAL") && line.contains("_MIN_SEQ_QUALITY=")) {
							intl.seq_quality = Integer.parseInt(line.substring(line.indexOf("=")+1));
						}
						
						//POSITION_PENALTY
						if(line.contains("PRIMER_LEFT") && line.contains("_POSITION_PENALTY=")) {
							left.position_penalty = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_RIGHT") && line.contains("_POSITION_PENALTY=")) {
							right.position_penalty = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						
						//END_STABILITY
						if(line.contains("PRIMER_LEFT") && line.contains("_END_STABILITY=")) {
							left.end_stability = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_RIGHT") && line.contains("_END_STABILITY=")) {
							right.end_stability = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						
						//template mispriming TODO check it
						if(line.contains("PRIMER_LEFT") && line.contains("_TEMPLATE_MISPRIMING=")) {
							left.template_mispriming = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_RIGHT") && line.contains("_TEMPLATE_MISPRIMING=")) {
							right.template_mispriming = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						
						//template mispriming thermo TODO check it
						if(line.contains("PRIMER_LEFT") && line.contains("_TEMPLATE_MISPRIMING_TH=")) {
							left.template_mispriming = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_RIGHT") && line.contains("_TEMPLATE_MISPRIMING_TH=")) {
							right.template_mispriming = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						
						if(line.contains("PRIMER_PAIR") && line.contains("_COMPL_ANY=")) {
							pp.compl_any = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_PAIR") && line.contains("_COMPL_ANY_TH=")) {
							pp.compl_any = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_PAIR") && line.contains("_COMPL_END=")) {
							pp.compl_any = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_PAIR") && line.contains("_COMPL_END_TH=")) {
							pp.compl_end = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_PAIR") && line.contains("_PRODUCT_SIZE=")) {
							pp.product_size = Integer.parseInt(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_PAIR") && line.contains("_PRODUCT_TM=")) {
							pp.product_tm = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_PAIR") && line.contains("_PRODUCT_TM_OLIGO_TM_DIFF=")) {
							pp.product_tm_oligo_tm_diff = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_PAIR") && line.contains("_T_OPT_A=")) {
							pp.t_opt_a = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_PAIR") && line.contains("_TEMPLATE_MISPRIMING=")) {
							pp.template_mispriming = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						if(line.contains("PRIMER_PAIR") && line.contains("_TEMPLATE_MISPRIMING_TH=")) {
							pp.template_mispriming = Double.parseDouble(line.substring(line.indexOf("=")+1));
						}
						
					}
				
					//storing pairs in their place
					if (new_retval.output_type == P3OutputType.primer_list) {
						if ((pa.get(indi).isPickLeftPrimer())) {
							new_retval.fwd.add_oligo_to_oligo_array(left);
						}
						if ((pa.get(indi).isPickRightPrimer())) {
							new_retval.rev.add_oligo_to_oligo_array(right);
						}
						if ((pa.get(indi).isPickInternalOligo())) {
							new_retval.intl.add_oligo_to_oligo_array(intl);
						}
						
					}else {
						pp.left = left;
						pp.right = right;
						if (pa.get(indi).isPickInternalOligo() ) {
							pp.intl = intl;
						}
						new_retval.best_pairs.add_pair(pp);
					}
				}// end if null
			}//end primer data loop
			result.add(new_retval);
			indi++; 
		}//end list of results
	
		
		return result;
	}
	
	public static OligoStats oligostats_from_explain_string(String s) {
		String aux;
		OligoStats ps = new OligoStats();
		List<String> valuesList = Arrays.asList(s.split(","));
		
		for (int i = 0; i < valuesList.size(); i++) {
			if (valuesList.get(i).contains("considered")) {
				aux = valuesList.get(i);
				aux = aux.replace("considered",""); 
				aux = aux.replace(" ","");
				ps.considered = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("would not amplify any of the ORF")) {
				aux = valuesList.get(i);
				aux = aux.replace("would not amplify any of the ORF",""); 
				aux = aux.replace(" ","");
				ps.no_orf = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("too many Ns")) {
				aux = valuesList.get(i);
				aux = aux.replace("too many Ns",""); 
				aux = aux.replace(" ","");
				ps.ns = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("overlap target")) {
				aux = valuesList.get(i);
				aux = aux.replace("overlap target",""); 
				aux = aux.replace(" ","");
				ps.target = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("overlap excluded region")) {
				aux = valuesList.get(i);
				aux = aux.replace("overlap excluded region",""); 
				aux = aux.replace(" ","");
				ps.excluded = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("GC content failed")) {
				aux = valuesList.get(i);
				aux = aux.replace("GC content failed",""); 
				aux = aux.replace(" ","");
				ps.gc = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("GC clamp failed")) {
				aux = valuesList.get(i);
				aux = aux.replace("GC clamp failed",""); 
				aux = aux.replace(" ","");
				ps.gc_clamp = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("low tm")) {
				aux = valuesList.get(i);
				aux = aux.replace("low tm",""); 
				aux = aux.replace(" ","");
				ps.temp_min = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("high tm")) {
				aux = valuesList.get(i);
				aux = aux.replace("high tm",""); 
				aux = aux.replace(" ","");
				ps.temp_max = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("high any compl")) {
				aux = valuesList.get(i);
				aux = aux.replace("high any compl",""); 
				aux = aux.replace(" ","");
				ps.compl_any = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("high end compl")) {
				aux = valuesList.get(i);
				aux = aux.replace("high end compl",""); 
				aux = aux.replace(" ","");
				ps.compl_end = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("high hairpin stability")) {
				aux = valuesList.get(i);
				aux = aux.replace("high hairpin stability",""); 
				aux = aux.replace(" ","");
				ps.hairpin_th = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("high repeat similarity")) {
				aux = valuesList.get(i);
				aux = aux.replace("high repeat similarity",""); 
				aux = aux.replace(" ","");
				ps.repeat_score = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("long poly-x seq")) {
				aux = valuesList.get(i);
				aux = aux.replace("long poly-x seq",""); 
				aux = aux.replace(" ","");
				ps.poly_x = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("low sequence quality")) {
				aux = valuesList.get(i);
				aux = aux.replace("low sequence quality",""); 
				aux = aux.replace(" ","");
				ps.seq_quality = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("high 3' stability")) {
				aux = valuesList.get(i);
				aux = aux.replace("high 3' stability",""); 
				aux = aux.replace(" ","");
				ps.stability = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("high template mispriming score")) {
				aux = valuesList.get(i);
				aux = aux.replace("high template mispriming score",""); 
				aux = aux.replace(" ","");
				ps.template_mispriming = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("lowercase masking of 3' end")) {
				aux = valuesList.get(i);
				aux = aux.replace("lowercase masking of 3' end",""); 
				aux = aux.replace(" ","");
				ps.gmasked = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("failed must_match requirements")) {
				aux = valuesList.get(i);
				aux = aux.replace("failed must_match requirements",""); 
				aux = aux.replace(" ","");
				ps.must_match_fail = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("not in any ok left region")) {
				aux = valuesList.get(i);
				aux = aux.replace("failed must_match requirements",""); 
				aux = aux.replace(" ","");
				ps.not_in_any_left_ok_region = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("not in any ok right region")) {
				aux = valuesList.get(i);
				aux = aux.replace("not in any ok right region",""); 
				aux = aux.replace(" ","");
				ps.not_in_any_right_ok_region = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("ok")) {
				aux = valuesList.get(i);
				aux = aux.replace("ok",""); 
				aux = aux.replace(" ","");
				ps.ok = Integer.parseInt(aux);
			}
		}
	
		
		return ps;
	}
	
	public static PairStats pairstats_from_explain_string(String s) {
		String aux;
		PairStats ps = new PairStats();
		List<String> valuesList = Arrays.asList(s.split(","));
		
		for (int i = 0; i < valuesList.size(); i++) {
			if (valuesList.get(i).contains("considered")) {
				aux = valuesList.get(i);
				aux = aux.replace("considered",""); 
				aux = aux.replace(" ","");
				ps.considered = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("no target")) {
				aux = valuesList.get(i);
				aux = aux.replace("no target",""); 
				aux = aux.replace(" ","");
				ps.target = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("unacceptable product size")) {
				aux = valuesList.get(i);
				aux = aux.replace("unacceptable product size",""); 
				aux = aux.replace(" ","");
				ps.product = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("low product Tm")) {
				aux = valuesList.get(i);
				aux = aux.replace("low product Tm",""); 
				aux = aux.replace(" ","");
				ps.low_tm = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("high product Tm")) {
				aux = valuesList.get(i);
				aux = aux.replace("high product Tm",""); 
				aux = aux.replace(" ","");
				ps.high_tm = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("tm diff too large")) {
				aux = valuesList.get(i);
				aux = aux.replace("tm diff too large",""); 
				aux = aux.replace(" ","");
				ps.temp_diff = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("high any compl")) {
				aux = valuesList.get(i);
				aux = aux.replace("high any compl",""); 
				aux = aux.replace(" ","");
				ps.compl_any = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("high end compl")) {
				aux = valuesList.get(i);
				aux = aux.replace("high end compl",""); 
				aux = aux.replace(" ","");
				ps.target = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("no internal oligo")) {
				aux = valuesList.get(i);
				aux = aux.replace("no internal oligo",""); 
				aux = aux.replace(" ","");
				ps.internal = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("high mispriming library similarity")) {
				aux = valuesList.get(i);
				aux = aux.replace("high mispriming library similarity",""); 
				aux = aux.replace(" ","");
				ps.repeat_sim = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("no overlap of required point")) {
				aux = valuesList.get(i);
				aux = aux.replace("no overlap of required point",""); 
				aux = aux.replace(" ","");
				ps.does_not_overlap_a_required_point = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("primer in pair overlaps a primer in a better pair")) {
				aux = valuesList.get(i);
				aux = aux.replace("primer in pair overlaps a primer in a better pair",""); 
				aux = aux.replace(" ","");
				ps.overlaps_oligo_in_better_pair = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("high template mispriming score")) {
				aux = valuesList.get(i);
				aux = aux.replace("high template mispriming score",""); 
				aux = aux.replace(" ","");
				ps.template_mispriming = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("not in any ok region")) {
				aux = valuesList.get(i);
				aux = aux.replace("not in any ok region",""); 
				aux = aux.replace(" ","");
				ps.not_in_any_ok_region = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("left primer to right of right primer")) {
				aux = valuesList.get(i);
				aux = aux.replace("left primer to right of right primer",""); 
				aux = aux.replace(" ","");
				ps.reversed = Integer.parseInt(aux);

			}else
			if (valuesList.get(i).contains("ok")) {
				aux = valuesList.get(i);
				aux = aux.replace("ok",""); 
				aux = aux.replace(" ","");
				ps.ok = Integer.parseInt(aux);
			}
		}
	
		
		return ps;
	}
}
