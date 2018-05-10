package com.primer3;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.List;
import java.util.Scanner;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import com.primer3.libprimer3.LibPrimer3;
import com.primer3.libprimer3.p3_global_settings;
import com.primer3.libprimer3.p3retval;
import com.primer3.libprimer3.seq_args;
import com.primer3.libprimer3.task;
import com.primer3.masker.formula_parameters;
import com.primer3.masker.masker;
import com.primer3.masker.masking_direction;
import com.primer3.thal.thallib;

public class Primer3Main {


	static String pr_release = LibPrimer3.libprimer3_release();
	static String pr_program_name = "" ;

	static CommandLine setArgs(String[] args)
	{
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();

		Option  argOption   = Option.builder("format_output").argName( "format_output" )
				.desc("format_output" )
				.build();
		options.addOption(argOption);

		argOption = Option.builder("d")
				.hasArg()
				.argName( "default_version" )
				.longOpt("default_version")    
				.desc("default_version" )
				.build();
		options.addOption(argOption);

		argOption = Option.builder("i")
				.hasArg()
				.argName( "io_version" )
				.longOpt("io_version")    
				.desc("io_version" )
				.build();
		options.addOption(argOption);


		argOption = Option.builder("p")
				.hasArg()
				.argName( "file_path" )
				.longOpt("p3_settings_file")    
				.desc("p3_settings_file" )
				.build();
		options.addOption(argOption);

		argOption = Option.builder()
				.longOpt("echo_settings_file")    
				.desc("echo_settings_file" )
				.build();
		options.addOption(argOption);

		argOption = Option.builder("s")
				.longOpt("strict_tags")    
				.desc("strict_tags" )
				.build();
		options.addOption(argOption);

		argOption = Option.builder("o")
				.hasArg()
				.argName( "output" )
				.longOpt("output")    
				.desc("output" )
				.build();
		options.addOption(argOption);

		argOption = Option.builder("e")
				.hasArg()
				.argName( "error" )
				.longOpt("error")    
				.desc("error" )
				.build();
		options.addOption(argOption);

		argOption = Option.builder("D")
				.desc("dump args" )
				.build();
		options.addOption(argOption);

		Option helpOption = new Option( "help", "Print Usage" );
		options.addOption(helpOption);

		// automatically generate the help statement
		HelpFormatter formatter = new HelpFormatter();
		formatter.setOptionComparator(null);
		formatter.setLongOptSeparator("=");
		try {
			// parse the command line arguments
			CommandLine line = parser.parse( options, args );

			if(line.getOptions().length == 0 )

			{
				formatter.printHelp(pr_program_name + " [options] <inputFile>", "head", options, "");
				//				formatter.printHelp( pr_program_name, options );
			}

			return line;		    
		}
		catch( ParseException exp ) {
			// oops, something went wrong
			System.err.println( "Parsing failed.  Reason: " + exp.getMessage() );
			//			formatter.printHelp( pr_program_name, options );
		}

		return null;
	}

	static public void main(String[] args)
	{

		args = "/data/softwares/primer3-primer3/test/primer_rat_input".split(" ");

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

		String p3_settings_path = null;
		String output_path = null;
		String error_path = null;
		StringBuilder fatal_parse_err = new StringBuilder();
		StringBuilder nonfatal_parse_err = new StringBuilder();
		StringBuilder warnings = new StringBuilder();



		p3_global_settings global_pa;
		seq_args sarg;
		// TODO :: Missing ??
		read_boulder_record_results read_boulder_record_res ;




		p3retval retval = null;
		int input_found=0;

		/* Get the program name for correct error messages */
		// pr_release = libprimer3_release();
		pr_program_name = args[0];
		LibPrimer3.p3_set_program_name(pr_program_name);
		CommandLine line = setArgs(args);

		if(line == null)
		{
			return;
		}

		List<String> leftArgs= line.getArgList();




		if(line.hasOption("a"))
			about=1;
		if(line.hasOption("p"))
			p3_settings_path = line.getOptionValue("p");

		if(line.hasOption("o"))
			output_path = line.getOptionValue("o");

		if(line.hasOption("e"))
			error_path = line.getOptionValue("e");

		if(line.getOptionValue("i","-1").equals("4"))
			io_version = 4;
		else
			io_version = 4;
		if(line.getOptionValue("d","").equals("1"))
			default_version = 4;
		else if(line.getOptionValue("d","").equals("2"))
			default_version = 2;
		else
			default_version = 2;
		if(line.hasOption("D"))
			dump_args = true ;

		// set redirfile
		if(error_path != null && !error_path.isEmpty())
		{
			try
			{
				PrintStream errStream = outputFile(error_path);
				System.setErr(errStream);
			}
			catch(FileNotFoundException fnEx)
			{
				System.err.println( "Error creating file : " + error_path );
				System.exit(-1);
			}
		}

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




		if(leftArgs.size() > 0)
		{
			String inputFileName = leftArgs.get(0);
			try {
				System.setIn(new FileInputStream(new File(inputFileName)));

			} catch (FileNotFoundException e) {
				System.err.println("Error opening file " + inputFileName);
				System.exit(-1);
			}
		}

		Scanner scan  = new Scanner(System.in);

		{
			//			InputStreamReader ir = new InputStreamReader( System.in);
			//			while(scan.hasNextLine())
			//			{
			//				System.out.println(scan.nextLine());
			//			}
		}



		/* Allocate the space for global settings and fill in default parameters */
		global_pa = p3_global_settings.p3_create_global_settings(default_version);

		if(global_pa == null)
		{
			print_usage();
			System.exit(-1);
		}

		global_pa.setDump(dump_args) ;


		/* Settings files have to be read in just below, and
	    	the functions need a temporary sarg */
		sarg = new seq_args();
		read_boulder_record_res = new read_boulder_record_results();
		/* Read data from the settings file until a "=" line occurs.  Assign
	     parameter values for primer picking to pa and sa. */
		if(p3_settings_path != null && !p3_settings_path.isEmpty())
		{
			boulder.read_p3_file(p3_settings_path, p3_file_type.settings,
					echo_settings && !format_output, 
					strict_tags, 
					global_pa, sarg, 
					fatal_parse_err,
					nonfatal_parse_err, warnings, read_boulder_record_res);

			/* Check if any thermodynamical alignment flag was given */
			if ((global_pa.thermodynamic_oligo_alignment == 1) || 
					(global_pa.thermodynamic_template_alignment == 1))
				read_thermodynamic_parameters();
			/* Check if masking template flag was given */
			if (global_pa.mask_template)
				validate_kmer_lists_path();


		}

		/* We also need to print out errors here because the loop erases all
	     errors at start. If there are fatal errors, write the proper
	     message and exit */
		if (!fatal_parse_err.toString().isEmpty()) {
			if (format_output) {
				boulder.format_error(sarg.sequence_name, fatal_parse_err.toString());
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
				boulder.format_error( sarg.sequence_name, nonfatal_parse_err.toString());
			} else {
				boulder.print_boulder_error(nonfatal_parse_err.toString());
			}
		}
		/* The temporary sarg is not needed any more */
		//	    destroy_seq_args(sarg);
		sarg = null ;


		while(true)
		{
			/* Create and initialize a seq_args data structure. sa (seq_args *) is
			 * initialized here because Values are _not_ retained across different
			 * input records. */
			sarg = new seq_args();
			read_boulder_record_res = new read_boulder_record_results();
			/* Reset all errors handlers and the return structure */
			fatal_parse_err = new StringBuilder();
			nonfatal_parse_err= new StringBuilder();
			warnings = new StringBuilder();
			retval = null;

			/* There were no more boulder records */
			// TODO :: read_boulder_record
			if(! boulder.read_boulder_record(scan,
					strict_tags,
					io_version,
					!format_output,
					p3_file_type.all_parameters,
					global_pa,
					sarg,
					fatal_parse_err,
					nonfatal_parse_err,
					warnings,
					read_boulder_record_res))
			{
				break; /* There were no more boulder records */
			}
			if(global_pa.mask_template){
				global_pa.lowercase_masking=global_pa.mask_template;
			}
			/* Check if any thermodynamical alignment flag was given and the
	        path to the parameter files changed - we need to reread them */
			if (((global_pa.thermodynamic_oligo_alignment == 1) ||
					(global_pa.thermodynamic_template_alignment == 1))
					&& (thermodynamic_path_changed ))
				read_thermodynamic_parameters();

			/* Check if template masking flag was given */
			if (global_pa.mask_template )
				validate_kmer_lists_path();


			/* Check that we found the thermodynamic parameters in case any thermodynamic flag was set to 1. */
			if (((global_pa.thermodynamic_oligo_alignment == 1) ||
					(global_pa.thermodynamic_template_alignment == 1))
					&& (thermodynamic_params_path == null)) {
				/* no parameter directory found, error */
				System.out.println("PRIMER_ERROR=thermodynamic approach chosen, but path to thermodynamic parameters not specified\n=\n");
				System.exit(-1);
			}

			/* Check that we found the kmer lists in case masking flag was set to 1. */
			if (global_pa.mask_template  && kmer_lists_path == null){
				System.out.println("PRIMER_ERROR=masking template chosen, but path to kmer lists not specified\n=\n");
				System.exit(-1);
			}

			/* Set up some masking parameters */
			/* edited by M. Lepamets */
			if (global_pa.mask_template ) {
				global_pa.mp.window_size = masker.DEFAULT_WORD_LEN_2;

				if (global_pa.pick_right_primer ) global_pa.mp.mdir = masking_direction.fwd;
				else if (global_pa.pick_left_primer ) global_pa.mp.mdir = masking_direction.rev;
				/* Check if masking parameters (k-mer list usage) have changed */
				if (global_pa.masking_parameters_changed) {
					//	            masker.delete_formula_parameters (global_pa.mp.fp, global_pa.mp.nlists);
					global_pa.mp.fp = formula_parameters.create_default_formula_parameters (global_pa.mp.list_prefix, kmer_lists_path);
					global_pa.masking_parameters_changed = false;
				}
			}

			input_found = 1;
			if ((global_pa.primer_task == task.generic)
					&& (global_pa.pick_internal_oligo )){
				PR_ASSERT(global_pa.pick_internal_oligo);
			}

			/* TODO :: If there are fatal errors, write the proper message and exit */
			//	       if (fatal_parse_err.data != NULL) {
			//	         if (format_output) {
			//	           format_error(stdout, sarg->sequence_name, fatal_parse_err.data);
			//	         } else {
			//	           print_boulder_error(fatal_parse_err.data);
			//	         }
			//	         fprintf(stderr, "%s: %s\n",
			//	                 pr_program_name, fatal_parse_err.data);
			//	         destroy_p3retval(retval);
			//	         destroy_seq_args(sarg);
			//	         exit(-4);
			//	       }

			if (!nonfatal_parse_err.toString().isEmpty()) {
				if (format_output) {
					boulder.format_error( sarg.sequence_name, nonfatal_parse_err.toString());
				} else {
					boulder.print_boulder_error(nonfatal_parse_err.toString());
				}
				//		        goto loop_wrap_up;
			}
			/* Print any warnings and continue processing */
			if (!warnings.toString().isEmpty()) {
				if (format_output) {
					boulder.format_warning( sarg.sequence_name, warnings.toString());
				} else {
					boulder.print_boulder_warning(warnings.toString());
				}
			}

			if (read_boulder_record_res.file_flag == 1 && sarg.sequence_name == null) {
				/* We will not have a base name for the files */
				if (format_output) {
					boulder.format_error((String)null,
							"Need PRIMER_SEQUENCE_ID if PRIMER_FILE_FLAG is not 0");
				} else {
					boulder.print_boulder_error("Need PRIMER_SEQUENCE_ID if PRIMER_FILE_FLAG is not 0");
				}
				//	    	      goto loop_wrap_up;
			}

			/* Pick the primers - the central function */
			LibPrimer3.p3_set_gs_primer_file_flag(global_pa,
					read_boulder_record_res.file_flag);
			retval = LibPrimer3.choose_primers(global_pa, sarg);
			//	       if (null == retval) exit(-2); /* Out of memory. */
			/* If it was necessary to use a left_input, right_input,
	       or internal_oligo_input primer that was
	       unacceptable, then add warnings. */
			if (global_pa.pick_anyway && format_output) {
				if (sarg.left_input != null) {
					retval.add_must_use_warnings( "Left primer", retval.fwd);
				}
				if (sarg.right_input != null) {
					retval.add_must_use_warnings("Right primer", retval.rev);
				}
				if (sarg.internal_input != null) {
					retval.add_must_use_warnings("Hybridization probe", retval.intl);
				}
			}
			// CONT :: Line 449


			if (retval.glob_err. length() == 0 && retval.per_sequence_err.length() == 0  ) {
				/* We need to test for errors before we call
		         p3_print_oligo_lists. This function only works on retval as
		         returned when there were no errors. */
				if (read_boulder_record_res.file_flag != 0) {
					/* Create files with left, right, and internal oligos. */
					LibPrimer3.p3_print_oligo_lists(retval, sarg, global_pa,
							retval.per_sequence_err,
							sarg.sequence_name);
				}
			}
			if (format_output) {
			    boulder. print_format_output(io_version, global_pa,
			                          sarg, retval, pr_release,
			                          read_boulder_record_res.explain_flag);
			    } else {
			      /* Use boulder output */
			    	retval.print_boulder(io_version, global_pa, sarg, 
			                    read_boulder_record_res.explain_flag != 0);
			    }

			
		} // end while(true)
	}



	private static void PR_ASSERT(boolean pick_internal_oligo) {
		// TODO Auto-generated method stub

	}



	static boolean   thermodynamic_path_changed = true;/* if this is set to 1, we need to re-read the thermodynamic parameters from new path */
	static String thermodynamic_params_path; /* path to thermodynamic parameter files */
	static String kmer_lists_path = null; /* path to kmer lists files */
	private static void validate_kmer_lists_path() {
		// TODO Auto-generated method stub

	}
	private static void read_thermodynamic_parameters() {

		thallib.get_thermodynamic_values();

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
}
