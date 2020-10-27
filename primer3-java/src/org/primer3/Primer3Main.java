/*
    This file is part of primer3 porting to java


	Original file are part of https://github.com/primer3-org/primer3
	Whitehead Institute for Biomedical Research, Steve Rozen
	(http://purl.com/STEVEROZEN/), Andreas Untergasser and Helen Skaletsky
	All rights reserved to Primer3 authors.

    Primer3 and the libprimer3 library are free software;
    you can redistribute them and/or modify them under the terms
    of the GNU General Public License as published by the Free
    Software Foundation; either version 2 of the License, or (at
    your option) any later version.

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this software (file gpl-2.0.txt in the source
    distribution); if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
package org.primer3;

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
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.primer3.libprimer3.LibPrimer3;
import org.primer3.libprimer3.P3GlobalSettings;
import org.primer3.libprimer3.P3Task;
import org.primer3.libprimer3.P3RetVal;
import org.primer3.libprimer3.SeqArgs;
import org.primer3.masker.formula_parameters;
import org.primer3.masker.masker;
import org.primer3.masker.masking_direction;
import org.primer3.p3_seq_lib.seq_lib;
import org.primer3.thal.ThAl;

public class Primer3Main {


	static String pr_release = LibPrimer3.libprimer3_release();
	static String pr_program_name = "" ;

	static CommandLine setArgs(String[] args)
	{
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();

		Option  argOption   = Option.builder("formated_output").argName( "formated_output" ).longOpt("formated_output")
				.desc("formated_output" )
				.build();
		options.addOption(argOption);





		argOption = Option.builder("fasta")
				.hasArg()
				.argName( "fasta" )
				.longOpt("fasta")    
				.desc("Fasta file contains input sequenes" )
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

			//			if(line.getOptions().length == 0 )
			//
			//			{
			//				formatter.printHelp(pr_program_name + " [options] <inputFile>", "head", options, "");
			//				//				formatter.printHelp( pr_program_name, options );
			//			}

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

		//		args = "/data/softwares/primer3-primer3/test/primer_rat_input".split(" ");

		//		args = "/data/softwares/primer3-primer3/example".split(" ");

		//		args = "/data/softwares/primer3-primer3/test/primer_not_ok_regions_input".split(" ");
		//		args = "/data/softwares/primer3-primer3/test/primer_first_base_index_input ".split( " " );

		//		args = "/data/softwares/primer3-primer3/test/primer_high_gc_load_set_input --p3_settings_file /data/softwares/primer3-primer3/test/primer_high_gc_load_set.set".split( " " );


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


		String p3_settings_path = null;
		String output_path = null;
		String error_path = null;
		StringBuilder fatal_parse_err = new StringBuilder();
		StringBuilder nonfatal_parse_err = new StringBuilder();
		StringBuilder warnings = new StringBuilder();



		P3GlobalSettings global_pa;
		SeqArgs sarg;
		// TODO :: Missing ??
//		read_boulder_record_results read_boulder_record_res ;




		P3RetVal retval = null;
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

		List<String> leftArgs = line.getArgList();




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
		if(line.hasOption("formated_output"))
			format_output = true;


		if(line.hasOption("fasta"))
			fastaInputFile = line.getOptionValue("fasta");


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


		try {
			while(true)
			{
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
				// TODO :: read_boulder_record
				if(! boulder.read_boulder_record(scan,
						strict_tags,
						io_version,
						!format_output,
						P3FileType.all_parameters,
						global_pa,
						sarg,
						fatal_parse_err,
						nonfatal_parse_err,
						warnings))
				{
					break; /* There were no more boulder records */
				}
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
						boulder.format_error( sarg.getSequenceName(), nonfatal_parse_err.toString());
					} else {
						boulder.print_boulder_error(nonfatal_parse_err.toString());
					}
					//		        goto loop_wrap_up;
				}
				/* Print any warnings and continue processing */
				if (!warnings.toString().isEmpty()) {
					if (format_output) {
						boulder.format_warning( sarg.getSequenceName(), warnings.toString());
					} else {
						boulder.print_boulder_warning(warnings.toString());
					}
				}

				if (global_pa.getFileFlag() == 1 && sarg.getSequenceName() == null) {
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
//				LibPrimer3.p3_set_gs_primer_file_flag(global_pa,
//						read_boulder_record_res.file_flag);
				retval = LibPrimer3.choose_primers(global_pa, sarg);

				if(sarg.isMultiplex() )
					lateResult.add(retval);

				//	       if (null == retval) exit(-2); /* Out of memory. */
				/* If it was necessary to use a left_input, right_input,
			   or internal_oligo_input primer that was
			   unacceptable, then add warnings. */
				if (global_pa.isPickAnyway() && format_output) {
					if (sarg.getLeftInput() != null) {
						retval.add_must_use_warnings( "Left primer", retval.fwd);
					}
					if (sarg.getRightInput() != null) {
						retval.add_must_use_warnings("Right primer", retval.rev);
					}
					if (sarg.getInternalInput() != null) {
						retval.add_must_use_warnings("Hybridization probe", retval.intl);
					}
				}
				// CONT :: Line 449


				if (retval.glob_err. length() == 0 && retval.per_sequence_err.length() == 0  ) {
					/* We need to test for errors before we call
			         p3_print_oligo_lists. This function only works on retval as
			         returned when there were no errors. */
					if (global_pa.getFileFlag() != 0) {
						/* Create files with left, right, and internal oligos. */
						LibPrimer3.p3_print_oligo_lists(retval, sarg, global_pa,
								retval.per_sequence_err,
								sarg.getSequenceName());
					}
				}
				if(! sarg.isMultiplex() )
				if (format_output) {
					boulder. print_format_output(io_version, global_pa,
							sarg, retval, pr_release,
							global_pa.getExplainFlag());
				} else {
					/* Use boulder output */
					retval.print_boulder(io_version, 
							global_pa.getExplainFlag() != 0);
				}


			} // end while(true)


			boolean multiplexHasResult = LibPrimer3.multiplexSearch.search();
			// this part is for multiplex after finishing each seq, we wait to all seqs within the same well are ready to search

			for(P3RetVal lateRetval : lateResult)
			{
				if (format_output) {
					boulder. print_format_output(io_version, global_pa,
							sarg, lateRetval, pr_release,
							global_pa.getExplainFlag());
				} else {
					/* Use boulder output */
					lateRetval.print_boulder(io_version, 
							global_pa.getExplainFlag() != 0);
				}
			}
			if(multiplexHasResult ) {
				LibPrimer3.multiplexSearch.print_boulder(io_version);
			}


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
}
