package primer3Test;

import static org.junit.Assert.*;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileDescriptor;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOError;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;

import org.apache.commons.io.IOUtils;
import org.junit.FixMethodOrder;
import org.junit.Ignore;
import org.junit.Test;
import org.junit.runners.MethodSorters;
import org.primer3.Primer3Main;

@FixMethodOrder(MethodSorters.JVM)
public class TestPrimer3Main {
	String message = "Hello World";	
	MessageUtil messageUtil = new MessageUtil(message);

	
	static String resourceFolder = "src/test/resources/";
	
	void testCase(String inputFile,String outputFile) throws FileNotFoundException, IOException
	{
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		System.setOut(new PrintStream(baos));
		Primer3Main.main((resourceFolder +inputFile).split(" "));
		
		List<String> resultlines =  IOUtils.readLines(new FileReader(new File(resourceFolder+outputFile)));
		
		String result = "";
		for(String line : resultlines)
			result += line + "\n";
		System.setOut(new PrintStream(new FileOutputStream(FileDescriptor.out)));
		String outputResult = baos.toString();
		
		assertEquals(result,outputResult);
	}
	
	
	void testCase(String msg, String inputFile,String outputFile) throws FileNotFoundException, IOException
	{
		testCase(msg,inputFile,outputFile, "");
	}
	
	
	void testCase(String msg, String inputFile,String outputFile, String addArgs) throws FileNotFoundException, IOException
	{
		if(addArgs == null )
			addArgs = "";
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		System.setOut(new PrintStream(baos));
		Primer3Main.main((resourceFolder +inputFile + " " + addArgs ).split(" "));
		
		List<String> resultlines =  IOUtils.readLines(new FileReader(new File(resourceFolder+outputFile)));
		
		String result = "";
		for(String line : resultlines)
			result += line + "\n";
		System.setOut(new PrintStream(new FileOutputStream(FileDescriptor.out)));
		String outputResult = baos.toString();
		
		assertEquals(msg ,result,outputResult);
	}
	
	
	@Test
	public void testPrimer3MainNoArgs() throws FileNotFoundException, IOException {	  
		
		
		
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		System.setOut(new PrintStream(baos));
		Primer3Main.main("src/test/resources/example".split(" "));
		
		List<String> resultlines =  IOUtils.readLines(new FileReader(new File("src/test/resources/example_result.txt")));
		
		String result = "";
		for(String line : resultlines)
			result += line + "\n";
		
		
		System.setOut(new PrintStream(new FileOutputStream(FileDescriptor.out)));
		String outputResult = baos.toString();
		
		assertEquals(result,outputResult);
	}

	
	@Ignore("Long test")
	@Test
	public void testP3Case1() throws FileNotFoundException, IOException {	  
		
		testCase("primer1_input","primer1_output");
		
	}
	
	
	
	
	@Test
	public void testP3CheckInput() throws FileNotFoundException, IOException {	  
		
		testCase("primer_check_input","primer_check_output");
	}
	
	@Test
	public void testP3EndPathology() throws FileNotFoundException, IOException {	  
		
		testCase("primer_end_pathology_input","primer_end_pathology_input","primer_end_pathology_output");
	}
	
	@Test
	public void testP3GCEnd() throws FileNotFoundException, IOException {	  
		
		testCase("primer_gc_end_input","primer_gc_end_input","primer_gc_end_output");
	}
	@Test
	public void testP3PrimerFirstBaseIndex() throws FileNotFoundException, IOException {	  
		
		testCase("primer_first_base_index_input","primer_first_base_index_input","primer_first_base_index_output");
	}
	
	@Test
	public void testP3HighGCLoadSet() throws FileNotFoundException, IOException {	  
		
		testCase("primer_high_gc_load_set_input",
				"primer_high_gc_load_set_input",
				"primer_high_gc_load_set_output",
				"--p3_settings_file " + resourceFolder + "primer_high_gc_load_set.set");
	}
	
	@Test
	public void testP3HighTMLoadSet() throws FileNotFoundException, IOException {	  
		
		testCase("primer_high_tm_load_set_input",
				"primer_high_tm_load_set_input",
				"primer_high_tm_load_set_output",
				"--p3_settings_file " + resourceFolder + "primer_high_tm_load_set.set");
	}
	
	@Test
	public void testP3InternalInput() throws FileNotFoundException, IOException {	  
		
		testCase("primer_internal_input",
				"primer_internal_input",
				"primer_internal_output");
	}
	
	
	@Test
	public void testP3OkRegionsInput1() throws FileNotFoundException, IOException {	  
		
		testCase("primer_ok_regions_input",
				"primer_ok_regions_input",
				"primer_ok_regions_output");
	}
	
	
//	@Ignore("Ignore for now")
	@Test
	public void testP3OkRegionsInput2() throws FileNotFoundException, IOException {	  
		
		testCase("primer_ok_regions2_input",
				"primer_ok_regions2_input",
				"primer_ok_regions2_output");
	}
	
	@Ignore("Ignore Not Implmented Now")
	@Test
	public void testP3LibAmbCodesInput() throws FileNotFoundException, IOException {	  
		
		testCase("primer_lib_amb_codes_input",
				"primer_lib_amb_codes_input",
				"primer_lib_amb_codes_output");
	}
	
	@Ignore("Ignore for Now")
	@Test
	public void testP3MustMatchInput() throws FileNotFoundException, IOException {	  
		
		testCase("primer_must_match_input",
				"primer_must_match_input",
				"primer_must_match_output");
	}
	
	
	@Test
	public void testP3MustOverlapPointInput() throws FileNotFoundException, IOException {	  
		
		testCase("primer_must_overlap_point_input",
				"primer_must_overlap_point_input",
				"primer_must_overlap_point_output");
	}
	
	
	@Test
	public void testP3NumBestInput() throws FileNotFoundException, IOException {	  
		
		testCase("primer_num_best_input",
				"primer_num_best_input",
				"primer_num_best_output");
	}
	
	@Ignore("Ignore Not Implmented Now")
	@Test
	public void testP3ObjFnInput() throws FileNotFoundException, IOException {	  
		
		testCase("primer_obj_fn_input",
				"primer_obj_fn_input",
				"primer_obj_fn_output");
	}
	
	
	@Test
	public void testP3OverlapJunctionInput() throws FileNotFoundException, IOException {	  
		
		testCase("primer_overlap_junction_input",
				"primer_overlap_junction_input",
				"primer_overlap_junction_output");
	}
	
	
	@Test
	public void testP3QualityBoundaryInput() throws FileNotFoundException, IOException {	  
		
		testCase("primer_quality_boundary_input",
				"primer_quality_boundary_input",
				"primer_quality_boundary_output");
	}
	
	
	@Test
	public void testP3RatInput() throws FileNotFoundException, IOException {	  
		
		testCase("primer_rat_input",
				"primer_rat_input",
				"primer_rat_output");
	}
	
	
	@Test
	public void testP3StartCodonInput() throws FileNotFoundException, IOException {	  
		
		testCase("primer_start_codon_input",
				"primer_start_codon_input",
				"primer_start_codon_output");
	}
	@Ignore("Very Long test")
	@Test
	public void testP3HumanInput() throws FileNotFoundException, IOException {	  
		
		testCase("primer_human_input",
				"primer_human_input",
				"primer_human_output");
	}
	
	@Test
	public void testP3NewTaskInput() throws FileNotFoundException, IOException {	  
		
		testCase("primer_new_tasks_input",
				"primer_new_tasks_input",
				"primer_new_tasks_output");
	}	
	
	
	
	
	
	
//	@Test
//	public void testP3Cases() throws FileNotFoundException, IOException {	  
//		
//		
//		String[] inputs = new String [] {
//				"src/test/resources/example",	
//				"src/test/resources/primer1_input"};
//		String[] outputs = new String [] {
//				"src/test/resources/example_result.txt",
//				"src/test/resources/primer1_output"};
//		
//		String[] expected = new String[outputs.length];
//		
//		for(int testCase=0;testCase < inputs.length; testCase++)
//		{
//			System.out.println("Input " + testCase);
//			String input = inputs[testCase];
//			String output = outputs[testCase];
//			ByteArrayOutputStream baos = new ByteArrayOutputStream();
//			System.setOut(new PrintStream(baos));
//			Primer3Main.main(input.split(" "));
//			List<String> resultlines =  IOUtils.readLines(new FileReader(new File(output)));
//			String result = "";
//			for(String line : resultlines)
//				result += line + "\n";
//			System.setOut(new PrintStream(new FileOutputStream(FileDescriptor.out)));
//			String outputResult = baos.toString();
//			assertEquals(result,outputResult);
//		}
//	}
	
}
