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
public class TestMultiplexInput {
	String message = "Hello World";	
	MessageUtil messageUtil = new MessageUtil(message);

	
	static String resourceFolder = "src/test/resources/multiplex_test/";
	

	
	
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
	public void testP3SpecificInput() throws FileNotFoundException, IOException {	  
		
		testCase("primer_task_specific_input",
				"primer_task_specific_input",
				"primer_task_specific_output",
				"-fasta src/test/resources/multiplex_test/primer_task_specific.fasta");
	}
	
	@Test
	public void testP3SpecificInput5T() throws FileNotFoundException, IOException {	  
		
		testCase("primer_task_specific_input_5T",
				"primer_task_specific_input_5T",
				"primer_task_specific_output",
				"-fasta src/test/resources/multiplex_test/targets_new.fasta");
	}	
	
	@Test
	public void testP3SpecificInput5Tm() throws FileNotFoundException, IOException {	  
		
		testCase("primer_task_specific_input_5T_m",
				"primer_task_specific_input_5T_m",
				"primer_task_specific_output",
				"-fasta src/test/resources/multiplex_test/targets_new.fasta");
	}	
	
	@Test
	public void testP3SpecificInput_7T() throws FileNotFoundException, IOException {	  
		
		testCase("primer_task_specific_input_7T",
				"primer_task_specific_input_7T",
				"primer_task_specific_output",
				"-fasta src/test/resources/multiplex_test/targets_new_fake.fasta");
	}	
	@Test
	public void testP3SpecificInput3() throws FileNotFoundException, IOException {	  
		
		testCase("primer_task_specific_input_3T",
				"primer_task_specific_input_3T",
				"primer_task_specific_output",
				"-fasta src/test/resources/multiplex_test/targets_new_fake2.fasta");
	}
	
	@Test
	public void testP3SpecificAlbInput() throws FileNotFoundException, IOException {	  
		
		testCase("primer_task_specific_alb_input",
				"primer_task_specific_alb_input",
				"primer_task_specific_output",
				"-fasta src/test/resources/multiplex_test/primer_task_specific.fasta");
	}	
	
}
