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


// mispairing test

@FixMethodOrder(MethodSorters.JVM)
public class TestGroupB {
	String message = "Hello World";	
	MessageUtil messageUtil = new MessageUtil(message);

	
	static String resourceFolder = "src/test/resources/groupB_test/";
	

	
	
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
	public void testP3Mispriming() throws FileNotFoundException, IOException {	  
		
		testCase("primer_mispriming_input",
				"primer_mispriming_input",
				"primer_mispriming_output");
	}	
	
	@Test
	public void testP3MisprimingTH() throws FileNotFoundException, IOException {	  
		
		testCase("primer_mispriming_th_input",
				"primer_mispriming_th_input",
				"primer_mispriming_th_output");
	}	
	@Test
	public void testP3MisprimingBoundary1() throws FileNotFoundException, IOException {	  
		
		testCase("primer_mispriming_boundary1_input",
				"primer_mispriming_boundary1_input",
				"primer_mispriming_boundary1_output");
	}
	@Test
	public void testP3MisprimingBoundary2() throws FileNotFoundException, IOException {	  
		
		testCase("primer_mispriming_boundary2_input",
				"primer_mispriming_boundary2_input",
				"primer_mispriming_boundary2_output");
	}	
	@Test
	public void testP3MisprimingLongLib() throws FileNotFoundException, IOException {	  
		
		testCase("primer_mispriming_long_lib_input",
				"primer_mispriming_long_lib_input",
				"primer_mispriming_long_lib_output");
	}	
	
	//@Ignore("Ignore Not Implmented Now")
	@Test
	public void testP3ObjFnInput() throws FileNotFoundException, IOException {	  
		
		testCase("primer_obj_fn_input",
				"primer_obj_fn_input",
				"primer_obj_fn_output");
	}
	@Ignore("Ignore Not Implmented Now")
	@Test
	public void testP3LibAmbCodesInput() throws FileNotFoundException, IOException {	  
		
		testCase("primer_lib_amb_codes_input",
				"primer_lib_amb_codes_input",
				"primer_lib_amb_codes_output");
	}
}
