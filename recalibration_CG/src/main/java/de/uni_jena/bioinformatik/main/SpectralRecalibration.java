/**
 *
 */

package de.uni_jena.bioinformatik.main;

import de.uni_jena.bioinformatik.input.InputData;
import org.rosuda.REngine.REXPMismatchException;
import org.rosuda.REngine.REngineException;

import java.io.IOException;

/**
 * SpectralRecalibration.
 * <h3>Usage</h3>
 * <ol>
 * <li>Class containing main method.</li>
 * <li>Creates an InputData object {@code InputData} and calls the userInput {@code userInput}
 * to ask for user input</li>
 * </ol>
 *
 * @author Purva Kulkarni
 * @Date 21 May 2014
 *
 */


public class SpectralRecalibration {

	// main method begins execution of java application
	public static void main(String[] args) throws IOException,REXPMismatchException, REngineException {

		//Start time to see how much time does the program take for a specific data set
		long startTime = System.currentTimeMillis();

		InputData f1 = new InputData();

		if(args.length == 0)
		{
			f1.userInput();
			long endTime = System.currentTimeMillis(); //end the counting of time elapsed
		}
		else if(args.length == 1) {

			String help = args[0];
			if(new String("--help").equals(help))
			{
				f1.welcomeMessage();
				f1.printhelp();
			}
			else
			{
				System.out.println("Wrong argument passed to MSICorrect");
			}
		}

		System.out.println("------------- * MSICorrect Exits * -------------");
	} //end main method

} // end class