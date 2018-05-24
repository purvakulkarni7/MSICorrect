/**
 * 
 */

package de.uni_jena.bioinformatik.input;

import java.io.*;
import java.util.ArrayList;

/**
 * ReadFile.
 * <h3>Usage</h3>
 * <ol>
 * <li>Reads file and returns data in the file as an ArrayList {@code fileInput}.</li>
 *  
 * </ol>
 * 
 * @author Purva Kulkarni
 * @Date 22 May 2014
 *
 */

public class ReadFile {

	/**
	 * Method which takes filename as input and returns mz and intensity values as ArrayList
	 * @param String fileName
	 * @return ArrayList<String>
	 */	
	
	public ArrayList<String> fileInput(String fileName) {

		// Initialize all variables
		String ch = null;

		ArrayList<String> fileData = new ArrayList<String>();
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(fileName)); // Take file name using BufferedReader
			while((ch = br.readLine()) != null) {
				if(ch.startsWith("0"))
					continue;
				else
					fileData.add(ch);
				}
			br.close();
				
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
			
		return fileData;
	}
}