/**
 * 
 */
package de.uni_jena.bioinformatik.input;

import java.io.File;
import java.util.Arrays;
import java.util.Comparator;

/**
 * SortFile.
 * <h3>Usage</h3>
 * <ol>
 * <li>Sorts the file according to the filenames (numeric sort).</li>
 *  
 * </ol>
 * 
 * @author Purva Kulkarni
 * @Date 21 May 2014
 *
 */

public class SortFile {

	/**
	 * Method which takes list of files as a fileArray and sorts them 
	 * 
	 * @param File[] files
	 * @return File[] files
	 */	
	
	public File[] sortByNumber(File[] files) {

		Arrays.sort(files, new Comparator<File>() {
			//@Override
			public int compare(File o1, File o2) {
				int n1 = extractNumber(o1.getName());
				int n2 = extractNumber(o2.getName());
				return n1 - n2;
			}

			private int extractNumber(String name) {
				int i = 0;
				try {
					int s = name.indexOf('_')+1;
					int e = name.lastIndexOf('.');
					String number = name.substring(s, e);
					i = Integer.parseInt(number);
				} catch(Exception e) {
					i = 0; // if filename does not match the format
					// then default to 0
				}
				return i;
			}
		});

		return files;

	}

}