/**
 * 
 */
package de.uni_jena.bioinformatik.input;

/**
 * RawDataStorage.
 * <h3>Usage</h3>
 * <ol>
 * <li>Class acting like a structure with 4 fields -</li>
 * <li>String coordinates,float[] mzArray, float[] intensityArray, String fileName </li>
 * </ol>
 * 
 * @author Purva Kulkarni
 * @Date 21 May 2014
 *
 */

public class RawDataStorage {
	
	public String coordinates;
	public float[] mzArray;
	public float[] intensityArray;
	public String fileName;
	
	RawDataStorage(String coordinates, float[] mzData, float[] intData, String fileName ){
		this.coordinates = coordinates;
		this.mzArray = mzData;
		this.intensityArray = intData;
		this.fileName = fileName;
	}
}
