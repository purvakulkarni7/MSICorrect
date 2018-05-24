/**
 * 
 */
package de.uni_jena.bioinformatik.input;

/**
 * PixelValues.
 * <h3>Usage</h3>
 * <ol>
 * <li>Class acting like a structure with 4 fields -</li>
 * <li>String coordinates,float[] mzArray, float[] intensityArray, String fileName, float pixelScore </li>
 * </ol>
 * 
 * @author Purva Kulkarni
 * @Date 31 May 2014
 *
 */

public class PixelValues {
	
	public String coordinates;
	public float[] mzArray;
	public String fileName;
	public float pixelScore;
	
	public PixelValues(String coordinates, float[] mzData, String fileName, float pixelScore ){
		this.coordinates = coordinates;
		this.mzArray = mzData;
		this.fileName = fileName;
		this.pixelScore = pixelScore;
	}
}
