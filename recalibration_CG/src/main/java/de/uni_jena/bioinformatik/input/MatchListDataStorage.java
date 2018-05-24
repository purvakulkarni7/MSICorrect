/**
 * 
 */

package de.uni_jena.bioinformatik.input;

/**
 * MatchListDataStorage.
 * <h3>Usage</h3>
 * <ol>
 * <li>Class acting like a structure with 2 fields -</li>
 * <li>String twoFileCoordinates,float distanceValue </li>
 * </ol>
 * 
 * @author Purva Kulkarni
 * @Date 22 May 2014
 *
 */

public class MatchListDataStorage {
	public String twoFileCoordinates;
	public float distanceValue;
	
	public MatchListDataStorage(String twoFileCoordinates, float distanceValue){
		this.twoFileCoordinates = twoFileCoordinates;
		this.distanceValue = distanceValue;
	}
}