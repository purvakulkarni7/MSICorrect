/**
 *
 */

package de.uni_jena.bioinformatik.functions;

import de.uni_jena.bioinformatik.input.MatchListDataStorage;

import java.io.IOException;
import java.util.ArrayList;

/**
 * CalculateDistance.
 * <h3>Usage</h3>
 * <ol>
 * <li>Checks for empty m/z lists {@code emptySpectra} </li>
 * <li>Compute the distance between 2 mzList {@code specDistance}</li>
 * </ol>
 *
 * @author Purva Kulkarni
 * @Date 21 May 2014
 *
 */

public class CalculateDistance {

	/**
	 * Checks if the specific mz spectra is empty, 
	 * returns false if the spectrum contains a value other than zero; or else returns true
	 *
	 * @param mzList
	 * @return boolean specValue 
	 */

	public boolean emptySpectra(float[] mzList) throws IOException
	{
		// orginal value of specValue is true; meaning the spectra is empty
		boolean specValue = true;
		int count = 0;
		for(int temp = 0;temp < mzList.length; temp++)
		{
			if (mzList[temp]!= 0)
			{
				specValue = false; // this means that the mzList contains a non-zero values hence it is not empty
				count++;
				if(count==2)
				{
					specValue = false;
					break;
				}
				else
				{
					specValue = true;
				}

			}
			else
				specValue = true;
		}
		return specValue;
	}


	/**
	 * Calculates the distance between two mzLists and returns a MatchListDataStorage object {@link MatchListDataStorage}
	 * The MatchListDataStorage object contains matchListCoordinates, distance as fields 
	 *
	 * @param consensusI, referenceJ, String coordinateI, String coordinateJ, String fileNameOfI, String fileNameOfJ
	 * @return MatchListDataStorage
	 */

	public MatchListDataStorage specDistance(float[] consensusI, float[] referenceJ, String coordinateI, String coordinateJ, String fileNameOfI, String fileNameOfJ, float epsilon) throws IOException
	{
		ArrayList<Float> matchListI = new ArrayList<Float>();
		ArrayList<Float> matchListJ = new ArrayList<Float>();

		ArrayList<Integer> matchListindexI = new ArrayList<Integer>();
		ArrayList<Integer> matchListindexJ = new ArrayList<Integer>();


		// Loop to generate match list
		for(int i=0; i<consensusI.length;i++)
		{
			//Checks if value is non-zero
			if(consensusI[i] > 0)
			{
				for(int j=0; j<referenceJ.length;j++)
				{
					//Checks if value is non-zero
					if(referenceJ[j] > 0)
					{
						// Checks if difference between consensusI[i]-referenceJ[j] is zero is less than epsilon 
						if(((consensusI[i]-referenceJ[j]) == 0) || (Math.abs(consensusI[i]-referenceJ[j])<= epsilon))
						{
							//If yes,
							matchListJ.add(referenceJ[j]); // adds value of referenceJ[j] in matchListJ
							matchListindexJ.add(j); // corresponding index in matchListindexJ
							matchListI.add(consensusI[i]); // adds value of consensusI[i] in matchListI
							matchListindexI.add(i); // corresponding index in matchListindexI
						}
						else
							//If the difference if greater than epsilon
							if(Math.abs(consensusI[i]-referenceJ[j])> epsilon)
							{
								continue; //Iterates to the next j value
							}
					}
					//If value of j is zero or less
					else
						continue; //Iterates to the next j value
				}	// End of j for loop

			}
			//If value of i is zero or less
			else
				continue;  //Iterates to the next i value
		} // End of i for loop


		String matchListCoordinates = coordinateI + "\t" + fileNameOfI +"\t" + coordinateJ + "\t" + fileNameOfJ;

		// Loop to iterate through matchListI to remove redundant values
		for(int k=0; k<(matchListI.size()-1); k++)
		{
			// Calculate the difference between matchList value k and k+1
			float sub1 = matchListI.get(k)-matchListI.get(k+1);

			// If difference is not equal to zero - means the two values are not equal
			if(sub1 != 0)
			{
				continue; // iterate to the next value of k

			}
			else
			// If the two values are equal
			{
				//Calculate the difference between corresponding index values in matchList J 
				float sub2 = matchListJ.get(k)-matchListJ.get(k+1);

				// If the difference is less than zero - means k+1 value is greater and hence has to be removed
				if(sub2 < 0)
				{
					matchListI.remove(k+1); // remove value at index k+1 and its index from matchList I 
					matchListindexI.remove(k+1);

					matchListJ.remove(k+1); // remove value at index k+1 and its index from matchList J
					matchListindexJ.remove(k+1);
				}

			}
		}

		int peakCount = matchListI.size(); // Count peaks
		float massRange=0;
		float distance = 0;

		// If peak Count is greater than 2
		if(peakCount > 2)
		{
			//Calculate the mass Range = first massValue - last massValue
			massRange = matchListI.get(peakCount-1)-matchListI.get(0);
		}
		else
		// If peak Count is less than 2 i.e it is 2, 1 or 0
		{
			massRange = 0; // Assign the mass range as 0.
		}


		// If peak Count is greater than 2
		if(peakCount > 2)
		{
			//Calculate the distance
			distance = (1/(peakCount*massRange));
			MatchListDataStorage mds = new MatchListDataStorage(matchListCoordinates, distance);
			return mds;
		}
		else
			distance = 1;
		MatchListDataStorage mds = new MatchListDataStorage(matchListCoordinates, distance);
		return mds;
	}
}