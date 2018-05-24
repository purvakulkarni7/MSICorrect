/**
 * 
 */

package de.uni_jena.bioinformatik.functions;

import de.uni_jena.bioinformatik.input.MatchListDataStorage;
import de.uni_jena.bioinformatik.input.PixelValues;

import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeMap;


/**
 * TopologicalGreedyOrdering.
 * <h3>Usage</h3>
 * <ol>
 * <li>Contains a method called greedyOrdering</li>
 * <li>Finds the largest pixelScore and uses it to build the crystal iteratively</li>
 * <li>When there is no other pixel un-traversed a consensus pixel, a new crystal is build up</li>
 *  <li>Generation of RankingList goes on simultaneously</li>
 *  <li>The method returns a RankedList</li>
 *
 * </ol>
 * 
 * @author Purva Kulkarni
 * @Date 2 June 2014
 *
 */


public class TopologicalGreedyOrdering {

	public TreeMap<Integer, String> greedyOrdering(ArrayList<String> fileNameList,
                                                   ArrayList<String> coordinateList, ArrayList<Float> pixelScoreList,
                                                   ArrayList<Object> spectrumArrayList, final int xMaximum,
                                                   final int yMaximum, float epsilon) throws IOException {
		float largestScore = 0;
		int startingIndex = 0;

		// Find the pixel with the largest score
		for(int i=0; i<pixelScoreList.size();i++)
		{
			if(largestScore < pixelScoreList.get(i))
			{
				largestScore = pixelScoreList.get(i);
				startingIndex = i;
			}
		}
//        System.out.println("Starting index: " + startingIndex);
//		System.out.println("Largest value in pixelScore ArrayList: " + largestScore + " at position " + coordinateList.get(startingIndex) + " Filename: " + fileNameList.get(startingIndex));
//		System.out.println("-----------------------------------------------------------");

        // Create a pvMax object which contains information on the starting pixel which will be the temporary consensus
		PixelValues pvMax = new PixelValues(coordinateList.get(startingIndex), (float[])spectrumArrayList.get(startingIndex), fileNameList.get(startingIndex), largestScore);


		// Step 1: Create method for crystal growth to start traversing
		// Step 2: Generate a ranking list based on the order of traversing
		// Step 3: Provide this ranking list to traverseSortedSet to generate a consensus spectrum

		// Create a TreeMap called rankList which will store the rank and the name of the file
		TreeMap<Integer, String> rankList = new TreeMap<Integer, String>();

        // Coordinates of the starting pixel assigned to consensusIcoordinate
		String consensusIcoordinate = pvMax.coordinates;
		String recalibrateJcoordinate;

        // Create an arraylist which store the traversed coordinates
		ArrayList<String> traversedCoordinates = new ArrayList<String>();
		float topPixelScore;
		float rightPixelScore;
		float bottomPixelScore;
		float leftPixelScore;
		float[] adjacentPixelArray = new float[4]; //Array in which the adjacent pixel score will be stored

		int iteration = 1;
		int crystalNumber = 0;

		// for loop to traverse the pixels and generate a ranking list
		for(int k  = coordinateList.size(); k > 0 ; k--)
		{
			// Split the coordinate string into individual x and y positions
            MatchListDataStorage mdsObject;
            de.uni_jena.bioinformatik.functions.CalculateDistance cd1 = new de.uni_jena.bioinformatik.functions.CalculateDistance();
            int consensusICoordinateIndex = coordinateList.indexOf(consensusIcoordinate); //get index of the coordinate position

            Object mzObject = spectrumArrayList.get(consensusICoordinateIndex);
            float[] consensusI = (float[])mzObject;
            String coordinateOfI = coordinateList.get(consensusICoordinateIndex);
            String FileNameOfI = fileNameList.get(consensusICoordinateIndex);

            String temp1[] = consensusIcoordinate.split("\t");
			int consensusxPos = Integer.valueOf(temp1[0]);
			int consensusyPos = Integer.valueOf(temp1[1]);

			//get pixel scores for adjoining pixels
			//1. Top pixel x-1,y
			String topPixel = (consensusxPos - 1) + "\t" + (consensusyPos);

			// check if there is a pixel which does not exist
			if(((consensusxPos-1) == 0) || (consensusyPos == 0) || ((consensusxPos-1) > xMaximum) || (consensusyPos > yMaximum) || (!(coordinateList.contains(topPixel))))
			{
				topPixelScore = 0;
				adjacentPixelArray[0] = topPixelScore;
			}					
			else
			{
				if(!(traversedCoordinates.contains(topPixel)))
				{
             		int topPixelIndex = coordinateList.indexOf(topPixel); //get index of the coordinate position

                    mzObject = spectrumArrayList.get(topPixelIndex);
                    float[] topPixelArray = (float[])mzObject;
                    String topPixelCoordinate = coordinateList.get(topPixelIndex);
                    String topPixelFileName = fileNameList.get(topPixelIndex);

                    mdsObject= cd1.specDistance(consensusI, topPixelArray, coordinateOfI, topPixelCoordinate, FileNameOfI, topPixelFileName, epsilon);
                    topPixelScore = mdsObject.distanceValue;

               //  topPixelScore = pixelScoreList.get(topPixelIndex); //get PixelScore //commented on 13 May 2015
					adjacentPixelArray[0] = topPixelScore;
				}
				else
				{
					topPixelScore = 0;
					adjacentPixelArray[0] = topPixelScore;
				}
			}

			//2. Right pixel x,y+1
			String rightPixel = (consensusxPos) + "\t" + (consensusyPos + 1);

			// check if there is a pixel which does not exist
			if((consensusxPos == 0) || ((consensusyPos + 1) == 0) || (consensusxPos > xMaximum) || ((consensusyPos+1) > yMaximum) || (!(coordinateList.contains(rightPixel))))
			{
				rightPixelScore = 0;
				adjacentPixelArray[1] = rightPixelScore;
			}					
			else
			{
				if(!(traversedCoordinates.contains(rightPixel)))
				{
					int rightPixelIndex = coordinateList.indexOf(rightPixel); //get index of the coordinate position

                    mzObject = spectrumArrayList.get(rightPixelIndex);
                    float[] rightPixelArray = (float[])mzObject;
                    String rightPixelCoordinate = coordinateList.get(rightPixelIndex);
                    String rightPixelFileName = fileNameList.get(rightPixelIndex);

                    mdsObject= cd1.specDistance(consensusI, rightPixelArray, coordinateOfI, rightPixelCoordinate, FileNameOfI, rightPixelFileName, epsilon);
                    rightPixelScore = mdsObject.distanceValue;

					adjacentPixelArray[1] = rightPixelScore;
				}
				else
				{
					rightPixelScore = 0;
					adjacentPixelArray[1] = rightPixelScore;
				}
			}

			//3. Bottom pixel x+1,y
			String bottomPixel = (consensusxPos + 1) + "\t" + (consensusyPos);

			// check if there is a pixel which does not exist
			if(((consensusxPos+1) == 0) || (consensusyPos == 0) || ((consensusxPos+1) > xMaximum) || (consensusyPos > yMaximum) || (!(coordinateList.contains(bottomPixel))))
			{
				bottomPixelScore = 0;
				adjacentPixelArray[2] = bottomPixelScore;
			}					
			else
			{
				if(!(traversedCoordinates.contains(bottomPixel)))
				{
					int bottomPixelIndex = coordinateList.indexOf(bottomPixel); //get index of the coordinate position

                    mzObject = spectrumArrayList.get(bottomPixelIndex);
                    float[] bottomPixelArray = (float[])mzObject;
                    String bottomPixelCoordinate = coordinateList.get(bottomPixelIndex);
                    String bottomPixelFileName = fileNameList.get(bottomPixelIndex);

                    mdsObject= cd1.specDistance(consensusI, bottomPixelArray, coordinateOfI, bottomPixelCoordinate, FileNameOfI, bottomPixelFileName, epsilon);
                    bottomPixelScore = mdsObject.distanceValue;

					adjacentPixelArray[2] = bottomPixelScore;
				}
				else
				{
					bottomPixelScore = 0;
					adjacentPixelArray[2] = bottomPixelScore;
				}
			}

			//4. Left pixel x,y-1
			String leftPixel = (consensusxPos) + "\t" + (consensusyPos-1);

			if((consensusxPos == 0) || ((consensusyPos-1) == 0) || (consensusxPos > xMaximum) || ((consensusyPos-1) > yMaximum) || (!(coordinateList.contains(leftPixel))))
			{
				leftPixelScore = 0;
				adjacentPixelArray[3] = leftPixelScore;
			}					
			else
			{
				if(!(traversedCoordinates.contains(leftPixel)))
				{
					int leftPixelIndex = coordinateList.indexOf(leftPixel); //get index of the coordinate position

                    mzObject = spectrumArrayList.get(leftPixelIndex);
                    float[] leftPixelArray = (float[])mzObject;
                    String leftPixelCoordinate = coordinateList.get(leftPixelIndex);
                    String leftPixelFileName = fileNameList.get(leftPixelIndex);

                    mdsObject= cd1.specDistance(consensusI, leftPixelArray, coordinateOfI, leftPixelCoordinate, FileNameOfI, leftPixelFileName, epsilon);
                    leftPixelScore = mdsObject.distanceValue;

					adjacentPixelArray[3] = leftPixelScore;
				}
				else
				{
					leftPixelScore = 0;
					adjacentPixelArray[3] = leftPixelScore;
				}
			}

			float maxAdjacentPixelScore = 0;
			int i;
			int count=0;

			for(i=0;i<adjacentPixelArray.length;i++)
			{
				if(adjacentPixelArray[i] > maxAdjacentPixelScore)
				{
					maxAdjacentPixelScore = adjacentPixelArray[i];
					count = i;
				}
			}

			// In case if all adjacent pixel scores are zero
			// In this case a new pixel with highest score has to be found
			float largestScoreforNewCrystal = 0;
			int startingIndexTemp = 0;
			String nextPixelTraversed = null;

			if (topPixelScore == 0 && rightPixelScore == 0 && bottomPixelScore == 0 && leftPixelScore == 0)
			{
				maxAdjacentPixelScore = 0;
				String lastPixel = (consensusxPos) + "\t" + (consensusyPos);
				//			System.out.println(lastPixel);
				nextPixelTraversed = lastPixel;

             	rankList.put(iteration, fileNameList.get(coordinateList.indexOf(nextPixelTraversed)));
				traversedCoordinates.add(nextPixelTraversed);
				iteration = iteration+1;
				crystalNumber = crystalNumber + 1;


				// Find the pixel with the largest score
				for(int l=0; l<coordinateList.size()-1;l++)
				{
					if(traversedCoordinates.contains(coordinateList.get(l)))
						continue;
					else
					{
						if(largestScoreforNewCrystal < pixelScoreList.get(l))
						{
							largestScoreforNewCrystal = pixelScoreList.get(l);
							startingIndexTemp = l;
						}
					}					
				}
                pvMax = new PixelValues(coordinateList.get(startingIndexTemp), (float[])spectrumArrayList.get(startingIndexTemp), fileNameList.get(startingIndexTemp), largestScoreforNewCrystal);
				consensusIcoordinate = pvMax.coordinates;	
				continue;
			}

			// Find which is the pixel position
			// String nextPixelTraversed = null;
			if(count == 0)
				nextPixelTraversed = topPixel;
			else
				if (count==1)
					nextPixelTraversed = rightPixel;
				else
					if (count==2)
						nextPixelTraversed = bottomPixel;
					else
						nextPixelTraversed = leftPixel;

			rankList.put(iteration, fileNameList.get(coordinateList.indexOf(consensusIcoordinate)));

			//Assign the next pixel traversed to the recalibrateJcoordinate
            recalibrateJcoordinate = nextPixelTraversed;

			//add this to the rankinglist
         	rankList.put(iteration+1, fileNameList.get(coordinateList.indexOf(recalibrateJcoordinate)));

			//increase the iteration
			iteration++;

			// Add consensus to traversedlist
			traversedCoordinates.add(consensusIcoordinate);

			// Add spectra to recalibrate in the traversedlist
			traversedCoordinates.add(recalibrateJcoordinate);

			consensusIcoordinate = recalibrateJcoordinate;
		}
		return(rankList);
	}

}
