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
 * <li>Contains a method called generateCrystal</li>
 * <li>Finds the largest pixelScore for the complete crystal block and use it to build the crystal iteratively</li>
 * <li>When there is no other pixel untraversed a consensus pixel, a new crystal is build up</li>
 *  <li>Generation of RankingList goes on simultaneously</li>
 *  <li>The method returns a RankedList</li>
 *
 * </ol>
 *
 * @author Purva Kulkarni
 * @Date 26 May 2015
 *
 */


public class CrystalGrowthOrdering {

    public TreeMap<Integer, String> generateCrystal(ArrayList<String> fileNameList,
                                                    ArrayList<String> coordinateList, ArrayList<Float> pixelScoreList,
                                                    ArrayList<Object> spectrumArrayList, final int xMaximum,
                                                    final int yMaximum, float epsilon) throws IOException {
        float largestScore = 0;
        int startingIndex = 0;

        // Find the pixel with the largest score
        for (int i = 0; i < pixelScoreList.size(); i++) {
            if (largestScore < pixelScoreList.get(i)) {
                largestScore = pixelScoreList.get(i);
                startingIndex = i;
            }
        }

        // Create a pvMax object which contains information on the starting pixel which will be the temporary consensus
        PixelValues pvMax = new PixelValues(coordinateList.get(startingIndex), (float[]) spectrumArrayList.get(startingIndex), fileNameList.get(startingIndex), largestScore);


        // Step 1: Create method for crystal growth to start traversing
        // Step 2: Generate a ranking list based on the order of traversing
        // Step 3: Provide this ranking list to traverseSortedSet to generate a consensus spectrum

        // Create a TreeMap called rankList which will store the rank and the name of the file
        TreeMap<Integer, String> rankList = new TreeMap<Integer, String>();

        // Coordinates of the starting pixel assigned to consensusIcoordinate
        String firstCoordinate = pvMax.coordinates;
        String recalibrateJcoordinate;

        // Create an arraylist which store the traversed coordinates
        ArrayList<String> traversedCoordinates = new ArrayList<String>();
        float topPixelScore;
        float rightPixelScore;
        float bottomPixelScore;
        float leftPixelScore;
        ArrayList<Float> adjacentPixelArrayList = new ArrayList<Float>(); //ArrayList in which the adjacent pixel score will be stored
        ArrayList<String> adjacentCoordinateArrayList = new ArrayList<String>(); //ArrayList in which the adjacent pixel score will be stored
        ArrayList<String> notYetTraversed = new ArrayList<String>();
        de.uni_jena.bioinformatik.functions.CalculateDistance cd1 = new de.uni_jena.bioinformatik.functions.CalculateDistance();
        String nextPixelTraversed = null;
        int iteration = 1;
        int crystalNumber = 0;

        for(int g = 0; g<coordinateList.size(); g++)
            notYetTraversed.add(coordinateList.get(g));



        rankList.put(iteration, fileNameList.get(coordinateList.indexOf(firstCoordinate)));
        traversedCoordinates.add(firstCoordinate);
        notYetTraversed.remove(firstCoordinate);

        while (traversedCoordinates.size() != coordinateList.size())
        {
            adjacentPixelArrayList.clear();
            adjacentCoordinateArrayList.clear();
            ArrayList<Float> adjacentDistanceArrayList = new ArrayList<Float>();
            MatchListDataStorage mdsObjectTemp;

            for (int e = 0; e < traversedCoordinates.size(); e++) {
                //   System.out.println("Printing non-traversed neighbors for :" + traversedCoordinates.get(e));
                String temp[] = traversedCoordinates.get(e).split("\t");
                int pixelxPos = Integer.valueOf(temp[0]);
                int pixelyPos = Integer.valueOf(temp[1]);


                // Get the mzArray and fileName for this pixel which is a part of the crystal
                String currentCrystalPixelCoordinate = pixelxPos + "\t" + pixelyPos;
                float[] CrystalPixelArray = (float[]) spectrumArrayList.get(coordinateList.indexOf(currentCrystalPixelCoordinate));
                String currentCrystalPixelFileName = fileNameList.get(coordinateList.indexOf(currentCrystalPixelCoordinate));


                String topP = (pixelxPos - 1) + "\t" + (pixelyPos);
                if (!traversedCoordinates.contains(topP) && coordinateList.contains(topP)) {

                    adjacentCoordinateArrayList.add(topP);

                    // Get the mzArray and fileName for topP
                    float[] topPPixelArray = (float[]) spectrumArrayList.get(coordinateList.indexOf(topP));
                    String topPPixelFileName = fileNameList.get(coordinateList.indexOf(topP));

                    // Calculate distance of topP with its neighboring pixel (which is a part of the crystal) and save it in a distance array
                    mdsObjectTemp = cd1.specDistance(CrystalPixelArray, topPPixelArray, currentCrystalPixelCoordinate, topP, currentCrystalPixelFileName, topPPixelFileName, epsilon);
                    adjacentDistanceArrayList.add(mdsObjectTemp.distanceValue);
                }

                String bottomP = (pixelxPos + 1) + "\t" + (pixelyPos);
                if (!traversedCoordinates.contains(bottomP) && coordinateList.contains(bottomP)){

                    adjacentCoordinateArrayList.add(bottomP);

                    // Get the mzArray and fileName for topP
                    float[] bottomPPixelArray = (float[]) spectrumArrayList.get(coordinateList.indexOf(bottomP));
                    String bottomPPixelFileName = fileNameList.get(coordinateList.indexOf(bottomP));

                    // Calculate distance of topP with its neighboring pixel (which is a part of the crystal) and save it in a distance array
                    mdsObjectTemp = cd1.specDistance(CrystalPixelArray, bottomPPixelArray, currentCrystalPixelCoordinate, bottomP, currentCrystalPixelFileName, bottomPPixelFileName, epsilon);
                    adjacentDistanceArrayList.add(mdsObjectTemp.distanceValue);
                }

                String leftP = (pixelxPos) + "\t" + (pixelyPos - 1);
                if (!traversedCoordinates.contains(leftP) && coordinateList.contains(leftP)) {

                    adjacentCoordinateArrayList.add(leftP);

                    // Get the mzArray and fileName for topP
                    float[] leftPPixelArray = (float[]) spectrumArrayList.get(coordinateList.indexOf(leftP));
                    String leftPPixelFileName = fileNameList.get(coordinateList.indexOf(leftP));

                    // Calculate distance of topP with its neighboring pixel (which is a part of the crystal) and save it in a distance array
                    mdsObjectTemp = cd1.specDistance(CrystalPixelArray, leftPPixelArray, currentCrystalPixelCoordinate, leftP, currentCrystalPixelFileName, leftPPixelFileName, epsilon);
                    adjacentDistanceArrayList.add(mdsObjectTemp.distanceValue);
                }

                String rightP = (pixelxPos) + "\t" + (pixelyPos + 1);
                if (!traversedCoordinates.contains(rightP) && coordinateList.contains(rightP) ){

                    adjacentCoordinateArrayList.add(rightP);

                    // Get the mzArray and fileName for topP
                    float[] rightPPixelArray = (float[]) spectrumArrayList.get(coordinateList.indexOf(rightP));
                    String rightPPixelFileName = fileNameList.get(coordinateList.indexOf(rightP));

                    // Calculate distance of topP with its neighboring pixel (which is a part of the crystal) and save it in a distance array
                    mdsObjectTemp = cd1.specDistance(CrystalPixelArray, rightPPixelArray, currentCrystalPixelCoordinate, rightP, currentCrystalPixelFileName, rightPPixelFileName, epsilon);
                    adjacentDistanceArrayList.add(mdsObjectTemp.distanceValue);
                }
            }

            //In case when the size of distance array list is zero (happens when all pixels are traversed or there are some disconnected pixels)
            if(adjacentDistanceArrayList.size() == 0)
            {
                if(notYetTraversed.size() != 0 )
                {
                    ArrayList<Float> notTravsersedPixelScores = new ArrayList<Float>();

                    // Finding pixel scores for all the non travsersed pixels
                    for(int l = 0; l< notYetTraversed.size(); l++)
                    {
                        notTravsersedPixelScores.add(pixelScoreList.get(coordinateList.indexOf(notYetTraversed.get(l))));
                    }

                    // Find the maximum pixel score
                    float maxPixelScoreNonTraversed =0;
                    for(int q = 0; q< notTravsersedPixelScores.size(); q++)
                    {
                        if(notTravsersedPixelScores.get(q) > maxPixelScoreNonTraversed)
                            maxPixelScoreNonTraversed = notTravsersedPixelScores.get(q);
                    }

                    String newCrystalInitiationPixel = notYetTraversed.get(notTravsersedPixelScores.indexOf(maxPixelScoreNonTraversed));
                    iteration++;
                    rankList.put(iteration, fileNameList.get(coordinateList.indexOf(newCrystalInitiationPixel)));
                    traversedCoordinates.add(newCrystalInitiationPixel);
                    notYetTraversed.remove(newCrystalInitiationPixel);
                    continue;
                }
            }

            float maxDistanceValue = adjacentDistanceArrayList.get(0);

            //Find maximum distance from the adjacent distance array
            for (int f = 0; f < adjacentDistanceArrayList.size() - 1; f++) {
                if (adjacentDistanceArrayList.get(f) > maxDistanceValue)
                    maxDistanceValue = adjacentDistanceArrayList.get(f);
            }

            recalibrateJcoordinate = adjacentCoordinateArrayList.get(adjacentDistanceArrayList.indexOf(maxDistanceValue));
            iteration++;
            rankList.put(iteration, fileNameList.get(coordinateList.indexOf(recalibrateJcoordinate)));

            traversedCoordinates.add(recalibrateJcoordinate);
            notYetTraversed.remove(recalibrateJcoordinate);
        }
        return(rankList);
    }

}
