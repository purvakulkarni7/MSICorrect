/**
 *
 */
package de.uni_jena.bioinformatik.input;
import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Scanner;
//import java.util.Map;
//import java.util.Scanner;
//import java.util.Set;
import java.util.TreeMap;

// import de.uni_jena.bioinformatik.functions.CalculateDistance;
import de.uni_jena.bioinformatik.functions.CrystalGrowthOrdering;
import de.uni_jena.bioinformatik.functions.TopologicalGreedyOrdering;
import de.uni_jena.bioinformatik.functions.TraverseSparseSet;
//import de.uni_jena.bioinformatik.functions.PixelScore;
//import de.uni_jena.bio.informatik.functions.MatrixFunctions;
//import de.uni_jena.bio.informatik.functions.SparseMatrix;
//import de.uni_jena.bio.informatik.functions.TraverseSparseSet;
import de.unijena.bioinf.ChemistryBase.ms.MutableSpectrum;
import de.unijena.bioinf.ChemistryBase.ms.Peak;
import de.unijena.bioinf.ChemistryBase.ms.utils.SimpleSpectrum;

import org.rosuda.REngine.REXP;
import org.rosuda.REngine.REXPMismatchException;
import org.rosuda.REngine.REngineException;
import org.rosuda.REngine.Rserve.RConnection;
import org.rosuda.REngine.Rserve.RserveException;
/**
 * InputData.
 * <h3>Usage</h3>
 * <ol>
 * <li>Accepts user input in the form of directory path {@code userInput}.
 *  Sorts the files in the directory </li>
 * <li>Creates an InputData object {@code InputData} and calls the userInput {@code userInput}
 * to ask for user input</li>
 * <li>Extracts mz and intensity values as array {@code extractMzIntensity} and checks for empty spectra</li>
 *  <li>Extracts mz and intensity values as array {@code extractMzIntensity}</li>
 *  <li> Generates match list using 2 spectra as input </li>
 *  <li>Calculates distance between two mz lists {@code specDistance}</li>
 *  <li>Transposes the distance matrix {@code transposeFloat}</li>
 *  <li>Generates a sparse matrix of this distance matrix {@code generateSparseMatrix}</li>
 *  <li>Sparse matrix is traversed to perform recalibration
 *  starting from the two mz lists having the minimum distance {@code traverseSortedSet}</li>
 *
 * </ol>
 *
 * @author Purva Kulkarni
 * @Date 21 May 2014
 *
 */
public class InputData {
    static String dirPath;
    static String preprocessedDirPath;
    static String recalibratedDirPath;
    static String mzRange;
    static String userInput;
    static float mzStart;
    static float mzEnd;
    public int xMaximum = 0;
    public int yMaximum = 0;

    public void userInput() throws IOException,REXPMismatchException, REngineException
    {
        welcomeMessage();

        System.out.println("MSICorrect needs the following user input");
        System.out.println("---");
        //For user input
        Scanner scanner = new Scanner(System.in );
        System.out.print("File path for input mass spectra: ");
        dirPath = scanner.nextLine(); // Takes the directory path as the user input

        Scanner scanner1 = new Scanner(System.in );
        System.out.print("File path to store recalibrated spectra: ");
        String pathName = scanner1.nextLine();

        Scanner scanner2 = new Scanner(System.in );
        System.out.print("Mass threshold for distance calculation between two peaklists: ");
        String value = scanner2.nextLine();
        float epsilon = Float.valueOf(value);

        Scanner scanner5 = new Scanner(System.in );
        System.out.print("Peaklist ordering approach (use TG or CG): ");
        String scanner5Input = scanner5.nextLine();

        Scanner scanner3 = new Scanner(System.in );
        System.out.print("Line distance to perform line-pair stabbing recalibration: ");
        String scanner3Input = scanner3.nextLine();
        double epsilon2 = Double.parseDouble(scanner3Input);

        Scanner scanner4 = new Scanner(System.in );
        System.out.print("Mass threshold for line-pair stabbing recalibration: ");
        String scanner4Input = scanner4.nextLine();
        double threshold = Double.parseDouble(scanner4Input);

        System.out.println("---");

        System.out.println("\nSTEP 1 of 6: Reading directory path and parsing data files");

        File folder = new File(dirPath);
        //checks whether the specified path contains a folder

        // Read the preprocessed files from the directory path
        if(folder.isDirectory())
        {
            File[] fileList = folder.listFiles();
            SortFile sf = new SortFile();

            // Calls sortByNumber method in class SortFile to list the files number wise
            fileList = sf.sortByNumber(fileList);

            //Creation of ArrayList to store individual data fields read from the file
            ArrayList<String> fileNameList = new ArrayList<String>();
            ArrayList<String> coordinateList = new ArrayList<String>();
            ArrayList<Float> totalIntensityList = new ArrayList<Float>();
            ArrayList<Float> newTotalIntensityList = new ArrayList<Float>();
            ArrayList<Float> averageSimilarityList = new ArrayList<Float>();
            ArrayList<Float> pixelScoreList = new ArrayList<Float>();
            ArrayList<Integer> positionList = new ArrayList<Integer>();
            ArrayList<Object> spectrumArrayList = new ArrayList<Object>();
            ArrayList<String> fileDataI = new ArrayList<String>();
            RawDataStorage dsr;

            // Create a new calculate distance object cd1
            de.uni_jena.bioinformatik.functions.CalculateDistance cd1 = new de.uni_jena.bioinformatik.functions.CalculateDistance();
            int fileCount = 0;

            int iterationCounter = 0;
            // Create a for loop to traverse the complete file List
            for(int i=0;i<(fileList.length);i++)
            {
                File spectrumFile = fileList[i];
                String fileNameI = new String();

                // Get file name and read the file
                fileNameI = folder + "/" + spectrumFile.getName();

                ReadFile rf = new ReadFile();
                // call the fileInput method and passe fileName to it
                fileDataI = rf.fileInput(fileNameI); //returns ArrayList stored in fileData


                // Call the listTo Array Method which converts the list containing coordinates + m/z and int columns to string array
                String[] fileDataArrayI = listToArray(fileDataI);

                float[] mzArray = new float[fileDataArrayI.length];
                float[] intensityArray = new float[fileDataArrayI.length];
                float[] mzList = new float[fileDataArrayI.length];


                // Returns a rawDataStorage Object
                dsr = extractMzIntensity(fileDataArrayI, mzArray, intensityArray, spectrumFile.getName());


                // Access the mzArray field and intensity in the dsr object and store in  a float array called mzList and intList
                mzList = dsr.mzArray;

                // Add all intensity values for a a single coordinate position (sum of intensity)
                float intensitySum = 0;

                // for loop to add all the intensity values at a specific coordinate position and store in intensitySum
                for(int j=0;j<dsr.intensityArray.length;j++)
                {
                    intensitySum = intensitySum + dsr.intensityArray[j];
                }
                totalIntensityList.add(intensitySum); //Add the total intensity at specific coordinate to totatlIntensity List

                // Get the total image dimensions by reading the coordinated of the last file
                // This will be the xMaximum and the yMaximum
                if(i==fileList.length-1) {
                    String lastCoordinateValues = dsr.coordinates;
                    String temp[] = lastCoordinateValues.split("\t");
                    //    System.out.println("last coordinate value: " + lastCoordinateValues);

                    // Extract the xMaximum and the yMaximum from the last Coordinate value
                    xMaximum = Integer.valueOf(temp[0]);
                    yMaximum = Integer.valueOf(temp[1]);
                }

                //Check if any of the spectra is empty
                if((cd1.emptySpectra(mzList))== true) // true means that the spectra is empty
                {
                    // System.out.println (i + " is empty");
                    fileCount++;
                    iterationCounter++;
                    continue;

                }
                else
                {
                    // Add the mzLists and its corresponding coordinate values to the respective Array Lists
                    iterationCounter++;
                    positionList.add(iterationCounter);
                    spectrumArrayList.add(mzList);
                    coordinateList.add(dsr.coordinates);
                    fileNameList.add(dsr.fileName);
                }
            }

            //Print the maximum coordinate value of the imaging data irrespective of empty spectra
            System.out.println("Total data coordinates: " + "x = " + xMaximum + " y = " + yMaximum);

            // Fill the new Total intensity list which is devoid of the sum of intensity for empty spectra
            for(int e = 0; e<positionList.size();e++)
            {
                newTotalIntensityList.add(totalIntensityList.get(positionList.get(e)-1));
            }

            System.out.println("There are " + fileCount + " empty peaklists in the provided mass lists.");

            // If all the files are empty
            if(fileCount == fileList.length)
                System.out.println("All "+ fileCount + " files contain empty spectra. The program will now terminate!");

            // Create two arrays which could store the x and y coordinates; these have the size identical to the coordinateList
            int[] xPositionArray = new int[coordinateList.size()];
            int[] yPositionArray = new int[coordinateList.size()];

            //   System.out.println("Printing coordinate list");
            for(int i=0; i< coordinateList.size();i++)
            {
                // Split the coordinateList to individual x and y coordinates
                String var[] = coordinateList.get(i).split("\t");
                // Add values to the position arrays
                xPositionArray[i] = Integer.valueOf(var[0]);
                yPositionArray[i] = Integer.valueOf(var[1]);
            }

            // Generate a MatchListDataStorage data object
            MatchListDataStorage mdsObject;
            MatchListDataStorage mdsObject2;

            System.out.println("---");
            System.out.print("STEP 2 of 6: Now performing distance calculation.");

            averageSimilarityList = AverageSimilarityCalculation(coordinateList, spectrumArrayList, fileNameList, epsilon);

            //Print all the calculated values in individual columns adjacent to each other
            for(int i=0; i < coordinateList.size();i++)
            {
                float pixelScoreValue = calculatePixelScore(newTotalIntensityList.get(i), averageSimilarityList.get(i));
                pixelScoreList.add(pixelScoreValue);
            }


            System.out.println("\n---");
            System.out.println("STEP 3 of 6: Generated pixel scores for all co-ordinate positions. Now performing " +
                    "peaklist ordering.");

            TreeMap<Integer, String> rankedList = new TreeMap<Integer, String>();

            // Create a TopologicalGreedyOrdering Object tco and then call the greedyOrdering method to create crystals
            // Also generates ranking list

            if(new String("TG").equals(scanner5Input)) {
                TopologicalGreedyOrdering tgo = new TopologicalGreedyOrdering();
                rankedList = tgo.greedyOrdering(fileNameList, coordinateList, pixelScoreList,
                        spectrumArrayList, xMaximum, yMaximum, epsilon);
            }
            else if(new String("CG").equals(scanner5Input)) {

                // Create a CrystalGrowthOrdering Object cco and then call the generateCrystal method to create crystals
                // Also generates ranking list
                CrystalGrowthOrdering cco = new CrystalGrowthOrdering();

                rankedList = cco.generateCrystal(fileNameList, coordinateList, pixelScoreList,
                        spectrumArrayList, xMaximum, yMaximum, epsilon);
            }

            System.out.println("---");
            System.out.println("STEP 4 of 6: Ranked list generated (Size: " + rankedList.size() + "). Now generating " +
                    "final consensus spectrum.");
            TraverseSparseSet tss = new TraverseSparseSet();

            SimpleSpectrum finalConsensus = tss.traverseSortedSet(rankedList, epsilon2, threshold);

            System.out.println("---");
            System.out.println("STEP 5 of 6: Final consensus spectrum generated. Now recalibrating the input " +
                    "peaklists against the final consensus spectrum.");


            int counter = 0;

            // Traverse the fileNameList and read every filename
            for(String item: fileNameList)
            {
                // Extract mz and intensity value from every item and store in dsrObject
                RawDataStorage dsrObject = extractMzIntensity(item);
                float[] mzArrayRecalibrate = dsrObject.mzArray; // Store mzArray

                double[] mzDoubleArrayRecalibrate = new double[mzArrayRecalibrate.length]; // Create a new array to store double mz values

                // Copy each value of mzArrayRecalibrate in mzDoubleArrayRecalibrate
                for(int j = 0; j< mzDoubleArrayRecalibrate.length;j++)
                    mzDoubleArrayRecalibrate[j] = mzArrayRecalibrate[j];

                float[] intArrayRecalibrate = dsrObject.intensityArray; // Store intensityArray
                double[] intDoubleArrayRecalibrate = new double[intArrayRecalibrate.length]; // Create a new array to store double intensity values

                // Copy each value of intArrayRecalibrate in intDoubleArrayRecalibrate
                for(int j = 0; j< intDoubleArrayRecalibrate.length;j++)
                    intDoubleArrayRecalibrate[j] = intArrayRecalibrate[j];

                // Create a simple Spectrum object
                SimpleSpectrum ssToRecalibrate = new SimpleSpectrum(mzDoubleArrayRecalibrate,intDoubleArrayRecalibrate);

                // Create a TraverseSparseSet object
                TraverseSparseSet tssObject = new TraverseSparseSet();

                // Call the performRecalibration function which then returns the recalibrated spectrum
                MutableSpectrum<Peak> newRecalibratedSpectrum = tssObject.performRecalibration(finalConsensus,
                        mzDoubleArrayRecalibrate, intDoubleArrayRecalibrate,
                        ssToRecalibrate, epsilon2, threshold);
                double[] newRecalibratedSpectraMz = new double[newRecalibratedSpectrum.size()];
                double[] newRecalibratedSpectraInt = new double[newRecalibratedSpectrum.size()];
                float[] newRecalibratedSpectraMzFloat = new float[newRecalibratedSpectrum.size()];
                float[] newRecalibratedSpectraIntFloat = new float[newRecalibratedSpectrum.size()];
                for(int k = 0; k<newRecalibratedSpectrum.size();k++)
                {
                    newRecalibratedSpectraMz[k] = newRecalibratedSpectrum.getMzAt(k);
                    newRecalibratedSpectraMzFloat[k] = (float) newRecalibratedSpectraMz[k];
                    newRecalibratedSpectraInt[k] = newRecalibratedSpectrum.getIntensityAt(k);
                    newRecalibratedSpectraIntFloat[k] = (float) newRecalibratedSpectraInt[k];
                }

                File dir = new File(pathName);
                //dir.mkdir();

                recalibratedDirPath = dir.getAbsolutePath();

                PrintWriter out1 = new PrintWriter(new FileWriter( dir +"/"+ item));
                // PrintWriter out1 = new PrintWriter(new FileWriter( recalibratedDirPath +"/"+ item));
                out1.println(coordinateList.get(counter));
                if(counter<fileNameList.size())
                    counter++;

                DecimalFormat df = new DecimalFormat();
                df.setMaximumFractionDigits(1);

                // Writing the recalibrated spectra in individual files
                for(int d = 0; d < newRecalibratedSpectrum.size(); d++)
                {
                    // double tempVal = Math.round(newRecalibratedSpectrum.getMzAt(d)*100.0)/100.0;
                    // (Removed precision limitation)

                    out1.println(newRecalibratedSpectrum.getMzAt(d)+ "\t" + newRecalibratedSpectrum.getIntensityAt(d));
                }
                out1.close();
            }
            System.out.println("---");
            System.out.println("STEP 6 of 6: All new recalibrated files successfully generated!\n");
        }
        else
        {
            System.out.println("The provided data path is either not a directory or a wrong path!");
        }
    }

    /**
     * Call the R script to generate intensity maps (generateIntensityMaps.R)
     *
     * @param
     * @return Generates a pdf file containing the ion intensity maps
     */

    public void welcomeMessage()
    {
        System.out.println("************************************************************************");
        System.out.print(
                "  __  __   _____  _____  _____                               _   \n" +
                        " |  \\/  | / ____||_   _|/ ____|                             | |  \n" +
                        " | \\  / || (___    | | | |      ___   _ __  _ __  ___   ___ | |_ \n" +
                        " | |\\/| | \\___ \\   | | | |     / _ \\ | '__|| '__|/ _ \\ / __|| __|\n" +
                        " | |  | | ____) | _| |_| |____| (_) || |   | |  |  __/| (__ | |_ \n" +
                        " |_|  |_||_____/ |_____|\\_____|\\___/ |_|   |_|   \\___| \\___| \\__|\n" +
                        "                                                                 \n");
        System.out.println("version 1.0.0\t\tContact: Purva Kulkarni (purva.kulkarni@uni-jena.de)");
        System.out.println("************************************************************************\n");
        System.out.println("Commandline usage: java -jar MSICorrectV_1_0_0.jar");
        System.out.println("Access help: java -jar MSICorrectV_1_0_0.jar --help");
        System.out.println("=================================================================\n");
    }

    public void printhelp() {
        System.out.println("MSICorrect is a Java-based tool that performs recalibration of mass spectrometry imaging " +
                "data.");

        System.out.println("\nRead and cite the following publication when using MSICorrect");
        System.out.println("Kulkarni, P., Kaftan, F., Kynast, P., Svatoš, A., Böcker, S. (2015).\nCorrecting mass " +
                "shifts: A lock mass-free recalibration procedure for mass spectrometry imaging data.\nAnal Bioanal " +
                "Chem, 407(25): 7603-13. DOI: 10.1007/s00216-015-8935-4");

        System.out.println("\n=================================================================");

        System.out.println("INFORMATION ON ESSENTIAL INPUT PARAMETERS IS PROVIDED BELOW");

        System.out.println("=================================================================");

        System.out.println("\n1. File path for input mass spectra \n" + "Diretory path that contains a list of " +
                "peaklists as *.txt files that have to be recalibrated. \nThe peaklists should have coordinates " +
                "values " +
                "x and y " +
                "separated by tab in the first line and then the mass and \nintensity value pairs separated by a tab " +
                "from " +
                "the " +
                "next line and onwards.\n");

        System.out.println("------\n");

        System.out.println("2. File path to store recalibrated spectra\n" + "Diretory path that  will " +
                "store the recalibrated peaklists. Make sure that you have permissions to write the files\nin the " +
                "specified directory.\n");

        System.out.println("------\n");

        System.out.println("3. Mass threshold for distance calculation between two peaklists\n" + "Suitable mass" +
                " window (in Da) depending on the dataset, to calculate pairwise similarity between the two peaklists" +
                "\n(example: 0.5).\n");

        System.out.println("------\n");

        System.out.println("4. Peaklist ordering approach (use TG or CG)\n" + "Before performing recalibration, it is"+
                "important that all the peaklists are ordered in a" +
                "specific manner.\nThis especially is helpful to minimize the error within the growing consensus " +
                "spectrum. MSICorrect implements \ntwo ordering approaches: type TG for topological greedy approach " +
                "and CG for Crytal growth ordering approach.\n");

        System.out.println("------\n");

        System.out.println("5. Line distance to perform line-pair stabbing recalibration\n"+"The Maximum " +
                "Line-Pair Stabbing (MLS) algorithm performs outlier detection and recalibration of a pair of " +
                "peaklists.\nThis parameter is a distance measure that separates the computational geometry " +
                "interpretation of two parallel lines\nin geometrical space (example: 0.5).\n");

        System.out.println("------\n");

        System.out.println("6. Mass threshold for for line-pair stabbing recalibration\n"+ "Suitable mass" +
                " window (in Da) depending on the dataset, to perform recalibration of two peaklists" +
                " (example: 0.8).\n");

        System.out.println("=====================");
        System.out.println("RECOMMENDED READING");
        System.out.println("=====================");

        System.out.println("S. Böcker and V. Mäkinen (2008). \n" +
                "Combinatorial approaches for mass spectra recalibration.\n" +
                "IEEE/ACM Trans Comput Biology Bioinform, 5(1):91-100.\n");
    }

    private void generateMaps(Scanner scanner, RConnection RC2) throws RserveException, REXPMismatchException
    {
        System.out.println("Enter the m/z range (within " + mzStart + " - " + mzEnd +" for the intensity maps (Format: mz1,mz2): ");
        mzRange = scanner.nextLine();

        String[] parts = mzRange.split(",");
        float mz1 = Float.parseFloat(parts[0]);
        float mz2 = Float.parseFloat(parts[1]);

        System.out.println("Step 12: Generating Intensity Maps for m/z range " + mz1 + "-" + mz2 +") and exporting it in a pdf file");

        REXP fileLocation = RC2.eval("generateIntensityMaps(\""+dirPath+"\",\""+preprocessedDirPath+"\", \""+recalibratedDirPath+"\","+mz1+","+mz2+")");
        System.out.println("Location of pdf intensity map file: " + fileLocation.asString());
    }


    /**
     * Converts the Array List to array and returns a string array (mz and intensity)
     *
     * @param ArrayList<String> fileData
     * @return fileDataArray
     */
    private String[] listToArray(ArrayList<String> fileData) {
        //Creates a new Array fileDataArray
        String[] fileDataArray = new String[fileData.size()];
        fileDataArray = fileData.toArray(fileDataArray); //Converts ArrayList to Array
        return fileDataArray;
    }
    /**
     * Extracts m/z and intensity list - method overriding 1.
     * Returns object of class RawDataStorage {@link RawDataStorage} containing the following fields
     * coordinateString, newMzArray, newIntensityArray, FileName.
     *
     * @param String[] fileDataArray
     * @param float[] mzArray
     * @param float[] intensityArray
     * @param String FileName
     * @return RawDataStorage
     */
    private RawDataStorage extractMzIntensity(String[] fileDataArray, float[] mzArray,float[] intensityArray, String FileName)
    {

        // Create new MzArray and intArray with length - 1 since coordinate will be removed from the file data
        float[] newMzArray = new float[mzArray.length-1];
        float[] newIntensityArray = new float[intensityArray.length-1];
        String coordinateString = null;



        for(int temp=0; temp<fileDataArray.length;temp++)
        {
            // If the line count is 0, split the variables as x and y
            if((temp == 0)) {
                // Stores the coordinate values in the string coordinateString
                coordinateString = fileDataArray[temp];
            }
            // If the line count is greater than 0, print the lines
            else
            if(temp > 0)
            {

                // Store the mz and int values in individual arrays
                String var[] = fileDataArray[temp].split("\t"); // Single space in the file separates mz and intensity
                String var3 = var[0];
                String var4 = var[1];
                mzArray[temp] = Float.parseFloat(var3);
                intensityArray[temp] = Float.parseFloat(var4);
            }
        }
        int tempVal=1;
        for(int l=0; l < newMzArray.length;l++)
        {
            newMzArray[l]= mzArray[tempVal];
            newIntensityArray[l]= intensityArray[tempVal];
            tempVal++;
        }

        // Create a new RawDataStorage Object ds which stores the coordinate, mz, int and file name and return it
        RawDataStorage ds = new RawDataStorage(coordinateString, newMzArray, newIntensityArray, FileName);
        return ds;
    }
    /**
     * Extracts m/z and intensity list when only a single file name is provided - method overriding 2.
     * Returns object of class RawDataStorage {@link RawDataStorage} containing the following fields
     * fileDataArrayI, mzArray, intensityArray, fileName.
     *
     * @param fileName
     * @return RawDataStorage
     */
    public RawDataStorage extractMzIntensity(String fileName) {
        //String fullFilePath = preprocessedDirPath + "/" + fileName;
        String fullFilePath = dirPath + "/" + fileName;
        ArrayList<String> fileDataI = new ArrayList<String>();
        RawDataStorage dsr;
        ReadFile rfI = new ReadFile();
        // call the fileInput method and passes fileName to it
        fileDataI = rfI.fileInput(fullFilePath); //returned ArrayList stored in fileData
        // Call the listTo Array Method
        String[] fileDataArrayI = listToArray(fileDataI);
        float[] mzArray = new float[fileDataArrayI.length];
        float[] intensityArray = new float[fileDataArrayI.length];
        dsr = extractMzIntensity(fileDataArrayI, mzArray, intensityArray, fileName);
        return dsr;
    }
    /**
     * pixelScore
     * This method calculates the score for every individual coordinate position
     * It basically adds the total intensity of a spectra and the similarity score with adjacent pixels
     *
     * @param float[] consensusI, float[] referenceJ, String coordinateI, String coordinateJ, String fileNameOfI, String fileNameOfJ
     * @return MatchListDataStorage
     */
    public float calculatePixelScore(float totalIntensity, float similarityScore) throws IOException
    {
        float pixelScore = totalIntensity + similarityScore;
        return pixelScore;
    }
    /**
     * AverageSimilarityCalculation
     * This method calculates the average similarity for each pixel in the coordinate list with all its neighbors(top, bottom, left, right)
     * @return ArrayList<Float> averageSimilarityList
     * Created on: 16 May 2015
     */
    public ArrayList<Float> AverageSimilarityCalculation(ArrayList<String> coordinateList, ArrayList<Object> spectrumArrayList, ArrayList<String> fileNameList, float epsilon) throws IOException
    {
        // Create a similarity list containing float distance values
        // Index 0 - top pixel
        // Index 1 - bottom pixel
        // Index 2 - left pixel
        // Index 3 - right pixel
        ArrayList<Float> similarityList = new ArrayList<Float>();
        ArrayList<Float> averageSimilarityList = new ArrayList<Float>();
        de.uni_jena.bioinformatik.functions.CalculateDistance cd = new de.uni_jena.bioinformatik.functions.CalculateDistance();
        MatchListDataStorage mdsObject;
        for (int k = 0; k<coordinateList.size();k++)
        {

            // Take the coordinate on position k as the current pixel
            String currentPixelCoordinate = coordinateList.get(k);
            String currentPixelFileName = fileNameList.get(k);

            //Split the coordinate value sting in x,y
            String temp1[] = currentPixelCoordinate.split("\t");
            int CurrentPixelXPos = Integer.valueOf(temp1[0]);
            int CurrentPixelYPos = Integer.valueOf(temp1[1]);

            Object mzObject = spectrumArrayList.get(k);
            float[] currentPixelMzArray = (float[]) mzObject; // this variable stores the m/z list of the current pixel
            String FileNameOfCurrentPixel = fileNameList.get(k);

            // Put a condition to check the current pixel is not at the end of a single row (i.e in the last column) of the pixel cube
            // (CurrentPixelYPos%yMaximum) != 0
            // Also, a second condition to check if the pixel is not at the beginning of the row
            if((CurrentPixelYPos%yMaximum) != 0 && (CurrentPixelYPos%yMaximum) != 1)
            {
                // Check if top, bottom, left, right neighbors of the currentPixel

                //topPixel (only the value of x changes - (x-1,y))
                String topPixelXPos = Integer.toString(CurrentPixelXPos - 1);
                String topPixelYPos = Integer.toString(CurrentPixelYPos);
                String topPixelCoordinate = topPixelXPos + "\t" + topPixelYPos;

                if (coordinateList.contains(topPixelCoordinate))
                {

                    //  System.out.println("Top Pixel exists in coordinate list: " + topPixelCoordinate);

                    // Get index of topPixel in coordinate list
                    int topPixelCoordinateIndex = coordinateList.indexOf(topPixelCoordinate);
                    // Get fileName of topPixel
                    String topPixelFileName = fileNameList.get(topPixelCoordinateIndex);
                    // Create topPixel m/z array
                    mzObject = spectrumArrayList.get(topPixelCoordinateIndex);
                    float[] topPixelMzArray = (float[]) mzObject;

                    // Calculate distance between current pixel and top Pixel
                    mdsObject= cd.specDistance(currentPixelMzArray,topPixelMzArray, currentPixelCoordinate, topPixelCoordinate, currentPixelFileName,topPixelFileName, epsilon);
                    // Add the top Pixel distance value to the averageSimilarityList at index 0
                    similarityList.add(0, mdsObject.distanceValue);
                }
                else
                {
                    similarityList.add(0, (float) 0);
                }

                //bottomPixel (only the value of x changes - (x+1,y))
                String bottomPixelXPos = Integer.toString(CurrentPixelXPos + 1);
                String bottomPixelYPos = Integer.toString(CurrentPixelYPos);
                String bottomPixelCoordinate = bottomPixelXPos + "\t" + bottomPixelYPos;

                if (coordinateList.contains(bottomPixelCoordinate))
                {
                    // Get index of bottomPixel in coordinate list
                    int bottomPixelCoordinateIndex = coordinateList.indexOf(bottomPixelCoordinate);
                    // Get fileName of bottomPixel
                    String bottomPixelFileName = fileNameList.get(bottomPixelCoordinateIndex);
                    // Create bottomPixel m/z array
                    mzObject = spectrumArrayList.get(bottomPixelCoordinateIndex);
                    float[] bottomPixelMzArray = (float[]) mzObject;

                    // Calculate distance between current pixel and bottom Pixel
                    mdsObject= cd.specDistance(currentPixelMzArray,bottomPixelMzArray,currentPixelCoordinate, bottomPixelCoordinate, currentPixelFileName,bottomPixelFileName, epsilon);
                    // Add the bottom Pixel distance value to the averageSimilarityList at index 0
                    similarityList.add(1, mdsObject.distanceValue);
                }
                else
                {
                    similarityList.add(1, (float) 0);
                }


                //leftPixel (only the value of y changes - (x,y-1))
                String leftPixelXPos = Integer.toString(CurrentPixelXPos);
                String leftPixelYPos = Integer.toString(CurrentPixelYPos - 1);
                String leftPixelCoordinate = leftPixelXPos + "\t" + leftPixelYPos;

                if (coordinateList.contains(leftPixelCoordinate))
                {
                    // Get index of leftPixel in coordinate list
                    int leftPixelCoordinateIndex = coordinateList.indexOf(leftPixelCoordinate);
                    // Get fileName of leftPixel
                    String leftPixelFileName = fileNameList.get(leftPixelCoordinateIndex);
                    // Create leftPixel m/z array
                    mzObject = spectrumArrayList.get(leftPixelCoordinateIndex);
                    float[] leftPixelMzArray = (float[]) mzObject;

                    // Calculate distance between current pixel and left Pixel
                    mdsObject= cd.specDistance(currentPixelMzArray,leftPixelMzArray, currentPixelCoordinate, leftPixelCoordinate, currentPixelFileName,leftPixelFileName, epsilon);
                    // Add the left Pixel distance value to the averageSimilarityList at index 0
                    similarityList.add(2, mdsObject.distanceValue);
                }
                else
                {
                    similarityList.add(2, (float) 0);
                }


                //rightPixel (only the value of y changes - (x,y+1))
                String rightPixelXPos = Integer.toString(CurrentPixelXPos);
                String rightPixelYPos = Integer.toString(CurrentPixelYPos + 1);
                String rightPixelCoordinate = rightPixelXPos + "\t" + rightPixelYPos;

                if (coordinateList.contains(rightPixelCoordinate))
                {
                    // Get index of rightPixel in coordinate list
                    int rightPixelCoordinateIndex = coordinateList.indexOf(rightPixelCoordinate);
                    // Get fileName of rightPixel
                    String rightPixelFileName = fileNameList.get(rightPixelCoordinateIndex);
                    // Create rightPixel m/z array
                    mzObject = spectrumArrayList.get(rightPixelCoordinateIndex);
                    float[] rightPixelMzArray = (float[]) mzObject;

                    // Calculate distance between current pixel and right Pixel
                    mdsObject= cd.specDistance(currentPixelMzArray,rightPixelMzArray, currentPixelCoordinate, rightPixelCoordinate, currentPixelFileName,rightPixelFileName, epsilon);
                    // Add the right Pixel distance value to the averageSimilarityList at index 0
                    similarityList.add(3, mdsObject.distanceValue);
                }
                else
                {
                    similarityList.add(3, (float) 0);
                }
            }
            else
            if((CurrentPixelYPos%yMaximum) == 0)
            {
                // Check if top, bottom and  left neighbors of the currentPixel

                //topPixel (only the value of x changes - (x-1,y))
                String topPixelXPos = Integer.toString(CurrentPixelXPos - 1);
                String topPixelYPos = Integer.toString(CurrentPixelYPos);
                String topPixelCoordinate = topPixelXPos + "\t" + topPixelYPos;

                if (coordinateList.contains(topPixelCoordinate))
                {
                    // Get index of topPixel in coordinate list
                    int topPixelCoordinateIndex = coordinateList.indexOf(topPixelCoordinate);
                    // Get fileName of topPixel
                    String topPixelFileName = fileNameList.get(topPixelCoordinateIndex);
                    // Create topPixel m/z array
                    mzObject = spectrumArrayList.get(topPixelCoordinateIndex);
                    float[] topPixelMzArray = (float[]) mzObject;

                    // Calculate distance between current pixel and top Pixel
                    mdsObject= cd.specDistance(currentPixelMzArray,topPixelMzArray, currentPixelCoordinate, topPixelCoordinate, currentPixelFileName,topPixelFileName, epsilon);
                    // Add the top Pixel distance value to the averageSimilarityList at index 0
                    similarityList.add(0, mdsObject.distanceValue);
                }
                else
                {
                    similarityList.add(0, (float) 0);
                }

                //bottomPixel (only the value of x changes - (x+1,y))
                String bottomPixelXPos = Integer.toString(CurrentPixelXPos + 1);
                String bottomPixelYPos = Integer.toString(CurrentPixelYPos);
                String bottomPixelCoordinate = bottomPixelXPos + "\t" + bottomPixelYPos;

                if (coordinateList.contains(bottomPixelCoordinate))
                {
                    // Get index of bottomPixel in coordinate list
                    int bottomPixelCoordinateIndex = coordinateList.indexOf(bottomPixelCoordinate);
                    // Get fileName of bottomPixel
                    String bottomPixelFileName = fileNameList.get(bottomPixelCoordinateIndex);
                    // Create bottomPixel m/z array
                    mzObject = spectrumArrayList.get(bottomPixelCoordinateIndex);
                    float[] bottomPixelMzArray = (float[]) mzObject;

                    // Calculate distance between current pixel and bottom Pixel
                    mdsObject= cd.specDistance(currentPixelMzArray,bottomPixelMzArray,currentPixelCoordinate, bottomPixelCoordinate, currentPixelFileName,bottomPixelFileName, epsilon);
                    // Add the bottom Pixel distance value to the averageSimilarityList at index 0
                    similarityList.add(1, mdsObject.distanceValue);
                }
                else
                {
                    similarityList.add(1, (float) 0);
                }

                //leftPixel (only the value of y changes - (x,y-1))
                String leftPixelXPos = Integer.toString(CurrentPixelXPos);
                String leftPixelYPos = Integer.toString(CurrentPixelYPos - 1);
                String leftPixelCoordinate = leftPixelXPos + "\t" + leftPixelYPos;

                if (coordinateList.contains(leftPixelCoordinate))
                {
                    // Get index of leftPixel in coordinate list
                    int leftPixelCoordinateIndex = coordinateList.indexOf(leftPixelCoordinate);
                    // Get fileName of leftPixel
                    String leftPixelFileName = fileNameList.get(leftPixelCoordinateIndex);
                    // Create leftPixel m/z array
                    mzObject = spectrumArrayList.get(leftPixelCoordinateIndex);
                    float[] leftPixelMzArray = (float[]) mzObject;

                    // Calculate distance between current pixel and left Pixel
                    mdsObject= cd.specDistance(currentPixelMzArray,leftPixelMzArray, currentPixelCoordinate, leftPixelCoordinate, currentPixelFileName,leftPixelFileName, epsilon);
                    // Add the left Pixel distance value to the averageSimilarityList at index 0
                    similarityList.add(2, mdsObject.distanceValue);
                }
                else
                {
                    similarityList.add(2, (float) 0);
                }

                similarityList.add(3, (float) 0);

            }
            else
            if((CurrentPixelYPos%yMaximum) == 1)
            {
                // Check if top, bottom and right neighbors of the currentPixel

                //topPixel (only the value of x changes - (x-1,y))
                String topPixelXPos = Integer.toString(CurrentPixelXPos - 1);
                String topPixelYPos = Integer.toString(CurrentPixelYPos);
                String topPixelCoordinate = topPixelXPos + "\t" + topPixelYPos;

                if (coordinateList.contains(topPixelCoordinate))
                {
                    // Get index of topPixel in coordinate list
                    int topPixelCoordinateIndex = coordinateList.indexOf(topPixelCoordinate);
                    // Get fileName of topPixel
                    String topPixelFileName = fileNameList.get(topPixelCoordinateIndex);
                    // Create topPixel m/z array
                    mzObject = spectrumArrayList.get(topPixelCoordinateIndex);
                    float[] topPixelMzArray = (float[]) mzObject;

                    // Calculate distance between current pixel and top Pixel
                    mdsObject= cd.specDistance(currentPixelMzArray,topPixelMzArray,currentPixelCoordinate, topPixelCoordinate, currentPixelFileName,topPixelFileName, epsilon);
                    // Add the top Pixel distance value to the averageSimilarityList at index 0
                    similarityList.add(0, mdsObject.distanceValue);
                }
                else
                {
                    similarityList.add(0, (float) 0);
                }

                //bottomPixel (only the value of x changes - (x+1,y))
                String bottomPixelXPos = Integer.toString(CurrentPixelXPos + 1);
                String bottomPixelYPos = Integer.toString(CurrentPixelYPos);
                String bottomPixelCoordinate = bottomPixelXPos + "\t" + bottomPixelYPos;

                if (coordinateList.contains(bottomPixelCoordinate))
                {
                    // Get index of bottomPixel in coordinate list
                    int bottomPixelCoordinateIndex = coordinateList.indexOf(bottomPixelCoordinate);
                    // Get fileName of bottomPixel
                    String bottomPixelFileName = fileNameList.get(bottomPixelCoordinateIndex);
                    // Create bottomPixel m/z array
                    mzObject = spectrumArrayList.get(bottomPixelCoordinateIndex);
                    float[] bottomPixelMzArray = (float[]) mzObject;

                    // Calculate distance between current pixel and bottom Pixel
                    mdsObject= cd.specDistance(currentPixelMzArray,bottomPixelMzArray, currentPixelCoordinate, bottomPixelCoordinate, currentPixelFileName,bottomPixelFileName, epsilon);
                    // Add the bottom Pixel distance value to the averageSimilarityList at index 0
                    similarityList.add(1, mdsObject.distanceValue);
                }
                else
                {
                    similarityList.add(1, (float) 0);
                }

                similarityList.add(2, (float) 0);

                //rightPixel (only the value of y changes - (x,y+1))
                String rightPixelXPos = Integer.toString(CurrentPixelXPos);
                String rightPixelYPos = Integer.toString(CurrentPixelYPos + 1);
                String rightPixelCoordinate = rightPixelXPos + "\t" + rightPixelYPos;

                if (coordinateList.contains(rightPixelCoordinate))
                {
                    // Get index of rightPixel in coordinate list
                    int rightPixelCoordinateIndex = coordinateList.indexOf(rightPixelCoordinate);
                    // Get fileName of rightPixel
                    String rightPixelFileName = fileNameList.get(rightPixelCoordinateIndex);
                    // Create rightPixel m/z array
                    mzObject = spectrumArrayList.get(rightPixelCoordinateIndex);
                    float[] rightPixelMzArray = (float[]) mzObject;

                    // Calculate distance between current pixel and right Pixel
                    mdsObject= cd.specDistance(currentPixelMzArray,rightPixelMzArray,currentPixelCoordinate, rightPixelCoordinate, currentPixelFileName,rightPixelFileName, epsilon);
                    // Add the right Pixel distance value to the averageSimilarityList at index 0
                    similarityList.add(3, mdsObject.distanceValue);
                }
                else
                {
                    similarityList.add(3, (float) 0);
                }
            }
            averageSimilarityList.add((similarityList.get(0) + similarityList.get(1) + similarityList.get(2) + similarityList.get(3)) / 4);
        }
        return averageSimilarityList;
    }
}