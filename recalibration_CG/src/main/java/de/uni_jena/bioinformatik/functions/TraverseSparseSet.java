/**
 * 
 */

package de.uni_jena.bioinformatik.functions;


import de.uni_jena.bioinformatik.input.RawDataStorage;
import de.uni_jena.bioinformatik.recalibrationModule.MzRecalibration;
import de.unijena.bioinf.ChemistryBase.ms.MutableSpectrum;
import de.unijena.bioinf.ChemistryBase.ms.Peak;
import de.unijena.bioinf.ChemistryBase.ms.utils.SimpleMutableSpectrum;
import de.unijena.bioinf.ChemistryBase.ms.utils.SimpleSpectrum;
import org.apache.commons.math3.analysis.UnivariateFunction;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Scanner;
import java.util.TreeMap;

/**
 * TraverseSparseSet.
 * <h3>Usage</h3>
 * <ol>
 * <li>Traverse the sparse set containing coordinate pairs and distance values</li>
 * <li>Takes the two mass lists with minimum distance and performs recalibration 
 * using functions available in the recalibration module {@code maxLinePairStabbing}</li>
 * </ol>
 * 
 * @author Purva Kulkarni
 * @Date 2 June 2014
 *
 */

public class TraverseSparseSet{
	public de.uni_jena.bioinformatik.input.InputData id = new de.uni_jena.bioinformatik.input.InputData();
	SimpleSpectrum newConsensus;
		int counter;;
// Default Constructor
	public  TraverseSparseSet() {}

	/**
	 * Method to traverse the rankedList and perform recalibration using the recalibration module 
	 * 
	 * @param rankedList
	 */

    public SimpleSpectrum traverseSortedSet(TreeMap<Integer, String> rankedList, double epsilon, double threshold)
    {
        // for loop which will traverse the complete ranked list
        for(int i = 1; i < rankedList.size(); i++)
        {
            if(i == 1)
            {
                long startTime = System.currentTimeMillis();
                // Uses generateConsensus method (Case I)
                newConsensus = generateConsensus(rankedList.get(i+1), rankedList.get(i), epsilon, threshold);
                long endTime = System.currentTimeMillis();
                i++;
            }
            else
            {
                // Uses generateConsensus method (Case II)
                long startTime1 = System.currentTimeMillis();
                newConsensus = generateConsensus(newConsensus, rankedList.get(i+1), epsilon, threshold);
                long endTime1 = System.currentTimeMillis();
            }
        }
        return newConsensus;
    }


    /**
	 * Generates a consensus spectrum when two spectra fileNames are provided (Case I)
	 * 
	 * @param recalibrateJFileName
	 * @param consensusIFileName
	 * @return SimpleSpectrum
	 */
	public SimpleSpectrum generateConsensus(String recalibrateJFileName,
			String consensusIFileName, double epsilon, double threshold) {
		RawDataStorage dsrConsensus;
		RawDataStorage dsrRecalibrate;

		// Extract mz and intensity from file 1 (considered as consensus)
		dsrConsensus = id.extractMzIntensity(consensusIFileName);

		// Create a Double array to store the mz Values 	
		double[] mzDoubleListConsensus = new double[dsrConsensus.mzArray.length];
		mzDoubleListConsensus = extractDoubleMzList(dsrConsensus, mzDoubleListConsensus);

		// Create a Double array to store the intensity Values 	
		double[] intensityDoubleListConsensus = new double[dsrConsensus.intensityArray.length];
		intensityDoubleListConsensus = extractDoubleIntensityList(dsrConsensus, intensityDoubleListConsensus);

		// Create SimpleSpectrum Object ssConsensus
		SimpleSpectrum ssConsensus = new SimpleSpectrum(mzDoubleListConsensus,intensityDoubleListConsensus);

		// Extract mz and intensity from file 2 (considered as recalibrate)
		dsrRecalibrate = id.extractMzIntensity(recalibrateJFileName);

		// Create a Double array to store the mz Values 	
		double[] mzDoubleListRecalibrate = new double[dsrRecalibrate.mzArray.length];
		mzDoubleListRecalibrate = extractDoubleMzList(dsrRecalibrate, mzDoubleListRecalibrate);

		// Create a Double array to store the intensity Values 	
		double[] intensityDoubleListRecalibrate = new double[dsrRecalibrate.intensityArray.length];
		intensityDoubleListRecalibrate = extractDoubleIntensityList(dsrRecalibrate, intensityDoubleListRecalibrate);

		// Create SimpleSpectrum Object ssRecalibrate
		SimpleSpectrum ssRecalibrate = new SimpleSpectrum(mzDoubleListRecalibrate,intensityDoubleListRecalibrate);



		//send to recalibration routine
		MutableSpectrum<Peak> newSMutable = performRecalibration(ssConsensus,
				mzDoubleListRecalibrate, intensityDoubleListRecalibrate,
				ssRecalibrate, epsilon, threshold);

        SimpleSpectrum newConsensus;

		//merge spectra
		newConsensus = mergeSpectrum(ssConsensus, newSMutable, epsilon);

		return newConsensus;
	}


	/**
	 * Generates a consensus spectrum when a newly generated consensus and a spectra fileName are provided (Case II)
	 * 
	 * @param newConsensus
	 * @param recalibrateJFileName
	 * @return 
	 */
	public SimpleSpectrum generateConsensus(SimpleSpectrum newConsensus, String recalibrateJFileName, double epsilon, double threshold) {

		RawDataStorage dsrRecalibrate;

		// Extract mz and intensity from file 2 (considered as recalibrate)
		dsrRecalibrate = id.extractMzIntensity(recalibrateJFileName);

		// Create a Double array to store the mz Values 	
		double[] mzDoubleListRecalibrate = new double[dsrRecalibrate.mzArray.length];
		mzDoubleListRecalibrate = extractDoubleMzList(dsrRecalibrate, mzDoubleListRecalibrate);

		// Create a Double array to store the intensity Values 	
		double[] intensityDoubleListRecalibrate = new double[dsrRecalibrate.intensityArray.length];
		intensityDoubleListRecalibrate = extractDoubleIntensityList(dsrRecalibrate, intensityDoubleListRecalibrate);

		// Create SimpleSpectrum Object ssRecalibrate
		SimpleSpectrum ssRecalibrate = new SimpleSpectrum(mzDoubleListRecalibrate,intensityDoubleListRecalibrate);


		//send to recalibration routine
		MutableSpectrum<Peak> newSMutable = performRecalibration(newConsensus,
				mzDoubleListRecalibrate, intensityDoubleListRecalibrate,
				ssRecalibrate, epsilon, threshold);

		//merge spectra
		newConsensus = mergeSpectrum(newConsensus, newSMutable, epsilon);

		return newConsensus;
	}


	/**
	 * Performs recalibration steps when 2 spectra are provided
	 * 
	 * @param ssConsensus
	 * @param mzDoubleListRecalibrate
	 * @param intensityDoubleListRecalibrate
	 * @param ssRecalibrate
	 * 
	 * @return MutableSpectrum
	 */
	public MutableSpectrum<Peak> performRecalibration(
			SimpleSpectrum ssConsensus, double[] mzDoubleListRecalibrate,
			double[] intensityDoubleListRecalibrate,
			SimpleSpectrum ssRecalibrate, double epsilon, double threshold) {

        long startTimerecal = System.currentTimeMillis();

        Scanner scanner = new Scanner(System.in );

		double[][] subsetArray = MzRecalibration.maxLinePairStabbing(ssRecalibrate, ssConsensus, epsilon, threshold);

		double[] subsetX = new double[subsetArray[0].length];
		double[] subsetY = new double[subsetArray[0].length];

		for(int i=0;i<subsetArray.length;i++)
		{
			if(i==0)
			{
				for(int j=0;j<subsetArray[i].length;j++)
				{
					subsetX[j] = subsetArray[i][j];
				}
			}
			else
				if(i==1)
				{
					for(int j=0;j<subsetArray[i].length;j++)
					{
						subsetY[j] = subsetArray[i][j];
					}
				}
		}

		UnivariateFunction uf = MzRecalibration.getLinearRecalibration(subsetX,subsetY);
		SimpleMutableSpectrum sMutable = new SimpleMutableSpectrum(mzDoubleListRecalibrate,intensityDoubleListRecalibrate);

		MutableSpectrum<Peak> newSMutable =  MzRecalibration.recalibrate(sMutable, uf);

        long endTimerecal = System.currentTimeMillis();

		return newSMutable;
	}


	/**
	 * Merges two spectra (as of here: old consensus and new consensus) 
	 * 
	 * @param oldConsensus, MutableSpectrum<Peak> mutatedSpectrum
	 * @return SimpleSpectrum
	 */		

	public SimpleSpectrum mergeSpectrum(SimpleSpectrum oldConsensus, MutableSpectrum<Peak> mutatedSpectrum, double epsilon)
	{

        long startTimeMerge = System.currentTimeMillis();

		ArrayList<Float> mzListNewConsensus = new ArrayList<Float>(); 
		ArrayList<Float> intListNewConsensus = new ArrayList<Float>();

		TraverseSparseSet tss = new TraverseSparseSet();

		float[][] mzIntOldConsensus;
		float[][] mzIntmutatedSpectrum;

		// Fill the 2D array with m/z and int values from old and mutated spectra
		mzIntOldConsensus = tss.getMzIntensityFromObject(oldConsensus);
		mzIntmutatedSpectrum = tss.getMzIntensityFromObject(mutatedSpectrum);


		float[] mzArrayOld = new float[oldConsensus.size()];	
		float[] intensityArrayOld = new float[oldConsensus.size()];

		// loop to fill the m/z and intensity values in the array mzArrayOld and intensityArrayOld
		for(int i=0;i<mzIntOldConsensus.length;i++)
		{
			for(int j=0;j<mzIntOldConsensus[i].length;j++)
			{
				if (j==0)
				{
					mzArrayOld[i] = mzIntOldConsensus[i][j];
				}
				else
					if(j==1)
					{
						intensityArrayOld[i] = mzIntOldConsensus[i][j];
					}
			}
		}

		float[] mzArrayMutated = new float[mutatedSpectrum.size()];	
		float[] intensityArrayMutated = new float[mutatedSpectrum.size()];

		// loop to fill the m/z and intensity values in the array mzArrayMutated and intensityArrayMutated
		for(int i=0;i<mzIntmutatedSpectrum.length;i++)
		{
			for(int j=0;j<mzIntmutatedSpectrum[i].length;j++)
			{
				if (j==0)
				{
					mzArrayMutated[i] = mzIntmutatedSpectrum[i][j];
				}
				else
					if(j==1)
					{
						intensityArrayMutated[i] = mzIntmutatedSpectrum[i][j];
					}
			}
		}

		// merge the mzArrayOld and the mzArrayMutated in a sorted order	

		float[] mzMergedTempArray = new float[mzArrayOld.length + mzArrayMutated.length];
		float[] intMergedTempArray = new float[intensityArrayOld.length + intensityArrayMutated.length]; 

		int p = 0, q = 0, r = 0;

		// Merging of mzArrayOld and mzArrayMutated
		while (p < mzArrayOld.length && q < mzArrayMutated.length)
		{
			if (mzArrayOld[p] < mzArrayMutated[q])
			{
				int value=p++;
				int value1 = r++;
				mzMergedTempArray[value1] = mzArrayOld[value];          
				intMergedTempArray[value1] = intensityArrayOld[value];

			}

			else
			{
				int value2=q++;
				int value3 =r++;
				mzMergedTempArray[value3] = mzArrayMutated[value2];    
				intMergedTempArray[value3] = intensityArrayMutated[value2]; 
			}
		}

		while (p < mzArrayOld.length)
		{
			int value4 = p++;
			int value5=r++;
			mzMergedTempArray[value5] = mzArrayOld[value4];
			intMergedTempArray[value5] = intensityArrayOld[value4];
		}


		while (q < mzArrayMutated.length) 
		{
			int value6 = q++;
			int value7=r++;
			mzMergedTempArray[value7] = mzArrayMutated[value6];
			intMergedTempArray[value7] = intensityArrayMutated[value6];   
		}

		for(int i=0;i < (mzMergedTempArray.length-1);i=i+2)
		{
			if(Math.abs(mzMergedTempArray[i] - mzMergedTempArray[i+1]) <= epsilon)
			{
				mzListNewConsensus.add(((intMergedTempArray[i] * mzMergedTempArray[i])+ (intMergedTempArray[i+1] * mzMergedTempArray[i+1]))/(intMergedTempArray[i] +intMergedTempArray[i+1]));
				intListNewConsensus.add(intMergedTempArray[i] +intMergedTempArray[i+1]);
				//i++;

			}
			else
			{
				mzListNewConsensus.add(mzMergedTempArray[i]);
				mzListNewConsensus.add(mzMergedTempArray[i+1]);
				intListNewConsensus.add(intMergedTempArray[i]);
				intListNewConsensus.add(intMergedTempArray[i+1]);
				//i++;

			}
		}

		// Make a while loop to do the merging iteratively
		int instanceCounter = 0;

		do{
			instanceCounter = 0;
			for(int s=0; s<(mzListNewConsensus.size()-1) ; s++)
			{
				if((Math.abs(mzListNewConsensus.get(s)-mzListNewConsensus.get(s+1)) <= epsilon))
				{
					instanceCounter++;
					break;
				}
				s++;

			}

			if(instanceCounter > 0)
			{
				ArrayList<Float> mzTempStorageList = new ArrayList<Float>();
				ArrayList<Float> intTempStorageList = new ArrayList<Float>();

				//Check if the distance between all the m/z points is greater than epsilon
				for(int s=0; s<(mzListNewConsensus.size()-1) ; s++)
				{
					if((Math.abs(mzListNewConsensus.get(s)-mzListNewConsensus.get(s+1)) <= epsilon))
					{
						mzTempStorageList.add(((intListNewConsensus.get(s) * mzListNewConsensus.get(s))+ (intListNewConsensus.get(s+1) * mzListNewConsensus.get(s+1)))/(intListNewConsensus.get(s) +intListNewConsensus.get(s+1)));
						intTempStorageList.add(intListNewConsensus.get(s) + intListNewConsensus.get(s+1));
						s++;

					}
					else
					{
						mzTempStorageList.add(mzListNewConsensus.get(s));
						mzTempStorageList.add(mzListNewConsensus.get(s+1));
						intTempStorageList.add(intListNewConsensus.get(s));
						intTempStorageList.add(intListNewConsensus.get(s+1));
						s++;
					}
				}

				// Clearing mzListNewConsensus so that mzTempStorageList can be transfered to it
				mzListNewConsensus.clear();
				intListNewConsensus.clear();

				//copy back into mzListNewConsensus
				for(int l=0;l<mzTempStorageList.size();l++)
				{
					mzListNewConsensus.add(mzTempStorageList.get(l));
					intListNewConsensus.add(intTempStorageList.get(l));

				}

				mzTempStorageList.clear();
				intTempStorageList.clear();
			}
			else
				break;

		}while(instanceCounter > 0);
		

		//Check if the distance between all the m/z points is greater than epsilon
		for(int s=0; s<mzListNewConsensus.size()-1 ; s++)
		{
			
			if((Math.abs(mzListNewConsensus.get(s)-mzListNewConsensus.get(s+1)) <= epsilon))
			{
				System.out.println("Report " + mzListNewConsensus.get(s) + "\t" + mzListNewConsensus.get(s+1));
				s++;
			}			
			else
			{
				s++;
			}
		}

		double[] mzArrayNewConsensus = new double[mzListNewConsensus.toArray().length];	
		double[] intensityArrayNewConsensus = new double[mzListNewConsensus.toArray().length];

		for(int i =0; i< mzListNewConsensus.size(); i++)
		{
			mzArrayNewConsensus[i] = mzListNewConsensus.get(i);
			intensityArrayNewConsensus[i] = intListNewConsensus.get(i);
		}

		SimpleSpectrum newConsensus = new SimpleSpectrum(mzArrayNewConsensus, intensityArrayNewConsensus);

        long endTimeMerge = System.currentTimeMillis();
		return newConsensus;
	}

	/**
	 * Removes duplicate values in a String ArrayList
	 * 
	 * @param arlList
	 * @return ArrayList<String>
	 */	
	public ArrayList<String> removeDuplicate(ArrayList<String> arlList)
	{
		HashSet<String> h = new HashSet<String>(arlList);
		arlList.clear();
		arlList.addAll(h);
		return arlList;
	}

	/**
	 * Extracts the mz and intensity values from a MutableSpectrum<Peak> object  
	 * 
	 * @param ms
	 * @return double[][] mzIntensityArray
	 */	
	public float[][] getMzIntensityFromObject(MutableSpectrum<Peak> ms)
	{
		float[][] mzIntensityArray = new float[ms.size()][2];
		for(int i=0; i< ms.size();i++)
		{
			mzIntensityArray[i][0] = (float) ms.getMzAt(i);
			mzIntensityArray[i][1] = (float) ms.getIntensityAt(i);
		}
		return mzIntensityArray;
	}

	/**
	 * Extracts the mz and intensity values from a SimpleSpectrum object 
	 * 
	 * @param
	 * @return double[][] mzIntensityArray
	 */	
	public float[][] getMzIntensityFromObject(SimpleSpectrum ss)
	{
		float[][] mzIntensityArray = new float[ss.size()][2];
		for(int i=0; i< ss.size();i++)
		{
			mzIntensityArray[i][0] = (float) ss.getMzAt(i);
			mzIntensityArray[i][1] = (float) ss.getIntensityAt(i);
		}
		return mzIntensityArray;
	}

	/**
	 * Extracts the double mz array from the RawDataStorage object 
	 * 
	 * @param  dsr
	 * @return double[] mzDoubleList
	 */	

	public double[] extractDoubleMzList(de.uni_jena.bioinformatik.input.RawDataStorage dsr, double[] mzDoubleArray) {
		int d=0;
		for(float arrayItem : dsr.mzArray){
			//Type cast the float arrayItem to double and store it in mzDoubleArray
			mzDoubleArray[d] = (double)(arrayItem) ;
			d++;
		}	
		return mzDoubleArray;
	}

	/**
	 * Extracts the double intensity array from the RawDataStorage object 
	 * 
	 * @param  dsr
	 * @return double[] intensityDoubleList
	 */	

	public double[] extractDoubleIntensityList(de.uni_jena.bioinformatik.input.RawDataStorage dsr, double[] intensityDoubleArray) {
		int d=0;
		for(float arrayItem : dsr.intensityArray){
			//Type cast the float arrayItem to double and store it in mzDoubleArray
			intensityDoubleArray[d] = (double)(arrayItem) ;
			d++;
		}	
		return intensityDoubleArray;
	}
}