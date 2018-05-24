## Main.R
## Main R file which takes in an imaging data file as input and calls different functions for further processing of data
## Input: File path/Directory Path
## Output: Multiple outputs
## Author: Purva Kulkarni
## Date: 7 May 2014
## Last Edited: 9 September 2014

Main <- function(FILE_PATH)
{ 
  ## load libraries
  library(MALDIquant)
  library(MALDIquantForeign)
  library(dcemriS4)
  require(gridExtra) # also loads grid
  library(lattice)
  library(fields)
  library(matlab)
  
  #Call the source files of the function which this script will use
  source('/home/purva/Server/Lab_proceedings_Current_Folder/Work/Mass_Recalibration/CrystalGrowthRecalibration/recalibration_CG/importImagingFile.R', echo=TRUE)
  source('/home/purva/Server/Lab_proceedings_Current_Folder/Work/Mass_Recalibration/CrystalGrowthRecalibration/recalibration_CG/plotIntensityMap.R', echo=TRUE)
  source('/home/purva/Server/Lab_proceedings_Current_Folder/Work/Mass_Recalibration/CrystalGrowthRecalibration/recalibration_CG/dataPreprocessing.R', echo=TRUE)
  source('/home/purva/Server/Lab_proceedings_Current_Folder/Work/Mass_Recalibration/CrystalGrowthRecalibration/recalibration_CG/exportFile.R', echo=TRUE)
  source('/home/purva/Server/Lab_proceedings_Current_Folder/Work/Mass_Recalibration/CrystalGrowthRecalibration/recalibration_CG/readRecalibratedSpectra.R', echo=TRUE)
  source('/home/purva/Server/Lab_proceedings_Current_Folder/Work/Mass_Recalibration/CrystalGrowthRecalibration/recalibration_CG/plotIMSSlice.R', echo=TRUE)
  source('/home/purva/Server/Lab_proceedings_Current_Folder/Work/Mass_Recalibration/CrystalGrowthRecalibration/recalibration_CG/plotIMSSlicePlusMinus.R', echo=TRUE)

# Call the function importAnalyzeFormat.R and provide FILE_PATH as input
spectralDataObjects <- importImagingFile(FILE_PATH)

# Extract the mass range of the data
totalMzValues<-mass(spectralDataObjects[[1]])
totalLength <- length(t(totalMzValues))
mzStart<-totalMzValues[1]
mzEnd<-totalMzValues[totalLength]

#s1 <- transformIntensity(spectralDataObjects, method="sqrt")

# Smooth Intensity
s2 <- smoothIntensity(spectralDataObjects, method="MovingAverage", halfWindowSize=2)

# Baseline Correction
#s3 <- removeBaseline(s2, method="TopHat", halfWindowSize=250)
s3 <- removeBaseline(s2, method="TopHat", halfWindowSize=100)


# Smooth Intensity
s4<-smoothIntensity(s3, method="MovingAverage", halfWindowSize=2)

# Perform Peak Picking
p <- detectPeaks(s4, method="MAD", halfWindowSize=13, SNR=3)

# Assign the p to preprocessedDataObjects
preprocessedDataObjects<-p


# Call the exportFile function to export individual text files with spectral information
 dir.create("PreprocessedSpectra", showWarnings = FALSE)
  setwd("PreprocessedSpectra")
 
for(i in 1:length(preprocessedDataObjects))
{
 coordinateValue<-metaData(preprocessedDataObjects[[i]])
 coordinates<-coordinateValue$imaging$pos 
 mzValues<-mass(preprocessedDataObjects[[i]])
 intensityValues<-intensity(preprocessedDataObjects[[i]])
 exportFile(coordinates,mzValues,intensityValues, i)
}

print("Spectrum files exported. Program will now terminate")
print("############################################################")

newFilesPath <- getwd();

combinedReturn <- c(newFilesPath,mzStart,mzEnd)

# Return the newFilesPath,mzStart and mzEnd values to the java program
return(combinedReturn)
}
