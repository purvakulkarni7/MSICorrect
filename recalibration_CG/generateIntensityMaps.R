## generateIntensityMaps.R
## File which takes in the paths of recalibrated spectra files and preprocessed spectra files and generates intensity maps
## Input: File path/Directory Path
## Output: Multiple outputs
## Author: Purva Kulkarni
## Date: 7 May 2014
## Last Edited: 9 September 2014

generateIntensityMaps <- function(FILE_PATH,FILE_PATH_PREPROCESSED,FILE_PATH_RECALIBRATED,mz1,mz2)
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
  source('/home/purva/Server/Lab_proceedings_Current_Folder/Work/Mass_Recalibration/CrystalGrowthRecalibration/recalibration_CG/readPreprocessedSpectra.R', echo=TRUE)
  source('/home/purva/Server/Lab_proceedings_Current_Folder/Work/Mass_Recalibration/CrystalGrowthRecalibration/recalibration_CG/readRecalibratedSpectra.R', echo=TRUE)
  source('/home/purva/Server/Lab_proceedings_Current_Folder/Work/Mass_Recalibration/CrystalGrowthRecalibration/recalibration_CG/plotIMSSlice.R', echo=TRUE)
  source('/home/purva/Server/Lab_proceedings_Current_Folder/Work/Mass_Recalibration/CrystalGrowthRecalibration/recalibration_CG/plotIMSSlicePlusMinus.R', echo=TRUE)
  
  # Read spectral data
  spectralDataObjects <- importImagingFile(FILE_PATH)
  
  # Read preprocessed data
  preprocessedDataObjects <- readPreprocessedSpectra(FILE_PATH_PREPROCESSED)
    
  # Read recalibrated data
  recalibratedDataObjects <- readRecalibratedSpectra(FILE_PATH_RECALIBRATED)
  
 # par(mfrow=c(2,2))
  
  # Original Image
  
png(filename = "RawImage.png")
imgRaw <- plotIMSSlice(spectralDataObjects, range = c(mz1, mz2), main = "Raw")
plot(imgRaw)
dev.off()
 
png(filename = "PeakPickedImage.png")
  # Peak Picked Image
  imgPeakPicked<-plotIMSSlice(preprocessedDataObjects, range = c(mz1, mz2), main = "Peak Picked")
plot(imgPeakPicked)
dev.off()

  
  # MassWise Peak picked
png(filename = "MasswisePeakPickedImage.png")
  imgMassWisePeakPicked<-plotIMSSlicePlusMinus(preprocessedDataObjects, range = c(mz1, mz2), main = "Mass-wise peakpicked")
plot(imgMassWisePeakPicked)
dev.off()
  
  # MassWise Recalibrated
png(filename = "MassWiseRecalibrated.png")
  imgMassWiseRecalibrated<-plotIMSSlicePlusMinus(recalibratedDataObjects, range = c(mz1, mz2), main = "Mass-wise recalibrated")
plot(imgMassWiseRecalibrated)
dev.off()
  
  # Arrange images in the form of a grid and print in a pdf file in landscape mode (a4r)
 # pdf("Fly_data_maps.pdf", width=11.00 ,height=9.00, paper="a4r")
 # grid.arrange(img1Raw, img2PeakPicked, img3MassWisePeakPicked, img4MassWiseRecalibrated)
 # dev.off()
  
  newFilesPath <- getwd();
  
  return(newFilesPath)
}

  