#Function readRecalibratedSpectra: 
#Takes the file path of the directory where the recalibrated spectra (.txt format) are present and returns the data in the form of a list
#Return: spectraldata, xpixel, ypixel
#Author: Purva Kulkarni
#Date: 29 July 2014

#######################################################################################################

## load libraries
library(MALDIquant)
library(MALDIquantForeign)
library(dcemriS4)
library(lattice)

readRecalibratedSpectra <- function(FILE_PATH_RECALIBRATED)
{

# Import text from the recalibrated files
spectralDataRecalibrated <- importTxt(FILE_PATH_RECALIBRATED, verbose = TRUE, skip=1)

spectralDataRecalibrated <- lapply(spectralDataRecalibrated, function(s) {

  ## fetch file name from metadata slot and read first line
  pos <- read.table(metaData(s)$file, nrows = 1, col.names = c("x", "y"))
  
  ## add coordinates to metadata
  metaData(s) <- modifyList(metaData(s), list(imaging=list(pos=pos)))
  
  ## return modified MassPeaks object
 s
})

return(spectralDataRecalibrated)
}