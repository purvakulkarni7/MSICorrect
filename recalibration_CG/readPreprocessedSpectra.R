#Function readPreprocessedSpectra: 
#Takes the file path of the directory where the preprocessed spectra (.txt format) are present and returns the data in the form of a list
#Return: spectraldata, xpixel, ypixel
#Author: Purva Kulkarni
#Date: 9 September 2014

#######################################################################################################

## load libraries
library(MALDIquant)
library(MALDIquantForeign)
library(dcemriS4)
library(lattice)

readPreprocessedSpectra <- function(FILE_PATH_PREPROCESSED)
{
  
  # Import text from the recalibrated files
  spectralDataPreprocessed <- importTxt(FILE_PATH_PREPROCESSED, verbose = TRUE, skip=1)
  
  spectralDataPreprocessed <- lapply(spectralDataPreprocessed, function(k) {
    
    ## fetch file name from metadata slot and read first line
    pos <- read.table(metaData(k)$file, nrows = 1, col.names = c("x", "y"))
    
    ## add coordinates to metadata
    metaData(k) <- modifyList(metaData(k), list(imaging=list(pos=pos)))
    
    ## return modified MassPeaks object
    k
  })
  
  return(spectralDataPreprocessed)
}