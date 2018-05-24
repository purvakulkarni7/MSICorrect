## exportFile.R
## Writes the co-ordinate information, mass and intensity values for a specific 
## pixel to individual text files
## Input: 
## x = List of MassSpectrum Objects (Reference: MALDIQuant Package) 
## Output: text file with coordinate information and m/z and intensity values
## Author: Purva Kulkarni
## Date: 19 May 2014

## define exportFile function
exportFile <- function(coordinate, mz, intensity, i)
{
#   filename<-filename<-paste(c('/Preprocessed_Spectra/Spectrum_',i,'.txt'), collapse='')
  
  filename <-paste(c('Spectrum_', i, '.txt'), collapse = '')
  write(coordinate, filename, ncolumns=2)
  spectralData <- cbind(mz, intensity, deparse.level = 0)
  foo<-c('%2.1f', '%2.9f')
  cbar <- sapply(1:2,function(j) sprintf(foo[j],spectralData[,j]))
  write(t(cbar),filename,ncolumns=2, append=TRUE)
}