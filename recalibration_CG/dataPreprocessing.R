## dataPreprocessing.R
## Performs preprocessing steps on the input data. The steps include:
## Sqrt transform, Smoothing, Baseline Correction, Peak Picking
## Input: 
## x = List of MassSpectrum Objects (Reference: MALDIQuant Package) 
## Output: Preprocessed list of MassSpectrum Objects (if peakpicking not applied)
## or, list of MassPeaks Objects (if peakpicking applied)
## Author: Purva Kulkarni
## Date: 15 May 2014

## define plotIntensityMap function
dataPreprocessing <- function(x)
{
  ## sqrt transform (for variance stabilization)
  s1 <- transformIntensity(x, method="sqrt")
  
  ## 21 point Savitzky-Golay-Filter for smoothing spectra
  ## (maybe you have to adjust the halfWindowSize;
  ## you could use a simple moving average instead)
  ## see ?smoothIntensity
  s2 <- smoothIntensity(s1, method="MovingAverage", halfWindowSize=2)
  
  
  ## remove baseline
  ## (maybe you have to adjust iterations to your spectra; high resolution
  ## spectra need a much lower iteration number (halfWindowSize, for some other
  ## baseline estimation algorithms)
  ## see ?removeBaseline, ?estimateBaseline
  s3 <- removeBaseline(s2, method="ConvexHull", iterations=100)
  
  
#   ## run peak detection
#   ## (maybe you need to adjust halfWindowSize [decreasing it for high resolution
#   ## spectra] and SNR [a higher value increase the True-Positive-Rate but decrease
#   ## sensitivity])
#   ## see ?detectPeaks, ?estimateNoise
   p <- detectPeaks(s3, method="MAD", halfWindowSize=20, SNR=2)

return(p)
}
