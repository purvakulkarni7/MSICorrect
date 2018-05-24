## plotIntensityMap.R
## Generates color intensity maps for imaging mass spectrometry data
## Input: 
## x = List of MassPeaks/ MassSpectrum Objects (Reference: MALDIQuant Package)
## range = Range of m/z values for which the plot has to be generated 
## main = Name of the plot
## Output: Color Intensity Map
## Author: Purva Kulkarni
## Date: 7 May 2014

## Example to call this function
## plotImage(s, range=c(156.95, 157.45), main="urinary bladder")

## define plotIntensityMap function
plotIntensityMap <- function(x, range, main) 
{
  ## display only mass in range
  x <- trim(x, range=range)
  
  ## find x and y positions
  pos <- lapply(x, function(y)metaData(y)$imaging$pos)
  pos <- do.call(rbind, pos)
  
  ## max x/y to build image matrix
  nc <- max(pos[, "x"])
  nr <- max(pos[, "y"])
  
  ## init matrix
  m <- matrix(NA, nrow=nr, ncol=nc)
  
  ## fill matrix with intensity values
  for (i in seq(along=x)) {
    m[pos[i, "y"], pos[i, "x"]] <- sum(intensity(x[[i]]), na.rm=TRUE)
  }
  
  ## scale matrix (better contrast)
  m <- m/max(m)
  
  ## build title
  main <- paste(main, " (m/z: ", range[1], "-", range[2], ")", sep="")
  
  ## prepare plot area
  #plot(NA, type="n", xlim=c(1, nc), ylim=c(1, nr), axes=TRUE, asp=1, xlab="", ylab="", main=main)
  ## plot image
  #rasterImage(m, xleft=1, xright=nc, ybottom=1, ytop=nr, interpolate=FALSE)
  
  #image(t(m), col=rainbow(100, start = 1/6, end = 1), asp=1)
  
  levelplot(t(m), scales=list(tick.number=0), col.regions=rainbow(100, start = 1/6, end = 1), xlab=" ", ylab=" ", main=main)
}