#' This function is a basic function for MALDI-IMS data.
#'
#' @title Plot MALDI-IMS data
#'
#' @param x list of \code{\linkS4class{MassSpectrum}} objects.
#' @param range \code{double}, length two, mass range
#' @param main \code{character} plot title.
#'
#' @return Nothing, plots an image.
#'
#' @seealso \url{https://github.com/sgibb/ims-shiny/}
#'
#' @export
#'
#'

plotIMSSlicePlusMinus <- function(x, range, main) {
  ## display only mass in range
  suppressWarnings(x <- trim(x, range=range))
  
  ## find x and y positions
  pos <- lapply(x, function(y)metaData(y)$imaging$pos)
  pos <- do.call(rbind, pos)
  
  ## max x/y to build image matrix
  nc <- max(pos[, "x"])
  nc
  nr <- max(pos[, "y"])
  nr
  
  ## init matrix
  m <- matrix(NA, nrow=nr, ncol=nc)
  
  ## fill matrix with intensity values
  for (i in seq(along=x)) {
    ## use mean to be sure we select only one mass per pixel
    w=sum(mass(x[[i]])*intensity(x[[i]]))/sum(intensity(x[[i]]))
    #m[pos[i, "y"], pos[i, "x"]] <- round(w, digits=1)
    m[pos[i, "y"], pos[i, "x"]] <- w
    # m[pos[i, "y"], pos[i, "x"]] <- mean(mass(x[[i]]), na.rm=TRUE)
    
    # taking weighted mean 
    #w <-((intensity(x[[i]]))/(max(intensity(x[[i]]))))
   # m[pos[i, "y"], pos[i, "x"]] <- weighted.mean(sort(intensity(x[[i]])), na.rm=TRUE, w)
  }
  
  ## create mass diff matrix
  m <- m-mean(m, na.rm=TRUE)

  
  ## scale matrix (better contrast)
 # m <- m/max(m, na.rm=TRUE)
  
  ## build title
 # main <- paste(main, " (m/z: ", range[1], "-", range[2], ")", sep="")
  
  ## prepare plot area
 # plot(NA, type="n", xlim=c(1, nc), ylim=c(1, nr), axes=FALSE, asp=1,xlab="", ylab="", main=main)
  
  ## plot image
 # rasterImage(t(m), xleft=1, xright=nc, ybottom=1, ytop=nr, interpolate=FALSE)
  
  #levelplot(t(m), scales=list(tick.number=0), col.regions=rainbow(100, start = 1/6, end = 1), xlab=" ", ylab=" ", main=main)
 
# s <- seq(from=range[1], to=range[2], by=0.1)
 
# levelplot(t(m), scales=list(tick.number=0), xlab=" ", ylab=" ", colorkey = list(at=s, labels=as.character(s)),col.regions = two.colors(n=256, start='blue', end='red', middle='black'), main=main)
 
  
# levelplot(t(m), scales=list(tick.number=0), xlab=" ", ylab=" ", colorkey=list(at=as.numeric(factor(c(seq(from=range[1], to=range[2], by=.1)))),labels=as.character(c( "327.1", "327.2", "327.3", "327.4", "327.5", "327.6", "327.7", "327.8", "327.9"))),col.regions = two.colors(n=256, start='red', end='blue', middle='black'), main=main)
 
 #levelplot(t(m), scales=list(tick.number=0), xlab=" ", ylab=" ", colorkey=list(at=as.numeric(factor(c(seq(from=range[1], to=range[2], by=.1))))),col.regions = two.colors(n=256, start='red', end='blue', middle='black'), main=main)
}
