sizeCaller <-
function(x_abifile,r_abifile=.r_abif,rpks=length(rsizes),
                        rsizes=c(50,75,100,139,150,160,200,NA,300,NA,350,400,450,490,500),
                        ranges=matrix(rep(c(2.3,8.0),sum(npeaks.per.channel)),nrow=sum(npeaks.per.channel),byrow=TRUE),
                        ranges.bp,
                        thres_start=2.0,
                        channels=1:4,
                        npeaks.per.channel=rep(1,length(channels)),
                        locus.names=paste(rep("V",length(channels)),1:sum(npeaks.per.channel),sep=""),
                        strain.name=deparse(substitute(x_abifile))) {
# validation steps
if (length(channels)!=length(npeaks.per.channel)) stop("Lenghts of 'channel' and 'npeaks.per.channel' must be equal!")
if (length(locus.names)!=sum(npeaks.per.channel)) stop("Length of locus.names differs from number of loci!")
# ranges.x holds the final ranges per locus
# variable bp specifies if values are in basepairs and need to be converted (TRUE)
# or if values are in time and can be directly used in seqinr::peakabif (FALSE)
bp <- FALSE
if (missing(ranges.bp)) {
  ranges.x <- ranges
} else {
  # if ranges.bp is present
  # validate ranges.bp
  if (!all(dim(ranges.bp)==c(sum(npeaks.per.channel),2))) {
    warning("ranges.bp could not be used, since ranges were not specified for each locus!")
    ranges.x <- ranges
  } else {
    bp <- TRUE
    ranges.x <- ranges.bp
  }
}
# default is FALSE
standardFound <- FALSE
# if unsuccessful, a matrix of NAs is returned
noPks <- matrix(rep(NA,sum(npeaks.per.channel)),nrow=1,dimnames=list(strain.name,locus.names))
result <- noPks
# standard peaks in x_abifile are located; see function searchStandardPeaks
sPks <- searchStandardPeaks(x_abifile,reference=r_abifile,npks=rpks)
pk <- numeric()
# if all standard peaks have been located the script proceeds
if (all(!is.na(sPks))) {
  standardFound <- TRUE
  # only peaks that are not assigned NA in rsizes are used for the generation of a function
  # converting time to basepairs
  sPks <- sPks[!is.na(rsizes)]
  rsz <- rsizes[!is.na(rsizes)]
  # conversion function time2bp is created
  time2bp <- splinefun(sPks,rsz)
  # inverse conversion function bp2time is created (is only needed if ranges.bp is given)
  bp2time <- splinefun(rsz,sPks)
  # values in ranges.x are converted if ranges.bp is given
  if (bp) ranges.x <- matrix(bp2time(ranges.x)/1000,nrow=nrow(ranges.x),ncol=ncol(ranges.x))
  # channels are analyzed
  for (i in 1:length(channels)) {
    # threshold is sequentially reduced if no peak is found
    for (j in seq(thres_start,0.1,-0.1)) {
      # an error or warning causes pk to be NULL
      tryCatch({ pk <- peakabif(x_abifile,chanel=channels[i],npeak=npeaks.per.channel[i],tmin=min(ranges.x[i,]),tmax=max(ranges.x[i,]),thres=j,fig=FALSE)$maxis },
                error=function(ex){pk<-NULL}, warning=function(ex){pk<-NULL})

      # if correct number of peaks is found ...
      # (NULL has length 0)
      if (length(pk)==npeaks.per.channel[i]) {
        # ... time is converted to basepairs using previously created conversion function
        if (i==1) { r.start <- 1 } else { r.start <- sum(npeaks.per.channel[1:i-1])+1 }
        r.end <- sum(npeaks.per.channel[1:i])
        # the result matrix is filled according to npeaks.per.channel
        result[r.start:r.end] <- time2bp(pk)
        pk <- NULL
        break
      }
    }
  }
}
if (standardFound) {
  return(result)
} else {
  return(noPks)
}
}

