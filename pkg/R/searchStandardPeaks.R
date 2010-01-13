searchStandardPeaks <-
function(x_abifile,refchannel=5,reference=.r_abif,npks=15,x_tmin=2.0,x_tmax=3.0,x_thres=0.5,x_step=0.1,r_tmin=2.3,r_thres=0.5,dist_thres=0.01) {
# default result is FALSE
result_found <- FALSE
# in case of error or unsuccessful retrieval of standard peaks an array of NAs is returned
xPks_error <- rep(NA,npks)
# peaks of reference are read in using the function seqinr::peakabif
refPks <- peakabif(reference,chanel=refchannel,npeak=npks,tmin=r_tmin,thres=r_thres,fig=FALSE)$maxis
# location of size standard peaks in x_abifile is attempted
tryCatch({
  # threshold is decreased sequentially by x_step (see seqinr::peakabif)
  for (j in seq(x_thres,0.1,-x_step)) {
    if (result_found) break
    # furthermore, tmin is increased sequentially by x_step (see seqinr::peakabif)
    for (i in seq(x_tmin,x_tmax,x_step)) {
      if (result_found) break
      tryCatch({xPks <- peakabif(x_abifile,chanel=refchannel,npeak=npks,tmin=i,thres=j,fig=FALSE)$maxis},
                error=function(ex){xPks <- numeric()},warning=function(ex){xPks <- numeric()})
      # successful retrieval requires same number of peaks as specified in npks ...
      if ((length(xPks)==npks) && (all(!is.na(xPks)))) {
        # ... and sufficient similarity of peaks with reference
        if (dist(rbind(normPeaks(refPks),normPeaks(xPks)))<dist_thres) {
          result_found <- TRUE
        } else {
          result_found <- FALSE
        }
      }
    }
  }
  # any error causes an array of NAs to be returned
}, error=function(ex){ return(xPks_error) })
if (result_found) {
  return(xPks)
} else {
  return(xPks_error)
}
}

