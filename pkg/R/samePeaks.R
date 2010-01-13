samePeaks <-
function(r1,r2,ch=5,npks=15,r1_tmin=2.3,r1_thres=0.5,r2_tmin=r1_tmin,r2_thres=r1_thres,dist_thres=0.01,s_warn=-1) {
# default result is FALSE
result <- FALSE
# retrieve current options for warning messages
warn.old <- getOption("warn")
# set new setting for warning messages
options(warn=s_warn)
# peak locations are retrieved using the function seqinr::peakabif
r1.pks <- peakabif(r1,chanel=ch,npeak=npks,tmin=r1_tmin,thres=r1_thres,fig=FALSE)$maxis
# if required number of peaks cannot be retrieved, function is aborted
if ((length(r1.pks)!=npks) || !(all(!is.na(r1.pks)))) {
  stop(cat("Number of peaks (",npks,") cannot be retrieved from first reference. \n",sep=""))
}
r2.pks <- peakabif(r2,chanel=ch,npeak=npks,tmin=r2_tmin,thres=r2_thres,fig=FALSE)$maxis
if ((length(r2.pks)!=npks) || !(all(!is.na(r2.pks)))) {
  stop(cat("Number of peaks (",npks,") cannot be retrieved from second reference. \n",sep=""))
}
# result is TRUE if, euclidian distance between are smaller than dist_thres
if (dist(rbind(normPeaks(r1.pks),normPeaks(r2.pks)))<=dist_thres) result <- TRUE
# handling of warning messages is reset to previous value
options(warn=warn.old)
return(result)
}

