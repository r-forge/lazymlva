normPeaks <-
function(pks) { return(((pks-min(pks))/(max(pks)-min(pks)))[-c(1,length(pks))]) }

