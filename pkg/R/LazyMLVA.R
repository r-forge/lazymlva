# Functions for LazyMLVA
# Author: Johannes Elias, jelias@hygiene.uni-wuerzburg.de

#
# Function: normPeaks
#
# This function normalizes the peak vector and removes extreme peaks
# which are always 0 and 1, respectively
#

normPeaks <- function(pks) { 
if (length(pks) < 3) stop("pks has to contain 3 peaks at least!")
return(((pks-min(pks))/(max(pks)-min(pks)))[-c(1,length(pks))]) 
}

##         TEST VERSION      R E M O V E   F R O M   P A C K A G E   S O U R C E
#r_abif <- read.abif("abif/1087_S1.fsa")  
#r_npk <- normPeaks(peakabif(r_abif,chanel=5,npeak=15,tmin=2.3,thres=0.5,fig=FALSE)$maxis)
#lazy.map <- read.csv2("lazy_map.csv")
#r_sizes <- c(50,75,100,139,150,160,200,NA,300,NA,350,400,450,490,500)
#save(r_abif,lazy.map,r_npk,r_sizes,file="sysdata.rda")
##         TEST VERSION      R E M O V E   F R O M   P A C K A G E   S O U R C E

#         PACKAGE VERSION   R E P L A C E   TEST VERSION   B Y   T H I S
## Generate r_npk from r_abif (r_abif is loaded via sysdata.rda in source directory "R")
#r_npk <- normPeaks(peakabif(r_abif,chanel=5,npeak=15,tmin=2.3,thres=0.5,fig=FALSE)$maxis)
#         PACKAGE VERSION   R E P L A C E   TEST VERSION   B Y   T H I S

#
# Function: extractNormPeaks
#
# This function combines extraction of peaks with seqinr::peakabif and normalization
# with normPeak
#
# function arguments (see peakabif):
# abifdata        the result returned by �read.abif� representing the ABIF file
# chanel          the channel number
# npeak           number of peaks to be extracted
# tmin            scaled starting time for the time axis
# thres           scaled threshold value
# fig             logical: should localized peaks be plotted
# ...             arguments passed to seqinr::peakabif
#

extractNormPeaks <- function(abifdata,chanel=5,npeak=15,tmin=2.3,thres=0.5,fig=FALSE,...) {
result <-  normPeaks(peakabif(abifdata=abifdata,chanel=chanel,npeak=npeak,tmin=tmin,thres=thres,fig=fig,...)$maxis)
return(result)
}

#
# Function: samePeaks
#
# This function retrieves peak locations in two reference electropherograms
# r1 and r2 and determines whether they are the same. It returns boolean values
# TRUE or FALSE.
# This function is useful for the search of reference peaks within an electropherogram.
#
# function arguments:
# e1          electropherogram 1
# e2          electropherogram 2
# ch          channel: default is 5
# npks        number of peaks to be retrieved from each electropherogram
# e1_tmin     starting time for the time axis for e1 (see function seqinr::peakabif)
# e1_thres    threshold, above which peaks are looked for in e1 (see function seqinr::peakabif)
# e2_tmin     starting time for the time axis for e2 (see function seqinr::peakabif)
# e2_thres    threshold, above which peaks are looked for in e2 (see function seqinr::peakabif)
# dist_thres  euclidian distance, below which function returns TRUE
# s_warn      sets the handling of warning messages; default: warnings ignored (see options())
#

samePeaks <- function(e1,e2,ch=5,npks=15,e1_tmin=2.3,e1_thres=0.5,e2_tmin=e1_tmin,e2_thres=e1_thres,dist_thres=0.01,s_warn=-1) {
# default result is FALSE
result <- FALSE
# validate npks
if (npks < 3) stop("npks has to be greater than 2!")
# retrieve current options for warning messages
warn.old <- getOption("warn")
# set new setting for warning messages
options(warn=s_warn)
# peak locations are retrieved using the function seqinr::peakabif
e1.pks <- peakabif(e1,chanel=ch,npeak=npks,tmin=e1_tmin,thres=e1_thres,fig=FALSE)$maxis
# if required number of peaks cannot be retrieved, function is aborted
if ((length(e1.pks)!=npks) || !(all(!is.na(e1.pks)))) {
  stop(cat("Number of peaks (",npks,") cannot be retrieved from first reference. \n",sep=""))
}
e2.pks <- peakabif(e2,chanel=ch,npeak=npks,tmin=e2_tmin,thres=e2_thres,fig=FALSE)$maxis
if ((length(e2.pks)!=npks) || !(all(!is.na(e2.pks)))) {
  stop(cat("Number of peaks (",npks,") cannot be retrieved from second reference. \n",sep=""))
}
# result is TRUE if, euclidian distance between are smaller than dist_thres
if (dist(rbind(normPeaks(e1.pks),normPeaks(e2.pks)))<=dist_thres) result <- TRUE
# handling of warning messages is reset to previous value
options(warn=warn.old)
return(result)
}

#
# Function: searchStandardPeaks
#
# This function retrieves locations of the internal standard size ladder.
# The default internal size ladder is GS500LIZ, which is represented 
# in reference r_npk. Nevertheless, any size ladder can be chosen by specifying
# reference accordingly.
#
# function arguments:
# x_abifile   list generated by seqinr::read.abif representing ABI file in .fsa format
# reference   contains normalized peak positions extracted with seqinr::peakabif 
#             generate your own with normPeaks + seqinr::peakabif:
#             e.g. reference <- normPeaks(peakabif(read.abif("YOUR_REFERENCE_FILE.fsa"),chanel=5,npeak=15,tmin=2.3,thres=0.5,fig=FALSE)$maxis)
# npks        number of peaks to be retrieved from x_abifile
# rsizes      sizes of reference peaks in basepairs. peaks not to be counted are to be given as NA.
# x_tmin      starting time for the time axis for x_abifile (see function seqinr::peakabif)
# x_tmax      ending time for the time axis for x_abifile (see function seqinr::peakabif)
# x_thres     threshold, above which peaks are looked for in x_abifile (see function seqinr::peakabif)
# x_step      step, by which x_thres is decreased and x_tmin increased until standard peaks are retrieved
# dist_thres  euclidian distance of peak vectors, below which identity of references in x_abifile and reference is assumed (see function samePeaks)
# plotit      logical value stating whether a time->bp plot should be shown
# ...         arguments passed to plot
#

searchStandardPeaks <- function(x_abifile,
                                refchannel=5,
                                reference=r_npk,
                                npks=15,
                                rsizes=r_sizes,
                                x_tmin=2.0,
                                x_tmax=3.0,
                                x_thres=0.5,
                                x_step=0.1,
                                dist_thres=0.01,
                                plotit=FALSE,
                                ...) {
# default result is FALSE
result_found <- FALSE
# in case of error or unsuccessful retrieval of standard peaks an array of NAs is returned
xPks_error <- rep(NA,npks)
# location of size standard peaks in x_abifile is attempted
tryCatch({
  # threshold is decreased sequentially by x_step (see seqinr::peakabif)
  for (j in seq(x_thres,0.1,-x_step)) {
    if (result_found) break
    # furthermore, tmin is increased sequentially by x_step (see seqinr::peakabif) until x_tmax is reached
    for (i in seq(x_tmin,x_tmax,x_step)) {
      if (result_found) break
      xPks <- numeric()
#      try({xPks <- peakabif(x_abifile,chanel=refchannel,npeak=npks,tmin=i,thres=j,fig=FALSE)$maxis}, TRUE) # try only suppresses error messages, not warnings
      tryCatch({xPks <- peakabif(x_abifile,chanel=refchannel,npeak=npks,tmin=i,thres=j,fig=FALSE)$maxis},
                error=function(ex){},warning=function(ex){})
      # successful retrieval requires same number of peaks as specified in npks ...
      if ((length(xPks)==npks) && (all(!is.na(xPks)))) {
        # ... and sufficient similarity of peaks with reference
        if (dist(rbind(reference,normPeaks(xPks))) < dist_thres) {
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
  if (plotit) {
    bp <- rsizes[!is.na(rsizes)]
    time <- xPks[!is.na(rsizes)]
    plot(time,bp,main=paste("Conversion function time->bp",deparse(substitute(x_abifile))),...)
    time.x <- seq(min(time),max(time),length.out=100)
    lines(time.x,splinefun(time,bp)(time.x),col="red")
  }
  return(xPks)
} else {
  warning(paste("standard peaks for ",deparse(substitute(x_abifile))," could not be retrieved."))
  return(xPks_error)
}
}

#
# Function: sizeCaller
#
# This function determines the sizes of products in selected channels 
#
# function arguments:
# x_abifile           list generated by seqinr::read.abif representing ABI file in .fsa format
# rsizes              sizes of reference peaks in basepairs. peaks not to be counted are to be given as NA.
# ranges              time range across which peaks are looked for in channels 
# ranges.bp           bp ranges across which peaks are looked for (values are first converted to time ranges)
# thres_start         threshhold, below which peaks are looked for
# channels            number of channels containing products of unknown sizes
# npeaks.per.channel  array containing number of peaks per channel (default: one peak per channel)
# locus.names         array of locus names
# strain.name         name of strain to be analyzed
# ...                 parameters passed on to searchStandardPeaks
#

sizeCaller <- function(x_abifile,
                        rsizes=r_sizes,
                        ranges=matrix(rep(c(2.3,8.0),sum(npeaks.per.channel)),nrow=sum(npeaks.per.channel),byrow=TRUE),
                        ranges.bp,
                        thres_start=2.0,
                        channels=1:4,
                        npeaks.per.channel=rep(1,length(channels)),
                        locus.names=paste(rep("V",length(channels)),1:sum(npeaks.per.channel),sep=""),
                        strain.name=deparse(substitute(x_abifile)),
                        ...) {
# validation steps
if (length(channels)!=length(npeaks.per.channel)) stop("\n  Lenghts of 'channel' and 'npeaks.per.channel' must be equal!")
if (length(locus.names)!=sum(npeaks.per.channel)) stop("\n  Length of locus.names differs from number of loci!")
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
# catch warning
sPks <- NA
tryCatch({sPks <- searchStandardPeaks(x_abifile,rsizes=rsizes,...)},warning=function(ex){})
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
      pk <- NULL # accept no warnings or errors
      tryCatch({ pk <- peakabif(x_abifile,chanel=channels[i],npeak=npeaks.per.channel[i],tmin=min(ranges.x[i,]),tmax=max(ranges.x[i,]),thres=j,fig=FALSE)$maxis },
                error=function(ex){}, warning=function(ex){})
      # if correct number of peaks is found ...
      # (btw NULL has length 0)
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
# generate variables holding faulty loci
f.loxnx <- character()
# issue results
if (standardFound) {
  # check if product sizes are available for all loci
  if (any(is.na(result))) {
    f.loxnx <- colnames(result)[is.na(result)]
    if (length(f.loxnx) > 1) {
      warning(paste("Product sizes pertaining to loci ",paste(f.loxnx,collapse=", "),", in file ",strain.name," could not be determined.",sep=""))
    } else {
      warning(paste("Product size pertaining to locus ",f.loxnx," in file ",strain.name," could not be determined.",sep=""))
    }
  }
  return(result)
} else {
  # format warning
  warning(paste("Standard peaks of file ",strain.name,", containing loci ",paste(locus.names,collapse=", "),", could not be retrieved.",sep=""))
  # return result
  return(noPks)
}
}

#
# Function: vntrLoci
#
# This function calls sizeCaller for all found files of a series and converts product sizes to repeat numbers according to
# a lazy map. data.frame lazy.map is loaded via sysdata.rda in source directory "R".
#
# function arguments:
# p                   path to directory containing abif files to be analyzed
# lz.map              a lazy map containing locus names, channel numbers, lengths of possible products
#                     see lazy.map
# file.ending         default is ".fsa"
# size.only           defaults to FALSE; if TRUE, only product sizes in basepairs are given
# filename.sep        defaults to "_". The character string in the filename before this separator is assumed to be the strain name
# wide.table          logical variable stating whether conversion of results to a wide table should be attempted
# ...                 arguments passed to sizeCaller
#

vntrLoci <- function(p,lz.map=lazy.map,file.ending=".fsa",size.only=FALSE,filename.sep="_",wide.table=TRUE,...) {
# parameter lz.map is validated
#cat(paste("\nValidating lz.map '",lz.map,"'...\n",sep=""))
if (class(lz.map) != "data.frame") stop("lz.map needs to be a data.frame!")
if (all(names(lz.map) != c("locus","from","to","repeats","series","channel"))) {
  stop("lz.map column names are not correct!") }
#cat(paste("...lz.map '",lz.map,"' seems ok.\n",sep=""))

# create list holding result
result <- list()

# convert several columns to factor if necessary
if (!is.factor(lz.map$locus)) lz.map$locus <- as.factor(lz.map$locus)
if (!is.factor(lz.map$series)) lz.map$series <- as.factor(lz.map$series)
if (!is.factor(lz.map$channel)) lz.map$channel <- as.factor(lz.map$channel)

# extract series from lz.map
series <- names(table(lz.map$series))

# estimate length of procedure: n.files is the number of files to be analyzed
n.files <- 0 
for (s in series) n.files <- n.files + length(dir(path=p,pattern=paste(s, file.ending, sep="")))
# generate txtProgressBar
tpb <- txtProgressBar(style=3,char=">")

# let�s go through each series in turn
for (s in series) {
  #cat(paste("Analyzing series '",s,"'...\n",sep=""))
  # retrieve the filenames of the series in question
  # it is assumed that filenames are of the form "XXX_S1.fsa", where XXX represents a strain specific string,
  # "_" is filename.separator, "S1" is the string representation of a series (e.g. series 1), and ".fsa" is the file ending.
  seriesFiles <- dir(path=p,pattern=paste(s, file.ending, sep=""),full.names=TRUE)
  # create a data.frame that holds the results of this series
  seriesResults <- data.frame()
  # generate a numeric variable that holds required channel numbers (these are to be fed to seqinr::peakabif)
  ch <- numeric()
  # generate a character variable that holds the locus names of a channel
  l.names <- character()
  # convert channels to numbers (so they can safely be fed to seqinr::peakabif)
  # if conversion fails the script is stopped, as there is no sense in proceeding
  tryCatch({ ch <- as.numeric(names(table(lz.map$channel))) },
          error=function(ex){ stop(paste("Column 'channel' of table",deparse(substitute(lz.map)),"has to be convertible to numeric!")) },
          warning=function(ex){ stop(paste("Column 'channel' of table",deparse(substitute(lz.map)),"has to be convertible to numeric!")) })
  # generate a variable that holds the number of peaks per channel (see sizeCaller)
  npks.ch <- numeric(length(ch))
  # generate a matrix that holds the ranges of each locus, see sizeCaller
  locus.ranges <- data.frame()
  # let's iterate through each channel of the series
  for (i in 1:length(ch)) {
    #cat(paste("\tchannel ",ch[i],"\n",sep=""))
    # determine number of entries (=loci) per channel (this will be 1 in most instances)
    temp.tab <- table(lz.map$locus,lz.map$series,lz.map$channel)[,s,ch[i]]
    # select loci, which have more than one entry
    temp.l.names <- names(subset(temp.tab,temp.tab!=0))
    # generate a data.frame that holds information about the loci of the channel
    temp.loci.info <- data.frame()
    # temp.locus.range holds the size range of the locus, i.e. products of this locus have sizes within this range
    temp.locus.range <- numeric()
    # if more than 1 locus is found associated with this channel (this will rarely be the case) ...
    if (length(temp.l.names)>1) {
      # save number of peaks (=products) in this channel
      npks.ch[i] <- length(temp.l.names)
      # extract information for each locus
      for (ln in temp.l.names) {
        # compute range of each locus
        temp.locus.range <- range(subset(lz.map,lz.map$locus==ln)[,c("from","to")])
        temp.loci.info <- rbind(temp.loci.info,data.frame(locus=ln,from=temp.locus.range[1],to=temp.locus.range[2]))
      }
      # ... order locus names according to product size
      temp.loci.info <- temp.loci.info[order(temp.loci.info$from),]
      # issue a warning if ranges overlap
      for (k in 1:(nrow(temp.loci.info))-1) if (temp.loci.info[k,"to"] > temp.loci.info[k+1,"from"]) warning(paste("Ranges of loci in series",s,",channel",ch,"overlap!"))
      # save locus names in l.names in correct order
      l.names <- c(l.names,temp.loci.info$locus)
      # attach locus ranges to locus.ranges in correct order
      locus.ranges <- rbind(locus.ranges,temp.loci.info[,c("from","to")]) 
    } else {
      # number of peaks (=products) in this channel is 1
      npks.ch[i] <- 1
      # compute range for locus
      temp.locus.range <- range(subset(lz.map,lz.map$locus==temp.l.names)[,c("from","to")])
      # save locus name in l.names
      l.names <- c(l.names,temp.l.names)
      # attach locus range to locus.ranges
      locus.ranges <- rbind(locus.ranges,data.frame(temp.locus.range[1],temp.locus.range[2]))
    }
  }
  #convert locus.ranges to a matrix
  locus.ranges <- as.matrix(locus.ranges)
  # send all files of this series to sizeCaller
  for (f in seriesFiles) {
    sn <- strsplit(f,filename.sep)[[1]][1]
    if (grepl("/",sn)) sn <- substr(sn,max(gregexpr("/",sn)[[1]])+1,nchar(sn))
    seriesResults <- rbind(seriesResults,sizeCaller(read.abif(f),
                                          ranges.bp=locus.ranges,
                                          channels=ch,
                                          npeaks.per.channel=npks.ch,
                                          locus.names=l.names,
                                          strain.name=sn,
                                          ...))
    #update txtProgressBar
    setTxtProgressBar(tpb, getTxtProgressBar(tpb) + 1/n.files)
  }
  result[[s]] <- seriesResults
}
# convert sizes to number of repeats if size.only==FALSE
if (!size.only) {
  #cat("Converting sizes to repeat lengths...\n")
  # iterate through series
  for (i in 1:length(result)) {
    # iterate through loci of a series
    for (n in colnames(result[[i]])) {
      # copy data from lz.map pertaining to a locus into a temporary data.frame
      temp.df <- subset(lz.map,lz.map$locus==n)
      # iterate through rows of a series table
      for (r in 1:nrow(result[[i]])) {
        # detect whether a conversion from basepairs to repeat number was successful
        conversion.found <- FALSE
        # look up the number of repeats in temp.df
        for (j in 1:nrow(temp.df)) {
          if (!is.na(result[[i]][r,n])) {
            if ((result[[i]][r,n] >= temp.df[j,"from"]) && (result[[i]][r,n] <= temp.df[j,"to"])) {
              conversion.found <- TRUE
              result[[i]][r,n] <- temp.df[j,"repeats"]
            }
          }
        }
        # issue a message if conversion was not successful
        if (!conversion.found && !is.na(result[[i]][r,n])) {
          warning(paste("The lazy.map provided did not allow conversion to repeat number in file ",rownames(result[[i]])[r],", locus ",n,". Only rounded product size in basepairs is given.",sep=""))
        }
      }
    }
    # round 
    result[[i]] <- round(result[[i]])
  }
}

# convert result (which is a list at this point) to a nice table with strain numbers as rows and loci as columns
# this only works if the same strains have been analyzed in each of the series
if (wide.table) {
  # check whether conversion is possible
  # iterate through tables in result
  same <- TRUE
  result.w <- result[[1]]
  for (i in 1:length(result)) if (!identical(rownames(result[[1]]),rownames(result[[i]]))) same <- FALSE
  if (same && (length(result) > 1)) {
    for (i in 2:length(result)) {
      result.w <- cbind(result.w,result[[i]])  
    }
    result <- result.w
  } else {
    warning("Conversion to wide table was not possible since files of different strains were analyzed in the series. Sorry.")
  }
}
# close txtProgressBar
close(tpb)
return(result)
}

#
# Function: plotabif2
#
# This function draws a simple plot of the chromatogram file similar to 
# seqinr::plotabif, yet additionally allows overplotting of several channels
#
# function arguments:
# abifdata        result returned by seqinr::read.abif
# channels        the channel numbers (default is 1:5)
# ch.names        names, by which channels are represented in abifdata
# ch.plotcolors   colors for channels 1:5
# ylab            label of y-axis
# xlab            label of x-axis
# tmin            unscaled starting time  
# tmax            unscaled ending time 
# ylim            displayed range in y-axis (see plot())
# add.legend      logical variable stating whether a legend is to be added
#

plotabif2 <- function(abifdata,
                      channels=1:5,
                      ch.plotcolors=c("blue","green","black","red","orange"),
                      main=paste("Channels",paste(channels,collapse = ", "),"of",deparse(substitute(abifdata))),
                      ylab = "Relative Fluorescent Units (RFU)",
                      xlab = "Time",
                      tmin = 1,
                      tmax = abifdata$Data[["SCAN.1"]],
                      ylim,
                      add.legend=TRUE,
                      ...) {
 # validate channels
 if (!is.numeric(channels)) stop("Please enter a numeric variable for the channels!")
 if (length(channels) > 5 || min(channels) < 1 || max(channels) > 5) stop("Please enter up to 5 channels designated 1-5!")
 # create ch.names
 ch.names = c(1:4,105)
 # if ylim is missing choose according to maximal expansion in all channels
 if (missing(ylim)) {
   ylim.min <- numeric()
   ylim.max <- numeric()
   # iterate through channels
   for (i in 1:length(channels)) {
    ylim.min <- c(ylim.min,min(abifdata$Data[[paste("DATA",ch.names[channels[i]],sep=".")]]))
    ylim.max <- c(ylim.max,max(abifdata$Data[[paste("DATA",ch.names[channels[i]],sep=".")]]))
   }
   ylim <- c(min(ylim.min),max(ylim.max))
 }
 
 # create variable xlim
 xlim <- c(tmin,tmax)
 
 # retrieve dye names as stated in abifdata
 dye.names <- character(length(channels))
 for (i in 1:length(channels)) dye.names[i] <- paste(abifdata$Data[[paste("DyeN", channels[i], sep = ".")]]," [",channels[i],"]",sep="")
 
 for (i in 1:length(channels)) {
  if (i==1) {
    plot(abifdata$Data[[paste("DATA",ch.names[channels[i]],sep=".")]],type="l",col=ch.plotcolors[channels[i]],xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,main=main,...)
  } else {
    lines(abifdata$Data[[paste("DATA",ch.names[channels[i]],sep=".")]],col=ch.plotcolors[channels[i]])
  }
 }
 
 # add a nice legend
 if (add.legend) legend(x="topright",legend=dye.names,col=ch.plotcolors[channels],lty=1,bty="n")                     
}                      