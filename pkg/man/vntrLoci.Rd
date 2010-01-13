\name{vntrLoci}
\alias{vntrLoci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
vntrLoci(lz.map = lazy.map, file.ending = ".fsa", size.only = FALSE, filename.sep = "_", export.path, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lz.map}{
%%     ~~Describe \code{lz.map} here~~
}
  \item{file.ending}{
%%     ~~Describe \code{file.ending} here~~
}
  \item{size.only}{
%%     ~~Describe \code{size.only} here~~
}
  \item{filename.sep}{
%%     ~~Describe \code{filename.sep} here~~
}
  \item{export.path}{
%%     ~~Describe \code{export.path} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(lz.map=lazy.map,file.ending=".fsa",size.only=FALSE,filename.sep="_",export.path,...) {
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

# let´s go through each series in turn
for (s in series) {
  #cat(paste("Analyzing series '",s,"'...\n",sep=""))
  # retrieve the filenames of the series in question
  # it is assumed that filenames are of the form "XXX_S1.fsa", where XXX represents a strain specific string,
  # "_" is filename.separator, "S1" is the string representation of a series (e.g. series 1), and ".fsa" is the file ending.
  seriesFiles <- dir(...,pattern=paste(s, file.ending, sep=""))
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
      for (k in 1:(nrows(temp.loci.info))-1) if (temp.loci.info[k,"to"] > temp.loci.info[k+1,"from"]) warning(paste("Ranges of loci in series",s,",channel",ch,"overlap!"))
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
    seriesResults <- rbind(seriesResults,sizeCaller(read.abif(f),
                                          ranges.bp=locus.ranges,
                                          channels=ch,
                                          npeaks.per.channel=npks.ch,
                                          locus.names=l.names,
                                          strain.name=strsplit(f,filename.sep)[[1]][1]))
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
      temp.df <- subset(lz.map,locus==n)
      # iterate through rows of a series table
      for (r in 1:nrow(result[[i]])) {
        # look up the number of repeats in temp.df
        for (j in 1:nrow(temp.df)) {
          if ((result[[i]][r,n] >= temp.df[j,"from"]) && (result[[i]][r,n] <= temp.df[j,"to"])) {
            result[[i]][r,n] <- temp.df[j,"repeats"]
          }
        }
      }
    }
  }
  }
#export results if export.path was specified
if (!missing(export.path)) {
  #cat("Exporting result into files...\n")
  if (!file.info(export.path)[["isdir"]]) {
    warning(paste("Results could not be exported, since export.path",export.path,"does not point to a directory!"))
  } else {
    # iterate through series and write one csv file per series
    for (i in 1:length(result)) {
      write.csv2(result[[i]],file=paste(export.path,"/",names(result)[i],".csv",sep=""))
    }
  }
  }
return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
