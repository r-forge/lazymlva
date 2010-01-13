\name{samePeaks}
\alias{samePeaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
samePeaks(r1, r2, ch = 5, npks = 15, r1_tmin = 2.3, r1_thres = 0.5, r2_tmin = r1_tmin, r2_thres = r1_thres, dist_thres = 0.01, s_warn = -1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r1}{
%%     ~~Describe \code{r1} here~~
}
  \item{r2}{
%%     ~~Describe \code{r2} here~~
}
  \item{ch}{
%%     ~~Describe \code{ch} here~~
}
  \item{npks}{
%%     ~~Describe \code{npks} here~~
}
  \item{r1_tmin}{
%%     ~~Describe \code{r1_tmin} here~~
}
  \item{r1_thres}{
%%     ~~Describe \code{r1_thres} here~~
}
  \item{r2_tmin}{
%%     ~~Describe \code{r2_tmin} here~~
}
  \item{r2_thres}{
%%     ~~Describe \code{r2_thres} here~~
}
  \item{dist_thres}{
%%     ~~Describe \code{dist_thres} here~~
}
  \item{s_warn}{
%%     ~~Describe \code{s_warn} here~~
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
