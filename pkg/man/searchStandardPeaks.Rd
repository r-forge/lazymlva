\name{searchStandardPeaks}
\alias{searchStandardPeaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
searchStandardPeaks(x_abifile, refchannel = 5, reference = .r_abif, npks = 15, x_tmin = 2, x_tmax = 3, x_thres = 0.5, x_step = 0.1, r_tmin = 2.3, r_thres = 0.5, dist_thres = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x_abifile}{
%%     ~~Describe \code{x_abifile} here~~
}
  \item{refchannel}{
%%     ~~Describe \code{refchannel} here~~
}
  \item{reference}{
%%     ~~Describe \code{reference} here~~
}
  \item{npks}{
%%     ~~Describe \code{npks} here~~
}
  \item{x_tmin}{
%%     ~~Describe \code{x_tmin} here~~
}
  \item{x_tmax}{
%%     ~~Describe \code{x_tmax} here~~
}
  \item{x_thres}{
%%     ~~Describe \code{x_thres} here~~
}
  \item{x_step}{
%%     ~~Describe \code{x_step} here~~
}
  \item{r_tmin}{
%%     ~~Describe \code{r_tmin} here~~
}
  \item{r_thres}{
%%     ~~Describe \code{r_thres} here~~
}
  \item{dist_thres}{
%%     ~~Describe \code{dist_thres} here~~
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
