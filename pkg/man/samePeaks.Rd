\name{samePeaks}
\alias{samePeaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
samePeaks
}
\description{
The function determines whether peaks in two electropherogram files (ABIF files) represent the same size standard or not. 
ABIF files have to be imported using the function \code{seqinr::read.abif}.
}
\usage{
samePeaks(e1, e2, ch = 5, npks = 15, e1_tmin = 2.3, e1_thres = 0.5, e2_tmin = e1_tmin, e2_thres = e1_thres, dist_thres = 0.01, s_warn = -1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{e1}{
electropherogram 1
}
  \item{e2}{
electropherogram 1
}
  \item{ch}{
channel number; see (see function \code{seqinr::peakabif})
}
  \item{npks}{
expected number of peaks
}
  \item{e1_tmin}{
starting time for the time axis for e1 (see function \code{seqinr::peakabif})
}
  \item{e1_thres}{
threshold, above which peaks are looked for in e1 (see function \code{seqinr::peakabif})
}
  \item{e2_tmin}{
starting time for the time axis for e2 (see function \code{seqinr::peakabif})
}
  \item{e2_thres}{
threshold, above which peaks are looked for in e2 (see function \code{seqinr::peakabif})
}
  \item{dist_thres}{
euclidian distance, below which function returns TRUE
}
  \item{s_warn}{
sets the handling of warning messages; see \code{options()}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
logical() value TRUE or FALSE
}
\references{
package \code{seqinr} by Charif, D. and Lobry, J.R. (2007)
}
\author{
Johannes Elias
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
seqinr::peakabif
}
\examples{
# the following returns TRUE 
samePeaks(r_abif,r_abif)
}

