\name{extractNormPeaks}
\alias{extractNormPeaks}
\title{
extractNormPeaks
}
\description{
This function just combines extraction of peaks with \code{seqinr::peakabif} and 
normalization with \code{normPeaks}. It just makes the creation of normalized peaks
a bit easier.
}
\usage{
extractNormPeaks(abifdata, chanel = 5, npeak = 15, tmin = 2.3, thres = 0.5, fig = FALSE, ...)
}
\arguments{
  \item{abifdata}{
the result returned by \code{seqinr::read.abif} representing the ABIF file
}
  \item{chanel}{
the channel number of the size standard
}
  \item{npeak}{
number of peaks to extract
}
  \item{tmin}{
scaled starting time for the time axis (see \code{seqinr::peakabif})
}
  \item{thres}{
scaled threshold value (see \code{seqinr::peakabif})
}
  \item{fig}{
logical: should localized peaks be plotted (see \code{seqinr::peakabif})
}
  \item{\dots}{
arguments passed to \code{seqinr::peakabif}
}
}
\details{
}
\value{
a numeric variable containing \code{npeak - 2} values
}
\references{
package \code{seqinr} by Charif, D. and Lobry, J.R. (2007)
}
\author{
Johannes Elias
}
\note{
}

\seealso{
\code{seqinr::peakabif}, \code{normPeaks}
}
\examples{
extractNormPeaks(r_abif)
}
