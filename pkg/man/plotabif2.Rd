\name{plotabif2}
\alias{plotabif2}
plotabif2
\title{
plotabif2
}
\description{
This function draws a simple plot of the chromatogram file similar to 
\code{seqinr::plotabif}, yet additionally allows overplotting of several channels.
}
\usage{
plotabif2(abifdata, 
          channels = 1:5, 
          ch.plotcolors = c("blue", "green", "black", "red", "orange"), 
          main = paste("Channels", paste(channels, collapse = ", "), "of", deparse(substitute(abifdata))),
          ylab = "Relative Fluorescent Units (RFU)", 
          xlab = "Time", 
          tmin = 1, 
          tmax = abifdata$Data[["SCAN.1"]], 
          ylim, 
          add.legend = TRUE, 
          ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{abifdata}{
result returned by \code{seqinr::read.abif}
}
  \item{channels}{
the channel numbers (default is 1:5)
}
  \item{ch.plotcolors}{
colors for channels 1:5
}
  \item{main}{
main title, see \code{plot}
}
  \item{ylab}{
label of y-axis, see \code{plot}
}
  \item{xlab}{
label of x-axis, see \code{plot}
}
  \item{tmin}{
unscaled starting time
}
  \item{tmax}{
unscaled ending time
}
  \item{ylim}{
see \code{plot}
}
  \item{add.legend}{
logical variable stating whether a legend is to be added
}
  \item{\dots}{
arguments passed to \code{plot}
}
}
\details{
}
\value{
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
\code{seqinr::plotabif}
}
\examples{
# plot channels of reference abif file r_abif
plotabif2(r_abif)
}
