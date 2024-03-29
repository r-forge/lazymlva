\name{sizeCaller}
\alias{sizeCaller}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
sizeCaller
}
\description{
The function extracts positions of standard peaks from an ABIF file and calculates the sizes of products in specified channels.
}
\usage{
sizeCaller(x_abifile, 
            rsizes = r_sizes, 
            ranges = matrix(rep(c(2.3, 8), sum(npeaks.per.channel)), nrow = sum(npeaks.per.channel), byrow = TRUE), 
            ranges.bp, 
            thres_start = 2, 
            channels = 1:4, 
            npeaks.per.channel = rep(1, length(channels)), 
            locus.names = paste(rep("V", length(channels)), 1:sum(npeaks.per.channel), sep = ""), 
            strain.name = deparse(substitute(x_abifile)), 
            ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x_abifile}{
list generated by \code{seqinr::read.abif} representing the ABIF file
}
  \item{rsizes}{
sizes of reference peaks in basepairs. peaks not to be counted are given as NA.
}
  \item{ranges}{
matrix of time ranges across which peaks are looked for in channels specified 
in \code{channels}
}
  \item{ranges.bp}{
matrix of basepair ranges across which peaks are looked for in channels specified 
in \code{channels}. If given, this matrix overrides matrix \code{ranges}.
}
  \item{thres_start}{
threshhold, below which peaks are looked for (see function \code{seqinr::peakabif})
}
  \item{channels}{
channel numbers containing products of unknown sizes
}
  \item{npeaks.per.channel}{
array containing number of peaks per channel (defaults to one peak per channel)
}
  \item{locus.names}{
array of locus names
}
  \item{strain.name}{
name of strain to be analyzed
}
  \item{\dots}{
arguments passed on to \code{searchStandardPeaks(...)}
}
}
\details{
}
\value{
matrix returning the peak sizes in basepairs of each locus
}
\references{
package \code{seqinr} by Charif, D. and Lobry, J.R. (2007)
}
\author{
Johannes Elias
}
\note{
}
\section{Warning }{
Warnings are issued if sizes in one or several loci cannot be calculated or if
\code{searchStandardPeaks} is unable to extract the standard ladder from 
\code{x_abifile}.
}
\seealso{
\code{searchStandardPeaks}, \code{r_sizes}
}
\examples{
# returns the sizes of peaks in the standard ABIF file
sizeCaller(r_abif)
}