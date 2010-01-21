\name{normPeaks}
\alias{normPeaks}
\title{
normPeaks
}
\description{
The function normalizes a peak vector containing the locations of peaks to make it 
comparable to other peak vectors.
}
\usage{
normPeaks(pks)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pks}{
a numeric vector of length >= 3 containing the locations of peaks in an ABIF file.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
returns a numeric vector of length \code{length(pks)-2} containing values between 0 and 1.
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
\code{seqinr::peakabif}, \code{extractNormPeaks}
}
\examples{
# "normalize" any vector
normPeaks(c(3,5,6,9,10))
}