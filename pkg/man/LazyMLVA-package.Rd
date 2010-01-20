\name{LazyMLVA-package}
\alias{LazyMLVA-package}
\alias{LazyMLVA}
\docType{package}
\title{
LazyMLVA
}
\description{
The package automates assignment of MLVA types. MLVA stands for 'Multiple-Locus-VNTR-Analysis', which is a typing approach for bacterial pathogens. 
It extracts the size standard from ABIF electropherogram files and automatically looks for products of unknown sizes in the specified channels. 
Further to estimating the sizes of found products, it converts the sizes to repeat lengths according to data specified in a 'lazy map'. 
The package depends on seqinr.

Function \code{sizeCaller} extracts the size ladder from the electropherogram and handles conversion of time to basepairs (see \code{?sizeCaller} for details).

Function \code{vntrLoci} analyzes files within a specified directory according to details given in a 'lazy map' (see \code{?vntrLoci} and \code{?lazy.map} for details regarding the function and the 'lazy map', respectively). 
}
\details{
\tabular{ll}{
Package: \tab LazyMLVA\cr
Type: \tab Package\cr
Version: \tab 0.0-3\cr
Date: \tab 2010-01-20\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
}
\author{
Johannes Elias

Maintainer: Johannes Elias <jelias@hygiene.uni-wuerzburg.de>
}
\references{
\code{citation(package="seqinr")}

Applied Biosystems. (2009) Applied Biosystems Genetic Analysis Data File Format [Online].  Available from: http://www.appliedbiosystems.com/support/software_community/ABIF_File_Format.pdf (Accessed: 14 January 2010). 
}
\keyword{ MLVA, bacterial typing, VNTR }
\seealso{
\code{seqinr::peakabif}, \code{seqinr::read.abif} 
}
\examples{
}
