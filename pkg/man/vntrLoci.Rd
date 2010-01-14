\name{vntrLoci}
\alias{vntrLoci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
vntrLoci
}
\description{
This function automates the analysis of ABIF files within a directory according to specifications in a 'lazy map'.
}
\usage{
vntrLoci(lz.map = lazy.map, file.ending = ".fsa", size.only = FALSE, filename.sep = "_", export.path, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lz.map}{
a 'lazy map' is a table listing all loci that are looked for, their basepair ranges, corresponding repeat numbers, etc.
}
  \item{file.ending}{
file ending of ABIF files
}
  \item{size.only}{
if TRUE, only sizes in basepairs are given
}
  \item{filename.sep}{
ABIF filenames are assumed to be of the form NAME + filename.sep + SERIES + file.ending; e.g. 1087_S1.fsa specifies series 'S1' of strain '1087'
}
  \item{export.path}{
if given, result is additionally written to csv files in this directory
}
  \item{\dots}{
arguments passed on to \code{dir(...)} like \code{path="."}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a list of matrices (one matrix per series) containing the repeat numbers of each locus and strain 
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# process ABIF files in directory 'abif'
# RUN AFTER INSTALLATIOM
# load(system.file("extdata/lazyFiles.RData",package="LazyMLVA"))
# vntrLoci(path=system.file("abif",package="LazyMLVA"),full.names=TRUE)

# show sizes only
# RUN AFTER INSTALLATIOM
# vntrLoci(path=system.file("abif",package="LazyMLVA"),size.only=TRUE)