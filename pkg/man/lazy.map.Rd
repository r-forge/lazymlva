\name{lazy.map}
\alias{lazy.map}
\docType{data}
\title{
lazy.map
}
\description{
Table containing locus names, possible product sizes (from-to), corresponding repeat numbers, series designation, and channel number (in ABIF file). 
A series corresponds to a Multiplex PCR run, i.e. series "S1" refers to amplification of loci "V3-2", "V7-1", "V13-1", and "V21-2".
}
\usage{lazy.map}
\format{
  A data frame with 53 rows and 6 columns.
  \describe{
    \item{\code{locus}}{a factor representing loci \code{V13-1} \code{V21-2} \code{V3-2} \code{V4-5} \code{V6-1} \code{V7-1} \code{V7-2} \code{V9-1}}
    \item{\code{from}}{a numeric vector describing the lowest product size in basepairs corresponding to the repeat number in column \code{repeats}}
    \item{\code{to}}{a numeric vector describing the highest product size in basepairs corresponding to the repeat number in column \code{repeats}}
    \item{\code{repeats}}{a numeric vector describing the number of VNTR repeats that can be deduced from the actual product length}
    \item{\code{series}}{a factor with levels \code{S1} \code{S2} reflecting the number of PCR rounds; using an ABI sequencer, amplifications can be detected in four channels at the most. If e.g. you have one amplification per channel (you can specify more), you would amplify four products per round. In a design that amplifies 8 loci, you would therefore need two rounds at least.}
    \item{\code{channel}}{a numeric vector describing the channel number in the ABIF file}
  }
}
\details{
Loci in \code{lazy.map} follow nomenclature by Schouls (2006) and are applicable for typing of Neisseria meningitidis. 
Other bacteria obviously have different loci, VNTR-sizes etc, that can be specified in a lazy map. 
You can specify any number of products per channel, provided product ranges do not overlap. 
(\code{lazy.map} conveniently specifies only one locus=product per channel). 
If product ranges overlap, a warning is issued.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Schouls, L. M., van der Ende, A., Damen, M., van de Pol, I. (2006) 'Multiple-locus variable-number tandem repeat analysis of Neisseria meningitidis yields groupings similar to those obtained by multilocus sequence typing', Journal of Clinical Microbiology, 44 (4), pp. 1509-18.
}
\examples{
lazy.map
}
\keyword{datasets, MLVA, VNTR}
