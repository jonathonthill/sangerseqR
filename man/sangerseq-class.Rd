\docType{class}
\name{sangerseq-class}
\alias{sangerseq}
\alias{sangerseq,abif-method}
\alias{sangerseq,scf-method}
\alias{sangerseq-class}
\title{Sangerseq Class Objects}
\usage{
sangerseq(obj)

\S4method{sangerseq}{abif}(obj)

\S4method{sangerseq}{scf}(obj)
}
\arguments{
  \item{obj}{Can be either an \code{\link{abif}} or
  \code{\link{scf}} object.}
}
\description{
Sangerseq Class Objects contain data necessary for using
sangerseq package functions (e.g.
\code{\link{chromatogram}}, \code{\link{makeBaseCalls}}).
The exact content will depend on the source of the data
(for example, scf files do not have secondary Basecalls).
}
\section{Slots}{

\describe{
\item{\code{primarySeqID}}{Source of the primary basecalls. Functions that modify
 these calls, such as \code{\link{makeBaseCalls}} and
 \code{\link{setAllelePhase}} will also change this value.}

\item{\code{secondarySeqID}}{Source of the secondary basecalls. See above.}

\item{\code{primarySeq}}{The primary Basecalls formatted as a
 \code{\link[Biostrings]{DNAString}} object.}

\item{\code{secondarySeq}}{The secondary Basecalls formatted as a
 \code{\link[Biostrings]{DNAString}} object.}

\item{\code{traceMatrix}}{A numerical matrix containing 4 columns corresponding to the
 normalized signal values for the chromatogram traces. Column order =
 A,C,G,T.}

\item{\code{peakPosMatrix}}{A numerical matrix containing the position of the maximum
 peak values for each base within each Basecall window. If no peak was
 detected for a given base in a given window, then "NA". Column order =
 A,C,G,T.}

\item{\code{peakAmpMatrix}}{A numerical matrix containing the maximum peak amplitudes
 for each base within each Basecall window. If no peak was detected for a
 given base in a given window, then 0. Column order = A,C,G,T.}
}}
\section{Accessor methods}{
  \code{\link{primarySeqID}}, \code{\link{primarySeq}},
  \code{\link{secondarySeqID}}, \code{\link{secondarySeq}},
  \code{\link{traceMatrix}}, \code{\link{peakPosMatrix}},
  \code{\link{peakAmpMatrix}}
}
\examples{
#sample sangerseq object created from abif file
hetsangerseq <- readsangerseq(system.file("extdata",
                                          "heterozygous.ab1",
                                          package = "sangerseqR"))
str(hetsangerseq)
#same for scf files
homosangerseq <- readsangerseq(system.file("extdata",
                                           "homozygous.scf",
                                           package = "sangerseqR"))
str(homosangerseq)
}
\seealso{
\code{\link{abif}}, \code{\link{scf}}
}

