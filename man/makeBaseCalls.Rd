\docType{methods}
\name{makeBaseCalls}
\alias{makeBaseCalls}
\alias{makeBaseCalls,sangerseq-method}
\title{Make Basecalls}
\usage{
makeBaseCalls(obj, ratio = 0.33)

\S4method{makeBaseCalls}{sangerseq}(obj, ratio = 0.33)
}
\arguments{
  \item{obj}{\code{\link{sangerseq}} class object}

  \item{ratio}{cutoff ratio for separating signal and
  noise. Ratio is relative to maximum peak in basecall
  window.}
}
\value{
\code{\link{sangerseq}} s4 object with new basecalls in the
primarySeq and secondarySeq fields. Matrix values are also
updated to reflect newly called base positions and
amplitudes of all peaks within window.
}
\description{
Updates a \code{\link{sangerseq}} class object to contain
primary and secondary peak calls.
}
\details{
Scf files do not contain secondary basecalls and ABIF files
sometimes (but not always) contain secondary peak calls
that show the secondary peak even if clearly a negative
peak. This is problematic in sequence reads where
heterozygous sequence data is contained in the chromatogram
data. \code{makeBaseCalls} identifies basecall windows
containing more than one peak using the provided cutoff
ratio and then makes heterozygous calls for those windows.
The primarySeq will always contain the base corresponding
to the maximum peak amplitude within the window. The
secondaryPeak will have the same base if the peak was
classified as a homozygous peak, the base corresponding to
the second tallest peak if two peaks were above the cutoff,
or an ambiguous base if more than two peaks were identified
in the window that are higher than the cutoff ratio.
}
\examples{
hetsangerseq <- readsangerseq(system.file("extdata",
                                          "heterozygous.ab1",
                                          package = "sangerseqR"))
hetcalls <- makeBaseCalls(hetsangerseq, ratio = 0.33)
}
\seealso{
\code{\link{chromatogram}}, \code{\link{setAllelePhase}},
\code{\link{sangerseq}}
}

