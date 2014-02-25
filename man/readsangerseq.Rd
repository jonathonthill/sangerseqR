\name{readsangerseq}
\alias{readsangerseq}
\title{Read Scf or ABIF Files}
\usage{
readsangerseq(filename)
}
\arguments{
  \item{filename}{Location of the file.}
}
\value{
\code{\link{sangerseq}} s4 object
}
\description{
This is a convienience function for reading Scf or ABIF
files into a sangerseq object, which can be used by the
other sangerseq package functions. It is equivalent to
calling \code{\link{read.scf}} or \code{\link{read.abif}}
as appropriate and then calling \code{\link{sangerseq}}.
}
\examples{
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
\code{\link{read.abif}}, \code{\link{read.scf}},
\code{\link{abif}}, \code{\link{scf}},
\code{\link{sangerseq}}
}

