#'@name sangerseq-class
#'  
#'@aliases sangerseq
#'  
#'@title Sangerseq Class Objects
#'  
#'@description Sangerseq Class Objects contain data necessary for using
#'sangerseq package functions (e.g. \code{\link{chromatogram}},
#'\code{\link{makeBaseCalls}}). The exact content will depend on the source of
#'the data (for example, scf files do not have secondary Basecalls).
#'
#'@param obj Can be either an \code{\link{abif}} or \code{\link{scf}} object.
#'  
#'@slot primarySeqID Source of the primary basecalls. Functions that modify
#'  these calls, such as \code{\link{makeBaseCalls}} and
#'  \code{\link{setAllelePhase}} will also change this value.
#'@slot secondarySeqID Source of the secondary basecalls. See above.
#'@slot primarySeq The primary Basecalls formatted as a
#'  \code{\link[Biostrings]{DNAString}} object.
#'@slot secondarySeq The secondary Basecalls formatted as a
#'  \code{\link[Biostrings]{DNAString}} object.
#'@slot traceMatrix A numerical matrix containing 4 columns corresponding to the
#'  normalized signal values for the chromatogram traces. Column order =
#'  A,C,G,T.
#'@slot peakPosMatrix A numerical matrix containing the position of the maximum
#'  peak values for each base within each Basecall window. If no peak was
#'  detected for a given base in a given window, then "NA". Column order =
#'  A,C,G,T.
#'@slot peakAmpMatrix A numerical matrix containing the maximum peak amplitudes
#'  for each base within each Basecall window. If no peak was detected for a
#'  given base in a given window, then 0. Column order = A,C,G,T.
#'  
#'@section Accessor methods: \code{\link{primarySeqID}}, 
#'  \code{\link{primarySeq}}, \code{\link{secondarySeqID}}, 
#'  \code{\link{secondarySeq}}, \code{\link{traceMatrix}}, 
#'  \code{\link{peakPosMatrix}}, \code{\link{peakAmpMatrix}}
#'  
#'@seealso \code{\link{abif}}, \code{\link{scf}}
#'@examples
#' #sample sangerseq object created from abif file
#' hetsangerseq <- readsangerseq(system.file("extdata", 
#'                                           "heterozygous.ab1", 
#'                                           package = "sangerseqR"))
#' str(hetsangerseq)
#' #same for scf files
#' homosangerseq <- readsangerseq(system.file("extdata", 
#'                                            "homozygous.scf", 
#'                                            package = "sangerseqR"))
#' str(homosangerseq)
#'@export

setClass("sangerseq", 
         representation(
           primarySeqID="character",
           primarySeq="DNAString",
           secondarySeqID="character",
           secondarySeq="DNAString",
           traceMatrix="matrix",
           peakPosMatrix="matrix",
           peakAmpMatrix="matrix"
         )
)

#' @name abif-class
#'   
#' @aliases abif
#'   
#' @title ABIF Class Objects
#'   
#' @description S4 object returned by \code{\link{read.abif}} containing all
#' fields in the ABIF file format (see
#' \url{http://home.appliedbiosystems.com/support/software_community/
#' ABIF_File_Format.pdf}).
#' Data fields vary by machine and basecaller versions. Must be converted to
#' \code{\link{sangerseq}} to be used in other functions from this package.
#' 
#' @slot header Header information from the file.
#' @slot directory Directory information from file containing field names and
#'   information for reading binary data.
#' @slot data List object containing all data fields and values in file.
#'   Included fields vary by machine and basecaller versions.
#'   
#' @seealso \code{\link{read.abif}}, \code{\link{scf}}, \code{\link{sangerseq}}
#' @examples
#' hetab1 <- read.abif(system.file("extdata", 
#'                                 "heterozygous.ab1", 
#'                                 package = "sangerseqR")) 
#' str(hetab1)


setClass("abifHeader", 
         representation(
           abif="character",
           version="integer",
           name="raw",
           number="integer",
           elementtype="integer",
           elementsize="integer",
           numelements="integer",
           dataoffset="integer",
           datahandle="integer"
         )
)

setClass("abifDirectory", 
         representation(
           name="character",
           tagnumber="integer",
           elementtype="integer",
           elementsize="integer",
           numelements="integer",
           datasize="integer",
           dataoffset="integer"
         )
)

#' @export
setClass("abif", 
         representation(
           header="abifHeader",
           directory="abifDirectory",
           data="list"  #must be list because fields vary between file versions
         )
)

#' @name scf-class
#'   
#' @aliases scf
#'   
#' @title Scf Class Objects
#'   
#' @description S4 object returned by \code{\link{read.scf}} containing all
#' fields in the SCF file format (see
#' \url{http://staden.sourceforge.net/manual/formats_unix_2.html}). Must be
#' converted to \code{\link{sangerseq}} to be used in other functions from this
#' package.
#' 
#' @slot header Header information from the file.
#' @slot sample_points Trace data matrix (Order = A, C, G, T).
#' @slot sequence_probs Matrix of the relative probabilities for each base at
#'   each position (Order = A, C, G, T).
#' @slot basecall_positions Vector containing trace matrix indices for each
#'   basecall.
#' @slot basecalls \code{\link[Biostrings]{DNAString}} object containing the
#'   basecalls.
#' @slot comments String containing any comments in the file.
#' @slot private Raw binary data containing any private data in the file.
#'   Generally not used.
#'   
#' @seealso \code{\link{read.scf}}, \code{\link{abif}}, \code{\link{sangerseq}}
#' @examples
#' homoscf <- read.scf(system.file("extdata", 
#'                                 "homozygous.scf", 
#'                                 package = "sangerseqR")) 
#' str(homoscf)


setClass("scfHeader", 
         representation(
           scf="character",
           samples="integer",
           samples_offset="integer",
           bases="integer",
           bases_left_clip="integer",
           bases_right_clip="integer",
           bases_offset="integer",
           comments_size="integer",
           comments_offset="integer",
           version="numeric",
           sample_size="integer",
           code_set="integer",
           private_size="integer",
           private_offset="integer"
         )
)

#' @export
setClass("scf", 
         representation(
           header="scfHeader",
           sample_points="matrix",
           sequence_probs="matrix",
           basecall_positions="integer",
           basecalls="character",
           comments="character",
           private="raw"
         )
)

