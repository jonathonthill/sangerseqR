figheight <- function(obj, trim5=0, trim3=0, width=100, pixelsperrow=200, showtrim=FALSE) {
  traces <- obj@traceMatrix
  basecalls1 <- unlist(strsplit(toString(obj@primarySeq), ""))
  basecalls2 <- unlist(strsplit(toString(obj@secondarySeq), ""))
  aveposition <- rowMeans(obj@peakPosMatrix, na.rm=TRUE)
  basecalls1 <- basecalls1[1:length(aveposition)] #####
  basecalls2 <- basecalls2[1:length(aveposition)] ######
  if(showtrim == FALSE) {
    if(trim5+trim3 > length(basecalls1)) basecalls1 <- ""
    else basecalls1 <- basecalls1[(1 + trim5):(length(basecalls1) - trim3)]
    if(trim5+trim3 > length(basecalls2)) basecalls2 <- ""
    else basecalls2 <- basecalls2[(1 + trim5):(length(basecalls2) - trim3)]
    aveposition <- aveposition[(1 + trim5):(length(aveposition) - trim3)] 
  }
  indexes <- 1:length(basecalls1)
  trimmed <- indexes <= trim5 | indexes > (length(basecalls1) - trim3) # all false if not trimmed
  if (!is.null(trim3)) {
    traces <- traces[1:(min(max(aveposition, na.rm=TRUE) + 10, nrow(traces))),]
  }
  if (!is.null(trim5)) {
    offset <- max(c(1, aveposition[1] - 10))
    traces <- traces[offset:nrow(traces),]
    aveposition <- aveposition - (offset-1)
  }
 
  valuesperbase <- nrow(traces)/length(basecalls1)
  tracewidth <- width*valuesperbase
  breaks <- seq(1,nrow(traces), by=tracewidth) 
  
  numplots <- length(breaks)
  return(numplots*pixelsperrow)
}

wrap_fixed = function(string, width=80) {
  pattern = paste("(.{1,", width, "})", sep="")
  res = gsub(pattern, "\\1\n", string)
  return(res)
}

cleanstring <- function(string) {
  string <- gsub("^>.*?\n", "", string)
  string <- toupper(string)
  string <- gsub("[^ACGTRYSWKMBDHVN]", "", string, perl=TRUE)
  return(string)
}

alignchromatogram <- function(data, block.width=50, trim=FALSE, refseq, trim5, trim3) {
  if (is.null(data)) return(NULL)
  d <- setAllelePhase(data, refseq, trim5, trim3)
  altseq <- toString(d@secondarySeq)
  if (trim == TRUE) {
    altseq <- toString(d@secondarySeq[(trim5 + 1):(nchar(altseq) - trim3)])
  }
  names(altseq) <- "Alt Allele"
  names(refseq) <- d@primarySeqID
  if (trim == TRUE) {
    pa <- pairwiseAlignment(altseq, refseq, type="local", gapExtension=-2)
  } else {
    pa <- pairwiseAlignment(altseq, refseq, type="global", gapExtension=-2)
  }
  alignment <- paste(capture.output(writePairwiseAlignments(pa, block.width=block.width)), collapse="\n")
  results <- list(altseq=altseq, refseq=gsub("-", "", refseq), alignment=alignment)
  return(results)
}