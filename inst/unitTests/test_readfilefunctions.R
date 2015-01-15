test_readfilefunctions <- function() {
  #method for checking whether DNAStrings are equal. Provided by H. Pages
  all.equal.XString <<- function(target, current, ...) #.GlobalEnv
  {
    ok1 <- target == current
    if (!ok1)
      ok1 <- "sequence mismatch"
    ok2 <- identical(class(target), class(current))
    if (!ok2)
      ok2 <- "class mismatch"
    ok3 <- all.equal(target@metadata, current@metadata)
    ok4 <- all.equal(target@elementMetadata, current@elementMetadata)
    ans <- character(0)
    if (!isTRUE(ok1))
      ans <- c(ans, ok1)
    if (!isTRUE(ok2))
      ans <- c(ans, ok2)
    if (!isTRUE(ok3))
      ans <- c(ans, ok3)
    if (!isTRUE(ok4))
      ans <- c(ans, ok4)
    if (length(ans) == 0L)
      return(TRUE)
    ans
  }
  
  #data files contains sangerseq obj from each file type
  expectedabif <- readRDS(
    system.file("testData", "heterozygousabif.rds", package="sangerseqR"))
  expectedscf <- readRDS(
    system.file("testData", "heterozygousscf.rds", package="sangerseqR"))
  #contains unprintable characters
  expectednoprintabif <- readRDS(
    system.file("testData", "noprintabif.rds", package="sangerseqR"))
  
  abif <- readsangerseq(
    system.file("extdata", "heterozygous.ab1", package="sangerseqR"))
  scf <- readsangerseq(
    system.file("extdata", "heterozygous.scf", package="sangerseqR"))
  noprintabif <- readsangerseq(
    system.file("extdata", "noprintabif.ab1", package="sangerseqR"))
  
  #Check if objects are the same
  checkEquals(abif, expectedabif)
  checkEquals(scf, expectedscf)
  checkEquals(noprintabif, expectednoprintabif)
}


