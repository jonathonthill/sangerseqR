#method for checking whether DNAStrings are equal. Provided by H. Pages
all.equal.XString <- function(target, current, ...)
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