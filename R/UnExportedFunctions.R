getpeaks <- function(trace) {
  r <- rle(trace)
  indexes <- which(rep(diff(sign(diff(c(-Inf, r$values, -Inf)))) == -2, 
                       times = r$lengths))
  cbind(indexes, trace[indexes])
}

peakvalues <- function(x, pstart, pstop) {
  region <- x[x[,1] > pstart & x[,1] < pstop, ,drop=FALSE]
  if (length(region[,1]) == 0) return(c(0, NA))
  else return(c(max(region[,2], na.rm=TRUE), region[which.max(region[,2]),1]))
}

#Create data frame of basecalls
basecalldf <- function(obj) {
  primary <- strsplit(toupper(primarySeq(obj, string=TRUE)), "")[[1]]
  secondary <- strsplit(toupper(secondarySeq(obj, string=TRUE)), "")[[1]]
  basecalls <- data.frame(primary=primary, 
                          secondary=secondary, 
                          stringsAsFactors=FALSE)
  basecalls$primary <- unname(IUPAC_CODE_MAP[basecalls$primary])
  basecalls$secondary <- unname(IUPAC_CODE_MAP[basecalls$secondary])
  basecalls$consensus <- basecalls$primary
  basecalls$consensus[basecalls$primary != basecalls$secondary 
                      | nchar(basecalls$consensus) > 1] <- "N"
  basecalls$possibilities <- basecalls$consensus
  basecalls$possibilities[basecalls$possibilities == "N"] <- 
    paste0(basecalls$primary[basecalls$possibilities == "N"], 
           basecalls$secondary[basecalls$possibilities == "N"])
  return(basecalls)
}


#functions for converting binary data into numbers/text
RTC <- function(x, ...) {
  string <- suppressWarnings(rawToChar(x, ...))
  if(length(string) > 1) string <- paste(string, collapse="")
  #found that some ab1 files have unprinted characters at the end of the string
  #this is designed to remove them
  string <- gsub("[^ -~]", "", string)
  return(string)
}

SInt32 <- function(f, n=length(f)/4) readBin(f, what = "integer", 
                                             signed = TRUE, endian = "big", 
                                             size = 4, n=n)
SInt16 <- function(f, n=length(f)/2) readBin(f, what = "integer", 
                                             signed = TRUE, endian = "big", 
                                             size = 2, n=n)
SInt8 <- function(f, n=length(f)) readBin(f, what = "integer", signed = TRUE, 
                                  endian = "big", size = 1, n=n)
UInt32 <- function(f, n=length(f)/4) readBin(f, what = "integer", 
                                             signed = FALSE, endian = "big", 
                                             size = 4, n=n)
UInt16 <- function(f, n=length(f)/2) readBin(f, what = "integer",
                                             signed = FALSE, endian = "big", 
                                             size = 2, n=n)
UInt8 <- function(f, n=length(f)) readBin(f, what = "integer", signed = FALSE, 
                                  endian = "big", size = 1, n=n)
f32 <- function(f, n=length(f)/4) readBin(f, what = "numeric", size = 4, n=n)
f64 <- function(f, n=length(f)/8) readBin(f, what = "numeric", size = 8, n=n)


#read binary data to string
readBinaryData <- function(filename) {
  fc <- file(filename, open = "rb")
  rawdata <- readBin(fc, what = "raw", n = 1.2*file.info(filename)$size)
  close(fc)
  return(rawdata)
}

#convert offsets to actual coordinates
convertPoints <- function(rawdata) {
  for (j in 1:2)  { #have to do twice
    previous <- 0
    for (i in 1:length(rawdata)) {
      rawdata[i] <- rawdata[i]  + previous
      previous <- rawdata[i]
    }
  }
  return(rawdata)
}

