#' Read Scf Files
#' 
#' Reads Scf sanger sequencing data files. Scf files are an open source binary 
#' sanger sequencing chromatogram data file 
#' (see \url{http://staden.sourceforge.net/manual/formats_unix_2.html}). 
#' The file is read and parsed into an \code{\link{scf}} class object.
#' 
#' @param filename Location of the file.
#' 
#' @return 
#' \code{\link{scf}} s4 object
#' 
#' @seealso
#' \code{\link{read.abif}}, \code{\link{readsangerseq}},  \code{\link{scf}}
#' 
#' @examples
#' homoscf <- read.scf(system.file("extdata", 
#'                                 "homozygous.scf", 
#'                                 package = "sangerseqR")) 
#' str(homoscf)
#' 
#' @export

read.scf <- function (filename) 
{
  #Load file
  rawdata <- readBinaryData(filename)
  
  #Initialize new scf object
  res <- new("scf")
  
  #Get Header
  res@header@scf <- RTC(rawdata[1:4])
  
  if (res@header@scf != ".scf") 
    stop("file not in SCF format")
  
  res@header@samples <- SInt32(rawdata[5:8])
  res@header@samples_offset <- SInt32(rawdata[9:12])
  res@header@bases <- SInt32(rawdata[13:16])
  res@header@bases_left_clip <- SInt32(rawdata[17:20]) #obsolete
  res@header@bases_right_clip <- SInt32(rawdata[21:24]) #obsolete
  res@header@bases_offset <- SInt32(rawdata[25:28])
  res@header@comments_size <- SInt32(rawdata[29:32])
  res@header@comments_offset <- SInt32(rawdata[33:36])
  res@header@version <- as.numeric(RTC(rawdata[37:40]))
  res@header@sample_size <- SInt32(rawdata[41:44])
  res@header@code_set <- SInt32(rawdata[45:48])
  res@header@private_size <- SInt32(rawdata[49:52])
  res@header@private_offset <- SInt32(rawdata[53:56])
  #res@header@unused <- SInt32(rawdata[57:128], n = 18)
  #res@header@unused[1:length(res@header@unused)]
  
  #Set data type boundaries
  samples_start <- res@header@samples_offset + 1
  bases_start <- res@header@bases_offset + 1
  comments_start <- res@header@comments_offset + 1 
  private_start <- res@header@private_offset + 1 
  samples_end <- res@header@samples_offset + res@header@samples*8
  bases_end <- res@header@bases_offset + res@header@bases*12
  comments_end <- res@header@comments_offset + res@header@comments_size
  private_end <- res@header@private_offset + res@header@private_size
   
  #Get Trace Data
  if(res@header@version > 2.9) {
    rawSamples <- matrix(SInt16(rawdata[samples_start:samples_end]), ncol=4, 
                         byrow=FALSE)
    res@sample_points <- apply(rawSamples, 2, convertPoints)
  } else {
    res@sample_points <- matrix(SInt16(rawdata[samples_start:samples_end]), ncol=4, 
                                byrow=TRUE)
  }
  
  #Get Basecall data
  rawBases <- rawdata[bases_start:bases_end]
  numBases <- res@header@bases
  
  if(res@header@version > 2.9) {
    pos_end <- numBases*4
    A_end <- pos_end + numBases
    C_end <- A_end + numBases
    G_end <- C_end + numBases
    T_end <- G_end + numBases
    Call_end <- T_end + numBases
    res@basecall_positions <- SInt32(rawBases[1:pos_end])
    res@basecalls <- paste(sapply(rawBases[(T_end+1):Call_end], RTC), 
                           collapse="")
    
    #get probabilities
    res@sequence_probs <- cbind(
      SInt8(rawBases[(pos_end+1):A_end]),
      SInt8(rawBases[(A_end+1):C_end]),
      SInt8(rawBases[(C_end+1):G_end]),
      SInt8(rawBases[(G_end+1):T_end])
    )
  } else {
    rawBaseMat <- matrix(rawBases, ncol=12, byrow=TRUE)
    #last 3 columns are empty
    res@basecall_positions <- apply(rawBaseMat[,1:4], 1, SInt32)
    res@sequence_probs <- cbind(
      SInt8(rawBaseMat[,5]),
      SInt8(rawBaseMat[,6]),
      SInt8(rawBaseMat[,7]),
      SInt8(rawBaseMat[,8])
      )
    res@basecalls <- paste(RTC(rawBaseMat[,9], multiple=TRUE), collapse="")
    res@basecalls <- gsub("-", "N", res@basecalls)
  }
  
  #Get Comments
  res@comments <- RTC(rawdata[comments_start:comments_end])
  #Get Private, leave as raw because don't know what will be put here
  res@private <- rawdata[private_start:
                           private_end]
  return(res)
}

#' Read ABIF Files
#' 
#' Reads ABIF sanger sequencing data files. ABIF files are a proprietary binary 
#' sanger sequencing chromatogram data file created by Applied Biosystems (see
#' \url{http://home.appliedbiosystems.com/support/software_community/
#' ABIF_File_Format.pdf}). The file is read and parsed into an 
#' \code{\link{abif}} class object. This method is based on the read.abif 
#' function in the seqinr package available on CRAN.
#' 
#' @param filename Location of the file.
#'   
#' @return \code{\link{abif}} s4 object
#' 
#' @references Charif, D. and Lobry, J.R. (2007) SeqinR 1.0-2: a contributed
#' package to teh R project for statistical computing devoted to biological
#' sequences retrieval and analysis. Structural approches to sequenc eevolution:
#' Molecules, networks, populations. pp. 207-232.
#' 
#' @seealso \code{\link{read.scf}}, \code{\link{readsangerseq}}, 
#' \code{\link{abif}}
#' 
#' @examples
#' hetab1 <- read.abif(system.file("extdata", 
#'                                 "heterozygous.ab1", 
#'                                 package = "sangerseqR")) 
#' str(hetab1)
#' 
#' @export

read.abif <- function (filename) {
  
  #initialize object
  res <- new("abif")
  
  #get data
  rawdata <- readBinaryData(filename)
  
  res@header@abif <- RTC(rawdata[1:4])
  
  if (res@header@abif != "ABIF") 
    stop("file not in ABIF format")
  
  res@header@version <- SInt16(rawdata[5:6])
  res@header@name <- rawdata[7:10]
  res@header@number <- SInt32(rawdata[11:14])
  res@header@elementtype <- SInt16(rawdata[15:16])
  res@header@elementsize <- SInt16(rawdata[17:18])
  res@header@numelements <- SInt32(rawdata[19:22])
  res@header@dataoffset <- SInt32(rawdata[27:30])
  dataoffset <- res@header@dataoffset + 1
  res@header@datahandle <- SInt32(rawdata[31:34])
  #res@header@unused <- SInt16(rawdata[35:128], n = 47)
  #res@header@unused[1:length(res@header@unused)] <- 0
  
  #get directory
  for (i in seq_len(res@header@numelements)) {
    deb <- (i - 1) * res@header@elementsize + dataoffset
    direntry <- rawdata[deb:(deb + res@header@elementsize)]
    res@directory@name <- c(res@directory@name, RTC(direntry[1:4]))
    res@directory@tagnumber <- c(res@directory@tagnumber, 
                                 SInt32(direntry[5:8]))
    res@directory@elementtype <- c(res@directory@elementtype, 
                                   SInt16(direntry[9:10]))
    res@directory@elementsize <- c(res@directory@elementsize, 
                                   SInt16(direntry[11:12]))
    res@directory@numelements <- c(res@directory@numelements, 
                                   SInt32(direntry[13:16]))
    res@directory@datasize <- c(res@directory@datasize, 
                                SInt32(direntry[17:20]))
    res@directory@dataoffset <- c(res@directory@dataoffset, 
                                  SInt32(direntry[21:24]))
  }
  #fix for error in some .ab1 files that have the wrong data type for the 
  #PCON fields. Usually is 2 ("character") but should be 1 ("Uint8")
  res@directory@elementtype[res@directory@name == "PCON"] <- as.integer(1)

  #get data list
  res@data <- vector("list", length(res@directory@name))
  names(res@data) <- paste(res@directory@name, 
                           res@directory@tagnumber, 
                           sep = ".")
  for (i in seq_len(res@header@numelements)) {
    deb <- (i - 1) * res@header@elementsize + dataoffset
    if (res@directory@datasize[i] > 4) {
      debinraw <- res@directory@dataoffset[i] + 1
    }
    else {
      debinraw <- deb + 20
    }
    elementtype <- res@directory@elementtype[i]
    numelements <- res@directory@numelements[i]
    elementsize <- res@directory@elementsize[i]
    data <- rawdata[debinraw:(debinraw + numelements * elementsize)]
    if (elementtype == 1) 
      res@data[[i]] <- UInt8(data, n = numelements)
    if (elementtype == 2) {
      res@data[[i]] <- tryCatch(RTC(data), finally = paste(rawToChar(data, 
                  multiple = TRUE), collapse = ""), 
                  error = function(er) {
                  cat(paste("an error was detected with the following message:",
                            er, 
                            " but this error was fixed\n", 
                            sep = " "))
                            })
    }
    if (elementtype == 3) 
      res@data[[i]] <- UInt16(data, n = numelements)
    if (elementtype == 4) 
      res@data[[i]] <- SInt16(data, n = numelements)
    if (elementtype == 5) 
      res@data[[i]] <- SInt32(data, n = numelements)
    if (elementtype == 7) 
      res@data[[i]] <- f32(data, n = numelements)
    if (elementtype == 8) 
      res@data[[i]] <- f64(data, n = numelements)
    if (elementtype == 10) 
      res@data[[i]] <- list(year = SInt16(data, n = 1), 
                            month = UInt8(data[-(1:2)], n = 1), 
                            day = UInt8(data[-(1:3)], 
                                        n = 1)
                            )
    if (elementtype == 11) 
      res@data[[i]] <- list(hour = UInt8(data, n = 1), 
                            minute = UInt8(data[-1], n = 1), 
                            second = UInt8(data[-(1:2)], n = 1), 
                            hsecond = UInt8(data[-(1:3)], n = 1)
                            )
    if (elementtype == 18) {
      n <- SInt8(rawdata[debinraw])
      pString <- RTC(rawdata[(debinraw + 1):(debinraw + 
                                               n)])
      res@data[[i]] <- pString
    }
    if (elementtype == 19) 
      res@data[[i]] <- RTC(data[1:(length(data) - 1)])
    if (elementtype >= 1024) 
      res@data[[i]] <- data
    if (elementtype %in% c(12, 13)) 
      warning("unimplemented legacy type found in file")
    if (elementtype %in% c(6, 9, 14, 15, 16, 17, 20, 128, 
                           256, 384)) 
      warning("unsupported legacy type found in file")
  }
  return(res)
}


#' Read Scf or ABIF Files
#' 
#' This is a convienience function for reading Scf or ABIF files into a
#' sangerseq object, which can be used by the other sangerseq package functions.
#' It is equivalent to calling \code{\link{read.scf}} or \code{\link{read.abif}}
#' as appropriate and then calling \code{\link{sangerseq}}.
#' 
#' @param filename Location of the file.
#'   
#' @return \code{\link{sangerseq}} s4 object
#' 
#' @seealso \code{\link{read.abif}}, \code{\link{read.scf}}, \code{\link{abif}},
#' \code{\link{scf}}, \code{\link{sangerseq}}
#' 
#' @examples
#' hetsangerseq <- readsangerseq(system.file("extdata", 
#'                                           "heterozygous.ab1", 
#'                                           package = "sangerseqR"))
#' str(hetsangerseq)
#' #same for scf files
#' homosangerseq <- readsangerseq(system.file("extdata", 
#'                                            "homozygous.scf", 
#'                                            package = "sangerseqR"))
#' str(homosangerseq)
#' 
#' @export

readsangerseq <- function (filename)  {
  #check file type
  fc <- file(filename, open = "rb")
  rawdata <- readBin(fc, what = "raw", n = 1.2*file.info(filename)$size)
  close(fc)
  
  #Get magic number
  filetype <- suppressWarnings(rawToChar(rawdata[1:4]))

  #execute appropriate function based on magic number
  if (filetype == ".scf") {
    seq <- read.scf(filename)
  }
  else if (filetype == "ABIF") {
    seq <- read.abif(filename)
  }
  else stop("Invalid File.")
  
  return(sangerseq(seq))
}