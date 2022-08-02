#' Write IBD probabilities to file.
#'
#' Writes IBD probabilities to a plain text, tab-delimited file.  Information
#' about the file format can be found in the vignette (
#' \code{vignette("IBDFileFormat", package = "statgenIBD")}).
#'
#' @param ibdProb An object of class \code{IBDprob} containing the IBD
#' probabilities.
#' @param file A character string specifying the path of the output file.
#' @param decimals An integer value specifying the number of decimals to include
#' in writing the output file. When negative, no decimal cut off is applied.
#' @param minProb A number between zero and one, specifying the minimum
#' probability cut off value. Probabilities below this cutoff are set to zero
#' and other probabilities are rescaled to make sure that the probabilities sum
#' up to one.
#'
#' @return No output. The output file is created as a result of calling this
#' function.
#'
#' @examples
#' ## Compute IBD probabilities for Steptoe Morex.
#' SxMIBD <- calcIBD(popType = "DH",
#'              markerFile = system.file("extdata/SxM", "SxM_geno.txt",
#'                                    package = "statgenIBD"),
#'              mapFile = system.file("extdata/SxM", "SxM_map.txt",
#'                                    package = "statgenIBD"))
#'
#' ## Write IBDs to file.
#' writeIBDs(SxMIBD, "SxM_IBDs.ibd")
#'
#' ## Write IBDs to file, set values <0.05 to zero and only print 3 decimals.
#' writeIBDs(SxMIBD, "SxM_IBDs.ibd", decimals=3, minProb=0.05)
#'
#' @export
writeIBDs <- function(ibdProb,
                      file,
                      decimals=-1,
                      minProb=0) {
  #if (!inherits(ibdProb, "IBDprob")) {
  #  stop("ibdProb should be an object of class IBDprob.\n")
  #}
  chkFile(file, fileType = "ibd")
  ibds <- ibdProb$markers
  markerNames <- dimnames(ibds)[[2]]
  genoNames <- dimnames(ibds)[[1]]
  numFounders <- dim(ibds)[[3]]
  numGenotypes <- dim(ibds)[[1]]
  fmt <- ifelse(decimals > 0, paste("%#.", decimals, "f", sep = ""), "%f")
  fileConn <- file(file, open = "w")
  writeLines(
    paste(c("Marker", "Genotype", dimnames(ibds)[[3]]), sep = "\t", collapse = "\t"),
    con = fileConn,
    sep = "\n"
  )
  for (mrkIndex in 1:dim(ibds)[[2]]) {
    curMarker <- ibds[,mrkIndex,]
    if (minProb > 0) {
      curMarker[curMarker < minProb] <- 0
      mrkProbs <- rowSums(curMarker, na.rm = FALSE)
      curMarker <- curMarker / mrkProbs
    }
    for (genoIndex in 1:dim(ibds)[[1]]) {
      writeLines(
        paste(
          c(
            markerNames[mrkIndex], ## Marker
            genoNames[genoIndex],  ## Genotype
            sub("\\.$", "",
              sub("0+$", "",
                  sprintf(curMarker[genoIndex,], fmt=fmt)
              )
            )
          ),
          sep = "\t",
          collapse = "\t"
        ),
        con = fileConn,
        sep = "\n"
      )
    }
  }
  close(fileConn)
}
