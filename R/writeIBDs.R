#' Write IBD probabilities to file.
#'
#' Writes IBD probabilities to a plain text, tab-delimited file. Information
#' about the file format can be found in the vignette (
#' \code{vignette("IBDFileFormat", package = "statgenIBD")}).
#'
#' @param IBDprob An object of class \code{IBDprob} containing the IBD
#' probabilities.
#' @param outFile A character string specifying the path of the output file.
#' @param decimals An integer value specifying the number of decimals to include
#' in writing the output file. When negative, no decimal cut off is applied.
#' @param minProb A numerical value between zero and 1 / number of parents,
#' specifying the minimum probability cutoff value. Probabilities below this
#' cutoff are set to zero and other probabilities are rescaled to make sure that
#' the probabilities sum up to one.
#'
#' @return No output. The output file is created as a result of calling this
#' function.
#'
#' @examples
#' ## Compute IBD probabilities for Steptoe Morex.
#' SxMIBD <- calcIBD(popType = "DH",
#'                  markerFile = system.file("extdata/SxM", "SxM_geno.txt",
#'                                          package = "statgenIBD"),
#'                  mapFile = system.file("extdata/SxM", "SxM_map.txt",
#'                                       package = "statgenIBD"))
#'
#' ## Write IBDs to file.
#' writeIBDs(IBDprob = SxMIBD, outFile = "SxM_IBDs.ibd")
#'
#' ## Write IBDs to file, set values <0.05 to zero and only print 3 decimals.
#' writeIBDs(IBDprob = SxMIBD, outFile = "SxM_IBDs2.ibd",
#'          decimals = 3, minProb = 0.05)
#'
#' @export
writeIBDs <- function(IBDprob,
                      outFile,
                      decimals = -1,
                      minProb = 0) {
  if (!inherits(IBDprob, "IBDprob")) {
    stop("IBDprob should be an object of class IBDprob.\n")
  }
  chkFile(outFile, fileType = "ibd")
  if (!is.numeric(decimals) || length(decimals) > 1) {
    stop("decimals should be a single numerical value.\n")
  }
  if (!is.numeric(minProb) || length(minProb) > 1 || minProb < 0 ||
      minProb >= 1 / length(IBDprob$parents)) {
    stop("minProb should be a numerical value between 0 and ",
         "1 / number of parents.\n")
  }
  markers <- IBDprob$markers
  markerNames <- rownames(markers)
  genoNames <- colnames(markers)
  parents <- IBDprob$parents
  fmt <- if (decimals > 0) paste0("%#.", decimals, "f") else "%f"
  if (minProb > 0) {
    ## Set values < minProb to zero and rescale.
    markers[markers < minProb] <- 0
    markers <- simplify2array(apply(X = markers, MARGIN = 2, FUN = function(x) {
      x / rowSums(x)
    }, simplify = FALSE))
    markers <- aperm(markers, c(1, 3, 2))
  }
  ## Create base data.frame for storing data.
  markersLongBase <- expand.grid(Genotype = genoNames, Marker = markerNames)
  markersLongBase <- markersLongBase[c("Marker", "Genotype")]
  for (parent in parents) {
    ## Format output.
    ## Trailing zeros are removed and number of decimals is adjusted.
    markersLongBase[[parent]] <-
      sub("\\.$", "",
          sub("0+$", "",
              sprintf(t(markers[, , parent]), fmt = fmt)
              )
    )
  }
  write.table(markersLongBase, file = outFile,
              quote = FALSE, sep = "\t", na = "-", row.names = FALSE,
              col.names = TRUE)
}
