#' Read IBD probabilities from file
#'
#' Reads IBD probabilities from a plain text, tab-delimited file. Information
#' about the file format can be found in the vignette (
#' \code{vignette("IBDFileFormat", package = "statgenIBD")}).
#'
#' @param infile A character string specifying the path of the file.
#'
#' @return An object of class \code{IBDprob} containing data a data.frame with
#' the IBD probabilities.
#'
#' @examples
#' ## Read IBD probabilities for Steptoe Morex.
#' SxMIBD <- readIBDs(system.file("extdata/SxM", "SxM_IBDs.ibd",
#'                                package = "statgenIBD"))
#'
#' ## Print summary.
#' summary(SxMIBD)
#'
#' @importFrom utils hasName read.table
#' @export
readIBDs <- function(infile) {
  if (!is.character(infile) || length(infile) > 1) {
    stop("infile path should be a single character string.\n")
  }
  if (!file.exists(infile)) {
    stop("infile not found.\n")
  }
  ## Read file.
  inDat <- read.table(infile, sep = "\t", header = TRUE)
  ## Check that data has required columns.
  if (!all(colnames(inDat)[1:2] == c("Marker", "Genotype"))) {
    stop("First two columns in infile should be named Marker and Genotype.\n")
  }
  if (ncol(inDat) < 4) {
    stop("At least two parent columns should be present in input.\n")
  }
  inDat[, -c(1,2)] <- apply(inDat[, -c(1,2)], 2, function(x) as.numeric(x))
  genoNamesIn <- unique(inDat[["Genotype"]])
  markerNamesIn <- unique(inDat[["Marker"]])
  ## Sort input data to get everything in expected order.
  inDat <- inDat[order(inDat[["Marker"]], inDat[["Genotype"]]), ]
  genoNames <- unique(inDat[["Genotype"]])
  markerNames <- unique(inDat[["Marker"]])
  parents <- colnames(inDat)[-c(1, 2)]
  nGeno <- length(genoNames)
  nMarkers <- length(markerNames)
  nPar <- length(parents)
  markers <- array(
    dim = c(nMarkers, nGeno, nPar),
    dimnames = list(markerNames, genoNames, parents)
  )
  for (parent in parents) {
    for (geno in genoNames) {
     markers[, geno, parent] <- inDat[inDat$Genotype == geno, parent]
    }
  }
  markers <- markers[markerNamesIn, genoNamesIn, ]
  res <- structure(list(map = NULL,
                        markers = markers,
                        popType = NULL,
                        parents = parents),
                   class = c("IBDprob", "list"))
  return(res)
}
