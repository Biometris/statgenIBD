#' Read IBD probabilities from file
#'
#' Reads IBD probabilities from a plain text, tab-delimited file. Information
#' about the file format can be found in the vignette (
#' \code{vignette("IBDFileFormat", package = "statgenIBD")}).
#'
#' @param inFile A character string specifying the path of the file.
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
readIBDs <- function(inFile) {
  if (!is.character(inFile) || length(inFile) > 1) {
    stop("inFile path should be a single character string.\n")
  }
  if (!file.exists(inFile)) {
    stop("inFile not found.\n")
  }
  ## Read file.
  inDat <- read.table(inFile, sep = "\t", header = TRUE)
  ## Check that data has required columns.
  if (!all(colnames(inDat)[1:2] == c("Marker", "Genotype"))) {
    stop("First to columns in inFile should be named Marker and Genotype.\n")
  }
  if (ncol(inDat) < 4) {
    stop("At least 2 parent columns should be present in input.\n")
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
    data = as.matrix(inDat[, -c(1, 2)]),
    dim = c(nGeno, nMarkers, nPar),
    dimnames = list(genoNames, markerNames, parents)
  )
  markers <- markers[genoNamesIn, markerNamesIn, ]
  res <- structure(list(map = NULL,
                        markers = markers,
                        popType = NULL,
                        parents = parents,
                        multiCross = NULL),
                   class = "IBDprob")
  return(res)
}
