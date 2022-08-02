#' Read IBD probabilities from file
#'
#' Reads IBD probabilities from a plain text, tab-delimited file. Information
#' about the file format can be found in the vignette (
#' \code{vignette("IBDFileFormat", package = "statgenIBD")}).
#'
#' @param file A character string specifying the path of the file.
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
#' @importFrom utils read.table
#' @export
readIBDs <- function(file) {
  if (!is.character(file) || length(file) > 1) {
    stop("File path should be a single character string.")
  }
  if (!file.exists(file)) {
    stop("File not found.")
  }

  data <- read.table(file, sep="\t", header=TRUE)
  data[, -c(1,2)] <- apply(data[, -c(1,2)], 2, function(x) as.numeric(x))
  genotypes <- unique(data$Genotype)
  markers <- unique(data$Marker)
  founders <- colnames(data)[-c(1,2)]
  numGenotypes <- length(genotypes)
  numMarkers <- length(markers)
  numFounders <- length(founders)

  if (length(unique(data[seq_len(numGenotypes), 1])) == 1) {
    # Grouped by marker
    mat <- array(
       data = as.matrix(data[,-c(1,2)]),
       dim = c(numGenotypes, numMarkers, numFounders),
       dimnames = list(genotypes, markers, founders)
    )
  } else if (length(unique(data[seq_len(numMarkers), 2])) == 1) {
    # Grouped by genotype
    stop("Reading by genotype order not implemented.")
  } else {
    # Not grouped by genotype or marker; lets read it the hard way
    mat <- array(rep(NA, numGenotypes * numMarkers * numFounders), dim=c(numGenotypes, numMarkers, numFounders))
    dimnames(mat) <- list(genotypes, markers, founders)
    for (ix in seq_len(nrow(data))) {
      mat[data[ix, 2], data[ix, 1],] <- as.numeric(data[ix, -c(1,2)])
    }
  }
  res <- structure(list(map = NULL,
                        markers = mat,
                        popType = NULL,
                        parents = founders,
                        multiCross = NULL),
                   class = "IBDprob")
  return(res)
}
