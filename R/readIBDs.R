#' Read IBD probabilities from file
#'
#' Reads IBD probabilities from a plain text, tab-delimited file. Information
#' about the file format can be found in the vignette (
#' \code{vignette("IBDFileFormat", package = "statgenIBD")}).
#'
#' @param infile A character string specifying the path of the input .txt file.
#' Compressed .txt files with extension ".gz" or ".bz2" are supported as well.
#'
#' @return An object of class \code{IBDprob}.
#'
#' @examples
#' ## Read IBD probabilities for Steptoe Morex.
#' SxMIBD <- readIBDs(system.file("extdata/SxM", "SxM_IBDs.txt",
#'                                package = "statgenIBD"))
#'
#' ## Print summary.
#' summary(SxMIBD)
#'
#' @importFrom utils hasName read.table
#' @export
readIBDs <- function(infile) {
  if (missing(infile) || !is.character(infile) || length(infile) > 1 ||
      file.access(infile, mode = 4) == -1 ||
      ## Compressed .csv files can be read by fread and should be
      ## allowed as inputs as well.
      !(tools::file_ext(infile) == "txt" ||
        (tools::file_ext(infile) %in% c("gz", "bz2") &&
         tools::file_ext(tools::file_path_sans_ext(infile)) == "txt"))) {
    stop("infile should be a character string indicating a readable .txt file")
  }
  fileExt <- tools::file_ext(infile)
  if (fileExt %in% c("gz", "bz2")) {
    decompFile <- tempfile(tmpdir = tempdir())
    R.utils::decompressFile(filename = infile, destname = decompFile,
                            ext = NULL,
                            FUN = if (fileExt == "gz") gzfile else bzfile,
                            remove = FALSE)
    infile <- decompFile
    on.exit(unlink(decompFile), add = TRUE)
  }
  ## Read file.
  inDat <- data.table::fread(infile, sep = "\t", header = TRUE)
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
  markers <- array(NA_real_,
    dim = c(nGeno, nMarkers, nPar),
    dimnames = list(genoNames, markerNames, parents)
  )
  for (parent in parents) {
    for (geno in genoNames) {
     markers[geno, , parent] <- inDat[inDat$Genotype == geno, parent]
    }
  }
  markers <- markers[genoNamesIn, markerNamesIn, ]
  res <- structure(list(map = NULL,
                        markers = markers,
                        popType = NULL,
                        parents = parents),
                   class = c("IBDprob", "list"))
  return(res)
}
