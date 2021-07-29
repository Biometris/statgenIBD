#' Extract Probabilities for markers
#'
#' Extract IBD probabilities one or more markers from an object of class
#' calcIBD.
#'
#' @param IBD An object of class calcIBD.
#' @param markers A character vector of markers that should be extracted.
#'
#' @return A data.frame with IBD probabilities for the extracted markers in the
#' column and genotypes in the rows.
#'
#'  @examples
#' ## Compute IBD probabilities for Steptoe Morex.
#' SxMIBD <- calcIBD(poptype = "DH",
#'                   locfile = system.file("extdata/SxM", "SxM_geno.txt",
#'                                         package = "statgenIBD"),
#'                   mapfile = system.file("extdata/SxM", "SxM_map.txt",
#'                                         package = "statgenIBD"))
#'
#' ## Get probabilities for a single marker.
#' probOne <- getProbs(IBD = SxMIBD,
#'                     markers = "plc")
#' head(probOne)
#'
#' ## Get probabilities for a multiple markers.
#' probMult <- getProbs(IBD = SxMIBD,
#'                      markers = c("plc", "tuba1"))
#' head(probMult)
#'
#' @export
getProbs <- function(IBD,
                     markers) {
  ## Checks.
  if (!inherits(IBD, "calcIBD")) {
    stop(deparse(substitute(IBD)), " should be an object of class calcIBD\n")
  }
  if (is.null(markers) || !is.character(markers)) {
    stop("markers should be a character vector\n")
  }
  missQTL <- markers[!markers %in% rownames(IBD$markers)]
  if (length(missQTL) > 0) {
    stop("The following markers are not in ", deparse(substitute(IBD)), ": ",
         paste(missQTL, collapse = ", "), "\n")
  }
  probs <- lapply(X = markers, FUN = function(marker) {
    prob <- IBD$markers[marker, , ]
    colnames(prob) <- paste0(marker, "_", colnames(prob))
    return(prob)
  })
  probs <- as.data.frame(do.call(cbind, probs))
  genoCross <- attr(x = IBD, which = "genoCross")
  if (!is.null(genoCross)) {
    probs <- merge(genoCross, probs, by.x = "geno", by.y = "row.names")
    probs <- probs[c("cross", "geno",
                     setdiff(colnames(probs), c("cross", "geno")))]
  } else {
    probs[["geno"]] <- rownames(probs)
    rownames(probs) <- NULL
    probs <- probs[c("geno", setdiff(colnames(probs), "geno"))]
  }
  return(probs)
}
