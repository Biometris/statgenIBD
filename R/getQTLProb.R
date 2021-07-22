#' Extract Probabilities for QTLS
#'
#' Extract IBD probabilities one or more QTLS from an object of class calcIBD.
#'
#' @param IBD An object of class calcIBD.
#' @param QTLS A character vector of QTLS that should be extracted.
#'
#' @return A data.frame with IBD probabilities for the extracted QTLS in the
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
#' ## Get probabilities for a single QTL.
#' probOne <- getQTLProb(IBD = SxMIBD,
#'                       QTLS = "plc")
#' head(probOne)
#'
#' ## Get probabilities for a multiple QTLs.
#' probMult <- getQTLProb(IBD = SxMIBD,
#'                        QTLS = c("plc", "tuba1"))
#' head(probMult)
#'
#' @export
getQTLProb <- function(IBD,
                       QTLS) {
  ## Checks.
  if (!inherits(IBD, "calcIBD")) {
    stop(deparse(substitute(IBD)), " should be an object of class calcIBD\n")
  }
  if (is.null(QTLS) || !is.character(QTLS)) {
    stop("QTLS should be a character vector\n")
  }
  missQTL <- QTLS[!QTLS %in% rownames(IBD$markers)]
  if (length(missQTL) > 0) {
    stop("The following QTLs are not in ", deparse(substitute(IBD)), ": ",
         paste(missQTL, collapse = ", "), "\n")
  }
  markerQTLS <- lapply(X = QTLS, FUN = function(QTL) {
    markerQTL <- IBD$markers[QTL, , ]
    colnames(markerQTL) <- paste0(QTL, "_", colnames(markerQTL))
    return(markerQTL)
  })
  markerQTLS <- as.data.frame(do.call(cbind, markerQTLS))
  genoCross <- attr(x = IBD, which = "genoCross")
  if (!is.null(genoCross)) {
    markerQTLS <- merge(genoCross, markerQTLS, by.x = "geno", by.y = "row.names")
    markerQTLS <- markerQTLS[c("cross", "geno",
                             setdiff(colnames(markerQTLS), c("cross", "geno")))]
  } else {
    markerQTLS[["geno"]] <- rownames(markerQTLS)
    rownames(markerQTLS) <- NULL
    markerQTLS <- markerQTLS[c("geno", setdiff(colnames(markerQTLS), "geno"))]
  }
  return(markerQTLS)
}
