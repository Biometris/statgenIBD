#' @export
getQTL <- function(IBD,
                   QTLS) {
  ## Checks.
  if (!inherits(IBD, "statgenIBD")) {
    stop(deparse(substitute(IBD)), "should be an object of class statgenIBD\n")
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
  markerQTLS <- do.call(cbind, markerQTLS)
  return(markerQTLS)
}
