#' Write to Flapjack format
#'
#' Export the results of an IBD calculation to flapjack format so it can be
#' visualalized there.
#'
#' @param x An object of class \code{\link{IBD}}.
#' @param outFileMap A character string, the full path to the output map file.
#' @param outFileGeno A character string, the full path to the output genotype
#' file.
#'
#' @export
writeFlapjack <- function(x,
                          outFileMap = "ibd_map.txt",
                          outFileGeno = "ibd_geno.txt") {
  map <- x$map
  markers <- x$markers
  parents <- x$parents
  nPar <- length(parents)
  nGeno <- dim(markers)[[2]]
  nMarkers <- dim(markers)[[1]]
  ## Convert to long format.
  markersLong <- markers3DtoLong(x)
  markersWide <- reshape(data = markersLong, idvar = c("genotype", "snp"),
                         timevar = "parent", direction = "wide")
  for (i in 1:nrow(markersWide)) {
     parentsI <- parents[markersWide[i , 3:(2 + nPar)] > (0.85 / nPar)]
     markersWide[i, "res"] <- paste0(parentsI, collapse = "/")
  }
  res <- matrix(data = markersWide[["res"]], nrow = nGeno, byrow = TRUE,
                dimnames = dimnames(markers)[2:1])
  ## Write map file.
  cat(file = outFileMap, "# fjFile = MAP\n")
  write.table(map, file = outFileMap,
              quote = FALSE, sep = "\t", na = "-", row.names = TRUE,
              col.names = FALSE, append = TRUE)
  ## Write marker file.
  cat(file = outFileGeno, "# fjFile = GENOTYPE\n\t")
  write.table(res, file = outFileGeno,
              quote = FALSE, sep = "\t", na = "-", row.names = TRUE,
              col.names = TRUE, append = TRUE)
}