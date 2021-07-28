#' Summary function for objects of class calcIBD
#'
#' Prints a short summary for objects of class \code{\link{calcIBD}}. The
#' summary consists of the population type, number of evaluation points,
#' number of individuals and names of the parents in the object.
#'
#' @param object An object of class \code{\link{calcIBD}}.
#' @param ... Not used.
#'
#' @return No return value, a summary is printed.
#'
#' @examples
#' ## Compute IBD probabilities for Steptoe Morex.
#' SxMIBD <- calcIBD(poptype = "DH",
#'                   locfile = system.file("extdata/SxM", "SxM_geno.txt",
#'                                         package = "statgenIBD"),
#'                   mapfile = system.file("extdata/SxM", "SxM_map.txt",
#'                                         package = "statgenIBD"))
#'
#' ## Print summary
#' summary(SxMIBD)
#'
#' @export
summary.calcIBD <- function(object,
                            ...) {
  cat("population type: ", object$poptype, "\n")
  cat("Number of evaluation points: ", nrow(object$markers), "\n")
  cat("Number of individuals: ", ncol(object$markers),"\n")
  cat("Parents: ", substring(dimnames(object$markers)[[3]], first = 2), "\n")
}

#' Concatenate function for objects of class calcIBD
#'
#' Concatenates objects of class \code{\link{calcIBD}}. All objects that are
#' concatenated  should have the same population type and the same map. The
#' function is mainly meant for combining information for multiple crosses
#' with overlapping parents.
#'
#' @param ... Objects of class \code{\link{calcIBD}}
#'
#' @return An object of class \code{\link{calcIBD}} containing data for all
#' concatenated objects.
#'
#' @export
c.calcIBD <- function(...) {
  args <- list(...)
  if (length(args) == 1) {
    return(args[[1]])
  }
  for (arg in args) {
    if (!inherits(arg, "calcIBD")) {
      stop("All inputs should be of class calcIBD.\n")
    }
  }
  pops <- unique(sapply(X = args, FUN = `[[`, "poptype"))
  if (!length(pops) == 1) {
    stop("All inputs should have the same population type.\n")
  }
  maps <- unique(lapply(X = args, FUN = `[[`, "map"))
  if (!length(maps) == 1) {
    stop("All inputs should have the same map.\n")
  } else {
    map <- maps[[1]]
  }
  markerLst <- lapply(X = args, FUN = `[[`, "markers")
  parentsNw <- unique(unlist(sapply(X = markerLst, FUN = function(mrk) {
    dimnames(mrk)[3]
  })))
  genoNw <- unlist(sapply(X = markerLst, FUN = colnames))
  nGeno <- sapply(X = markerLst, FUN = ncol)
  genoCross <- data.frame(cross = paste0("cross",
                                         rep(seq_along(nGeno), times = nGeno)),
                          geno = genoNw)
  markersNw <- array(dim = c(nrow(markerLst[[1]]), length(genoNw),
                             length(parentsNw)),
                     dimnames = list(rownames(markerLst[[1]]), genoNw, parentsNw))
  for (i in 1:nrow(markerLst[[1]])) {
    markersNw[i, , ] <- as.matrix(dfBind(lapply(X = markerLst, FUN = function(mrk) {
      as.data.frame(mrk[i, , ])
    })))
  }
  res <- structure(list(map = map,
                        markers = markersNw,
                        poptype = pops,
                        multicross = TRUE),
                   class = "calcIBD",
                   genoCross = genoCross)
  return(res)
}

#' Plot function for objects of class calcIBD
#'
#' Creates a plot for an object of class calcIBD.
#'
#' @param x An object of class \code{GWAS}.
#' @param ... Further arguments. Unused.
#' @param title A character string, the title of the plot.
#' @param output Should the plot be output to the current device? If
#' \code{FALSE}, only a ggplot object is invisibly returned.
#'
#' @export
plot.calcIBD <- function(x,
                         ...,
                         genotype,
                         title = NULL,
                         output = TRUE) {
  map <- x$map
  markers <- x$markers
  ## Convert to long format for plotting.
  markersLong <- expand.grid(snp = dimnames(markers)[[1]],
                             genotype = dimnames(markers)[[2]],
                             parent = dimnames(markers)[[3]])
  markersLong[["prob"]] <- c(markers)
  ## Merge map info to probabilities.
  plotDat <- merge(markersLong, map, by.x = "snp", by.y = "row.names")
  ## Restrict to selected genotype.
  plotDat <- plotDat[plotDat[["genotype"]] == genotype, ]
  ggplot2::ggplot(plotDat,
                  ggplot2::aes_string(x = "pos", y = "parent",
                                      fill = "prob")) +
    ggplot2::geom_tile(width = 3) +
    ggplot2::facet_grid(". ~ chr", scales = "free", space = "free",
                        switch = "both") +
    ggplot2::scale_fill_binned(type = "viridis") +
    ggplot2::labs(title = genotype) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}

