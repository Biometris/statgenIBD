#' Summary function for objects of class IBDprob
#'
#' Prints a short summary for objects of class \code{IBDprob}. The
#' summary consists of the population type, number of evaluation points,
#' number of individuals and names of the parents in the object.
#'
#' @param object An object of class \code{IBDprob}.
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
summary.IBDprob <- function(object,
                            ...) {
  cat("population type: ", object$poptype, "\n")
  cat("Number of evaluation points: ", nrow(object$markers), "\n")
  cat("Number of individuals: ", ncol(object$markers),"\n")
  cat("Parents: ", object$parents, "\n")
}

#' Concatenate function for objects of class IBDprob
#'
#' Concatenates objects of class \code{IBDprob}. All objects that are
#' concatenated  should have the same population type and the same map. The
#' function is mainly meant for combining information for multiple crosses
#' with overlapping parents.
#'
#' @param ... Objects of class \code{IBDprob}.
#'
#' @return An object of class \code{IBDprob} containing data for all
#' concatenated objects.
#'
#' @export
c.IBDprob <- function(...) {
  args <- list(...)
  if (length(args) == 1) {
    return(args[[1]])
  }
  for (arg in args) {
    if (!inherits(arg, "IBDprob")) {
      stop("All inputs should be of class IBDprob.\n")
    }
  }
  pops <- unique(sapply(X = args, FUN = `[[`, "poptype"))
  if (!length(pops) == 1) {
    stop("All inputs should have the same population type.\n")
  }
  parents <- unique(c(sapply(X = args, FUN = `[[`, "parents")))
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
                        parents = parents,
                        multicross = TRUE),
                   class = "IBDprob",
                   genoCross = genoCross)
  return(res)
}

#' Plot function for objects of class IBDprob
#'
#' Creates a plot for an object of class \code{IBDprob}.
#'
#' @param x An object of class \code{IBDprob}.
#' @param ... Further arguments. Unused.
#' @param genotype A character string indicating the genotype for which the
#' plot should be made.
#' @param title A character string, the title of the plot.
#' @param output Should the plot be output to the current device? If
#' \code{FALSE}, only a ggplot object is invisibly returned.
#'
#' @return A ggplot object is invisibly returned.
#'
#' @export
plot.IBDprob <- function(x,
                         ...,
                         genotype,
                         title = genotype,
                         output = TRUE) {
  map <- x$map
  markers <- x$markers
  ## Input checks.
  if (!inherits(genotype, "character") || length(genotype) > 1) {
    stop("genotype should be a character string.\n")
  }
  if (!inherits(title, "character") || length(title) > 1) {
    stop("title should be a character string.\n")
  }
  if (!inherits(genotype, "character") || length(genotype) > 1) {
    stop("genotype should be a character string.\n")
  }
  if (!genotype %in% dimnames(markers)[[2]]) {
    stop("genotype should be in markers.\n")
  }
  ## Convert to long format for plotting.
  markersLong <- markers3DtoLong(x)
  ## Merge map info to probabilities.
  plotDat <- merge(markersLong, map, by.x = "snp", by.y = "row.names")
  ## Restrict to selected genotype.
  plotDat <- plotDat[plotDat[["genotype"]] == genotype, ]
  p <- ggplot2::ggplot(plotDat,
                  ggplot2::aes_string(x = "pos", y = "parent",
                                      fill = "prob")) +
    ggplot2::geom_tile(width = 3) +
    ggplot2::facet_grid(". ~ chr", scales = "free", space = "free",
                        switch = "both") +
    ggplot2::scale_fill_gradient(low = "white", high = "black") +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::labs(title = genotype) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.background = ggplot2::element_blank(),
      panel.spacing.x = ggplot2::unit(4, "mm"))
  if (output) {
    plot(p)
  }
  invisible(p)
}

#' Helper function for converting 3D probability matrix to df.
#'
#' Helper function for converting 3D probability matrix to df.
#'
#' @noRd
#' @keywords internal
markers3DtoLong <- function(x) {
  markers <- x$markers
  parents <- x$parents
  markerCols <- dimnames(markers)[[3]]
  ## Create base data.frame for storing long format data.
  markersLongBase <- expand.grid(snp = dimnames(markers)[[1]],
                                 genotype = dimnames(markers)[[2]])
  markersLong <- NULL
  for (parent in parents) {
    ## Construct parent column.
    parentCol <- paste0("p", parent)
    ## Get other columns containing parent.
    parentSubCols <- markerCols[grep(pattern = parent, x = markerCols)]
    parentSubCols <- parentSubCols[-which(parentSubCols == parentCol)]
    ## Add values for parent to base.
    markersParent <- markersLongBase
    markersParent[["parent"]] <- parent
    ## Compute probability for parent.
    ## (2 * pPar + psubPar) / 2
    markersParent[["prob"]] <- c(markers[, , parentCol] +
      apply(X = markers[, , parentSubCols], MARGIN = 1:2, FUN = sum) / 2)
    ## Add to markersLong
    markersLong <- rbind(markersLong, markersParent)
  }
  return(markersLong)
}


