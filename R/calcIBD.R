#' @export
summary.calcIBD <- function(object,
                            ...) {
  cat("population type: ", object$poptype, "\n")
  cat("Number of evaluation points: ", nrow(object$markers), "\n")
  cat("Number of individuals: ", ncol(object$markers),"\n")
  cat("Parents: ", substring(dimnames(object$markers)[[3]], first = 2), "\n")
}


#' @export
c.calcIBD <- function(...) {
  args <- list(...)
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
  res <- structure(list(map = maps,
                        markers = markersNw,
                        poptype = pops),
                   class = "calcIBD",
                   multicross = length(args) > 1,
                   genoCross = genoCross)
  return(res)
}
