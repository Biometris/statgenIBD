#' Helper function for plotting mean of probabilities.
#'
#' @noRd
#' @keywords internal
meanProbsPlot <- function(markers,
                          map,
                          parents,
                          chr = NULL,
                          title = NULL) {
  if (!all(chr %in% map[["chr"]])) {
    stop("chr not found in map.\n")
  }
  ## Compute means.
  plotDat <- apply(X = markers, MARGIN = c(1, 3), FUN = mean)
  colnames(plotDat) <- parents
  ## Restrict map to selected chromosomes.
  if (!is.null(chr)) {
    map <- map[map[["chr"]] %in% chr, ]
    plotDat <- plotDat[rownames(plotDat) %in% rownames(map), ]
  }
  ## convert to data.frame:
  plotDat <- as.data.frame.table(plotDat)
  ## Construct title.
  if (is.null(title)) {
    title <- "Coverage per parent"
  }
  p <- ggplot2::ggplot(plotDat,
                       ggplot2::aes_string(x = "Var1", y = "Freq",
                                           color = "Var2", group = "Var2")) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Position", y = "Representation", title = title) +
    ggplot2::scale_color_discrete(name = "parent") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank())
  return(p)
}
