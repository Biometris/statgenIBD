#' Helper function for plotting the pedigree
#'
#' @keywords internal
pedPlot <- function(x,
                    title) {
  if (x$multiCross) {
    stop("Plot not implemented yet for multicross.\n")
  }
  pedDat <-  x$pedigree
  popType <- x$popType
  ## Restrict to parents and first progeny.
  pedDat <- pedDat[1:which.max(pedDat[["type"]] == popType), ]
  pedDat[nrow(pedDat), "ID"] <- "F1"
  generation <- as.numeric(factor(pedDat[["type"]],
                                  levels = unique(pedDat[["type"]]))) - 1

  ## Determine the row and column numbers in the plot.
  plotCols <- max(table(generation))
  plotRows <- length(unique(generation)) + 1
  ## Determine x and y positions of all parents and individuals in the plot
  xPos <- yPos <- NULL
  for (g in unique(generation)) {
    nGen <- sum(generation == g)
    if (nGen == plotCols) {
      xPosGen <- seq(1, plotCols, length.out = nGen)
    } else {
      xPosGen <- seq(1, plotCols, by = (plotCols - 1) / (nGen + 1))
      xPosGen <- xPosGen[-c(1, length(xPosGen))]
    }
    yPosGen <- rep(plotRows - g, nGen)
    xPos <- c(xPos, xPosGen)
    yPos <- c(yPos, yPosGen)
  }
  pedDat[["xPos"]] <- xPos
  pedDat[["yPos"]] <- yPos
  ## Construct data for plotting arrows.
  arrowDat <- rbind(merge(pedDat,
                          pedDat[, !colnames(pedDat) %in% c("par1", "par2")],
                          by.x = "par1", by.y = "ID"),
                    merge(pedDat,
                          pedDat[, !colnames(pedDat) %in% c("par1", "par2")],
                          by.x = "par2", by.y = "ID"))
  arrowDat[["linetype"]] <- "solid"
  ## Add extra arrow to bottom of plot.
  extArrow <- tail(arrowDat, 1)
  extArrow[["linetype"]] <- "dotted"
  extArrow[["yPos.y"]] <- extArrow[["yPos.x"]]
  extArrow[["yPos.x"]] <- extArrow[["yPos.x"]] - 1
  extArrow[["xPos.y"]] <- extArrow[["xPos.x"]]
  arrowDat <- rbind(arrowDat, extArrow)
  arrowDat[["linetype"]] <- factor(arrowDat[["linetype"]],
                                   levels = unique(arrowDat[["linetype"]]))
  ## Construct data for labels.
  labDat <- arrowDat[colnames(arrowDat) != "ID"]
  labDat <- merge(labDat, pedDat[c("ID", "xPos", "yPos")],
                  by.x = c(c("xPos.y", "yPos.y")), by.y = c("xPos", "yPos"))
  extLab <- extArrow
  extLab[["yPos.y"]] <- extLab[["yPos.x"]]
  extLab[["ID"]] <- popType
  labDat <- rbind(labDat, extLab)
  ## Construct texts.
  textDat <- data.frame(xPos.y = c(0 ,0, 0, tail(labDat[["xPos.y"]], 1)),
                        yPos.y = c(plotRows, 1, 0, 0),
                        text = c("Parent:", "Population type:", "size:",
                                 dim(x$markers)[2]))
  extText <- tail(arrowDat, 1)
  extText[["text"]] <- "selfing"
  extText[["xPos.y"]] <- (extText[["xPos.x"]] + extText[["xPos.y"]]) / 2
  extText[["yPos.y"]] <- (extText[["yPos.x"]] + extText[["yPos.y"]]) / 2
  textDat <- rbind(textDat, extText[c("xPos.y", "yPos.y", "text")])
  ## Construct title.
  if (is.null(title)) {
    title <- "pedigree"
  }
  ## Plot pedigree.
  ## Segments, arrows, labels and text separately and an additional
  ## segment for the arrow in the explanation part.
  p <- ggplot2::ggplot(arrowDat,
                  ggplot2::aes_string(x = "xPos.y", y = "yPos.y")) +
    ggplot2::geom_segment(ggplot2::aes_string(xend = "xPos.x",
                                              yend = "yPos.x",
                                              linetype = "linetype"),
                          size = 1, color = "blue") +
    ggplot2::geom_segment(ggplot2::aes_string(xend = "(xPos.x + xPos.y) / 2",
                                              yend = "(yPos.x + yPos.y) / 2"),
                          size = 1, color = "blue",
                          data = arrowDat[-nrow(arrowDat), ],
                          arrow = ggplot2::arrow(length = ggplot2::unit(0.3, "cm"),
                                                 type = "closed")) +
    ggplot2::geom_label(ggplot2::aes_string(label = "ID"),
                        data = labDat, fill = "white") +
    ggplot2::geom_text(ggplot2::aes_string(label = "text"),
                       data = textDat) +
    ggplot2::geom_segment(x = 0, y = plotRows - 0.2, xend = 0, yend  = 1.2,
                          size = 1, color = "blue", arrow = ggplot2::arrow()) +
    ggplot2::xlim(-0.5, plotCols + 0.5) +
    ggplot2::ylim(-0.5, plotRows + 0.5) +
    ggplot2::labs(x = "", y = "", title = title) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
          plot.background = ggplot2::element_blank(),
          strip.background = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(fill = NA),
          axis.text = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(hjust = 0.5),
          legend.position = "none")
  return(p)
}
