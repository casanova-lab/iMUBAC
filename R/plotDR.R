#' @rdname plotDR
#' @title Plot reduced dimensions
#'
#' @description
#' Dimension reduction plot colored by expression values or column metadata.
#' To avoid plotting too many dots as a vector graphic, \code{geom_point_rast} from \pkg{ggrastr}
#' is used to internally replace the \code{geom_point} layer.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param dimred A character string specifying which dimension reduction to use.
#'   Should be one of \code{reducedDimNames(sce)}.
#' @param colour_by A character string corresponding to a \code{colData(sce)} column.
#'   Specifies the color coding. If factor, a discrete color scale will be used. Otherwise, a cividis color scale will be applied by default.
#' @param text_by A character string corresponding to a \code{colData(sce)} column.
#'   Specifies the text coding.
#' @param ... Additional arguments passed to \code{plotReducedDim} in \pkg{scater}.
#'
#' @return A \code{ggplot} object.
#'
#' @export

plotDR <- function(
  sce,
  dimred="UMAP",
  colour_by="condition",
  text_by=NULL,
  ...
){
  p <- scater::plotReducedDim(
    sce,
    dimred=dimred,
    colour_by=colour_by,
    text_by=text_by,
    ...
  )
  p <- ggedit::remove_geom(p, geom="point")
  p$layers <- c(
    ggrastr::geom_point_rast(aes(colour=colour_by), shape="."),
    p$layers
  )
  if(is.factor(sce[[colour_by]])){
    nk <- nlevels(sce[[colour_by]])
    if(nk > length(myCols)){
      cols <- colorRampPalette(myCols)(nk)
    }else{
      cols <- myCols[seq_len(nk)]
    }
    p <- p +
      scale_colour_manual(name=colour_by, values=cols) +
      guides(colour=guide_legend(override.aes=list(shape=19, size=3, alpha=1)))
  }else{
    p <- p +
      viridis::scale_color_viridis(option="cividis")
  }
  p <- p +
    ggpubr::theme_pubr(14) +
    theme(aspect.ratio=1) +
    ggpubr::rremove("legend.title")
  return(p)
}
