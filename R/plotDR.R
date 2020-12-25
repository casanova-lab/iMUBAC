#' @rdname plotDR
#' @title Plot reduced dimensions
#'
#' @description
#' Dimension reduction plot colored by expression values or column metadata.
#' To avoid plotting too many dots as a vector graphic, \code{geom_point_rast} from \pkg{ggrastr}
#' is used to internally replace the \code{geom_point} layer.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param sce_bg An optional \code{\link[SingleCellExperiment]{SingleCellExperiment}} object for background data. If provided, the data would be depicted in a light grey in a layer behind the main plot.
#' @param dimred A character string specifying which dimension reduction to use.
#'   Should be one of \code{reducedDimNames(sce)}.
#' @param colour_by A character string corresponding to a \code{colData(sce)} column.
#'   Specifies the color coding. If factor, a discrete color scale will be used. Otherwise, a cividis color scale will be applied by default.
#' @param text_by A character string corresponding to a \code{colData(sce)} column.
#'   Specifies the text coding.
#' @param point_shape,point_size,point_alpha Parameters for the plot. Note: setting the shape to "." allows users to nicely plot a large number of dots.
#' @param rasterize Logical. Whether the plot needs to be rasterized using \pkg{ggrastr}. Note: facetting may cause the Cairo device to fail.
#' @param ... Additional arguments passed to \code{plotReducedDim} in \pkg{scater}.
#'
#' @return A \code{ggplot} object.
#'
#' @export

plotDR <- function(
  sce,
  sce_bg=NULL,
  dimred="UMAP",
  colour_by="condition",
  text_by=NULL,
  point_shape=19,
  point_size=0.5,
  point_alpha=1,
  rasterize=T,
  ...
){
  p <- scater::plotReducedDim(
    sce,
    dimred=dimred,
    colour_by=colour_by,
    text_by=text_by,
    ...
  )
  if(rasterize){
    p <- ggedit::remove_geom(p, geom="point")
    p$layers <- c(
      ggrastr::geom_point_rast(
        aes(colour=colour_by),
        shape=point_shape,
        size=point_size,
        alpha=point_alpha,
        raster.dpi=300
      ),
      p$layers
    )
  }else{
    p <- ggedit::remove_geom(p, geom="point")
    p$layers <- c(
      geom_point(
        aes(colour=colour_by),
        shape=point_shape,
        size=point_size,
        alpha=point_alpha
      ),
      p$layers
    )
  }
  if(!is.null(sce_bg)){
    p_bg <- scater::plotReducedDim(
      sce_bg,
      dimred=dimred
    )
    if(rasterize){
      p$layers <- c(
        ggrastr::geom_point_rast(
          data=p_bg$"data",
          colour="grey90",
          shape=point_shape,
          size=point_size,
          alpha=point_alpha,
          raster.dpi=300
        ),
        p$layers
      )
    }else{
      p$layers <- c(
        geom_point(
          data=p_bg$"data",
          colour="grey90",
          shape=point_shape,
          size=point_size,
          alpha=point_alpha
        ),
        p$layers
      )
    }
  }
  if(is.factor(sce[[colour_by]])){
    nk <- nlevels(sce[[colour_by]])
    if(nk > length(myCols)){
      cols <- colorRampPalette(myCols)(nk)
    }else{
      cols <- myCols[seq_len(nk)]
    }
    p <- p +
      scale_colour_manual(name=colour_by, values=cols) +
      guides(colour=guide_legend(nrow=1, override.aes=list(shape=19, size=3, alpha=1)))
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
