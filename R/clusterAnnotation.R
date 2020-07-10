#' @rdname clusterAnnotation
#' @title Integration of manual cluster annotations
#'
#' @description
#' \code{clusterAnnotation} takes both an SCE object returned by
#' \code{clusterPropagation} and a data.frame of annotations as inputs.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#'   \code{colData(sce)} must contain the cluster_id column.
#' @param df_anno A data.frame which has the following columns: cluster_id, celltype,
#' celltype_detailed, and order. The order column will be used to define the order of
#' the factor levels for cell types.
#'
#' @return A \code{\link{SingleCellExperiment-class}} object.
#'
#' @export

clusterAnnotation <- function(
  sce,
  df_anno
){
  df_anno <- df_anno %>%
    dplyr::arrange(order)
  sce$"celltype" <- factor(
    sce$"cluster_id",
    levels=df_anno$"cluster_id",
    labels=df_anno$"celltype"
  )
  sce$"celltype_detailed" <- factor(
    sce$"cluster_id",
    levels=df_anno$"cluster_id",
    labels=df_anno$"celltype_detailed"
  )
  return(sce)
}

