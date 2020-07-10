#' @rdname runUMAP
#' @title Perform dimension reduction through UMAP
#'
#' @description
#' \code{runUMAP} takes an SCE object and computes UMAP coordinates using \pkg{uwot}.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param by_exprs_values A character string.
#'   Specifies which assay data to use for clustering.
#' @param name A character string.
#'   Specifies the name of the reduced dimensions.
#' @param n_components,n_neighbors,min_dist,scale,n_threads See \code{\link[uwot]{umap}}.
#' @param seed Numeric. Sets a random seed.
#'
#' @return A \code{\link{SingleCellExperiment-class}} object.
#'
#' @export

runUMAP <- function(
  sce,
  by_exprs_values="exprs",
  name="UMAP",
  n_components=2,
  n_neighbors=100,
  min_dist=0.5,
  scale=T,
  n_threads=parallel::detectCores(logical=F),
  seed=12345
){
  # Set the random seed
  set.seed(seed)

  # UMAP
  names <- paste0(name, 1:n_components)
  SingleCellExperiment::reducedDim(sce, type=name, withDimnames=T) <- uwot::umap(
    t(assay(sce, by_exprs_values)),
    n_components=n_components,
    n_neighbors=n_neighbors,
    min_dist=min_dist,
    scale=scale,
    n_threads=n_threads,
    verbose=T
  ) %>%
    magrittr::set_colnames(names)

  # Output
  return(sce)
}

