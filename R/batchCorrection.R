#' @rdname batchCorrection
#' @title Batch-correction for CyTOF data
#'
#' @description
#' \code{batchCorrection} takes an SCE object returned by \code{prepSCE} and
#' performs batch correction using \pkg{Harmony}.
#' Down-sampling is done prior to batch correction for improving speed.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#'   \code{colData(sce)} must contain file_name and batch columns.
#' @param maxN Numeric.
#'   Specifies the maximum number of unique cells per batch. For example,
#'   if one batch contains five samples, the function randomly takes \code{maxN/5}
#'   cells from each of the samples.
#' @param seed Numeric. Sets a random seed.
#'
#' @return A \code{\link{SingleCellExperiment-class}} object.
#'
#' @export

batchCorrection <- function(
  sce,
  maxN=100000,
  seed=12345
){
  # Set the random seed
  set.seed(seed)

  # Down-sampling
  sce <- downsampleSCE(
    sce,
    maxN=maxN,
    group_by=c("batch"),
    indiv_by="file_name",
    seed=seed
  )

  # Batch correction
  assay(sce, "normexprs") <- t(harmony::HarmonyMatrix(
    data_mat=t(assay(sce, "exprs")),
    meta_data=as.data.frame(colData(sce)),
    vars_use=c("file_name", "batch"),
    do_pca=F,
    plot_convergence=T,
    verbose=T
  ))

  # Output
  return(sce)
}

