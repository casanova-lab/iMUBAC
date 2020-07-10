#' @rdname downsampleSCE
#' @title Down-sampling an SCE object
#'
#' @description
#' \code{downsampleSCE} down-samples an SCE object based on the grouping variables
#' provided. For example, when grouped by a batch variable, the function takes
#' maxN cells per batch. Additionally, when \code{indiv_by} is provided (e.g.,
#' file_name), the function takes equal numbers of cells per individual, collecting
#' maxN cells in total per batch.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param maxN Numeric.
#'   Specifies the maximum number of unique cells per group.
#' @param group_by A character vector.
#'   Specifies the colData columns for grouping.
#' @param indiv_by A character string.
#'   Specifies the colData column for individual sample identifiers.
#' @param seed Numeric. Sets a random seed.
#'
#' @return A \code{\link{SingleCellExperiment-class}} object.
#'
#' @export

downsampleSCE <- function(
  sce,
  maxN=10000,
  group_by=c("batch","group"),
  indiv_by=NULL,
  seed=12345
){
  set.seed(seed)
  sce$"unique_cell_id" <- 1:ncol(sce)
  dt <- colData(sce) %>%
    as.data.frame() %>%
    data.table::as.data.table()
  if(!is.null(indiv_by)){
    dt$"indiv_by" <- dt[[indiv_by]]
    dt[,"n_individual":=length(unique(indiv_by)),by=group_by]
    dt <- dt[, sample(unique_cell_id, min(maxN/n_individual, .N)), by=indiv_by]
    colnames(dt)[length(indiv_by)+1] <- "unique_cell_id"
  }else{
    dt <- dt[, sample(unique_cell_id, min(maxN, .N)), by=group_by]
    colnames(dt)[length(group_by)+1] <- "unique_cell_id"
  }
  sce <- sce[,sce$"unique_cell_id" %in% dt$"unique_cell_id"]
  sce$"unique_cell_id" <- NULL
  return(sce)
}

