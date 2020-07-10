#' @rdname updateMetadataSCE
#' @title Updating metadata of an SCE object
#'
#' @description
#' \code{updateMetadataSCE} takes a new metadata object and updates the corresponding
#' columns of the colData of an SCE object.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param md A data.frame/data.table containing the following metadata:
#' \itemize{
#' \item{file_name: the file name of the underlying FCS file.}
#' }
#'   Variables that do not exist in the colData of the SCE object will be ignored.
#' @param merge_by A character vector.
#'   Specifies the names of the variables used for merging.
#'
#' @return A \code{\link{SingleCellExperiment-class}} object.
#'
#' @export

updateMetadataSCE <- function(
  sce,
  md,
  merge_by=c("file_name")
){
  dt_meta <- data.table::as.data.table(as.data.frame(colData(sce)))
  cols <- intersect(colnames(dt_meta), colnames(md))
  dt_meta <- merge(dt_meta[,merge_by,with=F], md, all.x=T, all.y=F, sort=F)
  for(var in setdiff(cols, merge_by)){
    v <- dt_meta[[var]]
    if(is.factor(v)) v <- droplevels(v)
    sce[[var]] <- v
  }
  return(sce)
}

