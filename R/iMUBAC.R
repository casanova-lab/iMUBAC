#' iMUBAC: Integration of Multi-Batch Cytometry datasets.
#'
#' @description High-dimentional cytometry, including spectral flow cytometry and mass cytometry (CyTOF), enables deep immunophenotyping at a single-cell resolution.
#'   Analysis of cytometry data from multiple batches of experiments remains to be challenging due to the high-dimentionality and batch effects.
#'   We designed iMUBAC (Integration of Multi-Batch Cytometry datasets) to enable a rational and streamlined inter-batch comparisons through i) preprocessing, ii) batch-correction, iii) unsupervised clustering, and iv) batch-specific cell-type identification.
#'
#' @name iMUBAC
#' @docType package
#' @aliases iMUBAC package-iMUBAC
#' @importFrom caret train trainControl
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus
#' @importFrom cydar poolCells dnaGate outlierGate
#' @importFrom data.table := data.table as.data.table rbindlist setcolorder
#' @importFrom dplyr %>% mutate mutate_all transmute filter bind_rows group_by summarize sample_n tally n
#' @importFrom extraTrees extraTrees
#' @importFrom flowCore flowFrame colnames colnames<- exprs exprs<- keyword Subset fsApply biexponentialTransform transformList transform estimateLogicle
#' @importFrom FlowSOM BuildSOM ReadInput
#' @importFrom ggedit remove_geom
#' @importFrom ggpubr compare_means
#' @importFrom ggrastr geom_point_rast
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom grid gpar
#' @importFrom harmony HarmonyMatrix
#' @importFrom igraph cluster_louvain
#' @importFrom janitor make_clean_names
#' @importFrom magrittr set_names set_colnames set_rownames
#' @importFrom Matrix rowMeans rowSums
#' @importFrom matrixStats rowMedians rowQuantiles colQuantiles
#' @importFrom methods is
#' @importFrom ncdfFlow read.ncdfFlowSet as.flowSet
#' @importFrom parallel detectCores
#' @importFrom purrr is_empty map_depth map
#' @importFrom RColorBrewer brewer.pal
#' @importFrom S4Vectors metadata metadata<- DataFrame
#' @importFrom scater plotReducedDim
#' @importFrom scran buildSNNGraph
#' @importFrom SummarizedExperiment assay assay<- colData colData<- rowData rowData<-
#' @importFrom stats dist hclust setNames predict
#' @importFrom stringr str_detect
#' @importFrom tidyr complete
#' @importFrom uwot umap
#' @importFrom viridis scale_color_viridis scale_colour_viridis scale_fill_viridis
#' @import SingleCellExperiment
#' @import ComplexHeatmap
#' @import ggplot2
NULL
