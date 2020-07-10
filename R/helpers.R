# ==============================================================================
# Good colors
# ------------------------------------------------------------------------------
myCols <- c(
  RColorBrewer::brewer.pal(7, "Dark2"),
  RColorBrewer::brewer.pal(9, "Set1")
)

# ==============================================================================
# scale expression to values b/w 0 and 1 using
# low (1%) and high (99%) quantiles as boundaries
# ------------------------------------------------------------------------------
scale_exprs <- function(x, margin=1) {
  if (!is.matrix(x)) x <- as.matrix(x)
  qs <- c(matrixStats::rowQuantiles, matrixStats::colQuantiles)[[margin]]
  qs <- qs(x, probs = c(.01, .99))
  x <- switch(margin,
              "1" = (x - qs[, 1]) / (qs[, 2] - qs[, 1]),
              "2" = t((t(x) - qs[, 1]) / (qs[, 2] - qs[, 1])))
  x[x < 0 | is.na(x)] <- 0
  x[x > 1] <- 1
  return(x)
}

# ==============================================================================
# split cell indices by cell metadata factor(s)
#   - x:   an SCE with rows = cells, columns = features
#   - by:  colData columns specifying factor(s) to aggregate by
# ------------------------------------------------------------------------------
split_cells <- function(x, by) {
  stopifnot(is.character(by), by %in% colnames(colData(x)))
  cd <- data.frame(colData(x))
  dt <- data.table::data.table(cd, i = seq_len(ncol(x)))
  dt_split <- split(dt, by = by, sorted = TRUE, flatten = FALSE)
  purrr::map_depth(dt_split, length(by), "i")
}

# ==============================================================================
# aggregation of single-cell to pseudobulk data;
# e.g., median expression by cluster- or cluster-sample
#   - x:     an SCE object with rows = cells, columns = features
#   - assay: assay type to be used for computation
#   - by:    colData columns specifying factor(s) to aggregate by
#   - fun:   aggregation function specifying the
#            summary statistic, e.g., sum, mean, median
# ------------------------------------------------------------------------------
aggregateData <- function(
  x,
  by_exprs_values="exprs",
  by="cluster_id",
  fun=c("median", "mean", "sum")
){
  fun <- switch(match.arg(fun),
                median = rowMedians, mean = rowMeans, sum = rowSums)
  cs <- split_cells(x, by)
  pb <- purrr::map_depth(cs, -1, function(i) {
    if (length(i) == 0) return(numeric(nrow(x)))
    fun(assay(x, by_exprs_values)[, i, drop = FALSE])
  })
  purrr::map_depth(pb, -2, function(u) as.matrix(data.frame(
    u, row.names = rownames(x), check.names = FALSE)))
}
