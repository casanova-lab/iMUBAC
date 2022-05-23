#' @rdname clusterPropagation
#' @title Batch-wise propagation of clusters through machine learning
#'
#' @description
#' \code{clusterPropagation} takes both an original SCE object returned by
#' \code{prepSCE} and down-sampled, batch-corrected, and clustered SCE object
#' returned by \code{clustering} as inputs.
#' This function trains batch-specific classifiers using a specified algorithm,
#' and predicts cluster IDs in a batch-wise manner.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#'   \code{colData(sce)} must contain file_name and batch columns.
#' @param sce_down A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#'   \code{colData(sce_down)} must contain columns named as file_name, batch, and
#'   cluster_id.
#' @param by_exprs_values A character string.
#'   Specifies which assay data to use for classifier training and prediction.
#'   Must be the non-batch-corrected one in this case.
#' @param maxN Numeric.
#'   Specifies the maximum number of unique cells per cluster per batch to be used
#'   for classifier training.
#' @param numThreads Numeric. The number of threads for training classifiers.
#' @param seed Numeric. Sets a random seed.
#'
#' @return A \code{\link{SingleCellExperiment-class}} object.
#'
#' @export

clusterPropagation <- function(
  sce,
  sce_down,
  by_exprs_values="exprs",
  maxN=100,
  numThreads=4,
  seed=12345
){
  # The main function
  clusterPropagation_SingleBatch <- function(sce, sce_down){

    ## Set the random seed
    set.seed(seed)

    ## Down-sample to reduce computational burden & deal with class imbalance
    sce_down <- downsampleSCE(
      sce_down,
      maxN=maxN,
      group_by=c("batch","cluster_id"),
      indiv_by=NULL,
      seed=seed
    )

    ## Train a classifier
    cat("-Training a classifier...\n", sep="")
    require(caret)
    require(doParallel)
    cl <- makePSOCKcluster(numThreads)
    registerDoParallel(cl)
    trCtrls <- caret::trainControl(
      method="repeatedcv",
      number=10,
      repeats=5,
      sampling="up",
      verboseIter=F
    )
    mat <- t(assay(sce_down, by_exprs_values))
    tunegrid <- expand.grid(.mtry=5:15)
    #tunegrid <- expand.grid(.mtry=5:15, .numRandomCuts=1:2) ## for extraTrees
    caret.model <- caret::train(
      x=mat,
      y=droplevels(sce_down$"cluster_id"),
      preProcess=c("center","scale"),
      method="rf",
      metric="Kappa",
      tuneGrid=tunegrid,
      trControl=trCtrls
    )
    stopCluster(cl)

    ## Predict cell types
    cat("-Predicting clusters...\n", sep="")
    mat <- t(assay(sce, by_exprs_values))
    clusterIDs <- predict(caret.model, newdata=mat)
    clusterIDs <- as.numeric(as.character(clusterIDs)) ### remove factor levels

    ## Output
    return(clusterIDs)
  }

  # Initialize
  sce$"cluster_id" <- 0

  # Batch-wise classifier training & prediction
  for(b in levels(sce$"batch")){
    cat("Batch: ", b, "\n", sep="")
    clusterIDs <- clusterPropagation_SingleBatch(
      sce[,sce$"batch"==b],
      sce_down[,sce_down$"batch"==b]
    )
    sce[,sce$"batch"==b]$"cluster_id" <- clusterIDs
    gc();gc()
  }
  sce$"cluster_id" <- factor(
    sce$"cluster_id",
    levels=unique(c(0,as.numeric(as.character(levels(sce_down$"cluster_id")))))
  )

  # Output
  return(sce)
}

