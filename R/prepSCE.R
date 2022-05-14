#' @rdname prepSCE
#' @title Batch-wise preprocessing of CyTOF FCS files
#'
#' @description
#' \code{prepSCE} imports FCS files from a single batch as an NCDF format,
#' format channel names, internally gates out dead cells and doublets based
#' on the DNA and/or dead-cell-exclusion dyes provided, performs hyperbolic
#' arcsin transformation, and returns a \code{SingleCellExperiment} object.
#' When multiple batches of FCS files are provided, the preprocessing steps
#' are performed in a batch-wise manner (in other words, in a loop), and the
#' final results are merged into a single SCE object. Markers available in
#' only a subset of batches will be discarded.
#'
#' @param md A data.frame/data.table containing the following metadata:
#' \itemize{
#' \item{full_path: the full path to the FCS file.}
#' \item{batch: the batch ID. A factor is preferred.}
#' }
#'   Other columns can also be inherited in the final SCE object.
#' @param pd A list of data.frames/data.tables containing the following metadata:
#' \itemize{
#' \item{chennel: the channel names.}
#' \item{antigen: the corresponding antigen names.}
#' }
#'   Note that exclusion dyes do not have to be included here.
#'   The list must be named with the corresponding batch IDs.
#' @param file_name_keyword A character string.
#'   Specifies the keyword of the FCS files to be matched with the file names.
#'   Note that file name matching will be performed using the base, not full, file names.
#' @param channel_length A character string.
#'   Specifies the channel name for the event length. Can also be NULL.
#' @param channel_DNA A character vector of length two.
#'   Specifies the channel names for DNA dyes. Can also be NULL.
#' @param channel_LD A character string.
#'   Specifies the channel name for the live/dead exclusion dye. Can also be NULL.
#' @param type A string indicating the type of the underlying data. Can be either "CyTOF" or "Flow".
#'
#' @return A \code{\link{SingleCellExperiment-class}} object.
#'
#' @export

prepSCE <- function(
  md,
  pd,
  file_name_keyword="FILENAME",
  channel_length="Event_length",
  channel_DNA=c("Ir191Di","Ir193Di"),
  channel_LD="Rh103Di",
  type="CyTOF"
){
  # Parse FCS files
  if(!is.factor(md$"batch")) md$"batch" <- factor(md$"batch", levels=sort(unique(md$"batch")))
  md$"batch" <- droplevels(md$"batch")
  batch_list <- as.character(levels(md$"batch"))
  if(identical(type, "CyTOF")){
    trans <- F ### Note: The automatic options in the flowCore package are optimized for flow cytometry instead of mass cytometry data, so these options should be disabled.
  }else if(identical(type, "Flow")){
    trans <- "linearize"
  }else{
    stop("'type' must be either 'CyTOF' or 'Flow'!")
  }
  sce_list <- lapply(batch_list, function(b){
    ## Load data
    cat("Loading Batch ", b, "...\n", sep="")
    md <- as.data.frame(md[which(as.character(md$"batch")==b)])
    pd <- pd[[b]]
    fs <- suppressMessages(ncdfFlow::read.ncdfFlowSet(
      files=md$"full_path",
      transformation=trans,
      truncate_max_range=F
    ))

    ## Check file names
    ids <- try(flowCore::keyword(fs, file_name_keyword), silent=T)
    if(identical(class(ids),"try-error")) stop("File name keyword not found!")
    ids <- match(basename(ids), basename(md$"full_path"))
    fs <- fs[ids]

    ## Format channel names
    cat("-Formatting channel names...\n", sep="")
    pd$"channel" <- janitor::make_clean_names(pd$"channel", case="parsed")
    flowCore::colnames(fs) <- janitor::make_clean_names(flowCore::colnames(fs), case="parsed")
    channel_length <- janitor::make_clean_names(channel_length, case="parsed")
    channel_DNA <- janitor::make_clean_names(channel_DNA, case="parsed")
    channel_LD <- janitor::make_clean_names(channel_LD, case="parsed")
    chs_ex <- c(channel_length,channel_DNA,channel_LD)
    pd <- pd[pd$"channel" %in% flowCore::colnames(fs),]
    chs0 <- c(pd$"channel",chs_ex)
    chs <- c(pd$"antigen",chs_ex)
    fs <- fs[,chs0]
    flowCore::colnames(fs) <- chs

    ## Automated gating to exclude dead cells and doublets
    if(!purrr::is_empty(chs_ex)){
      cat("-Excluding dead cells and doublets...\n", sep="")
      ### Design gates using pooled cells from all samples in the batch
      pool.ff <- cydar::poolCells(fs)
      trans <- flowCore::estimateLogicle(pool.ff, chs_ex)
      proc.ff <- flowCore::transform(pool.ff, trans)
      if(!purrr::is_empty(channel_DNA)){
        gate.dna <- cydar::dnaGate(proc.ff, channel_DNA[1], channel_DNA[2], type="both") ### removing cells with too high or too low values
        proc.ff <- flowCore::Subset(proc.ff, gate.dna)
      }
      if(!purrr::is_empty(channel_LD)){
        gate.dead <- cydar::outlierGate(proc.ff, channel_LD, type="upper") ### removing cells with too high values
        proc.ff <- flowCore::Subset(proc.ff, gate.dead)
      }
      if(!purrr::is_empty(channel_length)){
        gate.db <- cydar::outlierGate(proc.ff, channel_length, type="upper") ### removing cells with too high values
        proc.ff <- flowCore::Subset(proc.ff, gate.db)
      }
      ### Apply gates for individual samples
      fs <- suppressMessages(flowCore::transform(fs, trans))
      if(!purrr::is_empty(channel_DNA)) fs <- flowCore::Subset(fs, gate.dna)
      if(!purrr::is_empty(channel_LD)) fs <- flowCore::Subset(fs, gate.dead)
      if(!purrr::is_empty(channel_length)) fs <- flowCore::Subset(fs, gate.db)
    }

    ## Convert to FlowSet
    fs <- ncdfFlow::as.flowSet(fs)

    ## Transformation
    if(identical(type, "CyTOF")){
      cat("-Hyperbolic arcsin transformation...\n", sep="")
      fs <- flowCore::fsApply(fs, function(ff){
        flowCore::exprs(ff) <- asinh(flowCore::exprs(ff)/5)
        return(ff)
      })
    }else if(identical(type, "Flow")){
      cat("-Logicle transformation...\n", sep="")
      tmax <- max(flowCore::fsApply(fs, function(ff){max(ff@exprs)}))
      fs <- flowCore::fsApply(fs, function(ff){
        trans <- flowCore::transformList(
          from=setdiff(chs, chs_ex),
          tfun=flowCore::logicleTransform(w=0.5, t=tmax, m=4.5, a=0)
        )
        #trans <- flowCore::estimateLogicle(ff, channels=setdiff(chs, chs_ex), t=max(ff@exprs))
        ff <- flowCore::transform(ff, trans)
        return(ff)
      })
    }else{
      stop("'type' must be either 'CyTOF' or 'Flow'!")
    }

    ## Construct an SCE object
    cat("-Constructing an SCE object...\n", sep="")
    es <- matrix(flowCore::fsApply(fs, exprs), byrow=T, nrow=length(chs), dimnames=list(chs, NULL))
    md$"n_cells" <- as.numeric(flowCore::fsApply(fs, nrow))
    rd <- DataFrame(row.names=chs, channel_name=chs0, marker_name=chs)
    cd <- DataFrame(lapply(md[c("full_path","batch")], function(u){rep(u, md$"n_cells")}), row.names=NULL)
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays=list(exprs=es),
      rowData=rd,
      colData=cd
    )
    sce <- sce[which(!rownames(sce) %in% chs_ex),] ### discard channels only for dead cell/doublet exclusion

    ## Quality-check
    rm <- apply(assay(sce,"exprs"),2,function(x){all(x==0)})
    if(sum(rm)>=1){
      message(paste0(sum(rm), " events (", round(sum(rm)/length(rm)*100), "%) with zeros for all markers are removed."))
      sce <- sce[,which(!rm)]
    }

    ## Output
    return(sce)
  })
  names(sce_list) <- batch_list

  # Sort rownames for consistencies between batches
  ## Markers missing in some of the sce objects are discarded.
  rn_int <- sort(Reduce("intersect", lapply(sce_list, rownames)))
  rn_union <- sort(Reduce("union", lapply(sce_list, rownames)))
  for(i in 1:length(sce_list)){
    sce_list[[i]] <- sce_list[[i]][rn_int,]
  }
  if(!all(rn_union %in% rn_int)) message("Markers missing in some of the batches are discarded.")

  # Merge SCE objects from multiple batches
  sce <- Reduce("rbind", sce_list)
  # sce <- suppressMessages(sce_cbind(
  #   sce_list,
  #   method="intersect",
  #   cut_off_batch=0,
  #   cut_off_overall=0,
  #   exprs="exprs",
  #   colData_names=c("full_path","batch"),
  #   batch_names=batch_list
  # ))
  sce$"batch" <- factor(sce$"batch", levels=c(batch_list))

  # Integrate metadata
  x <- setdiff(colnames(md), c("full_path","batch"))
  if(!purrr::is_empty(x)){
    dt_meta <- data.table::data.table("full_path"=sce$"full_path", "batch"=sce$"batch")
    dt_meta <- merge(dt_meta, md, by=c("full_path","batch"), all.x=T, all.y=F, sort=F)
    for(y in x){
      sce[[y]] <- dt_meta[[y]]
      if(is.factor(sce[[y]])) sce[[y]] <- droplevels(sce[[y]])
    }
  }
  if(!"file_name" %in% colnames(colData(sce))) sce$"file_name" <- basename(sce$"full_path")
  sce$"full_path" <- NULL

  # Output
  return(sce)
}

