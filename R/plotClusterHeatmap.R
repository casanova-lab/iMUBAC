#' @rdname plotClusterHeatmap
#' @title Plot cluster heatmap
#'
#' @description
#' Heatmaps summarizing the clustering result.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param features A character vector.
#'   Specifies which antigens to use for clustering.
#' @param clusters A character vector.
#'   Specifies the cluster IDs.
#' @param by_exprs_values A character string.
#'   Specifies which assay data to use for plotting.
#' @param fun A character string.
#'   Specifies the function for computing the summary statistic.
#' @param cluster_rows Logical.
#'   Specifies whether rows should be reordered by hierarchical clustering.
#' @param cluster_anno Logical.
#'   Specifies whether clusters should be annotated.
#' @param split_by A character string.
#'   Must corresponds to a column name of \code{rowData(sce)}.
#'   If specified, the data will be subset according to this variable,
#'   and multiple heatmaps will be drawn.
#' @param scale Logical.
#'   Specifies whether scaled values should be plotted.
#'   Scaled values corresponds to cofactor arcsinh-transformed expression
#'   values scaled between 0 and 1 using 1% and 99% percentiles as boundaries.
#'   Note that hierarchical clustering is performed on the unscaled data.
#' @param draw_dend Logical.
#'   Specifies if the row dendrogram should be drawn.
#' @param draw_freqs Logical.
#'   Specifyies whether to display cell counts and proportions.
#' @param hm2 A character string. Specifies the right-hand side heatmap.
#'   One of: \itemize{
#'   \item{\code{"abundances"}: cluster frequencies across samples}
#'   \item{a character string/vector corresponding to one/multiple marker(s):
#'     median marker expressions across samples and clusters}
#'   }
#'   Depending on argument \code{hm2}, the side heatmap can contain one of:
#'   \itemize{
#'   \item{relataive cluster abundances by sample}
#'   \item{median (scaled, arcsinh-transformed) marker expressions by sample}
#'   }
#'
#' @return a \code{\link{HeatmapList-class}} object.
#'
#' @export

plotClusterHeatmap <- function(
    sce,
    features=rownames(sce),
    clusters=sce$"cluster_id",
    by_exprs_values="exprs",
    fun="median",
    scale=T,
    cluster_rows=T,
    cluster_anno=F,
    draw_dend=T,
    draw_freqs=T,
    split_by=NULL,
    hm2=NULL
){
    # check validity of input arguments
    u <- c("abundances", features)
    if (!is.null(hm2)) stopifnot(hm2 %in% u)

    # Assign cluster IDs into a column named "cluster_id"
    ## This is essential for subsequent codes to work!
    if(is.null(levels(clusters))) clusters <- factor(clusters)
    sce$"cluster_id" <- clusters

    # Hierarchical clustering on marker medians by cluster
    nk <- nlevels(sce$"cluster_id")
    ms_by_k <- t(aggregateData(sce, by_exprs_values=by_exprs_values, by="cluster_id", fun=fun))
    d <- dist(ms_by_k[, features])
    if(cluster_rows){
        row_clustering <- hclust(d, method="average")
    }else{
        row_clustering <- FALSE
    }

    # Clustering row annotation
    if(cluster_anno){
        anno <- levels(sce$"cluster_id")
        if(nk > length(myCols)){
            cols <- colorRampPalette(myCols)(nk)
        }else{
            cols <- myCols[seq_len(nk)]
        }
        cols <- setNames(cols, anno)
        cluster_anno <- Heatmap(
            matrix=anno,
            col=cols,
            name="cluster_id",
            rect_gp=grid::gpar(col="white"),
            width=unit(.4, "cm"),
            cluster_rows=row_clustering,
            cluster_columns=F,
            show_row_dend=draw_dend,
            row_dend_reorder=F
        )
    }

    # Split cell indices by colData factor
    many <- !is.null(split_by)
    cs <- seq_len(ncol(sce))
    if(many) groups <- split(cs, sce[[split_by]]) else groups <- list(cs)

    # Scale expression matrix (optional)
    if(scale) assay(sce, by_exprs_values) <- scale_exprs(assay(sce, by_exprs_values) , 1)
    pals <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
    hm_cols <- colorRampPalette(pals)(100)

    # generate main heatmaps
    hms <- lapply(seq_along(groups), function(i) {
        idx <- groups[[i]]
        cs_by_k <- split(idx, sce$"cluster_id"[idx])
        # left-hand side heatmap:
        # median cell-type marker expressions across clusters
        if(!many){
            if(scale){
                hm1_es <- t(aggregateData(sce, by_exprs_values=by_exprs_values, by="cluster_id", fun=fun))
            }else{
                hm1_es <- ms_by_k
            }
        }else{
            hm1_es <- t(aggregateData(sce[, idx], by_exprs_values=by_exprs_values, by="cluster_id", fun=fun))
        }
        hm1 <- Heatmap(
            matrix=hm1_es[, features],
            col=hm_cols, name="expression",
            column_names_gp=grid::gpar(fontsize=8),
            rect_gp=grid::gpar(col='white'), na_col="lightgrey",
            cluster_rows=row_clustering, cluster_columns=FALSE,
            show_row_dend=draw_dend, column_title=names(groups)[i][many])

        # cluster frequencies
        freq_bars <- freq_anno <- NULL
        if(draw_freqs){
            fq <- round(tabulate(sce$"cluster_id"[idx]) / length(idx) * 100, 2)
            freq_bars <- rowAnnotation(
                "Freq [%]"=row_anno_barplot(
                    fq, axis=TRUE, border=FALSE, bar_with=.8,
                    gp=grid::gpar(fill="grey50", col="white")), width=unit(2, "cm")
                )
            labs <- paste0(levels(sce$"cluster_id"), " (", fq, "%)")
            freq_anno <- rowAnnotation(
                text=row_anno_text(labs),
                width=max_text_width(labs)
            )
        }

        # combine row annotations, heatmap,
        # and frequency bars & labels
        p <- hm1 + freq_bars + freq_anno
        if(is(cluster_anno, "Heatmap"))
            p <- cluster_anno + p

        # right-hand side heatmap
        if(!is.null(hm2)){
            if(hm2=="abundances"){
                # cluster frequencies across samples
                cs <- table(sce$"cluster_id"[idx], sce$"sample_id"[idx])
                fq <- as.matrix(unclass(prop.table(cs, 2)))
                fq <- fq[, !is.na(colSums(fq)), drop = FALSE]
                p <- p + Heatmap(
                    matrix=fq, name="frequency",
                    na_col="lightgrey", rect_gp=grid::gpar(col="white"),
                    show_row_names=FALSE, column_names_gp=grid::gpar(fontsize=8),
                    cluster_rows=row_clustering, cluster_columns=FALSE
                )
            }else{
                for(ch in hm2){
                    # median marker expression across samples & clusters
                    ms <- aggregateData(sce[ch, idx], by_exprs_values=by_exprs_values, by=c("cluster_id", "sample_id"), fun=fun)
                    ms <- do.call("rbind", ms)
                    rownames(ms) <- levels(sce$"cluster_id")
                    p <- p + Heatmap(
                        matrix=ms, col=hm_cols,
                        na_col="lightgrey", rect_gp=grid::gpar(col='white'),
                        show_heatmap_legend=FALSE, show_row_names=FALSE,
                        cluster_rows=row_clustering, cluster_columns=FALSE,
                        column_title=ch, column_names_gp=grid::gpar(fontsize=8)
                    )
                }
            }
        }
        return(p)
    })
    hm_list <- NULL
    for (i in seq_along(hms)) hm_list <- hm_list + hms[[i]]
    draw(hm_list)
    invisible(hm_list)
}
