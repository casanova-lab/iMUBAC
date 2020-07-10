iMUBAC: Integration of Multi-Batch Cytometry Datasets
===============================================

The 'iMUBAC' package provides a structured framework for objective inter-batch comparisons and unbiased immunophenotyping of high-dimentional cytometry datasets.

System Requirements
------------------------
iMUBAC has been tested on R versions >= 4.0 on Windows platform. 

Installation
------------------------
Install the latest version as follows:
``` r
if(!require(devtools)) install.packages("devtools")
devtools::install_github("VPetukhov/ggrastr")
devtools::install_github("immunogenomics/harmony")
devtools::install_github("casanova-lab/iMUBAC")
devtools::install_github("masato-ogishi/plotUtility")
```
-  You may be prompted to install some packages before installing iMUBAC. Follow the messages.
-  This package depends on some packages in the [*Bioconductor*](https://www.bioconductor.org/) (e.g., flowCore) that may not be automatically installed. Please check the DESCRIPTION file for more details on required packages and install them manually as required.
-  You need an appropriate rJava setting beforehand.

Usage
------------------
0. Working environment
``` r
options(java.parameters="-Xmx60G")  ## allow JAVA to use large memory space
library(tidyverse)
library(data.table)
library(SingleCellExperiment)
library(iMUBAC)
```
1. Data
-  For demonstration, previously published CyTOF datasets from Krieg et al. are included in the package.
-  Pre-treatment datasets stained with the Panel 3 (myeloid cell panel) from two different batches are included.
``` r
# Sample metadata
md <- data.table::fread(system.file("Metadata.csv", package="iMUBAC"))
md[,full_path:=file.path(system.file(package="iMUBAC"),file_name)]
md[,batch:=factor(batch, levels=c("Data23","Data29"))]
md[,group:=factor(group, levels=c("HD","R","NR"))]

# Panel design
pd <- list(
  data.table::fread(system.file("PanelInfo_Data23_Panel3.csv", package="iMUBAC")),
  data.table::fread(system.file("PanelInfo_Data29_Panel3.csv", package="iMUBAC"))
) %>% magrittr::set_names(c("Data23","Data29"))

# Preprocessing
sce <- iMUBAC::prepSCE(
  md=md[panel=="Panel3"],
  pd=pd,
  channel_length=NULL,      ## Doublet detection by event length is disabled.
  channel_DNA=c("Ir191Di","Ir193Di"),
  channel_LD="Pt198Di",
  type="CyTOF"
)
colData(sce) <- colData(sce)[c("file_name","panel","batch","donor_id","group","treatment")]
saveRDS(sce, "CyTOF_SCE_Merged_Panel3.rds")
```
2. Batch-correction
-  Features can be calculated as follows. Note: This computation is time-consuming and resource-intensive. Computation can be resumed if temporary files are stored in the temporary directory provided.

``` r
# Load the preprocessed data
sce <- readRDS("CyTOF_SCE_Merged_Panel3.rds")

# Take the subset of data from healthy donors
sce_down <- sce[,sce$"group"=="HD"]

# Batch-correction
sce_down <- iMUBAC::batchCorrection(
  sce_down,
  maxN=50000, ## A maximum of 50000 cells are randomly selected from each batch.
  seed=12345  ## a random seed
)
```
-  Optionally, we can explore the impact of batch-correction through UMAP.
``` r
# UMAP
sce_down <- iMUBAC::runUMAP(
  sce_down,
  by_exprs_values="exprs",
  name="UMAP",
  n_neighbors=25,
  min_dist=0.4,
  scale=T,
  n_threads=4, ## the number of threads for parallel computing
  seed=12345
)
sce_down <- iMUBAC::runUMAP(
  sce_down,
  by_exprs_values="normexprs",
  name="UMAPnorm",
  n_neighbors=25,
  min_dist=0.4,
  scale=T,
  n_threads=4,
  seed=12345
)

# Plot
plt1 <- ggpubr::ggarrange(
  iMUBAC::plotDR(
    sce_down,
    dimred="UMAP",
    colour_by="batch"
  ) +
    scale_colour_brewer(palette="Dark2"),
  iMUBAC::plotDR(
    sce_down,
    dimred="UMAPnorm",
    colour_by="batch"
  ) +
    scale_colour_brewer(palette="Dark2"),
  ncol=2, nrow=1, common.legend=T, legend="right"
)
plotUtility::savePDF(plt1, o="Plt1.pdf", w=15, h=6)
```
3. Unsupervised clustering
-  Currently two methods are implemented. 
-  The first approach groups cells into \code{xdim}x\code{ydim} clusters using \pkg{FlowSOM}, and then performs metaclustering with \pkg{ConsensusClusterPlus} into \code{maxK} clusters. 
-  The second approach performs dimention reduction through UMAP, constructs shared nearest-neighbor graphs using \pkg{scran}, and then performs community detection using the Leuvein algorithm in \pkg{igraph}.
``` r
# FlowSOM-guided clustering & ConsensusClusterPlus-guided metaclustering
sce_down <- iMUBAC::clustering(
  sce_down,
  features=rownames(sce_down), ## Using all markers for clustering
  by_exprs_values="normexprs",
  method="FlowSOM",
  xdim=20,
  ydim=20,
  maxK=40, ## the number of metaclusters returned
  seed=12345
)

# Alternatively, SNN-graph-guided clustering
## Not used in this demonstration
sce_down_snn <- iMUBAC::clustering(
  sce_down,
  features=rownames(sce_down), ## Using all markers for clustering
  by_exprs_values="normexprs",
  method="SNNGraph",
  n_components=10, ## the number of reduced dimentions for constructing the SNN graph
  n_neighbors=25,  ## the parameters for UMAP-based dimention reduction
  min_dist=0.4,    ## the parameters for UMAP-based dimention reduction
  seed=12345
)
``` 
-  We can explore the clusters through the median expression heatmap and UMAP plots.
``` r
# Median expression heatmap
plt2 <- iMUBAC::plotClusterHeatmap(
  sce_down,
  features=rownames(sce_down),
  clusters=sce_down$"cluster_id",
  by_exprs_values="normexprs",
  fun="median",
  scale=T,
  cluster_rows=T,
  cluster_anno=T,
  draw_dend=T,
  draw_freqs=T
)
plotUtility::savePDF(plt2, o="Plt2.pdf", w=12, h=8)

# UMAP plot
plt3 <- iMUBAC::plotDR(
  sce_down,
  dimred="UMAPnorm",
  colour_by="cluster_id",
  text_by="cluster_id"    ## to overlay cluster ids on each of the clusters
) +
  ggpubr::rremove("legend")
plotUtility::savePDF(plt3, o="Plt3.pdf", w=10, h=6)
```

4. Batch-specific cluster propagation through machine learning
-  Classifiers are trained for each batch utilizing the Extreme Randomized Trees algorithm.
-  There are both abundant and rare clusters. To mitigate the class imbalance problem and also to reduce the computational burden, here a maximum of 100 cells are randomly selected for each of the clusters before classifier training.
``` r
sce <- iMUBAC::clusterPropagation(
  sce,           ## the original data containing cells from controls and patients (non-batch-corrected)
  sce_down,      ## down-sampled data containing cells from controls (batch-corrected and clustered)
  by_exprs_values="exprs", ## Non-batch-corrected expression values are used for classifier training.
  maxN=100,
  numThreads=4,  ## the number of threads for parallel computing
  seed=12345
)
```
5. Cell-type identification through manual annotation
-  Users can provide cell-type annotations for each cluster. The data.frame should contain four columns: "cluster_id", "celltype", "celltype_detailed", and "order". 
-  The "celltype" and "celltype_detailed" columns can be identical but may be useful to provide different layers of annotation (e.g., "CD4 T" and "CD4 TEMRA").
-  The "order" column will be used to define the order of factor levels for cell types.
``` r
# Annotation
df_celltype <- readxl::read_excel("CyTOF_ClusterCellTypes_Panel3.xlsx")
sce <- clusterAnnotation(sce, df_celltype)
sce_down <- clusterAnnotation(sce_down, df_celltype)
saveRDS(sce, "CyTOF_SCE_Cluster_Panel3.rds")
saveRDS(sce_down, "CyTOF_SCE_Cluster_Panel3_Down.rds")

# Median expression heatmap
## Removing unnecessary clusters
sce_down <- sce_down[,!sce_down$"celltype_detailed" %in% c("Basophil","Unidentified")]
sce_down$"celltype_detailed" <- droplevels(sce_down$"celltype_detailed")
plt4 <- iMUBAC::plotClusterHeatmap(
  sce_down,
  features=rownames(sce_down),
  clusters=sce_down$"celltype_detailed",
  by_exprs_values="normexprs",
  fun="median",
  scale=T,
  cluster_rows=F,
  cluster_anno=F,
  draw_dend=T,
  draw_freqs=T
)

# UMAP
## Recompute UMAP coordinates after removing unnecessary clusters
sce_down <- sce_down[,!sce_down$"celltype_detailed" %in% c("Basophil","Unidentified")]
sce_down$"celltype_detailed" <- droplevels(sce_down$"celltype_detailed")
sce_down <- iMUBAC::runUMAP(
  sce_down,
  by_exprs_values="normexprs",
  name="UMAP",
  n_neighbors=25,
  min_dist=0.4,
  scale=T,
  n_threads=4,
  seed=12345
)
plt5 <- iMUBAC::plotDR(
  sce_down,
  dimred="UMAP",
  colour_by="celltype_detailed",
  text_by="celltype_detailed",
  point_size=1.5
) +
  theme(legend.position="right",
        legend.direction="vertical") +
  ggpubr::rremove("legend.title")
```
6. Differential abundance analysis
- A number of workflows are available for differential abundance (DA) analysis. Here, we use the QLF test in edgeR to detect DA subsets between non-responders and responders to PD-1 blockade immunotherapy.
``` r
library(edgeR)
sce <- readRDS("CyTOF_SCE_Cluster_Panel3.rds")
dt_da <- colData(sce) %>%
  as.data.frame() %>%
  dplyr::group_by(donor_id, batch, group, celltype_detailed) %>%
  dplyr::tally() %>%
  dplyr::mutate(percent=n/sum(n)*100) %>%
  tidyr::complete(tidyr::nesting(donor_id, batch, group), celltype_detailed, fill=list(n=0, percent=0)) %>%
  data.table::as.data.table()
dt_da[,batch:=factor(batch, levels=c("Data23","Data29"), labels=c("Batch1","Batch2"))]
dt_da[,variable:=paste0(group,"_",batch)]
dt_da[,variable:=factor(variable, levels=c("HD_Batch1", "HD_Batch2", "R_Batch1", "R_Batch2", "NR_Batch1", "NR_Batch2"))]
dt_da[,label:=paste0(group,"_",batch,"_",donor_id)]
dt_da_wide <- data.table::dcast(dt_da, celltype_detailed~label, value.var="n")
celltypeLabels <- dt_da_wide$"celltype_detailed"
dt_da_wide[,celltype_detailed:=NULL]
mat <- as.matrix(dt_da_wide)
rownames(mat) <- celltypeLabels
groups <- factor(
  colnames(mat),
  levels=unique(dt_da, by="label")$"label",
  labels=unique(dt_da, by="label")$"group"
)
batches <- factor(
  colnames(mat),
  levels=unique(dt_da, by="label")$"label",
  labels=unique(dt_da, by="label")$"batch"
)
d <- edgeR::DGEList(counts=mat, lib.size=colSums(mat), group=groups, remove.zeros=T)
design <- model.matrix(~0+groups+batches, data=d$samples)
d <- edgeR::estimateDisp(d, design)
fit <- edgeR::glmQLFit(d, design, robust=T)
dt_da_qlf <- edgeR::glmQLFTest(
  fit, 
  contrast=makeContrasts(NRvsR=groupsNR-groupsR, levels=design)
) %>%
  edgeR::topTags(n=Inf) %>%
  as.data.frame() %>%
  dplyr::transmute(
    celltype_detailed=rownames(.),
    log2FC=logFC, 
    PValue=PValue,
    AdjPValue=FDR
  ) %>%
  data.table::as.data.table()
```

Reference
------------------------
Ogishi, M et al. (2020) "Multi-batch cytometry data integration toward unbiased and streamlined immunophenotyping." *bioRxiv*. 
Krieg, C et al. (2018) "High-dimensional Single-Cell Analysis Predicts Response to anti-PD-1 Immunotherapy." *Nature Medicine*. https://pubmed.ncbi.nlm.nih.gov/29309059/
