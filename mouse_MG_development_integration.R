### INTEGRATION OF EPITHELIAL MOUSE MG CELLS 
# R version 3.6.3

library(Seurat)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(cowplot)
library(future)
library(session) 
plan("multiprocess", workers=16) 
options(future.globals.maxSize = 10000 * 1024^2) 

w_dir <- "" # Set working directory HERE
if(!dir.exists(w_dir)) dir.create(w_dir)
metadata_columns <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent_mt", "doublets_score", "S.Score", "G2M.Score", "RNA_snn_res.0.4", "scMCA", "scMCA_simple")


# vectors with the barcodes of the epithelial cells for each sample
E16_epi_cells <- "E16_epi_cells.rds" 
E18_epi_cells <- "E18_epi_cells.rds"
P4_epi_cells <- "P4_epi_cells.rds"
WK5_epi_cells <- "WK5_epi_cells.rds"
WK12_epi_cells <- "WK12_epi_cells.rds"


### INTEGRATION
### ----------------------------------------------------------------------------
# load seurat object and subset epithelial cells
E16 <- readRDS("E16.rds")
E16_subset <- subset(E16, subset = cell_barcode %in% E16_epi_cells)
E16_subset@meta.data <- E16_subset@meta.data[,metadata_columns]

E18 <- readRDS("E18.rds")
E18_subset <- subset(E18, subset = cell_barcode %in% E18_epi_cells)
E18_subset@meta.data <- E18_subset@meta.data[,metadata_columns]

P4 <- readRDS("P4.rds")
P4_subset <- subset(P4, subset = cell_barcode %in% P4_epi_cells)
P4_subset@meta.data <- P4_subset@meta.data[,metadata_columns]

WK5 <- readRDS("WK5.rds")
WK5_subset <- subset(WK5, subset = cell_barcode %in% WK5_epi_cells)
WK5_subset@meta.data <- WK5_subset@meta.data[,metadata_columns]

WK12 <- readRDS("WK12.rds")
WK12_subset <- subset(WK12, subset = cell_barcode %in% WK12_epi_cells)
WK12_subset@meta.data <- WK12_subset@meta.data[,metadata_columns]

umi_list <- list(E16_subset, E18_subset, P4_subset, WK5_subset, WK12_subset)


# normalize and find variable features for each object in the list
umi_list_norm <- lapply(X=umi_list, FUN=function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method="vst", nfeatures=2000)
})


# integration
anchors <- FindIntegrationAnchors(object.list=umi_list_norm, dims=1:30)
umi_integrated <- IntegrateData(anchorset=anchors, dims=1:30)

DefaultAssay(object=umi_integrated) <- "integrated"
saveRDS(umi_integrated, file=file.path(w_dir, "mouseMG_epithelial_Integration.rds"))
### ----------------------------------------------------------------------------



### PREPROCESSING
### ----------------------------------------------------------------------------
pre_dir <- file.path(w_dir,"preprocessing/")
if(!dir.exists(pre_dir)) dir.create(pre_dir)

# scale data with cell cycle regression
regr_var <- c("S.Score", "G2M.Score")
umi_integrated_prepro <- ScaleData(umi_integrated, vars.to.regress=regr_var, verbose=FALSE)


# run PCA
umi_integrated_prepro <- RunPCA(umi_integrated_prepro, npcs=50, verbose=F)

pdf(paste0(pre_dir, "PCElbowPlot.pdf"), width=15, height=15, useDingbats=FALSE)
print(ElbowPlot(object=umi_integrated_prepro, ndims=50))
invisible(dev.off())


# find clusters
pc <- 17
res <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

umi_integrated_prepro <- FindNeighbors(umi_integrated_prepro, reduction="pca", dims=1:pc, force.recalc=T)
umi_integrated_prepro <- FindClusters(umi_integrated_prepro, resolution=res)


# dimensionality reduction
umi_integrated_prepro <- RunUMAP(umi_integrated_prepro, dims=1:pc)
### ----------------------------------------------------------------------------



### PLOTS
### ----------------------------------------------------------------------------
cluster_dir <- file.path(w_dir,"clustering/")
if(!dir.exists(cluster_dir)) dir.create(cluster_dir)
umi_prepro <- umi_integrated_prepro
DefaultAssay(umi_prepro) <- "integrated"

reduction <- c("UMAP", "PCA")
pt.s <- 0.5

# plot clusters
for (red in reduction){
  red_label <- ifelse(red == "PCA", "PC", red)
  plot <- NULL
  for (r in res){
    i <- match(r,res)
    res_col <- paste0("integrated_snn_res.",r)
    title <- paste("Clusters res",r)
    plot[[i]] <- DimPlot(umi_prepro, group.by=res_col, reduction=tolower(red), label=T, pt.size=pt.s, label.size = 10) + ggtitle(title) + xlab(paste(red_label,"1")) + ylab(paste(red_label,"2")) + NoLegend()
  }
  pdf(file.path(cluster_dir, paste0("mouseMG_integration_",red,"_clusters.pdf")), width=3*7, height=2*7, useDingbats=FALSE)
  print(plot_grid(plotlist=plot, ncol=3))
  invisible(dev.off())
}


# plot samples 
for (red in reduction){
  red_label <- ifelse(red == "PCA", "PC", red)
  pdf(file.path(cluster_dir, paste0("mouseMG_integration_",red,"_origident.pdf")), width=7, height=7, useDingbats=FALSE)
  print(DimPlot(umi_prepro, group.by="orig.ident", reduction=tolower(red), label=F, pt.size=pt.s) + ggtitle("Sample") + xlab(paste(red_label,"1")) + ylab(paste(red_label,"2")))
  invisible(dev.off())
}


# plot annotation
for (red in reduction){
  red_label <- ifelse(red == "PCA", "PC", red)
  print(red)
  p <- DimPlot(umi_prepro, reduction=tolower(red), label=F, pt.size=pt.s, group.by="scMCA_simple", cols = cols) + ggtitle("scMCA annotation") + 
    theme(plot.title=element_text(hjust=0.5)) + xlab(paste(red_label,"1")) + ylab(paste(red_label,"2"))
  legend <- get_legend(p)
  pdf(file = file.path(cluster_dir, paste0("mouseMG_integration_",red,"_Annotation.pdf")), width=2*7, height=7, useDingbats=F)
  print(plot_grid(p + NoLegend(), legend, ncol=2))
  invisible(dev.off())
}

saveRDS(umi_prepro, file=file.path(w_dir, "mouseMG_epithelial_Integration_prepro.rds"))



