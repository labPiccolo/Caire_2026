### INTEGRATION OF EPITHELIAL NORMAL BREAST CELLS 
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
metadata_columns <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent_mt", "doublets_score", "S.Score", "G2M.Score", "BpEn.sc.main.labels", "epi_class")

hMG1_epi_cells <- "hMG1_epi_cells.rds" # vector with the barcodes of the epithelial cells in hMG1, as specified in hMG
hMG2_epi_cells <- "hMG2_epi_cells.rds" # vector with the barcodes of the epithelial cells in hMG2, as specified in hMG
hMG3_epi_cells <- "hMG3_epi_cells.rds" # vector with the barcodes of the epithelial cells in hMG3, as specified in hMG
# epi_class column in metadata contains the annotation of the epithelial cells in hMG, as specified in hMG (L1, L2, Basal)

### INTEGRATION
### ----------------------------------------------------------------------------
# load object and subset epithelial cells
hMG1 <- readRDS("hMG1.rds") # load the hMG1 Seurat object - change path if needed
hMG1_subset <- subset(hMG1, subset = cell_barcode %in% hMG1_epi_cells)
hMG1_subset@meta.data <- hMG1_subset@meta.data[,metadata_columns]

hMG2 <- readRDS("hMG2.rds") # load the hMG2 Seurat object - change path if needed
hMG2_subset <- subset(hMG2, subset = cell_barcode %in% hMG2_epi_cells)
hMG2_subset@meta.data <- hMG2_subset@meta.data[,metadata_columns]

hMG3 <- readRDS("hMG3.rds") # load the hMG3 Seurat object - change path if needed
hMG3_subset <- subset(hMG3, subset = cell_barcode %in% hMG3_epi_cells)
hMG3_subset@meta.data <- hMG3_subset@meta.data[,metadata_columns]

umi_list <- list(hMG1_subset, hMG2_subset, hMG3_subset)


# normalize and find variable features for each object in the list
umi_list_norm <- lapply(X=umi_list, FUN=function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method="vst", nfeatures=2000)
})



# CCA integration
anchors_cca <- FindIntegrationAnchors(object.list=umi_list_norm, dims=1:30, anchor.features=2000, reduction="cca")
umi_integrated_cca <- IntegrateData(anchorset=anchors_cca, dims=1:30)

DefaultAssay(object=umi_integrated_cca) <- "integrated"
saveRDS(umi_integrated_cca, file=file.path(w_dir, "hMG_epithelialIntegration.rds"))
### ----------------------------------------------------------------------------



### PREPROCESSING
### ----------------------------------------------------------------------------
scaling_data <- function(umi, sample = "UMI", vars = NULL, pre_dir = getwd(), suffix = "noReg", mkdir=TRUE, get_umi=FALSE){
  reg_dir <- paste0(pre_dir,suffix,"/")
  if (mkdir==TRUE) {if (!file.exists(reg_dir)){dir.create(reg_dir)}} else {reg_dir <- pre_dir}
  print(reg_dir)
  
  DefaultAssay(umi) <- "integrated"
  all.genes <- rownames(umi)
  
  
  print("Scaling data, vars to regress:")
  print(vars)
  umi_scaled <- ScaleData(object=umi, features=all.genes, vars.to.regress=vars)
  print("Run PCA")
  umi_scaled <- RunPCA(umi_scaled, npcs = 30, verbose = FALSE)
  
  print("Fast clustering using 10 pcs and res 0.05 and 0.6")
  umi_scaled <- FindNeighbors(umi_scaled, reduction="pca", dims=1:10, verbose=F)
  umi_scaled <- FindClusters(umi_scaled, resolution=c(0.05,0.6), verbose=F)
  umi_scaled$integrated_snn_res.0.05 <- as.character(umi_scaled$integrated_snn_res.0.05)
  umi_scaled$integrated_snn_res.0.6 <- as.character(umi_scaled$integrated_snn_res.0.6)
  
  print("Run dimensional reduction")
  umi_scaled <- RunTSNE(umi_scaled, dims=1:10, seed.use=1, perplexity=30)
  umi_scaled <- RunUMAP(umi_scaled, dims=1:10, verbose=F)
  
  print("Generating plots")
  pt.s <- 0.5
  reduction <- c("PCA","tSNE","UMAP")
  to.plot <- c("nFeature_RNA","nCount_RNA","percent_mt","S.Score","G2M.Score","doublets_score","KRT5","SLPI","ANKRD30A","orig.ident","integrated_snn_res.0.05","integrated_snn_res.0.6")
  title <- c("nFeature_RNA","nCount_RNA","percent_mt","S.Score","G2M.Score","doublets_score","KRT5","SLPI","ANKRD30A","Sample","Fast clusters res 0.05","Fast clusters res 0.6")
  
  DefaultAssay(umi_scaled) <- "RNA"
  
  for (red in reduction){
    print(red)
    red_label <- ifelse(red == "PCA", "PC", red)
    plot <- NULL
    for (m in to.plot){
      i <- match(m,to.plot)
      if (m=="orig.ident") {
        plot[[i]] <- DimPlot(umi_scaled, reduction = tolower(red), group.by = m, pt.size = pt.s, label = FALSE) + ggtitle(title[i]) + xlab(paste0(red_label, " 1")) + ylab(paste0(red_label, " 1"))}
      if (m=="integrated_snn_res.0.05" || m=="integrated_snn_res.0.6") {
        plot[[i]] <- DimPlot(umi_scaled, reduction = tolower(red), group.by = m, pt.size = pt.s, label = TRUE, label.size = 6, repel = TRUE) + ggtitle(title[i]) + xlab(paste0(red_label, " 1")) + ylab(paste0(red_label, " 1"))}
      if (m=="S.Score" || m=="G2M.Score" || m=="doublets_score" || m=="nFeature_RNA" || m=="nCount_RNA" || m=="percent_mt" || m=="KRT5" || m=="SLPI" || m=="ANKRD30A") {
        plot[[i]] <- FeaturePlot(umi_scaled, reduction = tolower(red), features = m, pt.size = pt.s) + ggtitle(title[i]) + xlab(paste0(red_label, " 1")) + ylab(paste0(red_label, " 1"))}
    }
    pdf(paste0(reg_dir,sample,"_",red,"_",suffix,".pdf"), width=3*7, height=4*7, useDingbats=FALSE)
    print(plot_grid(plotlist=plot, ncol=3))
    invisible(dev.off())
  }
  
  DefaultAssay(umi_scaled) <- "integrated"
  
  pdf(paste0(reg_dir, sample,"_PCHeatmap_",suffix,".pdf"), width=9, height=30, useDingbats=FALSE)
  DimHeatmap(umi_scaled, dims=1:30, cells=500, balanced=T)
  invisible(dev.off())
  
  if (get_umi==TRUE) {return(umi_scaled)}
}

pre_dir <- file.path(w_dir,"preprocessing/")
if(!dir.exists(pre_dir)) dir.create(pre_dir)

# chosen to proceed with no REGRESSION
regr_var <- NULL
umi_integrated_cca_prepro <- scaling_data(umi_integrated_cca, sample="hMG_PD_integrated", vars=regr_var, pre_dir=pre_dir, suffix="noReg", mkdir=FALSE, get_umi=TRUE)


# 15 PCs chosen
pc <- 15

### find clusters
set.seed(123)
umi_integrated_cca_prepro <- FindNeighbors(umi_integrated_cca_prepro, reduction="pca", dims=1:pc, force.recalc=T)
umi_integrated_cca_prepro <- FindClusters(umi_integrated_cca_prepro, resolution=c(0.05,0.1,0.2,0.4,0.6,0.8))

# @res 0.05
# Number of communities: 3
# @res 0.1
# Number of communities: 5
# @res 0.2
# Number of communities: 7
# @res 0.4
# Number of communities: 11
# @res 0.6
# Number of communities: 14
# @res 0.8
# Number of communities: 17

# dimensionality reduction
umi_integrated_cca_prepro <- RunUMAP(umi_integrated_cca_prepro, dims=1:pc)
### ----------------------------------------------------------------------------



### PLOTS
### ----------------------------------------------------------------------------
cluster_dir <- file.path(w_dir,"clustering/")
if(!dir.exists(cluster_dir)) dir.create(cluster_dir)
umi_prepro <- umi_integrated_cca_prepro
DefaultAssay(umi_prepro) <- "integrated"

reduction <- c("UMAP", "PCA")
res <- c(0.05,0.1,0.2,0.4,0.6,0.8)
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
  pdf(file.path(cluster_dir, paste0("hMG_integration_",red,"_clusters.pdf")), width=3*7, height=2*7, useDingbats=FALSE)
  print(plot_grid(plotlist=plot, ncol=3))
  invisible(dev.off())
}


# samples distribution
for (red in reduction){
  red_label <- ifelse(red == "PCA", "PC", red)
  pdf(file.path(cluster_dir, paste0("hMG_integration_",red,"_origident.pdf")), width=7, height=7, useDingbats=FALSE)
  print(DimPlot(umi_prepro, group.by="orig.ident", reduction=tolower(red), label=F, pt.size=pt.s) + ggtitle("Sample") + xlab(paste(red_label,"1")) + ylab(paste(red_label,"2")))
  invisible(dev.off())
}


# plot markers
DefaultAssay(umi_prepro) <- "RNA"

markers <- c("CD3D","CD79A","CD14","PECAM1","COL1A1","ACTA2","MKI67","EPCAM", "KRT5","KRT18","SLPI","ANKRD30A")
markers_to_plot <- markers[markers %in% rownames(umi_prepro$RNA)]

for (red in reduction){
  red_label <- ifelse(red == "PCA", "PC", red)
  plot <- NULL
  for (m in markers_to_plot) {
    i <- match(m,markers_to_plot)
    plot[[i]] <- FeaturePlot(umi_prepro, features=m, reduction=tolower(red), pt.size=pt.s) + xlab(paste(red_label,"1")) + ylab(paste(red_label,"2")) + theme(plot.title = element_text(hjust = 0))
  }
  pdf(file.path(cluster_dir, paste0("hMG_integration_",red,"_markers.pdf")), width=3*7, height=4*7, useDingbats=FALSE)
  print(plot_grid(plotlist=plot, ncol=3, nrow=4))
  invisible(dev.off())
}


# Plot annotation
cols <- c("#F6766D","#00B838","#619AFF","#EEEEEE")
names(cols) <- c("Basal","L1","L2","Others")

for (red in reduction){
  red_label <- ifelse(red == "PCA", "PC", red)
  print(red)
  p <- DimPlot(umi_prepro, reduction=tolower(red), label=F, pt.size=pt.s, group.by="epi_class", cols = cols) + ggtitle("Manual Annotation") + 
    theme(plot.title=element_text(hjust=0.5)) + xlab(paste(red_label,"1")) + ylab(paste(red_label,"2"))
  legend <- get_legend(p)
  pdf(file = file.path(cluster_dir, paste0("hMG_integration_",red,"_Annotation.pdf")), width=2*7, height=7, useDingbats=F)
  print(plot_grid(p + NoLegend(), legend, ncol=2))
  invisible(dev.off())
}

saveRDS(umi_prepro, file=file.path(w_dir, "hMG_PD_epithelialIntegration_prepro.rds"))
