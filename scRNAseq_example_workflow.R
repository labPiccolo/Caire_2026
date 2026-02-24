# R version 3.6.3

library(popsicleR) # version 0.2.1
library(Seurat)    # version 3.1.5
library(SingleR)   # version 1.0.6
library(infercnv)  # version 1.2.1
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(future) 
plan("multiprocess", workers=16) #parallel with 16 cores
options(future.globals.maxSize = 1500 * 1024^2) #1.5 GB

set.seed(123)

library(reticulate)
Sys.setenv(RETICULATE_PYTHON = ".../python3.6")    # set the path to your python executable here - used for Scrublet (version 0.2.1)
scrublet_path <- ".../Scrublet/do_scrublet_V.2.py" # set the path to scrublet.py here - used for Scrublet (version 0.2.1)

w_dir <- ""               # set the working directory HERE
data_dir <- ""            # set the data directory (where Cellranger output is located) HERE
sample_name <- ""         # set the sample name HERE
gene_position_file <- ""  # set the path to the gene positions file (e.g. cellranger-GRCh38-3.0.0_gen_pos.txt) HERE - used for infercnv (version 1.8.0)


### 1 - LOAD DATA AND QUALITY CONTROLS (USING POPSICLER)
### ----------------------------------------------------------------------------------------------------------
genes <- c("CD3D","VIM","PTPRC","EPCAM","MALAT1","CD79A") # exploratory genes for QC
umi <- PrePlots(
    sample_name = sample_name,
    input_data = data_dir,
    genelist = genes,
    percentage = 0.1, gene_filter = 200,
    cellranger = TRUE,
    organism = "human", # set to "mouse" if working with mouse data
    out_folder = w_dir
    )

# filtering based on QC metrics (popsicleR)

nFeature_RNA_lo <- 1000   # change this value based on the distribution of nFeature_RNA in your data
nFeature_RNA_hi <- NULL
nCount_RNA_lo   <- NULL
nCount_RNA_hi   <- 120000 # change this value based on the distribution of nCount_RNA in your data
percent_mt_hi   <- 25     # change this value based on the distribution of percent.mt in your data

umi_manualFilter <- FilterPlots(
    UMI = umi,
    G_RNA_low = nFeature_RNA_lo,
    U_RNA_hi = nCount_RNA_hi,
    percent_mt_hi = percent_mt_hi,
    out_folder = w_dir
    )
### ----------------------------------------------------------------------------------------------------------



### 2 - DOUBLET DETECTION (USING SCRUBLET)
### ----------------------------------------------------------------------------------------------------------
source_python(scrublet_path)
	
mtx_transposed <- BiocGenerics::t(umi_manualFilter@assays$RNA@counts)
doublets <- do_doublets(mtx_transposed, 0.1, FALSE)
file.rename("doublets_distribution.png", paste0(w_dir, sample_id, "_doublets_distribution.png"))

thr <- 0.28 # set the threshold for doublet score based on the distribution of scores in your data

names(doublets) <- c("score", "predicted")
doublets$predicted <- ifelse(doublets$score > thr,"TRUE","FALSE")

umi_manualFilter@meta.data$doublets_score <- doublets$score
umi_manualFilter@meta.data$doublets <- doublets$predicted

pdf(paste0(w_dir, sample_id, "_histogram_doublets.pdf"), useDingbats=FALSE)
h <- hist(doublets$score, breaks="Scott", plot=FALSE)
plot(h$mids, h$density, log="y", type='h', lwd = 4)
abline(v=thr, col="red")
invisible(dev.off())

cells <- rownames(umi_manualFilter@meta.data[umi_manualFilter$doublets=="FALSE",])
umi_manualFilter <- subset(umi_manualFilter, cells=cells)
### ----------------------------------------------------------------------------------------------------------



### 3 - NORMALIZE AND SCALE DATA (USING SEURAT)
### ----------------------------------------------------------------------------------------------------------
pre_dir <- file.path(w_dir,"prepro/")
if(!dir.exists(pre_dir)) dir.create(pre_dir)

# data normalization
umi_manualFilter <- NormalizeData(object=umi_manualFilter, normalization.method="LogNormalize", scale.factor=1e4)

pdf(paste0(pre_dir, sample_id, "_total_expression_before+after_norm.pdf"), useDingbats=FALSE)
par(mfrow = c(2,1))
hist(colSums(as.matrix(umi_manualFilter@assays$RNA@counts)), breaks=100, main="Total expression before normalization", xlab="Sum of expression")
hist(colSums(as.matrix(umi_manualFilter@assays$RNA@data)), breaks=100, main="Total expression after normalization", xlab="Sum of expression")
invisible(dev.off())

# detection of highly variable genes   
umi_manualFilter <- FindVariableFeatures(umi_manualFilter, selection.method="vst", nfeatures=2000)
  
top10 <- head(VariableFeatures(umi_manualFilter), 10)
  
plot1 <- VariableFeaturePlot(umi_manualFilter)
plot2 <- LabelPoints(plot=plot1, points=top10, repel=TRUE)
pdf(paste0(pre_dir, sample_id, "_plot_FindVariableGenes.pdf"), width=7, height=2*7, useDingbats=FALSE)
plot_grid(plot1, plot2, ncol=1)
invisible(dev.off())
  
# Cell-Cycle Scoring 
umi_manualFilter <- CellCycleScoring(umi_manualFilter, s.features=cc.genes$s.genes, g2m.features=cc.genes$g2m.genes, set.ident=T) # for mouse data change the cell cycle genes

pdf(paste0(pre_dir, sample_id, "_plot_CellCycleGenes.pdf"), useDingbats=FALSE)
print(RidgePlot(umi_manualFilter, features=c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol=2))
invisible(dev.off())


# Data scaling and linear dimensional reduction function
scaling_data <- function(umi, sample = "UMI", vars = NULL, pre_dir = getwd(), suffix = "noReg", mkdir=TRUE, get_umi=FALSE){
  reg_dir <- paste0(pre_dir,suffix,"/")
  if (mkdir==TRUE) {if (!file.exists(reg_dir)){dir.create(reg_dir)}} else {reg_dir <- pre_dir}
  print(reg_dir)
  all.genes <- rownames(umi)
  
  print("Scaling data, vars to regress:")
  print(vars)
  umi_scaled <- ScaleData(object=umi, features=all.genes, vars.to.regress=vars)
  print("Run PCA")
  umi_scaled <- RunPCA(object=umi_scaled, features=umi_scaled[["RNA"]]@var.features, npcs=30, verbose=F)
  
  print("Fast clustering using 10 pcs and res 0.6")
  umi_scaled <- FindNeighbors(umi_scaled, reduction="pca", dims=1:10, verbose=F)
  umi_scaled <- FindClusters(umi_scaled, resolution=0.6, verbose=F)
  
  print("Run dimensional reduction")
  umi_scaled <- RunTSNE(umi_scaled, dims=1:10, seed.use=1, perplexity=30)
  umi_scaled <- RunUMAP(umi_scaled, dims=1:10, verbose=F)
  
  print("Generating plots")
  pt.s <- 1
  reduction <- c("PCA","tSNE","UMAP")
  to.plot <- c("seurat_clusters","Phase","doublets_score","nFeature_RNA","nCount_RNA","percent_mt","S.Score","G2M.Score")
  title <- c("Clusters @ res0.6","CC phase","Doublets","","","","","")
  for (red in reduction){
    print(red)
    plot <- NULL
    for (m in to.plot){
      i <- match(m,to.plot)
      if (m=="orig.ident" || m=="Phase" || m=="Cell_type") plot[[i]] <- DimPlot(umi_scaled, group.by=m, reduction=tolower(red),label=F, pt.size=pt.s) + ggtitle(title[i])
      if (m=="seurat_clusters") plot[[i]] <- DimPlot(umi_scaled, group.by=m, reduction=tolower(red),label=T, pt.size=pt.s) + ggtitle(title[i])
      if (m=="nFeature_RNA" || m=="nCount_RNA" || m=="percent_mt" || m=="S.Score" || m=="G2M.Score" || m=="doublets_score") plot[[i]] <- FeaturePlot(umi_scaled, features=m, reduction=tolower(red),pt.size=pt.s) + theme(plot.title = element_text(hjust = 0))
    }
    pdf(paste0(reg_dir,"03d_",sample,"_",red,"_",suffix,".pdf"), width=3*7, height=3*7, useDingbats=FALSE)
    print(plot_grid(plotlist=plot, ncol=3))
    invisible(dev.off())
  }
  
  pdf(paste0(reg_dir,sample,"_PCHeatmap_",suffix,".pdf"), width=9, height=30, useDingbats=FALSE)
  DimHeatmap(umi_scaled, dims=1:30, cells=500, balanced=T)
  invisible(dev.off())
  
  if (get_umi==TRUE) {return(umi_scaled)}
}

# try to scale with no regression, and then with cell cycle regression to see if it improves the clustering and dimensional reduction results
regr_var <- NULL
scaling_data(umi_manualFilter, sample=sample_id, vars=regr_var, pre_dir=pre_dir, suffix="noReg")

regr_var <- c("S.Score", "G2M.Score")
scaling_data(umi_manualFilter, sample=sample_id, vars=regr_var, pre_dir=pre_dir, suffix="ccReg")

# proceed either with no regression or with cell cycle regression
all.genes <- rownames(umi_manualFilter)
regr_var <- NULL # regr_var <- c("S.Score", "G2M.Score")

umi_prepro <- ScaleData(object=umi_manualFilter, features=all.genes, vars.to.regress=regr_var)
umi_prepro <- RunPCA(object=umi_prepro, features=umi_prepro[["RNA"]]@var.features, npcs=30, verbose=F)

pdf(paste0(pre_dir, sample_id, "_PCElbowPlot.pdf"), width=15, height=15, useDingbats=FALSE)
print(ElbowPlot(object=umi_prepro, ndims=30))
invisible(dev.off())
### ----------------------------------------------------------------------------------------------------------



### 4 - CLUSTERING AND UMAP (USING SEURAT)
### ----------------------------------------------------------------------------------------------------------
cluster_dir <- file.path(w_dir,"clustering/")
if(!dir.exists(cluster_dir)) dir.create(cluster_dir)
pc <- 15 # change this value based on the ElbowPlot - tipically between 10 and 20
	
# find clusters
umi_prepro <- FindNeighbors(umi_prepro, reduction="pca", dims=1:pc, force.recalc=T)
umi_prepro <- FindClusters(umi_prepro, resolution=c(0.05,0.1,0.2,0.4,0.6,0.8))

# UMAP
umi_prepro <- RunUMAP(umi_prepro, dims=1:pc)


reduction <- c("PCA","UMAP")
res <- c(0.05,0.1,0.2,0.4,0.6, 0.8)
pt.s <- 1

# plot clusters
for (red in reduction){
	plot <- NULL
	for (r in res){
		i <- match(r,res)
		res_col <- paste0("RNA_snn_res.",r)
		title <- paste("Clusters res",r)
		plot[[i]] <- DimPlot(umi_prepro, group.by=res_col, reduction=tolower(red), label=T, pt.size=pt.s, label.size = 10) + ggtitle(title) + xlab(paste(red,"1")) + ylab(paste(red,"2")) + NoLegend()
	}
	pdf(file.path(cluster_dir, paste0(sample_id, "_",red,"_clusters.pdf")), width=3*7, height=2*7, useDingbats=FALSE)
	print(plot_grid(plotlist=plot, ncol=3))
	invisible(dev.off())
}

# plot QC
features_to_plot <- c("nFeature_RNA", "nCount_RNA", "percent_mt", "S.Score", "G2M.Score", "doublets_score")

for (red in reduction){
	plot <- NULL
	for (f in features_to_plot) {
		i <- match(f,features_to_plot)
		plot[[i]] <- FeaturePlot(umi_prepro, features=f, reduction=tolower(red), pt.size=pt.s) + xlab(paste(red,"1")) + ylab(paste(red,"2")) + theme(plot.title = element_text(hjust = 0))
	}
	pdf(file.path(cluster_dir, paste0(sample_id, "_",red,"_QC.pdf")), width=3*7, height=2*7, useDingbats=FALSE)
	print(plot_grid(plotlist=plot, ncol=3, nrow=2))
	invisible(dev.off())
}

features_to_plot <- c("MALAT1","NEAT1","GAPDH","nFeature_RNA","nCount_RNA","percent_mt", "doublets_score","S.Score","G2M.Score")
pdf(file.path(cluster_dir, paste0(sample_id, "_Vln_Additional_QC.pdf")), width=3*7, height=3*4, useDingbats=F)
for (r in res){
  plot <- NULL
  res_col <- paste0("RNA_snn_res.", r)
  Idents(umi_prepro) <- umi_prepro[[res_col]]
  
  for (f in features_to_plot) {
    i <- match(f, features_to_plot)
    plot[[i]] <- VlnPlot(umi_prepro, feature=f, pt.size=F) + NoLegend() + theme(plot.title = element_text(hjust = 0))
  }
  
  p <- plot_grid(plotlist=plot, ncol=3)
  title <- ggdraw() + draw_label(paste0("Violin Plots with clusters res.",r), fontface='bold')
  print(plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)))
}
invisible(dev.off())


# plot markers
markers <- c("EPCAM","VIM","PTPRC","CD3D","KRT5","COL1A1","CD14","VWF","JCHAIN","EPCAM","KRT5","KRT18","SLPI","ANKRD30A","ESR1","ERBB2","PRLR")
markers_to_plot <- markers[markers %in% rownames(umi_prepro)]

for (red in reduction){
	plot <- NULL
	for (m in markers_to_plot) {
		i <- match(m,markers_to_plot)
		plot[[i]] <- FeaturePlot(umi_prepro, features=m, reduction=tolower(red), pt.size=pt.s) + xlab(paste(red,"1")) + ylab(paste(red,"2")) + theme(plot.title = element_text(hjust = 0))
	}
	pdf(file.path(cluster_dir, paste0(sample_id, "_",red,"_markers.pdf")), width=3*7, height=6*7, useDingbats=FALSE)
	print(plot_grid(plotlist=plot, ncol=3, nrow=6))
	invisible(dev.off())
}


# Annotation of cell types using SingleR and BlueprintEncode reference dataset for human samples
# For mouse data cell type annotation was performed using scMCA (see methods for details)
BpEn.se <- BlueprintEncodeData()
singler.sc <- SingleR(test=umi_prepro@assays$RNA@data, ref=BpEn.se, labels=BpEn.se$label.main, method="single")
	
umi_prepro$BpEn.sc.main.labels <- singler.sc$labels
	
for (red in reduction){
	print(red)
	p <- DimPlot(umi_prepro, reduction=tolower(red), label=F, pt.size=pt.s, group.by="BpEn.sc.main.labels") + ggtitle("BpEn single cell annotation") + 
		theme(plot.title=element_text(hjust=0.5)) + xlab(paste(red,"1")) + ylab(paste(red,"2"))
	legend <- get_legend(p)
	pdf(file = file.path(anno.dir, paste0(sample_id, "_",red,"_BpEn.sc.main.labels.pdf")), width=2*7, height=7, useDingbats=F)
	print(plot_grid(p + NoLegend(), legend, ncol=2))
	invisible(dev.off())
}


# Custom annotation function
chosen_res <- 0.8 # change this value based on the clustering results and the resolution that best captures the main cell types in your data

customBpEnLabels <- function(umi, add_labels = NULL){
  umi$custom_BpEn <- umi$BpEn.sc.main.labels
  umi$custom_BpEn <- gsub("CD8\\+ |CD4\\+ ","",umi$custom_BpEn)
  umi$custom_BpEn <- gsub("Macrophages|DC|HSC|Monocytes","Myeloid",umi$custom_BpEn)
  umi$custom_BpEn <- gsub("Adipocytes|Chondrocytes|Fibroblasts","Stromal cells",umi$custom_BpEn)
  
  if(is.null(add_labels)){
    umi$custom_BpEn <- ifelse(umi$custom_BpEn %in% c("Epithelial cells","Stromal cells","Endothelial cells","Myeloid","T-cells","B-cells","NK cells"), umi$custom_BpEn, "Other")
  } else {
    umi$custom_BpEn <- ifelse(umi$custom_BpEn %in% c("Epithelial cells","Stromal cells","Endothelial cells","Myeloid","T-cells","B-cells","NK cells", add_labels), umi$custom_BpEn, "Other")
  }
  umi$custom_BpEn <- factor(umi$custom_BpEn, levels=c("Epithelial cells","Stromal cells","Endothelial cells","Myeloid","T-cells","B-cells","NK cells", add_labels, "Other"))
  
  message("\n")
  message("Check if all custom labels are within expected labels")
  print(unique(umi$custom_BpEn) %in% c("Epithelial cells","Stromal cells","Endothelial cells","Myeloid","T-cells","B-cells","NK cells","Other", add_labels))
  message("\n")
  message("Custom labels summary")
  print(summary(umi$custom_BpEn))
  return(umi)
}

umi_prepro <- customBpEnLabels(umi = umi_prepro, add_labels = "Neurons")

cols <- c("red","#596A98","deepskyblue","forestgreen","orchid","orange","blue4","#777777","black")
pt.s <- 0.5
for (red in reduction){
	print(red)
	pdf(file.path(anno.dir, paste0(sample_id, "_",red,"_custom_BpEn-labels.pdf")), width=2*7+2, height=7, useDingbats=F)
	p1 <- DimPlot(umi_prepro, group.by="custom_BpEn", reduction=tolower(red), label=F, pt.size=pt.s) + 
		scale_color_manual(breaks=levels(umi_prepro$custom_BpEn), values=cols) + 
		ggtitle("Custom BpEn sc labels") + xlab(paste(red,"1")) + ylab(paste(red,"2"))
	legend1 <- get_legend(p1)
	p2 <- DimPlot(umi_prepro, group.by=paste0("RNA_snn_res.", chosen_res), reduction=tolower(red), label=T, label.size = 8, pt.size=pt.s, repel = T) + 
	  ggtitle(paste("Cluster res", chosen_res)) + xlab(paste(red,"1")) + ylab(paste(red,"2"))
	
	print(plot_grid(p1 + NoLegend(), legend1, p2 + NoLegend(), ncol=3, rel_widths = c(7, 2, 7)))
	invisible(dev.off())
}


# dotplots
p1 <- DotPlot(
    umi_prepro,
    features = rev(markers_to_plot),
    group.by=paste0("RNA_snn_res.", chosen_res),
    cols = c("white", "red")) +
    RotatedAxis() +
    xlab("marker") + ylab(paste0("Chosen cluster res: ", chosen_res)) +
    ylab(col_group) +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5))

p2 <- DotPlot(
    umi_prepro,
    features = rev(markers_to_plot),
    group.by="custom_BpEn",
    cols = c("white", "red")) +
    RotatedAxis() +
    xlab("marker") + ylab("Custom BpEn labels") +
    ylab(col_group) +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5))

pdf(file.path(cluster_dir, paste0(sample_id, "_DotPlot.pdf")), useDingbats=F)
print(p1)
print(p2)
invisible(dev.off())
### ----------------------------------------------------------------------------------------------------------


### PIPELINE CONTINUES ONLY FOR HUMAN SAMPLES


### 5 - INFERCNV ANALYSIS (USING INFERCNV)
### ----------------------------------------------------------------------------------------------------------
cnv_dir <- file.path(w_dir, "inferCNV/")
if(!dir.exists(cnv_dir)) dir.create(cnv_dir)

epithelial <- c()       # specify the epithelial clusters to be tested against the others (taken as diploid reference)
clusters_remove <- c()  # specify the clusters to be removed from the analysis (e.g. clusters with very low quality cells)

res_col <- umi_prepro@meta.data[, paste0("RNA_snn_res.", chosen_res)]
clusters_keep <- unique(res_col)[! unique(res_col) %in% clusters_remove]
umi_prepro$keep <- res_col
umi_prepro <- subset(umi_prepro, subset = keep %in% clusters_keep)
res_col <- umi_prepro@meta.data[, paste0("RNA_snn_res.", chosen_res)]

reference <- as.character(unique(res_col)[! unique(res_col) %in% c(clusters_remove, epithelial)])

count <- as.matrix(umi_prepro[["RNA"]]@counts)

annot <- umi_prepro@meta.data[, paste0("RNA_snn_res.", resolution), drop = F]
umi_prepro@meta.data[, paste0("RNA_snn_res.", resolution)] <- as.character(annot[, 1])

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=count, annotations_file=annot, gene_order_file=gene_position_file, ref_group_names=reference)
infercnv_obj <- infercnv::run(infercnv_obj, cutoff=0.1, out_dir=cnv_dir, cluster_by_groups=T, denoise=T, HMM=T, num_threads=10)
### ----------------------------------------------------------------------------------------------------------



### 5 - TUMOR SUBSET RECLUSTERING (USING SEURAT)
### ----------------------------------------------------------------------------------------------------------
sub_dir <- file.path(w_dir,"tumor_subset/")
if(!dir.exists(sub_dir)) dir.create(sub_dir)

tumor_clusters <- c() # specify the tumor clusters to be subsetted and reclustered (using chosen_res resolution)

umi_prepro$keep <- ifelse(umi_prepro@meta.data[,old_clusters] %in% tumor_clusters, TRUE, FALSE)
umi_sub <- subset(umi_prepro, subset = keep == TRUE)
umi_sub <- DietSeurat(umi_sub)


### detection of highly variable genes   
Idents(umi_sub) <- "orig.ident"
umi_sub <- FindVariableFeatures(umi_sub, selection.method="vst", nfeatures=2000)
  
# Identify the 10 most highly variable genes
top20 <- head(VariableFeatures(umi_sub), 20)
  
# Cell-Cycle Scoring 
print(sum(cc.genes$s.genes %in% umi_sub[["RNA"]]@var.features))
print(sum(cc.genes$g2m.genes %in% umi_sub[["RNA"]]@var.features)) 

# Visualize the distribution of cell cycle markers 
pdf(file.path(sub_dir, paste0(sample_id, "_plot_CellCycleGene_TumorSubset.pdf")), useDingbats=FALSE)
print(RidgePlot(umi_sub, features=c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol=2))
invisible(dev.off())

# scale data with no regression
all.genes <- rownames(umi_sub)
regr_var <- NULL
umi_sub <- ScaleData(object=umi_sub, features=all.genes, vars.to.regress=regr_var)
umi_sub <- RunPCA(object=umi_sub, features=umi_sub[["RNA"]]@var.features, npcs=30, verbose=F)

pdf(file.path(sub_dir, paste0(sample_id, "_PCElbowPlot_TumorSubset.pdf")), width=15, height=15, useDingbats=FALSE)
print(ElbowPlot(object=umi_sub, ndims=30))
invisible(dev.off())

# find clusters and UMAP
pc <- 15 # change this value based on the ElbowPlot - tipically between 10 and 20
umi_sub <- FindNeighbors(umi_sub, reduction="pca", dims=1:pc, force.recalc=T)
umi_sub <- FindClusters(umi_sub, resolution=res)

umi_sub <- RunUMAP(umi_sub, dims=1:pc)

# plot QC
features_to_plot <- c("nFeature_RNA", "nCount_RNA", "percent_mt", "S.Score", "G2M.Score", "doublets_score")

for (red in reduction){
	plot <- NULL
	for (f in features_to_plot) {
		i <- match(f,features_to_plot)
		plot[[i]] <- FeaturePlot(umi_sub, features=f, reduction=tolower(red), pt.size=pt.s) + xlab(paste(red,"1")) + ylab(paste(red,"2")) + theme(plot.title = element_text(hjust = 0))
	}
	pdf(file.path(sub_dir, paste0(sample_id, "_",red,"_QC_TumorSubset.pdf")), width=3*7, height=2*7, useDingbats=FALSE)
	print(plot_grid(plotlist=plot, ncol=3, nrow=2))
	invisible(dev.off())
}

features_to_plot <- c("MALAT1","NEAT1","GAPDH", "PTPRC", "VIM", "PECAM1","nFeature_RNA","nCount_RNA","percent_mt", "doublets_score","S.Score","G2M.Score")
pdf(file.path(sub_dir, paste0(sample_id, "_Vln_Additional_QC_TumorSubset.pdf")), width=3*7, height=4*4, useDingbats=F)
for (r in res){
  plot <- NULL
  res_col <- paste0("RNA_snn_res.", r)
  Idents(umi_sub) <- umi_sub[[res_col]]
  
  for (f in features_to_plot) {
    i <- match(f, features_to_plot)
    plot[[i]] <- VlnPlot(umi_sub, feature=f, pt.size=F) + NoLegend() + theme(plot.title = element_text(hjust = 0))
  }
  
  p <- plot_grid(plotlist=plot, ncol=3)
  title <- ggdraw() + draw_label(paste0("Violin Plots with clusters res.",r), fontface='bold')
  print(plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)))
}
invisible(dev.off())


# plot markers
markers <- c("EPCAM","VIM","PTPRC","CD3D","KRT5","COL1A1","CD14","VWF","JCHAIN","EPCAM","KRT5","KRT18","SLPI","ANKRD30A","ESR1","ERBB2","PRLR")
markers_to_plot <- markers[markers %in% rownames(umi_sub)]

for (red in reduction){
	plot <- NULL
	for (m in markers_to_plot) {
		i <- match(m,markers_to_plot)
		plot[[i]] <- FeaturePlot(umi_sub, features=m, reduction=tolower(red), pt.size=pt.s) + xlab(paste(red,"1")) + ylab(paste(red,"2")) + theme(plot.title = element_text(hjust = 0))
	}
	pdf(file.path(sub_dir, paste0(sample_id, "_",red,"_markers_TumorSubset.pdf")), width=3*7, height=6*7, useDingbats=FALSE)
	print(plot_grid(plotlist=plot, ncol=3, nrow=6))
	invisible(dev.off())
}
### ----------------------------------------------------------------------------------------------------------
