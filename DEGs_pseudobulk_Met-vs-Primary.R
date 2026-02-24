#R version 3.6.3

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(Matrix)
library(Matrix.utils)
library(edgeR)
library(pheatmap)
library(png)
library(RColorBrewer)
library(data.table)


# get data
###-----------------------------------------------------------------------------
sample_list <- readRDS("sample_list.rds")
# 'sample list' is an rds object containing a list of 24 Seurat objects, each corresponding to a sample.
# Each Seurat object contains only the tumor cells from that sample, and has been preprocessed using the example scRNA-seq worflow
# Each sample is annotated with a "custom_label" metadata column indicating whether the sample is from a primary tumor or a metastasis.

str(sample_list, max.level = 1)
#List of 24
#$ BBM0      :Formal class 'Seurat' [package "Seurat"] with 12 slots
#$ BSM1      :Formal class 'Seurat' [package "Seurat"] with 12 slots
#$ BBM5      :Formal class 'Seurat' [package "Seurat"] with 12 slots
#$ BLM285    :Formal class 'Seurat' [package "Seurat"] with 12 slots
#$ BLM963    :Formal class 'Seurat' [package "Seurat"] with 12 slots
#$ BBM8      :Formal class 'Seurat' [package "Seurat"] with 12 slots
#$ BBM7      :Formal class 'Seurat' [package "Seurat"] with 12 slots
#$ BBM10     :Formal class 'Seurat' [package "Seurat"] with 12 slots
#$ BBM9      :Formal class 'Seurat' [package "Seurat"] with 12 slots
#$ BBM13     :Formal class 'Seurat' [package "Seurat"] with 12 slots
#$ #9290   :Formal class 'Seurat' [package "Seurat"] with 12 slots
#$ #9294   :Formal class 'Seurat' [package "Seurat"] with 12 slots
#$ #9289   :Formal class 'Seurat' [package "Seurat"] with 12 slots
#$ #9293   :Formal class 'Seurat' [package "Seurat"] with 12 slots
#$ #9291   :Formal class 'Seurat' [package "Seurat"] with 12 slots
#$ #9299   :Formal class 'Seurat' [package "Seurat"] with 12 slots
#$ #9306   :Formal class 'Seurat' [package "Seurat"] with 12 slots
#$ #9315   :Formal class 'Seurat' [package "Seurat"] with 12 slots
#$ #9317   :Formal class 'Seurat' [package "Seurat"] with 12 slots
#$ BC4       :Formal class 'Seurat' [package "Seurat"] with 12 slots
#$ BC3       :Formal class 'Seurat' [package "Seurat"] with 12 slots
#$ BC16      :Formal class 'Seurat' [package "Seurat"] with 12 slots
#$ #9282   :Formal class 'Seurat' [package "Seurat"] with 12 slots
#$ #9281   :Formal class 'Seurat' [package "Seurat"] with 12 slots


merged <- sample_list[[1]]
merged <- RenameCells(merged, new.names = paste0(names(sample_list)[1],"__", WhichCells(merged)))

for(i in 2:length(sample_list)){
  message("Merging ", names(sample_list)[i])
  to_be_merged <- sample_list[[i]]
  to_be_merged <- RenameCells(to_be_merged, new.names = paste0(names(sample_list)[i],"__", WhichCells(to_be_merged)))
  merged <- merge(merged, to_be_merged)
}

DefaultAssay(merged) <- "RNA"
###-----------------------------------------------------------------------------



# aggregate pseudobulks
###-----------------------------------------------------------------------------
merged$group_id <- merged$custom_label
merged$sample_id <- paste(merged$group_id, gsub("_|^\\d*_", "", merged$orig.ident), sep = "-")
merged$cluster_id <- "Tumor"
table(merged$sample_id, merged$group_id)
#                   Metastasis Primary
#Metastasis-BBM5          4773       0
#Metastasis-BBM10         4726       0
#Metastasis-BBM13         2931       0
#Metastasis-BBM7          1149       0
#Metastasis-BBM8           675       0
#Metastasis-BBM9           487       0
#Metastasis-BLM285        2560       0
#Metastasis-BLM963         243       0
#Metastasis-BBM0          2485       0
#Metastasis-BSM1          3756       0
#Primary-BC16                0    1171
#Primary-BC3                 0     364
#Primary-BC4                 0    2831
#Primary-GSM4909281          0    1152
#Primary-GSM4909282          0   13925
#Primary-GSM4909289          0    2625
#Primary-GSM4909290          0    3020
#Primary-GSM4909291          0    1315
#Primary-GSM4909293          0    3372
#Primary-GSM4909294          0    2179
#Primary-GSM4909299          0    1674
#Primary-GSM4909306          0    6305
#Primary-GSM4909315          0    3353
#Primary-GSM4909317          0    2043

counts <- GetAssayData(object = merged, slot = "counts", assay="RNA")
metadata <- merged@meta.data

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- metadata[, "sample_id"]
groups <- factor(groups)

# Aggregate across cluster-sample groups
# Each row corresponds to aggregate counts for a cluster-sample combo
aggr_counts <- aggregate.Matrix(t(counts), groupings = groups, fun = "sum") 
aggr_counts <- t(aggr_counts)

# Turn class "table" into a named vector of cells per sample
n_cells <- table(metadata$sample_id) %>%  as.vector()
names(n_cells) <- names(table(metadata$sample_id))
m <- match(names(n_cells), metadata$sample_id)

# Create the sample level metadata by selecting specific columns
aggr_metadata <- data.frame(metadata[m, ], 
                 n_cells, row.names = NULL) %>% 
  dplyr::select("sample_id", "group_id", "n_cells")
aggr_metadata$group_id <- factor(aggr_metadata$group_id, levels = c("Primary", "Metastasis"))
###-----------------------------------------------------------------------------


design <- model.matrix(~ 0 + aggr_metadata$group_id)

#filter low-expressed genes using filterByExpr by edgeR
keep <- filterByExpr(aggr_counts, design, min.prop = 0.25)
sum(keep)
#[1] 16680
aggr_counts <- aggr_counts[keep, ]



### EdgeR
###-----------------------------------------------------------------------------
# Double-check that both lists have same names
all(names(aggr_counts) == names(aggr_metadata))

rownames(design) <- aggr_metadata$sample_id
colnames(design) <- levels(factor(aggr_metadata$group_id))

# A positive FC is increased expression in the ko compared to het
(contrast <- makeContrasts("Metastasis-Primary", levels = design))

#edgeR test
y <- DGEList(aggr_counts, remove.zeros = TRUE)
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fitlrt <- glmFit(y, design)
res <- glmLRT(fitlrt, coef=2, contrast = contrast)
res <- topTags(res, n=Inf, adjust.method = "BH")


res_tbl <- res %>%
  data.frame() %>%
  mutate("gene" = rownames(res)) %>%
  as_tibble() %>%
  arrange(FDR)

write.table(res_tbl, file = "edgeR-Metastasis_vs_Primary.txt", row.names = FALSE, sep = "\t")


### Differential genes
# Set thresholds
padj_cutoff <- 0.05

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, FDR < padj_cutoff) %>%
  dplyr::arrange(FDR)

sig_res_up <- dplyr::filter(res_tbl, FDR < padj_cutoff) %>%
  dplyr::filter(logFC > 0) %>%
  dplyr::arrange(FDR, logFC)
sig_res_up

sig_res_down <- dplyr::filter(res_tbl, FDR < padj_cutoff) %>%
  dplyr::filter(logFC < 0) %>%
  dplyr::arrange(FDR, desc(logFC))
sig_res_down

# Check significant genes output
sig_res
write.table(sig_res, file = "edgeR-Metastasis_vs_Primary-SIG.txt", row.names = FALSE)

edgeR_sig <- list("up" = sig_res_up, "down" = sig_res_down)


# Extract normalized counts from dds object
normalized_counts <- cpm(y)
edgeR_normalized_counts <- normalized_counts


# TOP up genes
top20_sig_genes <- sig_res %>%
  dplyr::filter(logFC >= log2fc_cutoff) %>% 
  dplyr::arrange(FDR) %>%
  dplyr::pull(gene)

top20_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top20_sig_genes, ]
top20_sig_counts

top20_sig_df <- data.frame(top20_sig_counts)
top20_sig_df$gene <- rownames(top20_sig_counts)

top20_sig_df <- melt(setDT(top20_sig_df), 
                     id.vars = c("gene"),
                     variable.name = "sample_id") %>% 
  data.frame()

top20_sig_df$sample_id <- gsub("\\.", "-", top20_sig_df$sample_id)

top20_sig_df <- plyr::join(top20_sig_df, aggr_metadata,
                           by = "sample_id")

top20_sd <- top20_sig_df %>% group_by(gene) %>% summarize("sd" = sd(value))
top20_mean <- top20_sig_df %>% group_by(gene) %>% summarize("mean" = mean(value))
top20_sig_df <- top20_sig_df %>% left_join(top20_sd, by = "gene") %>% left_join(top20_mean, by = "gene") %>%
  mutate("label" = ifelse(value > mean+2.5*sd | value < mean-2.5*sd, TRUE, FALSE))

p_up <- ggplot(top20_sig_df, aes(y = value, x = group_id, col = group_id)) +
  geom_violin(alpha = 0.7, scale = "width", color = "black", aes(fill = group_id)) +
  geom_jitter(height = 0, width = 0.1, color = "black", size = 0.25) +
  geom_text_repel(aes(label=sample_id, alpha=label), size = 3, color = "black") + scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0)) +
  scale_y_continuous(trans = 'log') +
  scale_fill_manual(values = c("Metastasis" = "tomato", "Primary" = "steelblue")) +
  ylab("log of normalized expression level") +
  xlab("condition") +
  ggtitle("Significant up DE Genes") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ gene, ncol = 4)


# TOP down genes
top20_sig_genes <- sig_res %>%
  dplyr::filter(logFC < -log2fc_cutoff) %>% 
  dplyr::arrange(FDR) %>%
  dplyr::pull(gene)

top20_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top20_sig_genes, ]

top20_sig_df <- data.frame(top20_sig_counts)
top20_sig_df$gene <- rownames(top20_sig_counts)

top20_sig_df <- melt(setDT(top20_sig_df), 
                     id.vars = c("gene"),
                     variable.name = "sample_id") %>% 
  data.frame()

top20_sig_df$sample_id <- gsub("\\.", "-", top20_sig_df$sample_id)

top20_sig_df <- plyr::join(top20_sig_df, aggr_metadata,
                           by = "sample_id")

top20_sd <- top20_sig_df %>% group_by(gene) %>% summarize("sd" = sd(value))
top20_mean <- top20_sig_df %>% group_by(gene) %>% summarize("mean" = mean(value))
top20_sig_df <- top20_sig_df %>% left_join(top20_sd, by = "gene") %>% left_join(top20_mean, by = "gene") %>%
  mutate("label" = ifelse(value > mean+2.5*sd | value < mean-2.5*sd, TRUE, FALSE))

p_down <- ggplot(top20_sig_df, aes(y = value, x = group_id, col = group_id)) +
  geom_violin(alpha = 0.7, scale = "width", color = "black", aes(fill = group_id)) +
  geom_jitter(height = 0, width = 0.1, color = "black", size = 0.25) +
  geom_text_repel(aes(label=sample_id, alpha=label), size = 3, color = "black") + scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0)) +
  scale_y_continuous(trans = 'log') +
  scale_fill_manual(values = c("Metastasis" = "tomato", "Primary" = "steelblue")) +
  ylab("log of normalized expression level") +
  xlab("condition") +
  ggtitle("Significant down DE Genes") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ gene, ncol = 4)

pdf("edgeR-Metastasis_vs_Primary_sig_genes_expression.pdf", width = 12, height = 30, useDingbats = FALSE)
print(p_up)
print(p_down)
dev.off()


# Heatmap of all significant genes
sig_counts <- normalized_counts[rownames(normalized_counts) %in% sig_res$gene, ]
heat_colors <- rev(brewer.pal(11, "PRGn"))

pdf("edgeR-Metastasis_vs_Primary_heatmap.pdf", width = 8, height = 12, useDingbats = FALSE)
pheatmap(sig_counts, 
         color = heat_colors, 
         cluster_rows = TRUE, 
         show_rownames = FALSE,
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)  
dev.off()

logfc_cutoff <- 0

# Volcano plot
res_table_thres <- res_tbl[!is.na(res_tbl$FDR), ] %>% 
  mutate(threshold = case_when(FDR < padj_cutoff & logFC <= -logfc_cutoff ~ "Down",
                               FDR < padj_cutoff & logFC >= logfc_cutoff ~ "Up",
                               TRUE ~ "Not significant"))
res_table_thres_up <- res_table_thres %>%
  filter(threshold == "Up") %>%
  arrange(FDR) %>%
  head(20)
res_table_thres_down <- res_table_thres %>%
  filter(threshold == "Down") %>%
  arrange(FDR) %>%
  head(20)

pdf("edgeR-Metastasis_vs_Primary_volcano.pdf", width = 12, height = 8, useDingbats = FALSE)
ggplot(res_table_thres) +
  geom_point(aes(x = logFC, y = -log10(FDR), colour = threshold), size = 0.2) +
  geom_text_repel(data = mutate(edgeR_sig$up, threshold = "Up"), aes(label = gene, x = logFC, y = -log10(FDR), colour = threshold)) +
  geom_text_repel(data = mutate(edgeR_sig$down, threshold = "Down"), aes(label = gene, x = logFC, y = -log10(FDR), colour = threshold)) +
  scale_color_manual(values = c("Not significant" = "grey60", "Up" = "red3", "Down" = "blue")) +
  theme_minimal_grid() +
  ggtitle("Volcano plot of Metastasis vs Primary cells") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.3), hjust = 0.5),
        axis.title = element_text(size = rel(1.15)))  
dev.off()
###-----------------------------------------------------------------------------