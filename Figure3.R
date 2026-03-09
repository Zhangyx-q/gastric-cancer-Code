# ==============================================================================
# Project: Gastric Cancer Myeloid Analysis (Figure 3)
# Libraries
# ==============================================================================
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(reshape2)
library(UpSetR)
library(nichenetr)
library(ggpubr)
library(cowplot)
library(DESeq2)
library(SCENIC)

rm(list = ls())

# Global Aesthetics
Mcolors <- c("#F5AB5E", "#AFC778", "#947959", "#EFDBB9", "#5F97C6", "#A3CDEA")
site_colors <- c("#007891", "#129981", "#3A5182", "#EA977C", "#DA4B35")

# Load Data
Myeloid <- readRDS("./Myeloid.rds")

# ==============================================================================
# Figure 3A: DimPlot
# ==============================================================================
pdf("./Fig3A.pdf", 6.3, 3.7)
DimPlot(Myeloid, cols = Mcolors, label = TRUE) + theme_dr()
dev.off()

# ==============================================================================
# Figure 3B: Tissue Distribution & ROE Heatmap
# ==============================================================================
# 1. Calculate ROE (Observed/Expected)
mapz <- Myeloid@meta.data %>%
  group_by(cluster, Site) %>%
  summarise(n = n(), .groups = "drop")

wide_df <- dcast(mapz, cluster ~ Site, value.var = "n", fill = 0)
rownames(wide_df) <- wide_df[, 1]
wide_df <- wide_df[, -1]

row_totals <- rowSums(wide_df)
col_totals <- colSums(wide_df)
grand_total <- sum(wide_df)
expected <- outer(row_totals, col_totals) / grand_total

roe <- round(as.matrix(wide_df) / expected, 2)
roe <- roe[, c("LN", "MLN", "MMLN", "MT", "PT")] # Ensure column order
roe[roe > 3.5] <- 3.5

# Significance Matrix
symbol_matrix <- matrix("", nrow = nrow(roe), ncol = ncol(roe))
symbol_matrix[roe > 3] <- "+++"
symbol_matrix[roe > 2 & roe <= 3] <- "++"
symbol_matrix[roe > 1 & roe <= 2] <- "+"

pdf("./Fig3B_Heatmap.pdf", 3.5, 2.5)
pheatmap::pheatmap(
  roe,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = symbol_matrix,
  color = colorRampPalette(c("#f7f7f7", "#fad9c7", "#e68661", "#a91e2b"))(100)
)
dev.off()

# 2. Composition Barplots
mapz$Site <- factor(mapz$Site, levels = c("LN", "MLN", "MMLN", "PT", "MT"))

p1 <- ggplot(mapz, aes(x = cluster, y = n, fill = Site)) +
  geom_bar(stat = 'identity', position = 'fill', color = "grey", size = 0.2) +
  scale_fill_manual(values = site_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p2 <- ggplot(mapz, aes(x = Site, y = n, fill = cluster)) +
  geom_bar(stat = 'identity', position = 'fill', color = "grey", size = 0.2) +
  scale_fill_manual(values = Mcolors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pdf("./Fig3B_Barplot.pdf", 7, 4)
print(p1 | p2)
dev.off()

# ==============================================================================
# Figure 3C: Module Scores (M1, M2, DC, and Hallmarks)
# ==============================================================================
# Load Markers
marker <- read.table("./markers_33545035.txt", sep = "\t", header = TRUE)
marker_list <- apply(marker, 2, function(x) setdiff(as.character(x), ""))
Myeloid <- AddModuleScore(Myeloid, marker_list, name = names(marker_list))

# Hallmark Scoring
genez <- read.gmt("./h.all.v2024.1.Hs.symbols.gmt")
genelist_hallmark <- split(genez$gene, genez$term)
Myeloid <- AddModuleScore(Myeloid, genelist_hallmark, name = names(genelist_hallmark))

# Selection of specific pathways for visualization
target_pathways <- c(
  "HALLMARK_ANGIOGENESIS4", "HALLMARK_APOPTOSIS7", "HALLMARK_COMPLEMENT11",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING24", "HALLMARK_INFLAMMATORY_RESPONSE25",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB45", "HALLMARK_XENOBIOTIC_METABOLISM50"
)

# Helper function for consistent boxplots
plot_module <- function(obj, feature, title_str) {
  ggboxplot(obj@meta.data, x = "cluster", y = feature, width = 0.4,
            fill = "cluster", palette = Mcolors, outlier.shape = NA,
            title = title_str, legend = "none") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = NULL, y = NULL)
}

p1 <- plot_module(Myeloid, "M11", "M1 Score")
p2 <- plot_module(Myeloid, "M22", "M2 Score") + 
      geom_hline(yintercept = 0.2, linetype = "dashed")

# DC Markers
dc_marker <- read.table("./dc_markers.txt", header = TRUE, sep = "\t")
dc_list <- apply(dc_marker, 2, function(x) setdiff(as.character(x), ""))
Myeloid <- AddModuleScore(Myeloid, dc_list, name = names(dc_list))

p9 <- plot_module(Myeloid, "ActivatedDC1", "Activated DC")
p_hallmark <- plot_module(Myeloid, target_pathways[5], "Inflammatory Response")

pdf("./Fig3C.pdf", 8.8, 7.8)
plot_grid(p1, p2, p9, p_hallmark, ncol = 2)
dev.off()

# ==============================================================================
# Figure 3D: Time-course Analysis (TAM_C1QA)
# ==============================================================================
mac_sub <- subset(Myeloid, cluster == "M03_TAM_C1QA" & Site %in% c("LN", "MLN", "MMLN"))

# Aggregate for Pseudo-bulk
pb <- AggregateExpression(mac_sub, group.by = c("orig.ident", "Site"), 
                          assays = "RNA", slot = "counts")
counts <- pb$RNA

# Metadata for DESeq2
meta_pb <- data.frame(
  sample = colnames(counts),
  site = sapply(strsplit(colnames(counts), "_"), `[`, 2)
)
meta_pb$time_num <- case_when(
  meta_pb$site == "LN" ~ 0,
  meta_pb$site == "MLN" ~ 1,
  meta_pb$site == "MMLN" ~ 2
)
rownames(meta_pb) <- meta_pb$sample

# DESeq2 LRT Test
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta_pb, design = ~ time_num)
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)
res_lrt <- results(dds_lrt)
sig_genes <- rownames(res_lrt[which(res_lrt$padj < 0.05), ])

# Heatmap Preparation
vsd <- varianceStabilizingTransformation(dds_lrt, blind = FALSE)
avg_exp <- t(apply(assay(vsd)[sig_genes, ], 1, function(x) tapply(x, meta_pb$time_num, mean)))
scaled_exp <- t(scale(t(avg_exp)))

selected_genes <- c("IL1B", "STAT3", "CCL24", "CD33", "ATF3", "KLF4", "PFKFB3", "RUNX1", "RELA")

row_anno <- rowAnnotation(
  mark = anno_mark(at = which(rownames(scaled_exp) %in% selected_genes), 
                   labels = rownames(scaled_exp)[rownames(scaled_exp) %in% selected_genes])
)

pdf("./Fig3D.pdf", 7, 7)
Heatmap(scaled_exp, name = "Z-score", col = colorRamp2(c(-2, 0, 2), c("#007891", "white", "#b2182b")),
        show_row_names = FALSE, cluster_columns = FALSE, row_km = 4, right_annotation = row_anno)
dev.off()

# ==============================================================================
# Figure 3G & 3H: SCENIC TF Analysis
# ==============================================================================
# [Note: Assuming SCENIC loom and RSS calculations are pre-computed]
# rss_df from calcRSS...
# genes1 <- c(...) # Defined TF list

# Upset Plot (3G)
x_upset <- list(LN = genes[1:50], MLN = genes[51:100], MMLN = genes[101:150])
pdf("./Fig3G.pdf", 3.2, 2.8)
upset(fromList(x_upset), order.by = "freq")
dev.off()

# TF Heatmap (3H)
pdf("./Fig3H.pdf", 2, 4)
pheatmap::pheatmap(matrix_tf, cluster_rows = FALSE, cluster_cols = FALSE,
                   color = colorRampPalette(c("#ece5dd", "#f0b849", "#dd7500", "#af0606"))(20))
dev.off()

# ==============================================================================
# Figure 3I: NicheNet (Ligand-Target)
# ==============================================================================
# Filter clusters with low cell counts
Idents(LN_MLN) <- "cluster"
counts_table <- table(Idents(LN_MLN), LN_MLN$Site)
keep_clusters <- rownames(counts_table)[apply(counts_table, 1, min) >= 20]
LN_MLN_sub <- subset(LN_MLN, idents = keep_clusters)

nichenet_output <- nichenet_seuratobj_aggregate(
  seurat_obj = LN_MLN_sub, 
  sender = "all",
  receiver = "M03_TAM_C1QA", 
  condition_colname = "Site",
  condition_oi = "MLN",
  condition_reference = "LN",
  geneset = "up",
  expression_pct = 0.05,
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks
)

pdf("./Fig3I.pdf", 26, 7.3)
print(nichenet_output$ligand_activity_target_heatmap)
dev.off()
