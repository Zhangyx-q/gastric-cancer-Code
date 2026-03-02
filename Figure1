library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(ggplot2)
library(pheatmap)

source("./color.R")

GCdata <- readRDS("./GCdata.rds")

# ===============================
# Fig 1C – UMAP
# ===============================

pdf("./Fig1C.pdf", width = 9, height = 4.1)
DimPlot(GCdata, group.by = "cluster", cols = colorA, raster = FALSE)
dev.off()


# ===============================
# Fig 1D – Marker Heatmap
# ===============================

cluster_order <- c(
  "B01_naive_IGHD","B02_memery_CD27","B03_memery_TNFRSF13B",
  "B04_HSP_CCR7","B05_Plasma_JCHAIN","B06_GCB_LMO2","B07_Cycling_MKI67",
  "T01_CD4_LEF1","T02_CD4_SELL","T03_CD8_GZMB","T04_CD8_GZMK",
  "T05_CD4_CTLA4","T06_CD8_IL7R","T07_CD4_FOXP3","T08_CD4_ANXA1",
  "T09_gdT_GNLY","T10_CD8_CCL5","T12_CD4_FOS","T11_ILC3_AREG",
  "M01_cDC3_CCL19","M02_Mono_CXCL8","M03_TAM_C1QA",
  "M04_Mast_TPSAB1","M05_cDC2_CLEC10A","M06_cDC1_CLEC9A",
  "S1_Pericyte","S2_Fib_CXCL14","S3_Fib_SFRP2","S5_myCAF_COL1A1",
  "E1_TAEC_ESM1","E2_lymphatic EC_TNFRSF10C","E3_Capillary_RBP7",
  "E4_Vein_ACKR1","E5_Artery_NEAT1","E6_Proliferationg ECs"
)

gene_cell_exp <- AverageExpression(
  GCdata,
  features = top50_genes,
  group.by = "cluster",
  slot = "data"
)$RNA

mat <- t(scale(t(gene_cell_exp)))
mat <- MinMax(mat, max = 2.5, min = -2.5)
plotData <- mat[, intersect(cluster_order, colnames(mat))]

cluster.color <- c(Bcolors, Tcolors, Mcolors, Scolors, Ecolors)
names(cluster.color) <- cluster_order

col_fun <- colorRamp2(
  seq(-2.5, 2.5, length = 100),
  colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
)

ha <- HeatmapAnnotation(
  Cluster = colnames(plotData),
  col = list(Cluster = cluster.color),
  show_annotation_name = FALSE
)

add_dashed_lines <- function(ht, mat, rows = NULL, cols = NULL,
                             lty = 2, col = "black", lwd = 1.5) {
  
  ht_drawn <- draw(ht)
  ht_name <- names(ht_drawn@ht_list)[1]
  
  decorate_heatmap_body(ht_name, {
    nr <- nrow(mat)
    nc <- ncol(mat)
    
    if (!is.null(cols)) {
      for (c in cols) {
        grid.lines(x = unit(c / nc, "npc"),
                   y = unit(c(0, 1), "npc"),
                   gp = gpar(lty = lty, col = col, lwd = lwd))
      }
    }
    
    if (!is.null(rows)) {
      for (r in rows) {
        grid.lines(x = unit(c(0, 1), "npc"),
                   y = unit(c(1 - r / nr, 1 - r / nr), "npc"),
                   gp = gpar(lty = lty, col = col, lwd = lwd))
      }
    }
  })
}

pdf("./Fig1D.pdf", width = 8.5, height = 8)

ht <- Heatmap(
  plotData,
  name = "z-score",
  top_annotation = ha,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = col_fun,
  border = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE
)

add_dashed_lines(ht, plotData,
                 rows = c(336,674,893,1067),
                 cols = c(7,19,25,29))

dev.off()


# ===============================
# Fig 1E – UMAP split by site
# ===============================

pdf("./Fig1E.pdf", width = 11.7, height = 5)

DimPlot(
  GCdata,
  reduction = "umap",
  group.by = "cluster",
  cols = colorA,
  split.by = "Site",
  label = TRUE,
  raster = FALSE
)

dev.off()


# ===============================
# Fig 1F – Sample composition
# ===============================

sample_order <- c(
  "1NLN","5NLN","12NLN","13NLN","16NLN",
  "1NT","5NT","12NT","13NT","16NT",
  "2MLN","6MLN","7MLN","8MLN","9MLN","10MLN","11MLN","14MLN",
  "2MMLN","6MMLN","9MMLN","11MMLN",
  "2MT","6MT","7MT","8MT","9MT","10MT","11MT","14MT"
)

mapz <- GCdata@meta.data %>%
  group_by(orig.ident, celltypeMN, ss) %>%
  summarise(n = n(), .groups = "drop")

pdf("./sample_celltype.pdf", width = 7, height = 3)

p1 <- ggplot(mapz,
             aes(x = factor(orig.ident, levels = sample_order),
                 y = n,
                 fill = factor(celltypeMN))) +
  geom_bar(stat = "identity", position = "fill",
           color = "grey", size = 0.2) +
  scale_fill_manual(values = col) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(. ~ ss, scales = "free", space = "free")

print(p1)
dev.off()


# ===============================
# Fig 1G – Ro/e enrichment heatmap
# ===============================

GCdata$celltypeMN <- gsub("CD8\\+ T", "Tcell", GCdata$celltypeMN)
GCdata$celltypeMN <- gsub("CD4\\+ T", "Tcell", GCdata$celltypeMN)

mapz <- GCdata@meta.data %>%
  group_by(celltypeMN, Site) %>%
  summarise(n = n(), .groups = "drop")

wide_df <- dcast(mapz, celltypeMN ~ Site)
wide_df[is.na(wide_df)] <- 0

rownames(wide_df) <- wide_df[,1]
wide_df <- wide_df[,-1]

row_totals <- rowSums(wide_df)
col_totals <- colSums(wide_df)
grand_total <- sum(wide_df)

expected <- outer(row_totals, col_totals) / grand_total
roe <- as.matrix(wide_df) / expected
roe <- round(roe, 2)

Idents(GCdata) <- GCdata$celltypeMN
levels(GCdata) <- c("B","plasma","Tcell","ILC3","Myeloid",
                    "DC","Mast","endothelial",
                    "fibroblast","Epithelial")

roe <- roe[levels(GCdata),]
roe <- roe[,c(1,2,3,5,4)]

symbol_matrix <- matrix("", nrow = nrow(roe), ncol = ncol(roe))
symbol_matrix[roe > 3] <- "+++"
symbol_matrix[roe > 2 & roe <= 3] <- "++"
symbol_matrix[roe > 1 & roe <= 2] <- "+"

rownames(symbol_matrix) <- rownames(roe)
colnames(symbol_matrix) <- colnames(roe)

roe[roe > 3.5] <- 3.5

pdf("./Fig1G.pdf", width = 3.5, height = 4)

pheatmap(
  roe,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = symbol_matrix,
  color = colorRampPalette(
    c("#f7f7f7","#fad9c7","#e68661","#a91e2b")
  )(100)
)

dev.off()
