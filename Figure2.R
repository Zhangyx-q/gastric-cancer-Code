rm(list = ls())

library(Seurat)
library(ggplot2)
library(ggpubr)
library(scales)
library(dplyr)
library(reshape2)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(cowplot)
library(Startrac)

# ===============================
# Load data
# ===============================

Tcell <- readRDS("./Tcell.rds")
TCR   <- readRDS("./TCR.rds")

# ===============================
# Fig 2A – Clonal UMAP
# ===============================

Tcell$clonal <- "zNan"

clonal_cells <- TCR$Cell_Name[TCR$clone.status == "Clonal"]
clonal_cells <- intersect(colnames(Tcell), clonal_cells)

Tcell@meta.data[clonal_cells, "clonal"] <- 
  Tcell@meta.data[clonal_cells, "cluster"]

Tcolors <- rev(c(
  "#F49568","#ED7A6A","#77DCDD","#A5D1B0","#CE8A8D","#FBC99A",
  "#FFF7C1","#ADD3F4","#F7C9CF","#FEE4E8",
  "#7CA3B8","#BFB8D6","#bdc1c8"
))

Idents(Tcell) <- Tcell$clonal
levels(Tcell) <- sort(levels(Tcell))

p_cluster <- DimPlot(Tcell, group.by = "cluster",
                     cols = rev(Tcolors)[-13])

p_clonal  <- DimPlot(Tcell, group.by = "clonal",
                     cols = Tcolors,
                     order = levels(Tcell)[-13])

pdf("./Fig2A.pdf", 5.7, 4)
p_cluster | p_clonal
dev.off()


# ===============================
# Fig 2B – Ro/e heatmap
# ===============================

mapz <- Tcell@meta.data %>%
  group_by(cluster, Site) %>%
  summarise(n = n(), .groups = "drop")

wide_df <- dcast(mapz, cluster ~ Site)
wide_df[is.na(wide_df)] <- 0

rownames(wide_df) <- wide_df[,1]
wide_df <- wide_df[,-1]
wide_df <- wide_df[c(1:10,12,11),]

row_totals <- rowSums(wide_df)
col_totals <- colSums(wide_df)
grand_total <- sum(wide_df)

expected <- outer(row_totals, col_totals) / grand_total
roe <- round(as.matrix(wide_df) / expected, 2)
roe <- roe[,c(1,2,3,5,4)]
roe[roe > 3.5] <- 3.5

symbol_matrix <- matrix("", nrow = nrow(roe), ncol = ncol(roe))
symbol_matrix[roe > 3] <- "+++"
symbol_matrix[roe > 2 & roe <= 3] <- "++"
symbol_matrix[roe > 1 & roe <= 2] <- "+"

rownames(symbol_matrix) <- rownames(roe)
colnames(symbol_matrix) <- colnames(roe)

pdf("./Fig2B.pdf", 3.5, 2.5)
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


# ===============================
# Fig 2C – Cluster proportions
# ===============================

prop_df <- Tcell@meta.data %>%
  group_by(orig.ident, cluster, Site) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(prop = cell_count / sum(cell_count)) %>%
  ungroup()

plot_list <- list()

for (cl in unique(prop_df$cluster)) {
  
  df_sub <- prop_df[prop_df$cluster == cl,]
  
  p <- ggplot(df_sub,
              aes(x = factor(Site,
                             levels=c("LN","MLN","MMLN","PT","MT")),
                  y = prop,
                  fill = Site)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.7) +
    scale_fill_manual(values = c(
      "#007891","#129981","#3A5182",
      "#DA4B35","#EA977C"
    )) +
    theme_classic() +
    scale_y_continuous(labels = percent_format()) +
    labs(title = cl, x = "Group", y = "Cell Proportion") +
    stat_compare_means(method = "wilcox.test")
  
  plot_list[[cl]] <- p
}

pdf("./Fig2C.pdf", 21, 11)
plot_grid(plot_list[[3]], plot_list[[4]],
          plot_list[[6]], plot_list[[10]],
          ncol = 2)
dev.off()


# ===============================
# Fig 2D – STARTRAC indices
# ===============================

TCR$patient <- paste0(TCR$patient,"_",TCR$loc)
out <- Startrac.run(TCR, proj="GC", verbose=FALSE)

cluster_data <- out@cluster.sig.data
cluster_data <- cluster_data[cluster_data$aid != "GC",]

Tcolors <- c(
  "#F49568","#ED7A6A","#77DCDD","#A5D1B0",
  "#CE8A8D","#FBC99A","#FFF7C1",
  "#ADD3F4","#F7C9CF","#FEE4E8",
  "#7CA3B8","#BFB8D6"
)

pdf("./Fig2D.pdf",3.5,3.6)
ggplot(cluster_data,
       aes(x = majorCluster, y = value,
           color = majorCluster)) +
  geom_boxplot() +
  scale_color_manual(values = Tcolors) +
  facet_wrap(~ index, scales = "free_y",
             ncol = 1, strip.position = "right") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")
dev.off()


# ===============================
# Utility: Shared clone fraction
# ===============================

shared_fraction <- function(target_cluster,
                            reference_cluster,
                            location){
  
  ref <- TCR %>%
    filter(majorCluster == reference_cluster,
           clone.status == "Clonal") %>%
    pull(clone.id) %>% unique()
  
  target <- TCR %>%
    filter(majorCluster == target_cluster,
           loc == location) %>%
    pull(clone.id) %>% unique()
  
  length(intersect(ref, target)) / length(ref)
}


# ===============================
# Fig 2H – Transition heatmap
# ===============================

dat.plot <- as.matrix(
  subset(out@pIndex.tran, aid == out@proj)[,c(-1,-2,-3)]
)

rownames(dat.plot) <- 
  subset(out@pIndex.tran, aid == out@proj)[,3]

dat.plot[is.na(dat.plot)] <- 0

CD8 <- dat.plot[
  grep("CD8", rownames(dat.plot)),
  grep("CD8", colnames(dat.plot))
]

col_fun <- colorRamp2(
  seq(0, max(CD8), length=10),
  colorRampPalette(
    c("#efecea","#fec600","#a11826")
  )(10)
)

p_heat <- Heatmap(
  CD8,
  name="pIndex.tran",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE
)

pdf("./Fig2H.pdf",5,5)
print(p_heat)
dev.off()


# ===============================
# Fig 2L – Clone distribution heatmaps
# ===============================

prepare_clone_matrix <- function(df){
  
  df <- df[df$clone.status == "Clonal",]
  df$ClusterSite <- paste0(df$majorCluster,"_",df$loc)
  
  mat <- as.data.frame.matrix(
    table(df$clone.id, df$ClusterSite)
  )
  
  mat <- as.matrix(mat)
  mat[mat > 8] <- 9
  
  mat
}

TCR_subset <- TCR[
  TCR$majorCluster %in%
    c("T03_CD8_GZMB",
      "T04_CD8_GZMK",
      "T10_CD8_CCL5"),]

TCR_pt <- TCR_subset[TCR_subset$loc %in% c("PT","LN"),]
TCR_mt <- TCR_subset[!TCR_subset$loc %in% c("PT","LN"),]

mat_pt <- prepare_clone_matrix(TCR_pt)
mat_mt <- prepare_clone_matrix(TCR_mt)

col_fun <- colorRamp2(
  c(0:8, max(mat_pt)),
  c("lightgrey","#ffffcc","#ffeda0","#fed976",
    "#feb24c","#fd8d3c","#fc4e2a",
    "#e31a1c","#b10026","#800026")
)

heat_pt <- Heatmap(
  mat_pt,
  name="Cells",
  col=col_fun,
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  show_row_names=FALSE,
  column_title="TCR PT clone distribution"
)

heat_mt <- Heatmap(
  mat_mt,
  name="Cells",
  col=col_fun,
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  show_row_names=FALSE,
  column_title="TCR MT clone distribution"
)

pdf("./Fig2L.pdf",3.4,5.3)
heat_pt | heat_mt
dev.off()
