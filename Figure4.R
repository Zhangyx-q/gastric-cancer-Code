rm(list = ls())

library(Seurat)

Bcell <- readRDS("./Bcell.rds")

Bcolors = c("#EF98A1","#FBE3C0","#E2F2CD","#DBAA77","#B6DAA7","#F9D5D5","#AFC8E2")

pdf("./Fig4A.pdf",6.3,3.7)
DimPlot(Bcell,cols = Bcolors,label = T)
dev.off()


# ==============================================================================
# Figure 3B: Tissue Distribution & ROE Heatmap
# ==============================================================================
# 1. Calculate ROE (Observed/Expected)
mapz <- Bcell@meta.data %>%
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

pdf("./Fig4B_Heatmap.pdf", 3.5, 2.5)
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

pdf("./Fig4B_Barplot.pdf", 7, 4)
print(p1 | p2)
dev.off()


rm(list = ls())



GCdata <- readRDS("./GCdata.rds")

GCdata@meta.data$ss="Tumor"

GCdata@meta.data[grep("LN$",GCdata@meta.data$orig.ident),"ss"]<-"Lymph"

tumor<-subset(GCdata,ss=="Tumor")
tumor<-subset(tumor,cluster!="Epithelial")
tumor<-subset(tumor,cluster!="EEC")

cluster<-c("T02_CD4_SELL","T05_CD4_CTLA4","B07_Cycling_MKI67",
           "T06_CD8_IL7R","T01_CD4_LEF1","T08_CD4_ANXA1","E2_lymphatic EC_TNFRSF10C",
           "E1_TAEC_ESM1","E5_Artery_NEAT1")

prop_df <- tumor@meta.data %>%
  dplyr::group_by(orig.ident, cluster, Site) %>%
  summarise(cell_count = n()) %>%
  dplyr::group_by(orig.ident) %>%
  mutate(total_cells = sum(cell_count),
         prop = cell_count / total_cells) %>%
  ungroup()
library(dplyr)
library(ggplot2)
library(scales)
library(ggpubr)



library(dplyr)

prop_df1 <- prop_df %>%
  filter(cluster %in% c("T02_CD4_SELL", "T05_CD4_CTLA4", 
                        "B07_Cycling_MKI67", "T06_CD8_IL7R"))



library(ggplot2)
library(ggridges)
library(dplyr)

# 假设你的数据框叫 df，包含：Cluster_Label、NI_group、Ratio

# 先计算每组的平均值
mean_df <- prop_df1 %>%
  group_by(cluster, Site) %>%
  summarise(mean_ratio = mean(prop), .groups = "drop")
mean_df <- mean_df %>%
  mutate(y_value = as.numeric(factor(cluster, levels = rev(unique(prop_df1$cluster)))))

# 主图：ridge plot + 平均值线
p1=ggplot(prop_df1, aes(x = prop, y = cluster, fill = Site)) +
  geom_density_ridges(alpha = 0.6, scale = 1.2, color = "black") +
  geom_segment(data = mean_df,
               aes(x = mean_ratio, xend = mean_ratio,
                   y = y_value - 0.25, yend = y_value + 0.25,  # 控制短线高度
                   color = Site),
               size = 0.8, inherit.aes = FALSE) +
  scale_fill_manual(values = c("PT" = "#264653", "MT" = "#DA4B35")) +
  scale_color_manual(values = c("PT" = "#264653", "MT" = "#DA4B35")) +
  theme_minimal(base_size = 14) +
  labs(x = "Ratio", y = NULL) +
  theme(
    legend.title = element_blank(),
    axis.text.y = element_text(face = "bold")
  )

pdf("./Fig4C.pdf",6.4,4.3)
print(p1)
dev.off()



rm(list = ls())
Bcell <- readRDS("./Bcell.rds")
B5=subset(Bcell,cluster=="B05_Plasma_JCHAIN")
table(B5@meta.data$Site)
B5=subset(B5,Site%in%c("PT","MT"))

#IGxx <- c("IGHM","IGHA1","IGHG1","IGHG3","JCHAIN")
Idents(B5)=factor(B5@meta.data$Site,levels = c("PT","MT"))
aa=FindAllMarkers(B5)
aa=aa[which(aa$cluster=="MT"),]

mm=aa

source("./dotzz.R")

AA=dotzz(B5,assay = "RNA",features = aa$gene,group.by = "Site",scale = F)
head(AA)


library(ggplot2)

# 确保 id 是因子，并且顺序是 PT 在前、MT 在后

AA1=AA[,c(3,4,5)]

library(tidyr)
library(dplyr)

AA_wide <- AA1 %>%
  dplyr::select(features.plot, id, avg.exp.scaled) %>%
  pivot_wider(
    names_from = id,
    values_from = avg.exp.scaled
  )



# 把 gene 设为行名（可选）
AA_wide <- as.data.frame(AA_wide)
rownames(AA_wide) <- AA_wide$features.plot
AA_wide=AA_wide[,-1]
AA_wide1=cbind(AA_wide,aa)


AA_wide1 <- 
  AA_wide1 %>% 
  mutate(change = as.factor(ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 1,'Significant','No Significant'))) 


library(ggplot2)

p1=ggplot(AA_wide1, aes(x = PT, y = MT)) +
  geom_point(color = "grey80", size = 1) +
  #geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Average expression of Plasma in PT",
    y = "Average expression of Plasma in MT"
  ) +
  theme_classic(base_size = 14)+
  geom_point(aes(color = change),
             size = 2, 
             alpha = 0.5) +
  scale_color_manual(values = c("#caccd1","#a71930"))  +
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank(),
        legend.position = 'right') +
  # 添加标签：
  geom_text_repel(data = filter(AA_wide1, abs(avg_log2FC) > 1.2 & -log10(p_val_adj) > 10),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label = gene),
                  #color = change),
                  size = 4)

pdf("./Fig4D.pdf",6.15,4.7)
print(p1)
dev.off()



rm(list = ls())
GCdata <- readRDS("./GCdata.rds")

pt<-subset(GCdata,Site=="PT")
pt@meta.data$samples=pt@meta.data$sampleId

MT<-subset(GCdata,Site=="MT")
MT@meta.data$samples=MT@meta.data$sampleId

# Tcell<-subset(Tcell,Site=="MLN")
# 
# genez<-c("CXCL10","CXCL9","STAT1","IRF3","IRF7","IRF9","ISG15")
# 
# DotPlot(Tcell,assay = "RNA",genez)
# 

pt<-subset(GCdata,Site=="MLN")

rm(GCdata)

# 假设 seurat_obj 是你的 Seurat 对象，cluster 信息在 seurat_obj$seurat_clusters
# 假设你想对所有 cluster 限制每簇最多 5000 个细胞

max_cells_per_cluster <- 5000
library(Seurat)
# 遍历所有 cluster，随机抽样细胞
selected_cells <- unlist(lapply(levels(Idents(pt)), function(cluster) {
  cells_in_cluster <- WhichCells(pt, idents = cluster)
  n_cells <- length(cells_in_cluster)
  
  if (n_cells > max_cells_per_cluster) {
    # 超过5000，随机选5000
    sample(cells_in_cluster, max_cells_per_cluster)
  } else {
    # 不超过5000，保留所有
    cells_in_cluster
  }
}))

# 用筛选后的细胞子集来做 CellChat
pt <- subset(pt, cells = selected_cells)

cellchat <- createCellChat(object=pt,group.by = "cluster")
cellchat
summary(cellchat)
str(cellchat)
levels(cellchat@idents)
#cellchat <- setIdent(cellchat, ident.use = "cell_type")
groupSize <- as.numeric(table(cellchat@idents))  


CellChatDB <- CellChatDB.human
#导入小鼠是CellChatDB <- CellChatDB.mouse
str(CellChatDB) #查看数据库信息
#包含interaction、complex、cofactor和geneInfo这4个dataframe

cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

cellchat <- identifyOverExpressedGenes(cellchat,thresh.p = 1)
cellchat <- identifyOverExpressedInteractions(cellchat)


cellchat <- computeCommunProb(cellchat, type = "truncatedMean",trim = 0.01)
computeAveExpr(cellchat, features = c("CCL2","CCR2"),type =  "truncatedMean",trim = 0.1)

cellchat1 <- filterCommunication(cellchat, min.cells = 10,  rare.keep =T)

df.net <- subsetCommunication(cellchat1)
which(df.net$receptor=="CCR2")
df.net1=df.net[grep("CCL2$",df.net$ligand),]

pdf("./Fig4G.pdf",10,7)

netVisual_individual(cellchat1, signaling = "CCL", pairLR.use = "CCL2_CCR2", 
                     color.use = colorA,
                     targets.use = "B05_Plasma_JCHAIN"
                     
)
dev.off()


library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

# 1. Load data
Fibroblast <- readRDS("./Fibroblast.rds")

# 2. Define colors
Scolors <- c(
  "S1" = "#d20962",
  "S2" = "#f47721",
  "S3" = "#7ac143",
  "S4" = "#00a78e"
)

Site_colors <- c(
  "LN" = "#007891",
  "MLN" = "#129981",
  "MMLN" = "#3A5182",
  "PT" = "#EA977C",
  "MT" = "#DA4B35"
)

# 3. UMAP
p_umap <- DimPlot(
  Fibroblast,
  reduction = "umap",
  group.by = "cluster",
  cols = Scolors,
  label = TRUE
) +
  theme_classic()

ggsave(
  "./Fig4H.pdf",
  p_umap,
  width = 6.3,
  height = 3.7
)

# 4. Calculate cell numbers
mapz <- Fibroblast@meta.data %>%
  count(cluster, Site, name = "n") %>%
  mutate(
    Site = factor(
      Site,
      levels = c("LN", "MLN", "MMLN", "PT", "MT")
    )
  )

# 5. Tissue composition in each cluster
p_cluster <- ggplot(
  mapz,
  aes(x = cluster, y = n, fill = Site)
) +
  geom_col(
    position = "fill",
    color = "grey",
    linewidth = 0.2
  ) +
  scale_fill_manual(values = Site_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  ) +
  labs(x = NULL, y = "Proportion", fill = "Site")

ggsave(
  "./Fig4H1.pdf",
  p_cluster,
  width = 5,
  height = 4
)

# 6. Fibroblast composition in each tissue
p_site <- ggplot(
  mapz,
  aes(x = Site, y = n, fill = cluster)
) +
  geom_col(
    position = "fill",
    color = "grey",
    linewidth = 0.2
  ) +
  scale_fill_manual(values = Scolors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  ) +
  labs(x = NULL, y = "Proportion", fill = "Cluster")

ggsave(
  "./Fig4H2.pdf",
  p_site,
  width = 5,
  height = 4
)



library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

# 1. Load data
Fibroblast <- readRDS("./Fibroblast.rds")

# 2. Define colors
Scolors <- c(
  "S1" = "#d20962",
  "S2" = "#f47721",
  "S3" = "#7ac143",
  "S4" = "#00a78e"
)

Site_colors <- c(
  "LN" = "#007891",
  "MLN" = "#129981",
  "MMLN" = "#3A5182",
  "PT" = "#EA977C",
  "MT" = "#DA4B35"
)

# 3. UMAP
p_umap <- DimPlot(
  Fibroblast,
  reduction = "umap",
  group.by = "cluster",
  cols = Scolors,
  label = TRUE
) +
  theme_classic()

ggsave(
  "./Fig4H.pdf",
  p_umap,
  width = 6.3,
  height = 3.7
)

# 4. Calculate cell numbers
mapz <- Fibroblast@meta.data %>%
  count(cluster, Site, name = "n") %>%
  mutate(
    Site = factor(
      Site,
      levels = c("LN", "MLN", "MMLN", "PT", "MT")
    )
  )

# 5. Tissue composition in each cluster
p_cluster <- ggplot(
  mapz,
  aes(x = cluster, y = n, fill = Site)
) +
  geom_col(
    position = "fill",
    color = "grey",
    linewidth = 0.2
  ) +
  scale_fill_manual(values = Site_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  ) +
  labs(x = NULL, y = "Proportion", fill = "Site")

ggsave(
  "./Fig4H1.pdf",
  p_cluster,
  width = 5,
  height = 4
)

# 6. Fibroblast composition in each tissue
p_site <- ggplot(
  mapz,
  aes(x = Site, y = n, fill = cluster)
) +
  geom_col(
    position = "fill",
    color = "grey",
    linewidth = 0.2
  ) +
  scale_fill_manual(values = Scolors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  ) +
  labs(x = NULL, y = "Proportion", fill = "Cluster")

ggsave(
  "./Fig4H2.pdf",
  p_site,
  width = 5,
  height = 4
)



rm(list = ls())

library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# ==========================================
# Figure 4 H&L
# ==========================================

Fibroblast <- readRDS(
  "./Fibroblast.rds"
)

# Only retain S2 and S3 Fibroblast
B5 <- subset(
  Fibroblast,
  subset = cluster %in% c(
    "S2_Fib_CXCL14",
    "S3_Fib_SFRP2"
  )
)

Idents(B5) <- "cluster"

aa <- FindMarkers(
  B5,
  ident.1 = "S2_Fib_CXCL14",
  ident.2 = "S3_Fib_SFRP2",
  only.pos = FALSE,
  logfc.threshold = 0,
  min.pct = 0,
  return.thresh = 1
)

aa$gene <- rownames(aa)

geneList <- aa$avg_log2FC
names(geneList) <- aa$gene

# Remove NA and duplicated genes
geneList <- geneList[
  !is.na(geneList) &
    !is.na(names(geneList)) &
    names(geneList) != ""
]

geneList <- geneList[
  !duplicated(names(geneList))
]

# Sort from high to low
geneList <- sort(
  geneList,
  decreasing = TRUE
)

gse_GO <- gseGO(
  geneList = geneList,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE
)


AA <- dotzz(
  B5,
  assay = "RNA",
  features = aa$gene,
  group.by = "cluster",
  scale = FALSE
)


AA_wide <- AA %>%
  select(features.plot, id, avg.exp.scaled) %>%
  pivot_wider(
    names_from = id,
    values_from = avg.exp.scaled
  )


AA_wide1 <- AA_wide %>%
  left_join(
    aa,
    by = c("features.plot" = "gene")
  ) %>%
  mutate(
    change = case_when(
      p_val_adj < 0.05 & abs(avg_log2FC) > 1 ~ "Significant",
      TRUE ~ "No Significant"
    ),
    change = factor(
      change,
      levels = c("No Significant", "Significant")
    )
  )


label_data <- AA_wide1 %>%
  filter(
    abs(avg_log2FC) > 1.5,
    !is.na(p_val_adj),
    p_val_adj > 0,
    -log10(p_val_adj) > 10
  )

p1 <- ggplot(
  AA_wide1,
  aes(
    x = S2_Fib_CXCL14,
    y = S3_Fib_SFRP2
  )
) +
  geom_point(
    color = "grey80",
    size = 1
  ) +
  geom_point(
    aes(color = change),
    size = 2,
    alpha = 0.5
  ) +
  geom_text_repel(
    data = label_data,
    aes(label = features.plot),
    size = 4,
    max.overlaps = 20
  ) +
  scale_color_manual(
    values = c(
      "No Significant" = "#caccd1",
      "Significant" = "#a71930"
    )
  ) +
  labs(
    x = "Average expression of S2_Fib_CXCL14",
    y = "Average expression of S3_Fib_SFRP2",
    color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

# ==========================================
# 6. Save
# ==========================================

ggsave(
  filename = "./Figure4H.pdf",
  plot = p1,
  width = 7.5,
  height = 6
)




mm <- as.data.frame(gse_GO)

pathw <- c(
  "basement membrane organization",
  "cell-cell adhesion",
  "regulation of wound healing",
  "blood vessel development",
  "immune system process",
  "humoral immune response",
  "complement activation"
)

mm1 <- mm %>%
  filter(Description %in% pathw)

plot_data <- mm1 %>%
  mutate(
    logP = -log10(pmax(p.adjust, 1e-300)),
    logP_signed = ifelse(
      NES > 0,
      logP,
      -logP
    ),
    Description = factor(
      Description,
      levels = Description[order(NES)]
    )
  )


p1 <- ggplot(
  plot_data,
  aes(
    x = logP_signed,
    y = Description,
    fill = NES > 0
  )
) +
  geom_col(width = 0.6) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.5
  ) +
  scale_fill_manual(
    values = c(
      "TRUE" = "#FA8400",
      "FALSE" = "#59BD48"
    ),
    labels = c(
      "TRUE" = "NES > 0",
      "FALSE" = "NES < 0"
    )
  ) +
  scale_x_continuous(
    labels = abs
  ) +
  labs(
    x = expression(-log[10]("adjusted P")),
    y = NULL,
    fill = "NES"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 10),
    legend.position = "top"
  )

# ==========================================
# 8. Save figure
# ==========================================

ggsave(
  filename = "./Figure4L.pdf",
  plot = p1,
  width = 5.6,
  height = 4.5
)







