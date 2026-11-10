rm(list = ls())



cols=c("#058d77","#ce4a3a","#007891","#dd7d28","#aa1d38","#e3aab3","#28a0b8",
       "#b4a679")


pdf("./Figure7A.pdf",4.5,3.5)
DimPlot(Epithelial_meta,group.by = "cluster",cols = cols)
dev.off()

###################################

library(dplyr)
mapz <- Epithelial_meta@meta.data %>%
  dplyr::group_by(Site,cluster,orig.ident) %>%
  dplyr::summarise(n = dplyr::n())
library(ggplot2)
pdf("./Figure7B.pdf",4,3)
p3=ggplot(data = mapz, mapping = aes(x = cluster, y = n, fill = Site, group = interaction(Site, orig.ident,lex.order = T)))  + 
  geom_bar(stat = 'identity', position = 'fill', color = "grey", size = 0.2)+
  scale_fill_manual(values =c("#129981","#3A5182","#DA4B35"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
print(p3)
dev.off()



DotPlot(Epithelial_meta,features = c("CD24","KRAS","MMP7","IFITM1","IFITM3"))

  source(./dotzz.R")

pdf("./Figure7d.pdf",3.5,2.5)
  plotDotz(Epithelial_meta,assay = "RNA",features = c("CD24","KRAS","CD47","MMP7","IFITM1","IFITM3","CD74"))+
  
  scale_fill_gradientn(colours = c("#3A84B1","#CCE195","#F9EDA7","#ED8156","#CD464E"))+ #颜色p
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
dev.off()






############################
## Figure 7E
############################

library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)


output_dir <- "results/Figure4"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Identify positive and negative markers for all clusters
markerz <- FindAllMarkers(
  Epithelial_meta,
  assay = "RNA",
  only.pos = FALSE,
  min.pct = 0,
  logfc.threshold = 0,
  return.thresh = 1
)

# Save marker results
write.csv(
  markerz,
  file = file.path(output_dir, "Epithelial_marker_genes.csv"),
  row.names = FALSE
)


target_cluster <- 2

marker2 <- markerz %>%
  filter(cluster == target_cluster) %>%
  filter(!is.na(gene))

# Compatible with different Seurat versions
fc_column <- if ("avg_log2FC" %in% colnames(marker2)) {
  "avg_log2FC"
} else {
  "avg_logFC"
}


# Convert SYMBOL to ENTREZID
gene_annotation <- bitr(
  unique(marker2$gene),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Merge gene annotation with marker results
marker2_entrez <- marker2 %>%
  inner_join(
    gene_annotation,
    by = c("gene" = "SYMBOL")
  ) %>%
  filter(
    !is.na(ENTREZID),
    !is.na(.data[[fc_column]])
  ) %>%
  group_by(ENTREZID) %>%
  # Keep one entry for each ENTREZID
  slice_max(
    order_by = abs(.data[[fc_column]]),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

# Construct ranked gene list
geneList <- marker2_entrez[[fc_column]]
names(geneList) <- marker2_entrez$ENTREZID

geneList <- sort(
  geneList,
  decreasing = TRUE
)



go_enrich <- gseGO(
  geneList = geneList,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)



kegg_enrich <- gseKEGG(
  geneList = geneList,
  organism = "hsa",
  keyType = "ncbi-geneid",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)



# Read Hallmark gene sets
hallmark_file <- "data/h.all.v2024.1.Hs.symbols.gmt"

hall <- read.gmt(hallmark_file)

# Prepare SYMBOL-based ranked gene list
geneList_symbol <- marker2 %>%
  filter(
    !is.na(gene),
    !is.na(.data[[fc_column]])
  ) %>%
  group_by(gene) %>%
  slice_max(
    order_by = abs(.data[[fc_column]]),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

geneList_symbol <- geneList_symbol[[fc_column]]
names(geneList_symbol) <- marker2 %>%
  filter(
    !is.na(gene),
    !is.na(.data[[fc_column]])
  ) %>%
  distinct(gene, .keep_all = TRUE) %>%
  pull(gene)

# Sort genes by log2FC
geneList_symbol <- sort(
  geneList_symbol,
  decreasing = TRUE
)

# Hallmark GSEA
hallmark_gsea <- GSEA(
  geneList = geneList_symbol,
  TERM2GENE = hall,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

hallmark_result <- as.data.frame(hallmark_gsea)

write.csv(
  hallmark_result,
  file = file.path(output_dir, "cluster2_Hallmark_GSEA.csv"),
  row.names = FALSE
)


saveRDS(
  list(
    markers = markerz,
    cluster2_markers = marker2,
    GO_GSEA = go_enrich,
    KEGG_GSEA = kegg_enrich,
    GO_ORA = go_ora,
    Hallmark_GSEA = hallmark_gsea
  ),
  file = file.path(output_dir, "Figure7e.rds")
)

message("All enrichment analyses completed!")

