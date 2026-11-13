
# ============================================================
# STAD scRNA-seq data preprocessing and integration
# Author: 
# Description:
#   1. Read and merge multiple 10X samples
#   2. Quality control (QC)
#   3. Doublet detection using DoubletFinder
#   4. Cell and gene filtering
#   5. Normalization and highly variable gene selection
#   6. PCA and Harmony integration
#   7. UMAP visualization and clustering
#   8. Marker gene identification
#
# ============================================================


# ============================================================
# 0. Load packages and set working directory
# ============================================================

rm(list = ls())

library(Seurat)
library(harmony)
library(DoubletFinder)
library(ggplot2)
library(patchwork)
library(cowplot)
library(RColorBrewer)
library(dplyr)
library(stringr)

# Set working directory
setwd("/data/home/zhangyuxi/data/STAD/rawdata/run")

# Create output directory
dir.create("results", showWarnings = FALSE)
dir.create("results/QC", showWarnings = FALSE)
dir.create("results/UMAP", showWarnings = FALSE)
dir.create("results/markers", showWarnings = FALSE)
dir.create("results/RDS", showWarnings = FALSE)


# ============================================================
# 1. Read 10X data and create Seurat objects
# ============================================================

# Find all sample_filtered_feature_bc_matrix directories
data_dirs <- list.dirs(
  path = ".",
  recursive = TRUE,
  full.names = TRUE
)

data_dirs <- data_dirs[
  grepl("sample_filtered_feature_bc_matrix$", data_dirs)
]

# Keep directories under per_sample_outs
data_dirs <- data_dirs[
  grepl("per_sample_outs", data_dirs)
]

# Check sample directories
print(data_dirs)

# Extract sample names
sample_names <- basename(dirname(dirname(data_dirs)))

# Check sample names
print(sample_names)

# Read each sample and create Seurat object
sceList <- lapply(seq_along(data_dirs), function(i) {

  sample_id <- sample_names[i]

  message("Reading sample: ", sample_id)

  data <- Read10X(data.dir = data_dirs[i])

  # Add sample ID to cell barcodes
  colnames(data) <- paste0(
    sample_id, "_", colnames(data)
  )

  # Create Seurat object
  obj <- CreateSeuratObject(
    counts = data,
    project = sample_id
  )

  # Calculate mitochondrial gene percentage
  obj <- PercentageFeatureSet(
    obj,
    pattern = "^MT-",
    col.name = "percent.mt"
  )

  return(obj)
})

names(sceList) <- sample_names


# ============================================================
# 2. Merge all samples
# ============================================================

sceList1 <- merge(
  x = sceList[[1]],
  y = sceList[-1]
)

# Check merged object
sceList1

table(sceList1$orig.ident)


# ============================================================
# 3. Quality control (QC)
# ============================================================

# Define sample order
sample_order <- sample_names

sceList1$orig.ident <- factor(
  sceList1$orig.ident,
  levels = sample_order
)

Idents(sceList1) <- "orig.ident"

# Define colors for samples
sample_colors <- colorRampPalette(
  brewer.pal(12, "Set3")
)(length(sample_order))

names(sample_colors) <- sample_order

# QC violin plots
p_qc <- VlnPlot(
  sceList1,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  group.by = "orig.ident",
  cols = sample_colors,
  ncol = 3,
  pt.size = 0
) &
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )

ggsave(
  "results/QC/QC_violin.pdf",
  plot = p_qc,
  width = 20,
  height = 5
)

# Save metadata
write.csv(
  sceList1@meta.data,
  "results/QC/QC_metadata.csv"
)


# ============================================================
# 4. Split samples for DoubletFinder
# ============================================================

sce.all.list <- SplitObject(
  sceList1,
  split.by = "orig.ident"
)

# DoubletFinder parameters
npcs <- 30
expected_doublet_rate <- 0.05

# Run DoubletFinder independently for each sample
phe_lt <- lapply(names(sce.all.list), function(x) {

  message("Processing DoubletFinder: ", x)

  sce.all.filt <- sce.all.list[[x]]

  # Highly variable features
  sce.all.filt <- FindVariableFeatures(
    sce.all.filt,
    selection.method = "vst",
    nfeatures = 3000,
    verbose = FALSE
  )

  # Scaling and PCA
  sce.all.filt <- ScaleData(
    sce.all.filt,
    verbose = FALSE
  )

  sce.all.filt <- RunPCA(
    sce.all.filt,
    npcs = npcs,
    verbose = FALSE
  )

  # Parameter sweep
  sweep.res.list <- paramSweep(
    sce.all.filt,
    PCs = 1:npcs,
    sct = FALSE
  )

  sweep.stats <- summarizeSweep(
    sweep.res.list,
    GT = FALSE
  )

  bcmvn <- find.pK(sweep.stats)

  # Select optimal pK
  pK.1 <- as.numeric(
    as.character(
      bcmvn$pK[which.max(bcmvn$BCmetric)]
    )
  )

  # Expected number of doublets
  nExp <- round(
    ncol(sce.all.filt) * expected_doublet_rate
  )

  # Run DoubletFinder
  sce.all.filt <- doubletFinder(
    sce.all.filt,
    PCs = 1:npcs,
    pN = 0.25,
    pK = pK.1,
    nExp = nExp,
    sct = FALSE
  )

  # Rename classification column
  df_col <- grep(
    "^DF.classifications",
    colnames(sce.all.filt@meta.data),
    value = TRUE
  )

  sce.all.filt$doublet_info <- sce.all.filt@meta.data[
    [tail(df_col, 1)]
  ]

  return(sce.all.filt)
})

names(phe_lt) <- names(sce.all.list)

# Save DoubletFinder results
saveRDS(
  phe_lt,
  "results/RDS/DoubletFinder_results.rds"
)


# ============================================================
# 5. Merge DoubletFinder results
# ============================================================

sceList1 <- merge(
  x = phe_lt[[1]],
  y = phe_lt[-1]
)

# Check doublet classification
table(sceList1$doublet_info)


# ============================================================
# 6. Cell and gene filtering
# ============================================================

# Filter low-quality cells
sceList1 <- subset(
  sceList1,
  subset = percent.mt < 15 & nFeature_RNA > 500
)

# Retain singlets only
sceList1 <- subset(
  sceList1,
  subset = doublet_info == "Singlet"
)

# Check retained cells
table(sceList1$orig.ident)

# Add patient-level sample ID
# Samples from the same patient share the same Sample ID
sceList1$sampleId <- paste0(
  "Sample",
  sub("^(\\d+).*", "\\1", as.character(sceList1$orig.ident))
)

# Check patient IDs
table(sceList1$sampleId)


# ============================================================
# 7. Gene filtering
# ============================================================

# Define genes to remove
exclude_patterns <- c(
  "^MT-"
  # "^RP[SL]",    # Ribosomal protein genes
  # "^MRP[SL]",   # Mitochondrial ribosomal protein genes
  # "^RNA[18]S",  # rRNA
  # "^RNA28S",
  # "^RNA5S"
)

# Identify genes to remove
genes_to_remove <- grep(
  pattern = paste(exclude_patterns, collapse = "|"),
  x = rownames(sceList1),
  value = TRUE,
  ignore.case = TRUE
)

# Retain genes
keep_genes <- setdiff(
  rownames(sceList1),
  genes_to_remove
)

sceList1 <- subset(
  sceList1,
  features = keep_genes
)

# Check filtered object
sceList1


# ============================================================
# 8. Normalization and highly variable gene selection
# ============================================================

sceList1 <- NormalizeData(
  sceList1,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

sceList1 <- FindVariableFeatures(
  sceList1,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

# Scaling
sceList1 <- ScaleData(
  sceList1,
  verbose = FALSE
)

# PCA
sceList1 <- RunPCA(
  sceList1,
  npcs = 30,
  verbose = FALSE
)


# ============================================================
# 9. Harmony integration
# ============================================================

sceList1 <- RunHarmony(
  sceList1,
  group.by.vars = "sampleId",
  assay.use = "RNA",
  max_iter = 20
)

# Save integrated object
saveRDS(
  sceList1,
  "results/RDS/STAD_Harmony_integrated.rds"
)


# ============================================================
# 10. UMAP visualization based on PCA
# ============================================================

sceList_pca <- RunUMAP(
  sceList1,
  reduction = "pca",
  dims = 1:30
)

p_pca <- DimPlot(
  sceList_pca,
  reduction = "umap",
  group.by = "orig.ident",
  cols = sample_colors,
  label = TRUE
)

ggsave(
  "results/UMAP/UMAP_PCA_by_sample.pdf",
  plot = p_pca,
  width = 16,
  height = 15
)


# ============================================================
# 11. UMAP and clustering based on Harmony
# ============================================================

sceList1 <- RunUMAP(
  sceList1,
  reduction = "harmony",
  dims = 1:30
)

sceList1 <- FindNeighbors(
  sceList1,
  reduction = "harmony",
  dims = 1:30
)

# Clustering
sceList1 <- FindClusters(
  sceList1,
  resolution = 0.6
)

# UMAP by cluster
p_cluster <- DimPlot(
  sceList1,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  raster = FALSE
)

ggsave(
  "results/UMAP/UMAP_Harmony_clusters.pdf",
  plot = p_cluster,
  width = 16,
  height = 15
)

# UMAP by original sample
p_sample <- DimPlot(
  sceList1,
  reduction = "umap",
  group.by = "orig.ident",
  cols = sample_colors,
  label = TRUE
)

ggsave(
  "results/UMAP/UMAP_Harmony_samples.pdf",
  plot = p_sample,
  width = 16,
  height = 15
)


# ============================================================
# 12. Marker gene identification
# ============================================================

# Seurat v5: join expression layers before differential expression
if (packageVersion("SeuratObject") >= "5.0.0") {
  sceList1 <- JoinLayers(sceList1)
}

# Identify marker genes for all clusters
markers_all <- FindAllMarkers(
  sceList1,
  only.pos = FALSE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)

# Export marker genes
write.csv(
  markers_all,
  "results/markers/All_cluster_markers.csv",
  row.names = FALSE
)


# ============================================================
# 13. Save final Seurat object
# ============================================================

saveRDS(
  sceList1,
  "results/RDS/STAD_scRNAseq_final.rds"
)


new.cluster.ids<-0:21
new.cluster.ids[which(new.cluster.ids%in%c(1,	2,3,4,8,11))]<-"T"
new.cluster.ids[which(new.cluster.ids%in%c(0,5,6,10,12))]<-"B"
new.cluster.ids[which(new.cluster.ids%in%c(7))]<-"Epithelial"
new.cluster.ids[which(new.cluster.ids%in%c(9,14,15,18))]<-"Myeloid"
new.cluster.ids[which(new.cluster.ids%in%c(20))]<-"mural"
new.cluster.ids[which(new.cluster.ids%in%c(13))]<-"endothelial"
new.cluster.ids[which(new.cluster.ids%in%c(16))]<-"fibroblast"
new.cluster.ids[which(new.cluster.ids%in%c(19))]<-"enteroendocrine"
new.cluster.ids[which(new.cluster.ids%in%c(17,21))]<-"delete"

names(new.cluster.ids)<-0:21
sceList1@meta.data$celltype<-new.cluster.ids[match(sceList1@meta.data$seurat_clusters,names(new.cluster.ids))]

Idents(stad)<-sceList1$celltype

col=c(c(brewer.pal(12, "Set3")[-c(6,9)]),
      (brewer.pal(12, "Set3")[6]))

levels(sceList1)<-sort(levels(sceList1))

pdf("umap_cellcluster.pdf",6,5)
DimPlot(stad,reduction ="umap",group.by = "seurat_clusters" , label = T)
dev.off()

#############Tcell recluster###########

library(Seurat)

# 1. Read T-cell subset
Tcell <- readRDS("/T.rds")



# 3. Define T-cell subclustering function
seuratz <- function(sub_obj) {

  # Re-normalization and feature selection
  sub_obj <- NormalizeData(sub_obj)
  sub_obj <- FindVariableFeatures(
    sub_obj,
    selection.method = "vst"
  )
  sub_obj <- ScaleData(sub_obj)

  # PCA
  sub_obj <- RunPCA(sub_obj)

  # Neighbors and clustering
  sub_obj <- FindNeighbors(
    sub_obj,
    reduction = "harmony",
    dims = 1:30
  )

  sub_obj <- FindClusters(
    sub_obj,
    resolution = 0.9
  )

  # UMAP
  sub_obj <- RunUMAP(
    sub_obj,
    reduction = "harmony",
    dims = 1:30
  )

  return(sub_obj)
}

# 4. Set original cluster identities
Idents(Tcell) <- Tcell$seurat_clusters

# 5. Define cluster annotations
namez1 <- c(
  "T01_CD4_LEF1",
  "T02_CD4_SELL",
  "T03_CD8_GZMB",
  "T04_CD8_GZMK",
  "T05_CD4_CTLA4",
  "T06_CD8_IL7R",
  "T07_CD4_TNFRSF4",
  "T08_DNT_ANXA1",
  "T09_gdT_GNLY",
  "T10_CD8_CCL5",
  "T11_ILC-like_AREG",
  "T12_CD8_FOS",
  "T01_CD4_LEF1"
)

names(namez1) <- as.character(0:12)

# 6. Retain clusters 0–12
Tcell <- subset(Tcell, idents = as.character(0:12))

# 7. Add cell-type annotations
Tcell@meta.data$cluster <- namez1[
  match(
    as.character(Tcell@meta.data$seurat_clusters),
    names(namez1)
  )
]

# 8. Set annotated identities
Idents(Tcell) <- "cluster"
levels(Tcell) <- sort(levels(Tcell))

# 9. Save annotated object
saveRDS(Tcell, "Tcell.rds")

#######B cell#######

Bcell <- readRDS("./B.rds")
Bcell@meta.data$delete<-0
Bcell@meta.data[names(which(Bcell@assays$RNA@counts["CD3D",]>0)),"delete"]<-1
Bcell@meta.data[names(which(Bcell@assays$RNA@counts["CD3E",]>0)),"delete"]<-1
Bcell@meta.data[names(which(Bcell@assays$RNA@counts["CD3G",]>0)),"delete"]<-1

table(Bcell@meta.data$delete)
# table(Bcell@meta.data[which(Bcell@meta.data$delete==1),"orig.ident"])
Bcell<-subset(Bcell,delete==0)

Bcell<-seuratz(Bcell)

Bcell <- RunUMAP(Bcell, reduction = "harmony", dims = 1:30)

DimPlot(Bcell,label = T)

BcellMarker<-FindAllMarkers(Bcell)
write.csv(BcellMarker," \\Bcell\\BcellMarker6.csv",row.names = F)

Bcell<-subset(Bcell,idents = c(0, 1,2, 3, 4,5,6,7))

Bcolors = c("#EF98A1","#FBE3C0","#E2F2CD","#DBAA77","#B6DAA7",
            "#F9D5D5","#AFC8E2","#DEECF6")

# 0: B_naive_IGHD
# 1: B_naive_SELL
# 2: B_Naive_CCR7
# 3: Plasma JCHAIN
# 4: B_Cycling_MKI67
# 5: GC MEF2B

namez1<-c("B01_naive_IGHD","B02_memery_CD27","B03_memery_TNFRSF13B","B04_HSP_CCR7","B05_Plasma_JCHAIN","B06_GCB_LMO2","B07_Cycling_MKI67","B05_Plasma_JCHAIN")

names(namez1)<-c(0, 1, 2, 3, 4,5,6,7)

Bcell@meta.data$cluster<-namez1[match(Bcell@meta.data$seurat_clusters,names(namez1))]
Idents(Bcell)<-Bcell$cluster
levels(Bcell)<-sort(levels(Bcell))
saveRDS(Bcell,"Bcell.rds")


##############Myeloid################

Myeloid <- readRDS("/Myeloid.rds")

Myeloid<-seuratz(Myeloid)

DimPlot(Myeloid,label = T)
table(Myeloid@meta.data$RNA_snn_res.0.7)
MyeloidMarker<-FindAllMarkers(Myeloid)
write.csv(MyeloidMarker,"\\Myeloid\\MyeloidMarker7.csv",row.names = F)

Mcolors = c("#F5AB5E","#AFC778","#947959","#EFDBB9","#5F97C6","#A3CDEA","#F09496","#F4DFDD")

Myeloid1<-subset(Myeloid,idents = c(0, 1, 2, 3, 4,6,7,8,9))


namez1<-c("M01_cDC3_CCL19","M02_Neutrophil_CXCL8","M03_TAM_C1QA","M04_Mast_TPSAB1","M02_Neutrophil_CXCL8",
          "M05_cDC2_CLEC10A","M03_TAM_C1QA","M06_cDC1_CLEC9A")

names(namez1)<-c(0, 1, 2, 3, 4,6,7,8)

Myeloid1@meta.data$cluster<-namez1[match(Myeloid1@meta.data$seurat_clusters,names(namez1))]
Idents(Myeloid1)<-Myeloid1$cluster
levels(Myeloid1)<-sort(levels(Myeloid1))
saveRDS(Myeloid1,"Myeloid.rds")

###################################################

Fibroblast <- readRDS("/Fibroblast.rds")

Fibroblast<-seuratz(Fibroblast)

FibroblastMarker<-FindAllMarkers(Fibroblast)
write.csv(FibroblastMarker,"\\FibroblastMarker.csv",row.names = F)

Fibroblast1<-Fibroblast

namez1<-c("S1_Pericyte","S1_Pericyte","S2_CAF_CXCL14","S3_CAF_SFRP2","S2_CAF_CXCL14","S4_myCAF_COL1A1")


names(namez1)<-c( 0,1, 2, 3,4, 5)

Fibroblast1@meta.data$cluster<-namez1[match(Fibroblast1@meta.data$seurat_clusters,names(namez1))]
Idents(Fibroblast1)<-Fibroblast1$cluster
levels(Fibroblast1)<-sort(levels(Fibroblast1))
saveRDS(Fibroblast1,"Fibroblast.rds")


