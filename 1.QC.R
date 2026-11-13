
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

