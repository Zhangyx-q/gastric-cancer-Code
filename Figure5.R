# ==============================================================================
# 1. Global Dependencies & Setup
# ==============================================================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(ggpubr)
library(ggridges)
library(nichenetr)
library(GSVA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(limma)
library(ggalluvial)

rm(list = ls())

# ==============================================================================
# 2. Figure 5B: T11_ILC3_AREG Proportion Analysis
# ==============================================================================
Tcell <- readRDS("./Tcell.rds")

prop_df <- Tcell@meta.data %>%
  group_by(orig.ident, cluster, Site) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(total_cells = sum(cell_count),
         prop = cell_count / total_cells) %>%
  ungroup()

target_cluster <- "T11_ILC3_AREG"
prop_sub <- prop_df %>% filter(cluster == target_cluster)

p5b <- ggplot(prop_sub, aes(x = factor(Site, levels=c("LN", "MLN", "MMLN", "PT", "MT")), 
                            y = prop, fill = Site)) +
  stat_boxplot(geom = "errorbar", width = 0.3) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  scale_fill_manual(values = c("#007891","#129981","#3A5182","#DA4B35","#EA977C")) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw(base_size = 14) +
  labs(x = "Group", y = "Cell Proportion", title = target_cluster)

pdf("./Fig5B.pdf", 3.7, 2.5)
print(p5b)
dev.off()

# ==============================================================================
# 3. Figure 5F: ILC3 Marker Expression (Ridge Plot)
# ==============================================================================
ILC3_sub <- subset(Tcell, cluster == "T11_ILC3_AREG")
ridge_genes <- c("KLRB1", "IL7R", "KIT", "HLA-DRA", "CCR6", "LTB", "AREG", "TNFSF11")

df_ridge <- FetchData(ILC3_sub, vars = ridge_genes) %>%
  pivot_longer(cols = everything(), names_to = "gene", values_to = "expression") %>%
  filter(expression > 0)

df_ridge$gene <- factor(df_ridge$gene, levels = rev(ridge_genes))

p5f <- ggplot(df_ridge, aes(x = expression, y = gene)) +
  geom_density_ridges(scale = 1.3, rel_min_height = 0.01, fill = "#ad5555", 
                      color = "black", size = 0.75) +
  theme_classic(base_size = 10) +
  labs(x = "Log-normalized expression", y = NULL)

pdf("./Fig5F.pdf")
print(p5f)
dev.off()

# ==============================================================================
# 4. Figure 5H: NicheNet Ligand-Target Network
# ==============================================================================
GCdata <- readRDS("./GCdata.rds") %>% subset(Site %in% c("LN", "MLN", "MMLN"))
Idents(GCdata) <- "cluster"

# Load NicheNet assets
lr_network <- readRDS("./lr_network_human_21122021.rds") %>% distinct(from, to)
ligand_target_matrix <- readRDS("./ligand_target_matrix_nsga2r_final.rds")

# Identify Differentially Expressed Genes for Receiver
markers_ilc3 <- FindMarkers(GCdata, ident.1 = 'T11_ILC3_AREG', min.pct = 0.5, logfc.threshold = 0.5)
geneset_oi <- rownames(markers_ilc3) %>% .[. %in% rownames(ligand_target_matrix)]
background_genes <- rownames(GCdata) %>% .[. %in% rownames(ligand_target_matrix)]

# Predict Ligand Activity
ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi, 
  background_expressed_genes = background_genes, 
  ligand_target_matrix = ligand_target_matrix, 
  potential_ligands = lr_network$from
)

best_ligands <- ligand_activities %>% top_n(20, aupr_corrected) %>% pull(test_ligand)

# Prepare Visualisation
active_links <- best_ligands %>% 
  lapply(get_weighted_ligand_target_links, geneset = geneset_oi, 
         ligand_target_matrix = ligand_target_matrix, n = 200) %>% 
  bind_rows() %>% drop_na()

vis_matrix <- prepare_ligand_target_visualization(active_links, ligand_target_matrix, cutoff = 0.25)

p5h <- vis_matrix %>% 
  make_heatmap_ggplot("Prioritized ligands", "Predicted target genes", 
                      color = "mediumvioletred", legend_title = "Regulatory potential") +
  theme(axis.text.x = element_text(face = "italic"))

ggsave('./Fig5H.pdf', p5h, width = 8, height = 6)



library(dplyr)
library(tidyr)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(survival)
library(survminer)
library(ggplot2)
library(GSVA)

# ==========================================
# 1. Load expression data
# ==========================================

data <- read.table(
  "E:/胃癌/20250319/TCGA-STAD.star_fpkm.tsv",
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

# Remove Ensembl version suffix
data$ENSEMBL <- sub("\\..*$", "", data$Ensembl_ID)

# ==========================================
# 2. Ensembl ID -> Gene Symbol
# ==========================================

gene_map <- bitr(
  unique(data$ENSEMBL),
  fromType = "ENSEMBL",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
) %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

# Merge gene annotation
data1 <- data %>%
  select(-Ensembl_ID) %>%
  inner_join(gene_map, by = "ENSEMBL") %>%
  select(-ENSEMBL) %>%
  relocate(SYMBOL)

# Average expression of duplicated gene symbols
expr <- avereps(
  as.matrix(data1[, -1]),
  ID = data1$SYMBOL
)

expr <- as.matrix(expr)

# ==========================================
# 3. Load clinical data
# ==========================================

clinical <- read.table(
  "E:/胃癌/20250319/TCGA_clinicalNEW.txt",
  sep = "\t",
  header = TRUE,
  quote = "",
  check.names = FALSE
)

clinical$sample <- paste0(clinical$sample, "A")
rownames(clinical) <- clinical$sample

# Match samples
common_samples <- intersect(
  rownames(clinical),
  colnames(expr)
)

clinical <- clinical[common_samples, , drop = FALSE]
expr <- expr[, common_samples, drop = FALSE]

# Check sample matching
stopifnot(
  identical(rownames(clinical), colnames(expr))
)

# ==========================================
# 4. Extract PDGFC expression
# ==========================================

stopifnot("PDGFC" %in% rownames(expr))

clinical$PDGFC <- as.numeric(expr["PDGFC", ])

# ==========================================
# 5. Calculate ILC3 ssGSEA score
# ==========================================

ILC3_genes <- c(
  "AREG", "SPINK2", "LINC00299", "KIT",
  "SCN1B", "LST1", "SOX4", "TYROBP",
  "TMIGD2", "IL7R"
)

# Retain genes present in expression matrix
ILC3_genes <- intersect(
  ILC3_genes,
  rownames(expr)
)

# Check signature coverage
print(ILC3_genes)

# GSVA: ssGSEA
# Compatible with recent GSVA versions
ssgsea_par <- ssgseaParam(
  exprData = expr,
  geneSets = list(ILC3 = ILC3_genes),
  normalize = TRUE
)

scorez <- gsva(ssgsea_par)

clinical$ILC3 <- as.numeric(
  scorez["ILC3", rownames(clinical)]
)

# ==========================================
# 6. Save integrated data
# ==========================================

write.csv(
  clinical,
  "./TCGA_STAD_PDGFC_ILC3_survival_data.csv",
  row.names = FALSE
)



# ==========================================
# 7. Prepare 5-year survival data
# ==========================================

survival_data <- clinical %>%
  filter(
    !is.na(OS.time),
    !is.na(OS),
    OS.time > 0,
    !is.na(PDGFC),
    !is.na(ILC3)
  ) %>%
  mutate(
    time_5y = pmin(OS.time, 1825),
    event_5y = ifelse(OS.time > 1825, 0, OS)
  )

# ==========================================
# 8. Survival analysis function
# ==========================================

survival_analysis <- function(
    data,
    variable,
    output,
    split = "median"
) {

  # Remove missing values
  df <- data %>%
    filter(
      !is.na(.data[[variable]]),
      !is.na(time_5y),
      !is.na(event_5y)
    )

  # Continuous Cox regression
  formula_cox <- as.formula(
    paste0(
      "Surv(time_5y, event_5y) ~ ",
      variable
    )
  )

  cox_model <- coxph(
    formula_cox,
    data = df
  )

  print(summary(cox_model))

  # Grouping
  if (split == "median") {

    cutoff <- median(
      df[[variable]],
      na.rm = TRUE
    )

  } else if (split == "optimal") {

    cut <- surv_cutpoint(
      df,
      time = "time_5y",
      event = "event_5y",
      variables = variable,
      minprop = 0.4
    )

    cutoff <- cut$cutpoint$cutpoint

  } else {
    stop("split must be median or optimal")
  }

  df$group <- ifelse(
    df[[variable]] > cutoff,
    "High",
    "Low"
  )

  df$group <- factor(
    df$group,
    levels = c("Low", "High")
  )

  # Grouped Cox regression
  group_cox <- coxph(
    Surv(time_5y, event_5y) ~ group,
    data = df
  )

  print(summary(group_cox))

  # Kaplan-Meier
  km_fit <- survfit(
    Surv(time_5y, event_5y) ~ group,
    data = df
  )

  p <- ggsurvplot(
    km_fit,
    data = df,
    pval = TRUE,
    risk.table = FALSE,
    conf.int = FALSE,
    legend.title = variable,
    legend.labs = c("Low", "High"),
    xlab = "Time (days)",
    ylab = "Overall survival probability",
    xlim = c(0, 1825),
    break.time.by = 365,
    palette = c("#70b29c", "#a71930")
  )

  ggsave(
    paste0(output, ".pdf"),
    p$plot,
    width = 4.6,
    height = 4.1
  )

  return(
    list(
      data = df,
      cox = cox_model,
      group_cox = group_cox,
      km = km_fit,
      cutoff = cutoff,
      plot = p
    )
  )
}

# ==========================================
# 9. PDGFC survival analysis
# ==========================================

PDGFC_result <- survival_analysis(
  data = survival_data,
  variable = "PDGFC",
  output = "Figure5I",
  split = "median"
)

# ==========================================
# 10. ILC3 survival analysis
# ==========================================

ILC3_result <- survival_analysis(
  data = survival_data,
  variable = "ILC3",
  output = "Figure5J",
  split = "median"
)
