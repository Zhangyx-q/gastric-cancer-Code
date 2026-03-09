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



# ==============================================================================
# 6. Figure 5I & 5J: Clinical Correlation (TCGA-STAD)
# ==============================================================================
# Assuming 'clinInfo11' (Stage) and 'clinInfo111' (Grade) are prepared
# Plotting Correlation with Grade (Fig 5J)
p5j <- ggplot(clinInfo111, aes(x = tumor_grade, y = score, fill = tumor_grade)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.1) +
  scale_fill_manual(values = c("#c59b95", "#9c5a57", "#782d3e")) +
  geom_smooth(aes(x = tumor_grade_num, y = score, group = 1), method = "lm", color = "black") +
  stat_compare_means(method = "t.test", label = "p.signif") +
  theme_classic() +
  labs(title = "Grade vs ILC3 Score")

pdf("./Fig5J.pdf", 3.9, 2.7)
print(p5j)
dev.off()

# Plotting Correlation with Stage (Fig 5I)
p5i <- ggplot(clinInfo11, aes(x = ajcc_pathologic_stage1, y = score, fill = ajcc_pathologic_stage1)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.1) +
  scale_fill_manual(values = c("#b0a696","#c59b95", "#9c5a57", "#782d3e")) +
  geom_smooth(aes(x = tumor_grade_num, y = score, group = 1), method = "lm", color = "black") +
  theme_classic() +
  labs(title = "Stage vs ILC3 Score")

pdf("./Fig5I.pdf", 5, 3.5)
print(p5i)
dev.off()

# ==============================================================================
# 5. Figure 5L: External Validation (GSE5081)
# ==============================================================================
ilc3_signature <- list(c("AREG", "SPINK2", "LINC00299", "KIT", "SCN1B", "LST1", "SOX4", "TYROBP", "TMIGD2", "IL7R"))

gse_data <- read.table("./GSE5081.txt", sep = "\t", header = TRUE, row.names = 1) %>% as.matrix() %>% scale()
gse_clin <- read.table("./GSE5081_clin.txt", sep = "\t", row.names = 1, header = TRUE) %>% t() %>% as.data.frame()

ssgsea_scores <- GSVA::gsva(gse_data, ilc3_signature, method = "ssgsea") %>% t() %>% as.data.frame()
ssgsea_scores$Group <- factor(gse_clin$`!Sample_characteristics_ch1`, levels = c("HP-", "HP+"))

p5l <- ggplot(ssgsea_scores, aes(x = Group, y = V1, fill = Group)) +
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  stat_compare_means(method = "t.test") +
  scale_fill_manual(values = c("HP-" = "#007891", "HP+" = "#b22746")) +
  labs(x = "H. pylori infection", y = "ILC3 Signature Score", title = "GSE5081") +
  theme_classic()

pdf("./Fig5L.pdf", 3.9, 2.7)
print(p5l)
dev.off()
