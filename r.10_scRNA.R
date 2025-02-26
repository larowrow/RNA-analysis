rm(list = ls())
setwd('')
# Load required packages -------------------------------------------------------
library(Seurat)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(dplyr)
library(tidyr)
library(CellChat)
library(ggrepel)
# 1. Read data and basic preprocessing ----------------------------------------
seu1 <- readRDS("GSE154396.rds")  
seu2 <- readRDS("GSE225948.rds")

# Here we take seu1 as an example object for demonstration.
seu <- seu1

# Ensure correct grouping factors if needed (example: "treatment" column)
# Adjust the factor levels based on your dataset.
seu$treatment <- factor(seu$treatment, levels = c("Sham","MCAO"))

# Standard Seurat workflow
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 30)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, dims = 1:20)

# 2. Violin + boxplot + statistical test --------------------------------------
show_genes <- c("Mcemp1", "Clec4d", "Cacna1e")

p_vln <- VlnPlot(
  seu, features = show_genes, group.by = "celltype", split.by = "treatment", pt.size = 0
) &
  geom_boxplot(
    width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA
  ) &
  stat_compare_means(aes(group = treatment), method = "t.test", label = "p.signif") &
  scale_fill_d3() &
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_vln)

# 3. Extract expression data manually for ggplot (optional) -------------------
plot_data <- FetchData(seu, vars = c(show_genes, "celltype", "treatment")) %>%
  pivot_longer(cols = all_of(show_genes), names_to = "gene", values_to = "expression")

p_vln_ggplot <- ggplot(plot_data, aes(x = celltype, y = expression, fill = treatment)) +
  geom_violin(scale = "width", trim = TRUE, position = position_dodge(width = 0.9)) +
  geom_boxplot(
    width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA
  ) +
  facet_wrap(~gene, scales = "free_y") +
  stat_compare_means(
    aes(group = treatment), method = "t.test", label = "p.signif",
    hide.ns = TRUE, size = 3, position = position_dodge(width = 0.9)
  ) +
  scale_fill_d3() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_vln_ggplot)
# ggsave("features_vln_with_stats.pdf", p_vln_ggplot, width = 18, height = 8)

# 4. Dimensional reduction plot -----------------------------------------------
DimPlot(seu, group.by = "celltype", label = TRUE, repel = TRUE) +
  scale_color_d3() +
  theme_classic()

# 5. Differential expression analysis (example) -------------------------------
# Example: comparing "MCAO" vs "Sham" in "treatment"
Idents(seu) <- seu$treatment
degs <- FindMarkers(
  seu, ident.1 = "MCAO", ident.2 = "Sham",
  test.use = "wilcox", logfc.threshold = 0.25, min.pct = 0.1
)

# Volcano plot
degs$gene <- rownames(degs)
degs$log10P <- -log10(degs$p_val_adj + 1e-300)
sig_genes <- subset(degs, p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)

p_volcano <- ggplot(degs, aes(x = avg_log2FC, y = log10P)) +
  geom_point(aes(color = avg_log2FC), alpha = 0.6) +
  scale_color_gradientn(colors = c("#141378", "#30AF95", "#F0EA07")) +
  geom_text_repel(
    data = sig_genes, aes(label = gene),
    size = 3, color = "red"
  ) +
  labs(x = "log2 Fold Change", y = "-log10 P.adjust") +
  theme_classic()

print(p_volcano)
# ggsave("volcano.pdf", p_volcano, width = 6, height = 5)

# 6. CellChat basic workflow (example) ----------------------------------------
# Create CellChat object
data_input <- GetAssayData(seu, assay = "RNA", slot = "data")
labels <- Idents(seu)
meta <- data.frame(group = labels, row.names = names(labels))
cellchat_obj <- createCellChat(object = data_input, meta = meta, group.by = "group")

# Set mouse database and subset data
CellChatDB.use <- subsetDB(CellChatDB.mouse, search = "Secreted Signaling")
cellchat_obj@DB <- CellChatDB.use
cellchat_obj <- subsetData(cellchat_obj)
cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
cellchat_obj <- projectData(cellchat_obj, PPI.mouse)

# Compute communication probability and aggregate
cellchat_obj <- computeCommunProb(cellchat_obj)
cellchat_obj <- filterCommunication(cellchat_obj, min.cells = 10)
cellchat_obj <- computeCommunProbPathway(cellchat_obj)
cellchat_obj <- netAnalysis_computeCentrality(cellchat_obj)
cellchat_obj <- aggregateNet(cellchat_obj)

# Visualization example: circle plot of interactions
group.colors <- c("#E63946","#2A9D8F","#F4A261","#264653","#E76F51")
netVisual_circle(
  cellchat_obj@net$count,
  color.use = group.colors,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of Interactions"
)

# End of script ----------------------------------------------------------------

