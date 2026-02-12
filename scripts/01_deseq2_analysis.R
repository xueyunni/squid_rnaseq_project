# ==============================================================================
# Project: Doryteuthis pealeii Statocyst vs Gill RNA-seq 
# Author: Max Ni
# Date: 2025-07
# ==============================================================================

# -----Reproducibility-------
set.seed(42) # Fixes randomness in PCA/plotting
options(stringsAsFactors = FALSE)

# Record environment info (nice for GitHub)
message("R version: ", R.version.string)

# ---- libraries ----
suppressPackageStartupMessages({
  library(DESeq2) # for DGE analysis using a negative binomial distribution
  library(tibble) # reimagining of a dataframe, also provided by the tidyverse package
  library(dplyr) # a core package in the tidyverse,
  library(tidyr)
  library(rtracklayer) # for importing, exporting, and manipulating genomic data, esp for working with genome annotation files like BED, GFF
  library(ggplot2) # for data visualization, part of the tidyverse. You build plots by layering components like data, aesthetics, and geometirc objects
  library(ggrepel) # an extension to ggplot2 that helps you add text labels to plots without overlap
  library(EnhancedVolcano) # for making publication-ready volcano plots
  library(pheatmap) # for making pretty heatmaps with clustering and annotations
  library(matrixStats)
  library(clusterProfiler) # for functional enrichment analysis and visualization of gene clusters or gene list, widely used after DE to interpret biological meaning
  library(enrichplot) # a companion package to clusterProfiler, focusing on visualization of enrichment results (GO, KEGG, GSEA) from packages like clusterProfiler
  library(org.Hs.eg.db) # a bioconductor annotation package for human, Dpe genome was annotated against the human genome
  library(GO.db) # another bioconductor annotation package, providing the GO database itself
  library(AnnotationDbi)
  library(forcats)
  library(ggiraph) # for clusterprofile
})

# ==============================================================================
# 1. Load raw counts and clean
# ==============================================================================
# Read in the count matrix
counts_raw <- read.table("results/counts/gene_id_trimmed.counts", header = TRUE, row.names = 1, comment.char = "#", check.names = FALSE)
# File path
# First row contains column names
# Use the first column as row names
# Ignore lines starting with "#"
# Don't modify column names to be syntactically valid

# View the raw count matrix
class(counts_raw)
head(counts_raw)
colnames(counts_raw)

# Data cleaning
# Keep only the count columns
counts <- counts_raw[, 6:ncol(counts_raw)] # Start from the 6th column to the last column
# ncol() returns the total number of columns

# view the refined count matrix
head(counts)
class(counts)
colnames(counts)

# Convert gene ids into gene symbols
rownames(counts) <- sub("_Dpe.*", "", rownames(counts))
rownames(counts)

# view the count matrix
head(counts)

# ==============================================================================
# 2. DESeq2 object setup
# ==============================================================================
# Create the metadata(or coldata)
coldata <- data.frame(
  row.names = colnames(counts), # Use sample/column names as row names of coldata
  condition = c("Statocyst", "Statocyst", "Gill", "Gill")
)

# check if row names of coldata same as column names of count matrix
all(rownames(coldata) == colnames(counts))

# DE analysis
# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition) # The formula that tells DESeq2 what to compare

# Run the DESeq analysis
# Fit the model and estimate dispersion, automatic normalization, DESeq2 uses size factor normalization
dds <- DESeq(dds)
# Performs the full DESeq2 workflow on your dataset
# Estimate size factor - normalize for sequencing depth
# Estimate dispersion - variance of counts across replicates
# Fit the negative binomial GLM - according to your design
# Perform wald test or LRT - to compute DE
# Returns the updated DESeqDataSet ready for results extraction

# save DESeq2 object

dir.create("results/rds", recursive = TRUE, showWarnings = FALSE)

saveRDS(dds, "results/rds/dds_statocyst_vs_gill.rds")

# Extract the DE results table
res <- results(dds)
res
# Independent filtering 
# DESeq2 checks the mean normalized counts for each gene
# Genes below a certain threshold (mean count < 1 by default) are filtered out before testing, so they don't get p-values

# ------------------------------------------------------------------------------
# 2.1 Visualize gene dispersion
# ------------------------------------------------------------------------------
# Visualize the results in MA plot
pdf("results/plots/MA_plot.pdf", width = 7, height = 6)

plotMA(res,
       cex = 0.6, # Point size
       main = "MA plot",
       xlab = "Mean of Normalized Counts",
       ylab = "Log2 Fold Change",
       alpha = 0.05)

dev.off()

print("MA plot saved successfully.")

# Check the reference level!
levels(coldata$condition) # NULL

# Convert condition to a factor with the correct reference - Statocyst vs Gill, reference level = Gill
coldata$condition <- factor(coldata$condition, levels = c("Gill", "Statocyst"))

# Recheck the reference level
levels(coldata$condition) # [1] "Gill"      "Statocyst"

# Update the dds object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

# Rerun the DE analysis
dds <- DESeq(dds)

# Extract the DE results table
res <- results(dds) 

# Check the results summary
summary(res)

# confirm the contrast
resultsNames(dds)

# Replot MA 
pdf("results/plots/MA_plot.pdf", width = 7, height = 6)

plotMA(res,
       cex = 0.6, # Point size
       main = "MA plot",
       xlab = "Mean of Normalized Counts",
       ylab = "Log2 Fold Change",
       alpha = 0.05)

dev.off()

print("MA plot saved successfully.")

# ==============================================================================
# 3. Normalization and transformation
# ==============================================================================
# vst() normalization
# VST - variance stabilizing transformation, a normalization method in DESeq2, log-like scale, suitable for PCA, heatmaps, clustering
vsd <- vst(dds, blind = FALSE) # blind = FALSE, tells the function to use the experimental design(your ~ condition) when stabilizing variance
head(assay(vsd))

# ------------------------------------------------------------------------------
# 3.1 PCA
# ------------------------------------------------------------------------------
# PCA captures global patterns in the data using all genes
# PCA to confirm sample group separation
pca_plot <- plotPCA(vsd, intgroup = "condition") + 
  geom_point(size = 2) +
  ggtitle("PCA") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )

ggsave(filename = "results/plots/PCA_plot.pdf",
       plot = pca_plot,
       width = 7,
       height = 4
         )

print("PCA plot saved successfully.")

# ------------------------------------------------------------------------------
# 3.2 Heatmap
# ------------------------------------------------------------------------------
# Extract the matrix of vst-normalized counts
vsd_mat <- assay(vsd) # rows = geneID, columns = samples
class(vsd_mat) # [1] "matrix" "array" 
head(vsd_mat)

# Check the matrix
summary(vsd_mat)
any(is.na(vsd_mat))       # FALSE
any(is.infinite(vsd_mat))  # FLASE

# Change column names (Sample names)
colnames(vsd_mat) <- c("Statocyst 1", "Statocyst 2", "Gill 1", "Gill 2")

# Save 
write.csv(vsd_mat, "results/tables/vst_mat_all_genes.csv")

# Sample annotation
annotation_col <- as.data.frame(colData(dds)[, "condition", drop = FALSE])
rownames(annotation_col) <- colnames(vsd_mat)

# Compute row variances
row_var <- rowVars(vsd_mat)

# Keep rows with variance > 0
vsd_mat_filtered <- vsd_mat[row_var > 0, ]

# plot heatmap
pheatmap(
  vsd_mat_filtered,
  scale = "row",
  annotation_col = annotation_col,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_row = 0.5,
  main = "Heatmap of All Genes",
  filename = "results/plots/heatmap_all_genes.pdf",
  width = 7,
  height =7,
  dpi = 300
)

print("heatmap for all genes saved successfully.")

# ==============================================================================
# 4. Differential Expression analysis
# ==============================================================================
# Filter significant genes with padj < 0.05 & log2FC > 1
sig_genes <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
sig_genes # 8068 genes 

# Export sig_genes
write.csv(as.data.frame(sig_genes), "results/tables/DEG.csv")

print("DEG.csv saved successfully.")

# ------------------------------------------------------------------------------
# 4.1 Volcano plot
# ------------------------------------------------------------------------------
# Volcano plot for DGEs with top 100 being labeled
# Check NAs
any(is.na(sig_genes$log2FoldChange))  # FALSE
any(is.na(sig_genes$padj)) # FALSE

# Convert sig_genes to a dataframe
sig_genes_df <- as.data.frame(sig_genes)

# Highlight upregulated vs downregulated genes in different colors
sig_genes_df_up_down <- sig_genes_df %>%
  mutate(regulation = case_when(
    padj < 0.05 & log2FoldChange > 1 ~ "Up",
    padj < 0.05 & log2FoldChange < -1 ~ "Down",
    TRUE ~ "Not Significant")
  )

head(sig_genes_df_up_down)

# Filter the most significant genes
# top up-regulated (log2FC > 1) with padj < 0.05
top_up <- sig_genes_df %>% filter(log2FoldChange > 1, padj < 0.05) %>% arrange(desc(log2FoldChange))
top_up

# save top_up for GO analysis
write.csv(top_up, "results/tables/sig_gene_up.csv")

# top down-regulated
top_down <- sig_genes_df %>% filter(log2FoldChange < 1, padj < 0.05) %>% arrange(log2FoldChange)

# save top_down for GO analysis
write.csv(top_down, "results/tables/sig_gene_down.csv")

# Pick top 100 genes
top100_sig_genes <- sig_genes_df_up_down %>% 
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice_head(n = 100)

head(top100_sig_genes)

top100_sig_genes_volcano_plot <- ggplot(sig_genes_df_up_down, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = regulation), alpha = 0.7) +
  scale_color_manual(values = c("Up" = "firebrick", "Down" = "steelblue", "Not Significant" = "gray80")) +
  geom_text_repel(
    data = top100_sig_genes,
    aes(label = rownames(top100_sig_genes)),
    size = 1.5,
    max.overlaps = Inf) +
  theme_minimal() + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) +
  labs(
    title = "Volcano Plot of DEGs",
    x = "Log2 Fold Change",
    y = "-log10 padj",
    color = "Regulation"
  )

ggsave(filename = "results/plots/volcano_top100_DEGs.pdf", 
       plot = top100_sig_genes_volcano_plot,
       width = 7,
       height = 10)

print("Volcano plot for Top 100 sig genes saved successfully.")

# ------------------------------------------------------------------------------
# 4.2 EnhancedVolcano plot
# ------------------------------------------------------------------------------
# Enhanced volcano plot with top 100 genes labeled by significance (padj < 0.05) & log2FC > 1
# Remove NAs first to avoid errors
res_noNA <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]

# Order by padj first, then by absolute log2FC (descending)
res_ordered <- res_noNA[order(res_noNA$padj, -abs(res_noNA$log2FoldChange)), ]

# View the ordered results table
head(res_ordered, 5)

# Get top 100 genes 
top100_genes <- rownames(res_ordered)[1:100]

# Make Enhanced Volcano plot
enhancedvolcano_plot <- EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',         # column for fold change
                y = 'padj',                   # column for adjusted p-value
                xlab = bquote(~Log[2]~ 'FoldChange'),
                ylab = bquote(~-Log[10]~ 'padj'),
                pCutoff = 0.05,               # significance threshold
                FCcutoff = 1.0,               # log2FC threshold
                pointSize = 1.0,
                labSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.1,
                arrowheads = FALSE,
                colAlpha = 0.8,               # transparency
                legendLabels=c('NS','Log2FC','padj','padj & Log2FC'),
                legendPosition = 'right',
                title = 'Differential Expression Volcano',
                selectLab = top100_genes, # highlight specific genes
                # boxedLabels = TRUE,
                max.overlaps = Inf
) + 
  theme_minimal() + # use a clean background
  theme(
    panel.grid = element_blank(),   # remove grid lines
    # panel.border = element_blank(), # optional: remove border
    axis.ticks = element_line(),
    axis.line = element_line(), # keep axis lines visible
    plot.title = element_text(hjust = 0.5, face = "bold")   # centers title
  )

ggsave(filename = "results/plots/EnhancedVolcano.pdf",
       plot = enhancedvolcano_plot,
       width = 7, 
       height = 10)

print("EnhancedVolcano plot saved successfully.")

# ------------------------------------------------------------------------------
# 4.3 Heatmap for DEGs
# ------------------------------------------------------------------------------
# Only keep rows with matching geneIDs
keep_ids <- intersect(rownames(sig_genes_df), rownames(vsd_mat))
keep_ids
# Subset vst expression matrix
vsd_mat_sig <- vsd_mat[keep_ids, ]
head(vsd_mat_sig)

# Save 
write.csv(vsd_mat_sig, "results/tables/vst_mat_DEGs.csv")

# plot heatmap
pheatmap(
  vsd_mat_sig,
  scale = "row",
  annotation_col = annotation_col,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_row = 0.5,
  main = "Heatmap of DEGs",
  filename = "results/plots/heatmap_DEGs.pdf",
  width = 7,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------------------------
# 4.4 Heatmap for top 200 DEGs
# ------------------------------------------------------------------------------
top200_ids <- sig_genes_df %>%
  arrange(padj) %>%
  slice_head(n = 200) 

top200_vsd_mat_sig <- vsd_mat_sig[intersect(rownames(top200_ids), rownames(vsd_mat_sig)), ]

# Save
write.csv(top200_vsd_mat_sig, file = "results/tables/vst_DEGs_top200.csv")

pheatmap(
  top200_vsd_mat_sig, 
  scale = "row",
  annotation_col = annotation_col,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 1,
  main = "Heatmap of Top 200 DEGs",
  filename = "results/plots/heatmap_DEGs_top200.pdf",
  width = 7,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------------------------
# 4.5 Heatmap for top 100 DEGs
# ------------------------------------------------------------------------------
top100_ids <- sig_genes_df %>%
  arrange(padj) %>%
  slice_head(n = 100) 

top100_vsd_mat_sig <- vsd_mat_sig[intersect(rownames(top100_ids), rownames(vsd_mat_sig)), ]

# save
write.csv(top100_vsd_mat_sig, file = "results/tables/vst_DEGs_top100.csv")

pheatmap(
  top100_vsd_mat_sig, 
  scale = "row",
  annotation_col = annotation_col,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 2,
  main = "Heatmap of Top 100 DEGs",
  filename = "results/plots/heatmap_DEGs_top100.pdf",
  width = 7,
  height =7,
  dpi = 300
)

# ==============================================================================
# 5. Export results
# ==============================================================================
write.csv(as.data.frame(res), "results/tables/deseq2_results.csv")

# ==============================================================================
# 6. Save session info
# ==============================================================================
dir.create("results/session", showWarnings = FALSE, recursive = TRUE)
writeLines(capture.output(sessionInfo()), file.path("results/session/", "session_info_DE.txt"))
message("Done. Outputs in: ", normalizePath("results/session"))