# ==============================================================================
# Script: 02_enrichment_GO.R
# Project: Doryteuthis pealeii Statocyst vs Gill RNA-seq
# Purpose: GO enrichment (human ortholog mapping via org.Hs.eg.db) + plots + optional GO-term heatmaps
# Author: Max Ni
# Date: 2025-07 
# ==============================================================================

# ---- reproducibility ----
set.seed(42)
options(stringsAsFactors = FALSE)
message("R version: ", R.version.string)

# ---- libraries ----
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(forcats)
  
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(GO.db)
  library(AnnotationDbi)
  
  library(ggplot2)
  library(pheatmap)
  library(matrixStats)
})


# Gene ontology enrichment analysis, find what biological processes, functions, or cellular components are overrepresented
# Because Doryteuthis pealeii is not a model organism, standard GO/KEGG tools don't support it directly. The annotation provides: 
# JGI gene symbol and combined_symbols based on closest human matches
# This means you can map squid DEGs to human gene symbols, then perform enrichment using human databases

# ==============================================================================
# 1. Load all DEGs 
# ==============================================================================
deg_file <- "results/tables/DEG.csv"

if (!file.exists(deg_file)) stop("Missing: ", deg_file, "\nRun your DESeq2 script first.")
deg_df <- read.csv(deg_file, row.names = 1, check.names = FALSE)

uniq_degs <- unique(rownames(deg_df))
message("Loaded DEGs: ", length(uniq_degs))

# ==============================================================================
# 2. Convert gene symbols to Entrez IDs
# ==============================================================================
gene_df <- bitr(uniq_degs, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # 64.55% of input gene IDs failed to map
if (nrow(gene_df) == 0) stop("No genes mapped to Entrez IDs. Check DEG rownames are human-like symbols.")

entrez_ids <- gene_df$ENTREZID
message("Mapped to Entrez IDs: ", length(entrez_ids)) # Mapped to Entrez IDs: 2860

# ==============================================================================
# 3. Run GO Enrichment
# ==============================================================================
ego <- enrichGO(
  gene = entrez_ids, # a vector of Entrez Gene IDs not gene symbols
  OrgDb = org.Hs.eg.db, # organism-specific annotation package
  ont = "BP", #GO ontology BP = Biological Process
  pAdjustMethod = "BH", # multiple testing correction method, BH = Benjamini-Hochberg false discovery rate
  qvalueCutoff = 0.05, # threshold for padj(FDR), only terms with q-value < 0.05 will be considered significant
  readable = TRUE # converts Entrez IDs in results back to gene symbols for easier interpretation
)

ego_df <- as.data.frame(ego)
write.csv(ego_df, "results/tables/GO_enrichGO_BP.csv", row.names = FALSE)
message("GO terms found: ", nrow(ego_df))

# Dotplot
go_dotplot <- dotplot(ego, showCategory = 30, title = "GO Biological Processes") + 
  labs(
    color = "FDR", # rename color legend
    size = "Gene Count" # Rename size legend
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 8), # GO term labels
    axis.text.x = element_text(size = 6), # x-axis ticks
    axis.title = element_text(size = 8), # axis titles
    legend.title = element_text(size = 4, face = "bold"),
    legend.text = element_text(size = 4)
  )

ggsave(filename = "results/plots/GO_dotplot.pdf",
       plot = go_dotplot,
       width = 7,
       height = 10,
       dpi = 300
)

# x-axis: GeneRatio = number of your genes found in a GO term / total number of DEGs submitted
# Out of all the significant genes I tested, what proportion are annotated to this GO term?

# See the genes in the top enriched GO term
ego_df$Description[1] # name of the GO term
ego_df$geneID[1] # the genes involved

ego_df[1, c("Description", "geneID")]

# Find genes for specific GO terms
ego_df$Description

# ==============================================================================
# 4. Potential pathways involved in MET in squid statocyst
# ==============================================================================
vst_file <- "results/tables/vst_mat_all_genes.csv"
stopifnot(file.exists(vst_file))
vst_mat <- read.csv(vst_file, check.names = FALSE)
head(vst_mat[, 1])

# ------------------------------------------------------------------------------
# 4.1 inner ear development pathway
# ------------------------------------------------------------------------------
# Get the genes involved in the inner ear development pathway
# Retrieve gene symbols from that exact description
inner_ear_development <- ego_df[ego_df$Description == "inner ear development", "geneID"]
# This is a character vector of length 1 

# Split the string into a list by / and convert it into a simple character vector
inner_ear_development <- unlist(strsplit(inner_ear_development, "/"))

inner_ear_development

# Find common genes between the genes in inner ear development and the gene symbols in expression matrix
genes_to_vsd_inner <- intersect(inner_ear_development, rownames(vsd_mat))

# # Extract the expression values for only those genes from the expression matrix
expr_subset_inner <- vsd_mat[genes_to_vsd_inner, ]

# Heatmap for genes expressed in the inner ear development GO term
pheatmap(
  expr_subset_inner,
  scale = "row",
  annotation_col = annotation_col,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 4,
  fontsize_col = 6,
  main = "Genes expressd in Inner Ear Development",
  filename = "results/plots/heatmap_inner_ear.pdf",
  width = 7,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------------------------
# 4.2 sensory perception of mechanical stimulus pathway
# ------------------------------------------------------------------------------
mechanical_stimulus <- ego_df[ego_df$Description == "sensory perception of mechanical stimulus", "geneID"]
mechanical_stimulus <- unlist(strsplit(mechanical_stimulus, "/"))

genes_to_vsd_mechanical <- intersect(mechanical_stimulus, rownames(vsd_mat))

expr_subset_mechanical <- vsd_mat[genes_to_vsd_mechanical, ]

pheatmap(
  expr_subset_mechanical,
  scale = "row",
  annotation_col = annotation_col,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 4,
  fontsize_col = 6,
  main = "Genes expressd in Sensory Perception of Mechanical Stimulus",
  filename = "results/plots/heatmap_sensory_mechanical.pdf",
  width = 7,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------------------------
# 4.3 sensory perception of sound pathway
# ------------------------------------------------------------------------------
sound <- ego_df[ego_df$Description == "sensory perception of sound", "geneID"]
sound <- unlist(strsplit(sound, "/"))

genes_to_vsd_sound <- intersect(sound, rownames(vsd_mat))

expr_subset_sound <- vsd_mat[genes_to_vsd_sound, ]

pheatmap(
  expr_subset_sound,
  scale = "row",
  annotation_col = annotation_col,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 4,
  fontsize_col = 6,
  main = "Genes expressd in Sensory Perception of Sound",
  filename = "results/plots/heatmap_sensory_sound.pdf",
  width = 7,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------------------------
# 4.4 calcium ion transmembrane transport pathway
# ------------------------------------------------------------------------------
calcium_transport <- ego_df[ego_df$Description == "calcium ion transmembrane transport", "geneID"]
calcium_transport <- unlist(strsplit(calcium_transport, "/"))

genes_to_vsd_calcium <- intersect(calcium_transport, rownames(vsd_mat))

expr_subset_calcium <- vsd_mat[genes_to_vsd_calcium, ]

pheatmap(
  expr_subset_calcium,
  scale = "row",
  annotation_col = annotation_col,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 4,
  fontsize_col = 6,
  main = "Genes expressd in Calcium Ion Transmembrane Transport",
  filename = "results/plots/heatmap_Ca_transport.pdf",
  width = 7,
  height = 7,
  dpi =300
)

# ==============================================================================
# 5. Find which pathways a gene is involved in 
# ==============================================================================

# ------------------------------------------------------------------------------
# 5.1 Go terms involved in PKD2
# ------------------------------------------------------------------------------
# Convert gene symbol to Entrez ID
gene_symbol_pkd2 <- "PKD2"
entrez_id_pkd2 <- bitr(gene_symbol_pkd2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

# Find GO terms annotated to this gene
columns(org.Hs.eg.db)

go_annotations_pkd2 <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = entrez_id_pkd2,
    columns = c("GO", "ONTOLOGY", "SYMBOL"),
    keytype = "ENTREZID"
  )
) 

go_annotations_pkd2

# Get GO term descriptions, remove duplicated GO IDs

go_terms_pkd2 <- AnnotationDbi::select(
  GO.db,
  keys = unique(go_annotations_pkd2$GO),
  columns = "TERM",
  keytype = "GOID"
) %>%
  distinct(GOID, TERM)

go_terms_pkd2

# Merge GO annotations and term names
go_merged_unique_pkd2 <- go_annotations_pkd2 %>%
  distinct(ENTREZID, GO, ONTOLOGY) %>%
  left_join(go_terms_pkd2, by = c("GO" = "GOID"))

# Optional: count frequency if some terms repeat
# go_counts_unique <- go_merged_unique %>%
#  count(TERM, ONTOLOGY) %>%
#  arrange(desc(n))

go_plot_pkd2 <- go_merged_unique_pkd2 %>%
  filter(!is.na(TERM)) %>%
  mutate(dummy = 1)

# Add GO IDs
go_plot_pkd2 <- go_plot_pkd2 %>%
  mutate(label = paste0(GO, ": ", TERM))

# Barplot
go_pkd2 <- ggplot(go_plot_pkd2, aes(x = label, y = dummy, fill = ONTOLOGY)) +
  geom_col() +
  coord_flip() +
  labs(
    title = paste("GO Terms for", gene_symbol_pkd2),
    x = NULL,
    y = NULL,
    fill = "Ontology"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(size = 3),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.2, "cm"),
    plot.margin = margin(t = 5.5, r = 5.5, b = 25, l = 5.5)
    
  )

ggsave(filename = "results/plots/go_PKD2_refined.pdf",
       plot = go_pkd2,
       width = 7, 
       height = 7,
       dpi = 300)
# ------------------------------------------------------------------------------
# 5.2 Go terms involved in PKD1L2
# ------------------------------------------------------------------------------
# Convert gene symbol to Entrez ID
gene_symbol_pkd1l2 <- "PKD1L2"
entrez_id_pkd1l2 <- bitr(gene_symbol_pkd1l2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

# Find GO terms annotated to this gene
columns(org.Hs.eg.db)

go_annotations_pkd1l2 <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = entrez_id_pkd1l2,
    columns = c("GO", "ONTOLOGY", "SYMBOL"),
    keytype = "ENTREZID"
  )
) 

go_annotations_pkd1l2

# Get GO term descriptions, remove duplicated GO IDs

go_terms_pkd1l2 <- AnnotationDbi::select(
  GO.db,
  keys = unique(go_annotations_pkd1l2$GO),
  columns = "TERM",
  keytype = "GOID"
) %>%
  distinct(GOID, TERM)

go_terms_pkd1l2

# Merge GO annotations and term names
go_merged_unique_pkd1l2 <- go_annotations_pkd1l2 %>%
  distinct(ENTREZID, GO, ONTOLOGY) %>%
  left_join(go_terms_pkd1l2, by = c("GO" = "GOID"))

# Optional: count frequency if some terms repeat
# go_counts_unique_pkd1l2 <- go_merged_unique_pkd1l2 %>%
#   count(TERM, ONTOLOGY) %>%
#   arrange(desc(n))


go_plot_pkd1l2 <- go_merged_unique_pkd1l2 %>%
  filter(!is.na(TERM)) %>%
  mutate(dummy = 1)

# Add GO IDs
go_plot_pkd1l2 <- go_plot_pkd1l2 %>%
  mutate(label = paste0(GO, ": ", TERM))

# Barplot
go_pkd1l2 <- ggplot(go_plot_pkd1l2, aes(x = label, y = dummy, fill = ONTOLOGY)) +
  geom_col(width = 0.6) +
  coord_flip() +
  labs(
    title = paste("GO Terms for", gene_symbol_pkd1l2),
    x = NULL,
    y = NULL,
    fill = "Ontology"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(size = 4),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.2, "cm"),
    plot.margin = margin(t = 5.5, r = 5.5, b = 25, l = 5.5)
  )

ggsave(filename = "results/plots/go_PKD1L2_refined.pdf",
       plot = go_pkd1l2,
       width = 5, 
       height = 1.5,
       dpi = 300)


# ==============================================================================
# 6. Expression of Nompc (ANK1 and ANK1_2)
# ==============================================================================
# Check if ANK1 and ANK1_2 appear in vsd_mat
"ANK1" %in% rownames(vsd_mat) # TRUE
"ANK1_2" %in% rownames(vsd_mat) # TRUE
# This means that these two genes are DEGs

# Find its exact row
which(rownames(vsd_mat) == "ANK1") # 334
which(rownames(vsd_mat) == "ANK1_2") # 335

# Extract the row data
vsd_mat["ANK1", ]
# Statocyst 1 Statocyst 2      Gill 1      Gill 2 
#   6.785918    6.986949   10.901358   10.738620 
vsd_mat["ANK1_2", ]
# Statocyst 1 Statocyst 2      Gill 1      Gill 2 
#   6.785918    6.785918   10.816405   10.657601

# Check if ANK1 appear in GO
"ANK1" %in% gene_df$SYMBOL # returns TRUE, so ANK1 appears in GO term
"ANK1_2" %in% gene_df$SYMBOL # returns FALSE

# Find the exact row of ANK1
which(gene_df$SYMBOL == "ANK1") # 86

# Extract the row data of ANK1
gene_df[gene_df$SYMBOL == "ANK1", ] 
#     SYMBOL ENTREZID
# 149   ANK1      286

# Convert gene symbol to Entrez ID
gene_symbol_ank1 <- "ANK1"

entrez_id_ank1 <- bitr(gene_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

go_annotations_ank1 <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = entrez_id_ank1,
    columns = c("GO", "ONTOLOGY", "SYMBOL"),
    keytype = "ENTREZID"
  )
) 

go_annotations_ank1

# get GO term descriptions, remove duplicated GO IDs
go_terms_ank1 <- AnnotationDbi::select(
  GO.db,
  keys = unique(go_annotations_ank1$GO),
  columns = "TERM",
  keytype = "GOID"
) %>%
  distinct(GOID, TERM)

go_terms_ank1

# Merge GO annotations and term names
go_merged_unique_ank1 <- go_annotations_ank1 %>%
  distinct(ENTREZID, GO, ONTOLOGY) %>%
  left_join(go_terms_ank1, by = c("GO" = "GOID"))

# Optional: count frequency if some terms repeat
go_counts_unique_ank1 <- go_merged_unique_ank1 %>%
  count(TERM, ONTOLOGY) %>%
  arrange(desc(n))

# Barplot
go_ank1 <- ggplot(go_counts_unique_ank1, aes(x = reorder(TERM, n), y = n, fill = ONTOLOGY)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = paste("Unique GO Terms for", gene_symbol_ank1),
    x = "GO Term",
    y = "Count",
    fill = "Ontology"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 4),
    axis.text.y = element_text(size = 4),
    axis.title = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  )

ggsave(filename = "results/plots/go_ANK1.pdf",
       plot = go_ank1,
       width = 7,
       height = 7,
       dpi = 300)

# ==============================================================================
# 7. Additional analysis of genes of interest
# ==============================================================================

# ------------------------------------------------------------------------------
# 7.1 Compare ANK1 expression levels between conditions
# ------------------------------------------------------------------------------
# Extract ANK1 expression
gene <- "ANK1"
expr_ank1 <- vsd_mat[gene, ]

# Combine with condition labels
df_expr_ank1 <- data.frame(
  expression = as.numeric(expr_ank1),
  condition = coldata$condition
)

# Plot ANK1 expression by condition
expr_ank1 <- ggplot(df_expr_ank1, aes(x = condition, y = expression, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = paste("Expression of", gene), y = "Normalized Expression") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(filename = "results/plots/expression_ank1_condition.pdf",
       plot = expr_ank1,
       width = 6,
       height = 7,
       dpi =300)

# ------------------------------------------------------------------------------
# 7.2 Compare PKD2 expression between conditions
# ------------------------------------------------------------------------------
# Extract PKD2 expression
gene <- "PKD2"
expr_expr_pkd2 <- vsd_mat[gene, ]

# Combine with condition labels
df_expr_pkd2 <- data.frame(
  expression = as.numeric(expr_expr_pkd2),
  condition = coldata$condition)

# plot expression by condition
expr_pkd2 <- ggplot(df_expr_pkd2, aes(x = condition, y = expression, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  labs(title = paste("Expression of", gene), y = "vst normalized expression") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(filename = "results/plots/expression_pkd2_condition.pdf",
       plot = expr_pkd2,
       width = 6,
       height = 7,
       dpi = 300)


# ==============================================================================
# run GO for genes upregulated in statocyst
# ==============================================================================
# ----------------------
# 1. load up genes
# ----------------------
top_up_file <- "results/tables/sig_gene_up.csv"

if (!file.exists(top_up_file)) stop("Missing: ", top_up_file, "\nRun your DESeq2 script first.")
top_up_df <- read.csv(top_up_file, row.names = 1, check.names = FALSE)

uniq_top_up <- unique(rownames(top_up_df))
message("Loaded top up: ", length(uniq_top_up))

# ---------------------------------------
# 2. Convert gene symbols to Entrez IDs
# ---------------------------------------
top_up_df <- bitr(uniq_top_up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # 70.9% of input gene IDs are fail to map...
if (nrow(top_up_df) == 0) stop("No genes mapped to Entrez IDs. Check DEG rownames are human-like symbols.")

entrez_ids_up <- top_up_df$ENTREZID
message("Mapped to Entrez IDs: ", length(entrez_ids_up)) # Mapped to Entrez IDs: 1181

# ----------------------
# 3. Run GO Enrichment
# ----------------------
ego_up <- enrichGO(
  gene = entrez_ids_up, # a vector of Entrez Gene IDs not gene symbols
  OrgDb = org.Hs.eg.db, # organism-specific annotation package
  ont = "BP", #GO ontology BP = Biological Process
  pAdjustMethod = "BH", # multiple testing correction method, BH = Benjamini-Hochberg false discovery rate
  qvalueCutoff = 0.05, # threshold for padj(FDR), only terms with q-value < 0.05 will be considered significant
  readable = TRUE # converts Entrez IDs in results back to gene symbols for easier interpretation
)

class(ego_up)

ego_up_df <- as.data.frame(ego_up)

write.csv(ego_up_df, "results/tables/GO_enrichGO_up_BP.csv", row.names = FALSE)
message("GO terms found: ", nrow(ego_up_df)) # GO terms found: 227

# Dotplot
go_up_dotplot <- dotplot(ego_up, showCategory = 30, title = "GO Biological Processes in Statocyst") + 
  labs(
    color = "FDR", # rename color legend
    size = "Gene Count" # Rename size legend
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 8), # GO term labels
    axis.text.x = element_text(size = 6), # x-axis ticks
    axis.title = element_text(size = 8), # axis titles
    legend.title = element_text(size = 4, face = "bold"),
    legend.text = element_text(size = 4)
  )

ggsave(filename = "results/plots/GO_up_dotplot.pdf",
       plot = go_up_dotplot,
       width = 7,
       height = 10,
       dpi = 300
)

# x-axis: GeneRatio = number of your genes found in a GO term / total number of DEGs submitted
# Out of all the significant genes I tested, what proportion are annotated to this GO term?

# See the genes in the top enriched GO term
ego_up_df$Description[1] # name of the GO term
ego_up_df$geneID[1] # the genes involved

ego_up_df[1, c("Description", "geneID")]

# Find genes for specific GO terms
ego_up_df$Description


# ==============================================================================
# run GO for genes upregulated in gill
# ==============================================================================
# ----------------------
# 1. load down genes
# ----------------------
top_down_file <- "results/tables/sig_gene_down.csv"

if (!file.exists(top_down_file)) stop("Missing: ", top_down_file, "\nRun your DESeq2 script first.")
top_down_df <- read.csv(top_down_file, row.names = 1, check.names = FALSE)

uniq_top_down <- unique(rownames(top_down_df))
message("Loaded top down: ", length(uniq_top_down)) # Loaded top down: 4009

# ---------------------------------------
# 2. Convert gene symbols to Entrez IDs
# ---------------------------------------
top_down_df <- bitr(uniq_top_down, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # 58.12% of input gene IDs are fail to map...
if (nrow(top_down_df) == 0) stop("No genes mapped to Entrez IDs. Check DEG rownames are human-like symbols.")

entrez_ids_down <- top_down_df$ENTREZID
message("Mapped to Entrez IDs: ", length(entrez_ids_down)) # Mapped to Entrez IDs: 1679

# ----------------------
# 3. Run GO Enrichment
# ----------------------
ego_down <- enrichGO(
  gene = entrez_ids_down, # a vector of Entrez Gene IDs not gene symbols
  OrgDb = org.Hs.eg.db, # organism-specific annotation package
  ont = "BP", #GO ontology BP = Biological Process
  pAdjustMethod = "BH", # multiple testing correction method, BH = Benjamini-Hochberg false discovery rate
  qvalueCutoff = 0.05, # threshold for padj(FDR), only terms with q-value < 0.05 will be considered significant
  readable = TRUE # converts Entrez IDs in results back to gene symbols for easier interpretation
)

class(ego_down)

ego_down_df <- as.data.frame(ego_down)

write.csv(ego_down_df, "results/tables/GO_enrichGO_down_BP.csv", row.names = FALSE)
message("GO terms found: ", nrow(ego_down_df)) # GO terms found: 460

# Dotplot
go_down_dotplot <- dotplot(ego_down, showCategory = 30, title = "GO Biological Processes in Gill") + 
  labs(
    color = "FDR", # rename color legend
    size = "Gene Count" # Rename size legend
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 8), # GO term labels
    axis.text.x = element_text(size = 6), # x-axis ticks
    axis.title = element_text(size = 8), # axis titles
    legend.title = element_text(size = 4, face = "bold"),
    legend.text = element_text(size = 4)
  )

ggsave(filename = "results/plots/GO_down_dotplot.pdf",
       plot = go_down_dotplot,
       width = 7,
       height = 10,
       dpi = 300
)

# ==============================================================================
# 8. Save session info
# ==============================================================================
writeLines(capture.output(sessionInfo()), file.path("results/session/", "session_info_GO.txt"))
message("Done. Outputs in: ", normalizePath("results/session"))
