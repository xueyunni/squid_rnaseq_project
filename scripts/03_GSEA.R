# ==============================================================================
# Script: 03_GSEA.R
# Project: Statocyst vs Gill RNA-seq (DESeq2)
# Purpose: GSEA using ranked DESeq2 statistics
# Author: Max Ni
# Date: 2026-01
# ==============================================================================

# ---- reproducibility ----
set.seed(42)
options(stringsAsFactors = FALSE)
message("R version: ", R.version.string)

# ---- libraries ----

install.packages("readr")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  
  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)   # change if using another organism DB
})

# ---- paths ----
dds_rds <- "results/rds/dds_statocyst_vs_gill.rds"   
out_dir_go <- "results/GSEA_GO"
out_dir_kegg <- "results/GSEA_KEGG"
dir.create(out_dir_go, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_kegg, recursive = TRUE, showWarnings = FALSE)


# ---- load DESeq2 object ----
stopifnot(file.exists(dds_rds))
dds <- readRDS(dds_rds)

# ---- DESeq2 results (Wald stats) ----
res <- results(dds)

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  filter(!is.na(stat)) %>%          # stat is required for ranking
  filter(!is.na(pvalue))            # optional sanity filter
message("Genes with non-NA stat: ", nrow(res_df)) # Genes with non-NA stat: 21551

colnames(res_df)
# [1] "gene"           "baseMean"       "log2FoldChange" "lfcSE"          "stat"           "pvalue"         "padj"   

# convert gene symbols to entrez ids
symbol_entrezid <- bitr(
  res_df$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
# 64.13% of input gene IDs are fail to map...

# merge ids back into DESeq2 results
res_df2 <- merge(res_df, symbol_entrezid, by.x = "gene", by.y = "SYMBOL")
head(res_df2)

# ---- build ranked list for GSEA ----
# Recommended ranking for GSEA: Wald statistic (stat), not log2FC
gene_rank <- res_df2$stat
class(gene_rank) # [1] "numeric"

names(gene_rank) <- res_df2$ENTREZID

gene_rank <- sort(gene_rank, decreasing = TRUE)

# Optional: remove duplicated gene IDs (keep the one with largest |stat|)
# Useful if ortholog mapping created duplicates.
if (any(duplicated(names(gene_rank)))) {
  message("Detected duplicated gene IDs in ranking list. Collapsing duplicates by max |stat|...")
  gene_rank <- tapply(gene_rank, names(gene_rank), function(x) x[which.max(abs(x))])
  gene_rank <- sort(gene_rank, decreasing = TRUE)
}

summary(gene_rank)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -45.4484  -3.1322  -0.5977  -0.6741   1.7407  45.2861 

sum(duplicated(gene_rank)) # 98

# ---- run GSEA: GO Biological Process ----
# keyType: set to the ID type in names(gene_rank)
# Common: "SYMBOL" or "ENTREZID"
key_type_used <- "ENTREZID"   # <- CHANGE if needed (e.g., "ENTREZID")

ego_gsea <- gseGO(
  geneList     = gene_rank,
  OrgDb        = org.Hs.eg.db,
  keyType      = key_type_used,
  ont          = "BP",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.001,
  verbose      = FALSE
)

# ---- save result table ----
gsea_go_tbl <- as.data.frame(ego_gsea@result)
nrow(gsea_go_tbl) # 94
colnames(gsea_go_tbl)
#  [1] "ID"              "Description"     "setSize"         "enrichmentScore" "NES"             "pvalue"         
# [7] "p.adjust"        "qvalue"          "rank"            "leading_edge"    "core_enrichment"

write.csv(gsea_go_tbl, file.path(out_dir_go, "GSEA_GO_BP_table.csv"))
message("Saved: ", file.path(out_dir_go, "GSEA_GO_BP_table.csv"))

# ---- quick filtering helpers ----
# search for mechanosensory / ear / cilia terms
keyword_hits <- gsea_go_tbl %>%
  filter(grepl("ear|mechan|cilium|axoneme|sensory|synap", Description, ignore.case = TRUE)) %>%
  arrange(p.adjust)

write.csv(keyword_hits, file.path(out_dir_go, "GSEA_GO_BP_keywords.csv"))

# ---- dotplot (top terms) ----
p_dot <- dotplot(ego_gsea, showCategory = 30) +
  scale_y_discrete(labels = function(x) gsub("\n", " ", x)) +   # force one line per GO term
  theme(
    axis.text.y = element_text(size = 9),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 35)
  )

ggsave(
  filename = file.path(out_dir_go, "GSEA_GO_BP_dotplot.pdf"),
  plot     = p_dot,
  width    = 10,
  height   = 7
)

message("Saved: ", file.path(out_dir_go, "GSEA_GO_BP_dotplot.pdf"))

# ---- enrichment curve for a chosen term ----
# cilium organization
if (nrow(gsea_go_tbl) > 0) {
  top_id_1 <- gsea_go_tbl$ID[gsea_go_tbl$Description == "cilium organization"]
  desc_1 <- gsea_go_tbl$Description[gsea_go_tbl$Description == "cilium organization"]
  p_curve <- gseaplot2(ego_gsea, geneSetID = top_id_1, title = desc_1)
  
  ggsave(
    filename = file.path(out_dir_go, paste0("GSEA_curve_", top_id_1, ".pdf")),
    plot     = p_curve,
    width    = 9,
    height   = 5
  )
  message("Saved enrichment curve for: ", top_id_1)
}

# microtubule-based movement
if (nrow(gsea_go_tbl) > 0) {
  top_id_2 <- gsea_go_tbl$ID[gsea_go_tbl$Description == "microtubule-based movement"]
  desc_2 <- gsea_go_tbl$Description[gsea_go_tbl$Description == "microtubule-based movement"]
  p_curve <- gseaplot2(ego_gsea, geneSetID = top_id_2, title = desc_2)
  
  ggsave(
    filename = file.path(out_dir_go, paste0("GSEA_curve_", top_id_2, ".pdf")),
    plot     = p_curve,
    width    = 9,
    height   = 5
  )
  message("Saved enrichment curve for: ", top_id_2)
}

# cilium assembly
if (nrow(gsea_go_tbl) > 0) {
  top_id_3 <- gsea_go_tbl$ID[gsea_go_tbl$Description == "cilium assembly"]
  desc_3 <- gsea_go_tbl$Description[gsea_go_tbl$Description == "cilium assembly"]
  p_curve <- gseaplot2(ego_gsea, geneSetID = top_id_3, title = desc_3)
  
  ggsave(
    filename = file.path(out_dir_go, paste0("GSEA_curve_", top_id_3, ".pdf")),
    plot     = p_curve,
    width    = 9,
    height   = 5
  )
  message("Saved enrichment curve for: ", top_id_3)
}

# cilium movement 
if (nrow(gsea_go_tbl) > 0) {
  top_id_4 <- gsea_go_tbl$ID[gsea_go_tbl$Description == "cilium movement"]
  desc_4 <- gsea_go_tbl$Description[gsea_go_tbl$Description == "cilium movement"]
  p_curve <- gseaplot2(ego_gsea, geneSetID = top_id_4, title = desc_4)
  
  ggsave(
    filename = file.path(out_dir_go, paste0("GSEA_curve_", top_id_4, ".pdf")),
    plot     = p_curve,
    width    = 9,
    height   = 5
  )
  message("Saved enrichment curve for: ", top_id_4)
}

message("Done.")

# -------extract the top BPs in GO----------------
target_terms <- c(
  "cilium organization",
  "microtubule-based movement",
  "cilium assembly",
  "pattern specification process",
  "cilium movement",
  "exocytosis",
  "neurotransmitter transport",
  "vesicle-mediated transport in synapse",
  "cilium movement involved in cell motility",
  "cilium or flagellum-dependent cell motility",
  "cilium-dependent cell motility",
  "potassium ion transport",
  "regulation of microtubule−based process",
  "synaptic vesicle cycle",
  "regulated exocytosis",
  "microtubule−based transport",
  "microtubule bundle formation",
  "flagellated sperm motility",
  "sperm motility",
  "neurotransmitter secretion",
  "signal release from synapse",
  "regulation of exocytosis",
  "axoneme assembly",
  "synaptic vesicle exocytosis",
  "motile cilium assembly",
  "sperm flagellum assembly",
  "calcium−ion regulated exocytosis",
  "epithelial cilium movement involved in extracellular fluid movement",
  "synaptic vesicle priming",
  "sperm axoneme assembly"
)

top_go_terms <- gsea_go_tbl %>%
  dplyr::filter(Description %in% target_terms)

top_go_terms <- top_go_terms %>%
  dplyr::select(-last_col())

top_go_terms

# save the table
write.csv(top_go_terms,
          file = "results/GSEA_GO/top_go_terms.csv",
          row.names = FALSE)


# ---- run GSEA: KEGG----
gsea_kegg <- gseKEGG(
  geneList     = gene_rank,   # ranked vector (names = Entrez IDs)
  organism     = "hsa",       # human KEGG
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.001,
  verbose      = FALSE
)

# ---- save result table ----
gsea_kegg_tbl <- as.data.frame(gsea_kegg)
nrow(gsea_kegg_tbl) # 41
colnames(gsea_kegg_tbl)
#  [1] "ID"              "Description"     "setSize"         "enrichmentScore" "NES"             "pvalue"         
# [7] "p.adjust"        "qvalue"          "rank"            "leading_edge"    "core_enrichment"

write.csv(gsea_kegg_tbl, file.path(out_dir_kegg, "GSEA_kegg_table.csv"))
message("Saved: ", file.path(out_dir_kegg, "GSEA_kegg_table.csv"))



# --- Save session info----
writeLines(capture.output(sessionInfo()), file.path("results/session/", "session_info_GSEA.txt"))
message("Done. Outputs in: ", normalizePath("results/session"))



