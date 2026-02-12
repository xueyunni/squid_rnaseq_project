# ==============================================================================
# Script: 03_GSEA_GO.R
# Project: Statocyst vs Gill RNA-seq (DESeq2)
# Purpose: GSEA (GO Biological Process) using ranked DESeq2 statistics
# Author: Max Ni
# Date: 2026-02
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
dds_rds   <- "results/rds/dds_statocyst_vs_gill.rds"   
out_dir   <- "results/GSEA_GO"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

# ---- build ranked list for GSEA ----
# Recommended ranking for GSEA: Wald statistic (stat), not log2FC
gene_rank <- res_df$stat
names(gene_rank) <- res_df$gene

gene_rank <- sort(gene_rank, decreasing = TRUE)

head(gene_rank)

# Optional: remove duplicated gene IDs (keep the one with largest |stat|)
# Useful if ortholog mapping created duplicates.
if (any(duplicated(names(gene_rank)))) {
  message("Detected duplicated gene IDs in ranking list. Collapsing duplicates by max |stat|...")
  gene_rank <- tapply(gene_rank, names(gene_rank), function(x) x[which.max(abs(x))])
  gene_rank <- sort(gene_rank, decreasing = TRUE)
}

summary(gene_rank)
sum(duplicated(gene_rank))

# ---- run GSEA: GO Biological Process ----
# keyType: set to the ID type in names(gene_rank)
# Common: "SYMBOL" or "ENTREZID"
key_type_used <- "SYMBOL"   # <- CHANGE if needed (e.g., "ENTREZID")

ego_gsea <- gseGO(
  geneList     = gene_rank,
  OrgDb        = org.Hs.eg.db,
  keyType      = key_type_used,
  ont          = "BP",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

# ---- save result table ----
gsea_tbl <- as.data.frame(ego_gsea@result)

write.csv(gsea_tbl, file.path(out_dir, "GSEA_GO_BP_table.csv"))
message("Saved: ", file.path(out_dir, "GSEA_GO_BP_table.csv"))

# ---- quick filtering helpers ----
# search for mechanosensory / ear / cilia terms
keyword_hits <- gsea_tbl %>%
  filter(grepl("ear|mechan|cilium|axoneme|sensory|synap", Description, ignore.case = TRUE)) %>%
  arrange(p.adjust)

write_csv(keyword_hits, file.path(out_dir, "GSEA_GO_BP_keywords.csv"))

# ---- dotplot (top terms) ----
p_dot <- dotplot(ego_gsea, showCategory = 30) +
  scale_y_discrete(labels = function(x) gsub("\n", " ", x)) +   # force one line per GO term
  ggtitle("GSEA: GO Biological Process (Statocyst vs Gill)") +
  theme(
    axis.text.y = element_text(size = 9),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 35)
  )

ggsave(
  filename = file.path(out_dir, "GSEA_GO_BP_dotplot.pdf"),
  plot     = p_dot,
  width    = 10,
  height   = 7
)

message("Saved: ", file.path(out_dir, "GSEA_GO_BP_dotplot.pdf"))

# ---- enrichment curve for a chosen term ----
# Pick the top term or specify an ID from gsea_tbl$ID
if (nrow(gsea_tbl) > 0) {
  top_id <- gsea_tbl$ID[gsea_tbl$Description == "cilium movement"]
  desc <- gsea_tbl$Description[gsea_tbl$Description == "cilium movement"]
  p_curve <- gseaplot2(ego_gsea, geneSetID = top_id, title = desc)
  
  ggsave(
    filename = file.path(out_dir, paste0("GSEA_curve_", top_id, ".pdf")),
    plot     = p_curve,
    width    = 9,
    height   = 5
  )
  message("Saved enrichment curve for: ", top_id)
}

message("Done.")