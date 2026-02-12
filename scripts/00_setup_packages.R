# ==============================================================================
# 00_setup_packages.R
# Purpose: create/restore renv env and install packages
# ==============================================================================

# ---- renv bootstrap ----
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

# Initialize only if not already initialized
if (!file.exists("renv.lock")) {
  renv::init(bare = TRUE)  # bare avoids installing a bunch of suggested packages
} else {
  renv::restore(prompt = FALSE)
}

# ---- repos: keep Bioconductor + CRAN consistent ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
options(repos = BiocManager::repositories())

# ---- packages ----
bioc_pkgs <- c(
  "DESeq2",
  "rtracklayer",
  "EnhancedVolcano",
  "clusterProfiler",
  "enrichplot",
  "org.Hs.eg.db",
  "GO.db",
  "AnnotationDbi"
)

cran_pkgs <- c(
  "tibble",
  "dplyr",
  "tidyr",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "matrixStats",
  "forcats",
  "ggiraph"
)

# Install (renv will record versions; Bioc packages are ok here as long as repos are set)
renv::install(c(bioc_pkgs, cran_pkgs))

# Snapshot lockfile so others can restore exact versions
renv::snapshot(prompt = FALSE)

message("✅ Setup complete. Restart R, then run your analysis script.")
