###load libraries###

# load libraries
library(DESeq2) # for DGE analysis using a negative binomial distribution
library(tibble) # reimagining of a dataframe, also provided by the tidyverse package
library(dplyr) # a core package in the tidyverse, 
library(rtracklayer) # for importing, exporting, and manipulating genomic data, esp for working with genome annotation files like BED, GFF
library(ggplot2) # for data visualization, part of the tidyverse. You build plots by layering components like data, aesthetics, and geometirc objects
library(ggrepel) # an extension to ggplot2 that helps you add text labels to plots without overlap
library(EnhancedVolcano) # for making publication-ready volcano plots
library(pheatmap) # for making pretty heatmaps with clustering and annotations
library(clusterProfiler) # for functional enrichment analysis and visualization of gene clusters or gene list, widely used after DE to interpret biological meaning
library(enrichplot) # a companion package to clusterProfiler, focusing on visualization of enrichment results (GO, KEGG, GSEA) from packages like clusterProfiler
library(org.Hs.eg.db) # a bioconductor annotation package for human, Dpe genome was annotated against the human genome
library(GO.db) # another bioconductor annotation package, providing the GO database itself

###data cleaning###

# read in the count matrix and clean 
counts_raw <- read.table("results/counts/gene_id_trimmed.counts", header = TRUE, row.names = 1, comment.char = "#", check.names = FALSE)
# file path
# first row contains column names
# use the first column as row names
# ignore lines starting with "#"
# don't modify column names to be syntactically valid


class(counts_raw)

head(counts_raw)

# inspect column names
colnames(counts_raw)

# extract count data
counts <- counts_raw[, 6:ncol(counts_raw)] # start from the 6th column to the last column
# ncol() returns the total number of columns


# view the count matrix
head(counts)

class(counts)


# convert gene ids into gene symbols
rownames(counts) <- sub("_Dpe.*", "", rownames(counts))

rownames(counts)

# create sample conditions
coldata <- data.frame(
  row.names = colnames(counts), # use sample/column names as row names of coldata
  condition = c("Statocyst", "Statocyst", "Gill", "Gill")
)


# check if row names of coldata same as column names of count matrix
all(rownames(coldata) == colnames(counts))

###DEG analysis###

# create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition) # the formula that tells DESeq2 what to compare



# fit the model and estimate dispersion, automatic normalization, DESeq2 uses size factor normalization
dds <- DESeq(dds)
# performs the full DESeq2 workflow on your dataset
# estimate size factor - normalize for sequencing depth
# estimate dispersion - variance of counts across replicates
# fit the negative binomial GLM - according to your design
# perform wald test or LRT - to compute DE
# returns the updated DESeqDataSet ready for results extraction


# get normalized counts table - raw counts divided by size factor to correct raw library size difference, using linear scale
# norm_counts <- counts(dds, normalized = TRUE)

# extract differential expression results
res <- results(dds)
res



# visualize the results
plotMA(res,
       cex = 0.6, # point size
       main = "MA plot",
       xlab = "Mean of Normalized Counts",
       ylab = "Log2 Fold Change",
       alpha = 0.05)


# check the reference level!
coldata
levels(coldata$condition) # NULL

# convert condition to a factor with the correct reference - Statocyst vs Gill, reference level = Gill
coldata$condition <- factor(coldata$condition, levels = c("Gill", "Statocyst"))

# recheck the reference level
levels(coldata$condition) # [1] "Gill"      "Statocyst"

# update the dds object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

# rerun 
dds <- DESeq(dds)
# fits the negative binomial model and estimates size factors, dispersion, and coefficients 
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

# extract results
res <- results(dds) 
# independent filtering 
# deseq2 checks the mean normalized counts for each gene
# genes below a certain threshold (mean count < 1 by default) are filtered out before testing, so they don't get p-values

res

head(rownames(res))

res["PIEZO2", c("baseMean", "log2FoldChange", "pvalue", "padj")]

res["TRPC4", c("baseMean", "log2FoldChange", "pvalue", "padj")]


# check results summary
summary(res)

rownames(res)

# confirm the contrast
resultsNames(dds)

# replot 
plotMA(res,
       cex = 0.6, # point size
       main = "MA plot",
       xlab = "Mean of Normalized Counts",
       ylab = "Log2 Fold Change",
       alpha = 0.05)

##########################statistical testing for all genes (with top 100 genes labeled)#####################################

# Enhanced volcano plot to show the top 100 genes by significance (padj < 0.05) & log2FC > 1

# remove NAs first to avoid errors
res_noNA <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]

# order by padj first, then by absolute log2FC (descending)
res_ordered <- res_noNA[order(res_noNA$padj, -abs(res_noNA$log2FoldChange)), ]

# check
head(res_ordered, 20)

# get top 100 genes 
top100_genes <- rownames(res_ordered)[1:100]
top100_genes


EnhancedVolcano(res,
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
    plot.title = element_text(hjust = 0.5, face = "bold")   # <-- centers title
  )

##########################################################################################################################

####################################filter DEGs###########################################################################

# filter significant genes with padj < 0.05 & log2FC > 1
sig_genes <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
sig_genes # 8068 genes 
class(sig_genes)

# export sig_genes
write.csv(as.data.frame(sig_genes), "DEG_trimmed.csv")

##########################################################################################################################

##########################################volcano plot for top 100 DEGs#################################################

# check NAs
any(is.na(sig_genes$log2FoldChange))  # TRUE if there are NAs
any(is.na(sig_genes$padj))

# convert to a dataframe
sig_genes_df <- as.data.frame(sig_genes)
head(sig_genes_df)

# PKD2 expression 
sig_genes_df["PKD2", ]


# highlight upregulated vs downregulated genes in different colors
sig_genes_df_up_down <- sig_genes_df %>%
  mutate(regulation = case_when(
    padj < 0.05 & log2FoldChange > 1 ~ "Up",
    padj < 0.05 & log2FoldChange < -1 ~ "Down",
    TRUE ~ "Not Significant")
    )

head(sig_genes_df_up_down)

# filter the most significant genes
# top up-regulated (log2FC > 1) with padj < 0.05
top_up <- sig_genes_df %>% filter(log2FoldChange > 1, padj < 0.05) %>% arrange(desc(log2FoldChange))

# top down-regulated
top_down <- sig_genes_df %>% filter(log2FoldChange < 1, padj < 0.05) %>% arrange(log2FoldChange)

# pick top 100 genes

top100_sig_genes <- sig_genes_df_up_down %>% 
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice_head(n = 100)

head(top100_sig_genes)

ggplot(sig_genes_df_up_down, aes(x = log2FoldChange, y = -log10(padj))) +
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


#############################################################

sig_genes_df_up_down$padj[sig_genes_df_up_down$padj == 0] <- 1e-300
sig_genes_df_up_down$neg_log10_padj <- -log10(sig_genes_df_up_down$padj)

ggplot(sig_genes_df_up_down, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = regulation), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Up" = "firebrick", "Down" = "steelblue", "Not Significant" = "gray80")) +
  geom_text_repel(
    data = top100_sig_genes,
    aes(label = rownames(top100_sig_genes)),
    size = 2,
    max.overlaps = 30,
    box.padding = 0.4,
    segment.size = 0.2
  ) +
  coord_cartesian(clip = "off", ylim = c(0, max(sig_genes_df_up_down$neg_log10_padj) * 1.1)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) +
  labs(
    title = "Volcano Plot of DEGs",
    x = "Log2 Fold Change",
    y = "-log10 adjusted p-value",
    color = "Regulation"
  )









################################################VST Normalization#########################################################

# variance stabilizing transformation, a normalization method in DESeq2, log-like scale, suitable for PCA, heatmaps, clustering
vsd <- vst(dds, blind = FALSE) # blind = FALSE, tells the function to use the experimental design(your ~ condition) when stabilizing variance
head(assay(vsd))


#####################################################################################################

###################################################PCA############################################

# PCA to confirm sample group separation
q <- plotPCA(vsd, intgroup = "condition") 

q +
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
#PCA captures global patterns in the data, need to use all genes 

####################################Test below later############################
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
# Plot PCA
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCA of Samples",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = "Condition"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +
  scale_color_brewer(palette = "Set2")

##################################################################################





#####################check back later###############################
sample_info <- data.frame(condition = colData(dds)$condition)
rownames(sample_info) <- colnames(dds)

pheatmap(assay(vsd)[rownames(sig_genes), ], 
         cluster_rows = TRUE, 
         annotation_col = sample_info,
         show_rownames = FALSE)
#####################################################################


##########################################heatmap#############################


# extract the matrix of vst-normalized counts
vsd_mat <- assay(vsd) # rows = geneID, columns = samples
class(vsd_mat) # [1] "matrix" "array" 

head(vsd_mat)

# only keep rows with matching geneIDs
keep_ids <- intersect(rownames(sig_genes_df), rownames(vsd_mat))

# subset vst expression matrix
vsd_mat_sig <- vsd_mat[keep_ids, ]
class(vsd_mat_sig)

"PKD2" %in% rownames(vsd_mat_sig)

# vsd_mat_symbol[grep("^tmc", rownames(vsd_mat_symbol), ignore.case = TRUE), ]

# vsd_mat_symbol[grep("^piezo", rownames(vsd_mat_symbol), ignore.case = TRUE), ]

# vsd_mat_symbol[grep("^trp", rownames(vsd_mat_symbol), ignore.case = TRUE), ]

# replace rownames with gene symbol
# match gene symbols to rownames of the matrix
# symbol_map <- sig_genes_geneid %>% 
#  filter(geneID %in% keep_ids) %>%
#  distinct(geneID, symbol)

# reorder to match vsd matrix row order
# gene_symbols <- symbol_map$symbol[match(rownames(vsd_mat_symbol), symbol_map$geneID)]
# gene_symbols

# replace rownames
# rownames(vsd_mat_symbol) <- gene_symbols

# change column names
colnames(vsd_mat_sig) <- c("Statocyst 1", "Statocyst 2", "Gill 1", "Gill 2")

# save 
write.csv(vsd_mat_sig, "vst_mat_sig_genes.csv")

head(vsd_mat_sig)

# sample annotation
annotation_col <- as.data.frame(colData(dds)[, "condition", drop = FALSE])
rownames(annotation_col) <- colnames(vsd_mat_sig)

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
  main = "Heatmap of DEGs"
)

########################################heatmap of top 200 genes########################
# plot heatmap for top 200 genes
top200_ids <- sig_genes_df %>%
  arrange(padj) %>%
  slice_head(n = 200) 

top200_vsd_mat_sig <- vsd_mat_sig[intersect(rownames(top200_ids), rownames(vsd_mat_sig)), ]
# rownames(top200_vsd_mat_symbol) <- sig_genes_geneid$symbol[match(rownames(vsd_mat_symbol_top200), sig_genes_geneid$geneID)]
# colnames(vsd_mat_symbol_top200) <- c("Statocyst 1", "Statocyst 2", "Gill 1", "Gill 2")

# save
write.csv(top200_vsd_mat_sig, file = "vst_sig_genes_top200.csv")

pheatmap(
  top200_vsd_mat_sig, 
  scale = "row",
  annotation_col = annotation_col,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 2,
  main = "Heatmap of Top 200 DEGs"
)


##########################################################################


########################################heatmap of top 100 genes########################
# plot heatmap for top 100 genes
top100_ids <- sig_genes_df %>%
  arrange(padj) %>%
  slice_head(n = 100) 

top100_vsd_mat_sig <- vsd_mat_sig[intersect(rownames(top100_ids), rownames(vsd_mat_sig)), ]
# rownames(top200_vsd_mat_symbol) <- sig_genes_geneid$symbol[match(rownames(vsd_mat_symbol_top200), sig_genes_geneid$geneID)]
# colnames(vsd_mat_symbol_top200) <- c("Statocyst 1", "Statocyst 2", "Gill 1", "Gill 2")

# save
write.csv(top100_vsd_mat_sig, file = "vst_sig_genes_top100.csv")

pheatmap(
  top100_vsd_mat_sig, 
  scale = "row",
  annotation_col = annotation_col,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 2,
  main = "Heatmap of Top 100 DEGs"
)


##########################################################################

###Gene ontology enrichment analysis, find what biological processes, functions, or cellular components are overrepresented
## Because Doryteuthis pealeii is not a model organism, standard GO/KEGG tools don't support it directly. The annotation provides: 
## JGI gene symbol and combined_symbols based on closest human matches
## this means you can map squid DEGs TO human gene symbols, then perform enrichment using human databases

# step1 - get the gene symbols for squid DEGs, using sig_genes, now match them to gene symbols
# load gene id mapping
# id_map <- read.delim("pita_all_ids.txt", header = TRUE)
# head(id_map)

# inspect sig_genes row names
#sig_genes@rownames
# head(sig_genes) 
# the gene IDs in sig_genes don't match any columns in id_map. sig_genes has gene ids that contain gene symbols, 
# coordinates, and strand info

# extract the gene symbols from the beginning of each string and then match it to combined_symbols in id_map
# add column names to row names
# class(sig_genes)
# sig_genes_df <- as.data.frame(sig_genes)
# sig_genes_geneid <- rownames_to_column(sig_genes_df, var = "geneID")

# head(sig_genes_geneid)

# extract gene symbols
# sig_genes_geneid$symbol <- sub("_Dpe.*", "", sig_genes_geneid$geneID)
# head(sig_genes_geneid)

# merge significant genes with annotation by joining on gene symbols
# merged_sig <- left_join(sig_genes_geneid, id_map, by = c("symbol" = "combined_symbols"))
# head(merged_sig)

# load annotation file
gtf_file <- "annotation/pita_genes_v0.3.bed_12to.gtf"
gtf_data <- import(gtf_file)
gtf_df <- as.data.frame(gtf_data)

# gtf_data

### filter the most significant genes
# top up-regulated (log2FC > 1) with padj < 0.05
# top_up <- merged_sig %>% filter(log2FoldChange > 1, padj < 0.05) %>% arrange(desc(log2FoldChange))

# top down-regulated
# top_down <- merged_sig %>% filter(log2FoldChange < 1, padj < 0.05) %>% arrange(log2FoldChange)




# Save mergerd_sig 
# write.csv(merged_sig, "sig_genes_pita_map.csv")
# write.csv(top_up, "sig_genes_pita_up.csv")
# write.csv(top_down, "sig_genes_pita_down.csv")

# volcano plot to show 
# ggplot(merged_sig, aes(x = log2FoldChange, y = -log10(padj))) +
# geom_point(aes(color = padj < 0.05)) +
# theme_minimal() +
# labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10(FDR)")

# make volcano plot for top 100 genes
# highlight upregulated vs downregulated genes in different colors
# merged_sig_up_down <- merged_sig %>%
# mutate(regulation = case_when(
#   padj < 0.05 & log2FoldChange > 1 ~ "Up",
#   padj < 0.05 & log2FoldChange < -1 ~ "Down",
#   TRUE ~ "Not Significant"
#  ))


# pick top 100 genes
# top100 <- merged_sig %>% 
#  filter(!is.na(padj)) %>%
#  arrange(padj) %>%
#  slice_head(n = 100)

# ggplot(merged_sig_up_down, aes(x = log2FoldChange, y = -log10(padj))) +
#  geom_point(aes(color = regulation), alpha = 0.7) +
#  scale_color_manual(values = c("Up" = "firebrick", "Down" = "steelblue", "Not Significant" = "gray80")) +
# geom_text_repel(
#   data = top100,
#   aes(label = symbol),
#   size = 1.5,
#    max.overlaps = Inf
#  ) +
#  theme_minimal() +
#  theme(
#    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
#  ) +
#  labs(
#    title = "Volcano Plot of Top 100 Genes",
#    x = "Log2 Fold Change",
#    y = "-log10(FDR)",
#    color = "Regulation"
#  )

# plot heatmap

# extract the matrix of vst-normalized counts
vsd_mat <- assay(vsd) # rows = geneID, columns = samples

# only keep rows with matching geneIDs
keep_ids <- intersect(sig_genes_geneid$geneID, rownames(vsd_mat))

# subset vst expression matrix
vsd_mat_symbol <- vsd_mat[keep_ids, ]

"PKD2" %in% rownames(vsd_mat_symbol)

vsd_mat_symbol[grep("^tmc", rownames(vsd_mat_symbol), ignore.case = TRUE), ]

vsd_mat_symbol[grep("^piezo", rownames(vsd_mat_symbol), ignore.case = TRUE), ]

vsd_mat_symbol[grep("^trp", rownames(vsd_mat_symbol), ignore.case = TRUE), ]

# replace rownames with gene symbol
# match gene symbols to rownames of the matrix
symbol_map <- sig_genes_geneid %>% 
  filter(geneID %in% keep_ids) %>%
  distinct(geneID, symbol)

# reorder to match vsd matrix row order
gene_symbols <- symbol_map$symbol[match(rownames(vsd_mat_symbol), symbol_map$geneID)]
gene_symbols

# replace rownames
rownames(vsd_mat_symbol) <- gene_symbols

# change column names
colnames(vsd_mat_symbol) <- c("Statocyst 1", "Statocyst 2", "Gill 1", "Gill 2")

# save 
write.csv(vsd_mat_symbol, "sig_genes_symbol.csv")

# sample annotation
annotation_col <- as.data.frame(colData(dds)[, "condition", drop = FALSE])
rownames(annotation_col) <- colnames(vsd_mat_symbol)

# plot heatmap
pheatmap(
  vsd_mat_symbol,
  scale = "row",
  annotation_col = annotation_col,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 0.5,
  main = "Heatmap of DEGs"
)

# plot heatmap for top 200 genes
top200_ids <- sig_genes_geneid %>%
  arrange(padj) %>%
  slice_head(n = 200) %>%
  pull(geneID)

vsd_mat_symbol_top200 <- vsd_mat[intersect(top200_ids, rownames(vsd_mat)), ]
rownames(vsd_mat_symbol_top200) <- sig_genes_geneid$symbol[match(rownames(vsd_mat_symbol_top200), sig_genes_geneid$geneID)]
colnames(vsd_mat_symbol_top200) <- c("Statocyst 1", "Statocyst 2", "Gill 1", "Gill 2")

# save
write.csv(vsd_mat_symbol_top200, file = "vsd_Top200.csv")

pheatmap(
  vsd_mat_symbol_top200, 
  scale = "row",
  annotation_col = annotation_col,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 2,
  main = "Heatmap for Top 200 DEGs"
)

### GO enrichment analysis###

# step 1 - get unique gene symbols of significant genes
uniq_sig_genes <- merged_sig %>%
  filter(padj < 0.05, !is.na(symbol)) %>%
  pull(symbol) %>%
  unique()

uniq_sig_genes


# step 2 - convert gene symbols to Entrez IDs
gene_df <- bitr(uniq_sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # 64.43% of input gene IDs failed to map

entrez_ids <- gene_df$ENTREZID

# step 3 - run GO Enrichment
ego <- enrichGO(
  gene = entrez_ids, # a vector of Entrez Gene IDs not gene symbols
  OrgDb = org.Hs.eg.db, # organism-specific annotation package
  ont = "BP", #GO ontology BP = Biological Process
  pAdjustMethod = "BH", # multiple testing correction method, BH = Benjamini-Hochberg false discovery rate
  qvalueCutoff = 0.05, # threshold for padj(FDR), only terms with q-value < 0.05 will be considered significant
  readable = TRUE # converts Entrez IDs in results back to gene symbols for easier interpretation
)


# GO dotplot
p <- dotplot(ego, showCategory = 30, title = "GO Biological Processes")

# customize font sizes using theme()
p + 
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

# interpretation of GO enrichment plot
# x-axis: GeneRatio = number of your genes found in a GO term / total number of significant gene submitted
# out of all the significant genes I tested, what proportion are annotated to this GO term?

# view genes in a specific GO term
ego_df <- as.data.frame(ego)

# see the genes in the top enriched GO term

ego_df$Description[1] # name of the GO term
ego_df$geneID[1] # the genes involved

ego_df[1, c("Description", "geneID")]

# find genes for specific terms
ego_df$Description

# inner ear development pathway
inner_ear_development <- ego_df[ego_df$Description == "inner ear development", "geneID"]
inner_ear_development <- unlist(strsplit(inner_ear_development, "/"))

inner_ear_development

# get the genes present in vsd matrix
genes_to_vsd_inner <- intersect(inner_ear_development, rownames(vsd_mat_symbol))

# subset the matrix
expr_subset_inner <- vsd_mat_symbol[genes_to_vsd_inner, ]

"PKD2" %in% rownames(expr_subset_inner) # FALSE

rownames(expr_subset_inner)

# plot heatmap for genes expressed in the inner ear development GO term
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
  main = "Genes expressd in Inner Ear Development"
)



# sensory perception of mechanical stimulus
mechanical_stimulus <- ego_df[ego_df$Description == "sensory perception of mechanical stimulus", "geneID"]
mechanical_stimulus <- unlist(strsplit(mechanical_stimulus, "/"))

genes_to_vsd_mechanical <- intersect(mechanical_stimulus, rownames(vsd_mat_symbol))

expr_subset_mechanical <- vsd_mat_symbol[genes_to_vsd_mechanical, ]

rownames(expr_subset_mechanical)


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
  main = "Genes expressd in Sensory Perception of Mechanical Stimulus"
)



# sensory perception of sound
sound <- ego_df[ego_df$Description == "sensory perception of sound", "geneID"]
sound <- unlist(strsplit(sound, "/"))

genes_to_vsd_sound <- intersect(sound, rownames(vsd_mat_symbol))

expr_subset_sound <- vsd_mat_symbol[genes_to_vsd_sound, ]

rownames(expr_subset_sound)

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
  main = "Genes expressd in Sensory Perception of Sound"
)



# calcium ion transmembrane transport
# retrieve gene symbols from that exact description
calcium_transport <- ego_df[ego_df$Description == "calcium ion transmembrane transport", "geneID"]
# this is a character vector of length 1 

# split the string into a list by / and convert it into a simple character vector
calcium_transport <- unlist(strsplit(calcium_transport, "/"))

# find common genes between the genes in calcium_transport and the gene symbols in expression matrix
genes_to_vsd_calcium <- intersect(calcium_transport, rownames(vsd_mat_symbol))

# extract the expression values for only those genes from the expression matrix
expr_subset_calcium <- vsd_mat_symbol[genes_to_vsd_calcium, ]

rownames(expr_subset_calcium)

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
  main = "Genes expressd in Calcium Ion Transmembrane Transport"
)

# find which pathways a gene is involved in 

# convert gene symbol to Entrez ID
gene_symbol <- "PKD2"
entrez_id_pkd2 <- bitr(gene_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

# find GO terms annotated to this gene

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


# get GO term descriptions, remove duplicated GOIDs

go_terms <- AnnotationDbi::select(
  GO.db,
  keys = unique(go_annotations_pkd2$GO),
  columns = "TERM",
  keytype = "GOID"
) %>%
  distinct(GOID, TERM)

go_terms



# Merge GO annotations and term names
go_merged_unique <- go_annotations_pkd2 %>%
  distinct(ENTREZID, GO, ONTOLOGY) %>%
  left_join(go_terms, by = c("GO" = "GOID"))

# Optional: count frequency if some terms repeat
go_counts_unique <- go_merged_unique %>%
  count(TERM, ONTOLOGY) %>%
  arrange(desc(n))

# Barplot
ggplot(go_counts_unique, aes(x = reorder(TERM, n), y = n, fill = ONTOLOGY)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = paste("Unique GO Terms for", gene_symbol),
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


# check if ANK1 appear in vsd_mat_symbol object
"ANK1" %in% rownames(vsd_mat_symbol)
"ANK1_2" %in% rownames(vsd_mat_symbol)

# to find its exact row
which(rownames(vsd_mat_symbol) == "ANK1")
which(rownames(vsd_mat_symbol) == "ANK1_2")

# extract the row data
vsd_mat_symbol["ANK1", ]
vsd_mat_symbol["ANK1_2", ]

vsd_mat_symbol["PKD2", ]

# check if ANK1 appear in gene_df
"ANK1" %in% gene_df$SYMBOL # returns TRUE, so ANK1 appears in GO term
"ANK1_2" %in% gene_df$SYMBOL # returns FALSE

# to find its exact row
which(gene_df$SYMBOL == "ANK1")

# extract the row data
gene_df[gene_df$SYMBOL == "ANK1", ]




# ANK1
# convert gene symbol to Entrez ID
gene_symbol_ank1 <- "ANK1"

entrez_id_ank1 <- bitr(gene_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

# find GO terms annotated to this gene

columns(org.Hs.eg.db)

go_annotations_ank1 <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = entrez_id_ank1,
    columns = c("GO", "ONTOLOGY", "SYMBOL"),
    keytype = "ENTREZID"
  )
) 

go_annotations_ank1


# get GO term descriptions, remove duplicated GOIDs

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
ggplot(go_counts_unique_ank1, aes(x = reorder(TERM, n), y = n, fill = ONTOLOGY)) +
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



# compare ANK1 expression level between conditions
# extract ANK1 expression
gene <- "ANK1"
expr_value <- vsd_mat_symbol[gene, ]

# combine with condition labels
df <- data.frame(
  expression = as.numeric(expr_value),
  condition = coldata$condition
)

# plot expression by condition
ggplot(df, aes(x = condition, y = expression, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = paste("Expression of", gene), y = "Normalized Expression") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )





# compare PKD2 expression between conditions
# extract PKD2 expression
gene2 <- "PKD2"
expr_value2 <- vsd_mat_sig[gene2, ]

expr_value2

# combine with condition labels
df2 <- data.frame(
  expression = as.numeric(expr_value2),
  condition = coldata$condition)

# plot expression by condition
ggplot(df2, aes(x = condition, y = expression, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  labs(title = paste("Expression of", gene2), y = "vst normalized expression") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


vsd_mat_symbol["PKD2", ]

vsd_mat_symbol["PKD1L2_3", ]

vsd_mat_symbol["PKD1L2_4", ]

vsd_mat_symbol["TRPM2", ]

vsd_mat_symbol["TRPV5_2", ]

vsd_mat_symbol["TRPV5", ]

vsd_mat_symbol["TRPC4", ]

vsd_mat_symbol["TMC5", ]

vsd_mat_symbol["PIEZO2", ]


