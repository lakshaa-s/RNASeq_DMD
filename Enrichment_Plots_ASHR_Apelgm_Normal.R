
# Install and load necessary packages
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "enrichplot", "cowplot", "rtracklayer", "org.Hs.eg.db"))
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(cowplot)
library(rtracklayer)
library(org.Hs.eg.db)

# Set working directory
setwd("/Users/lakshaasrishankar/Downloads")

# Read the CSV file without setting row.names initially
count_data_ashr_raw <- read.csv("~/Downloads/mut_vs_ctrl_ashr_filtered.csv", header = TRUE)

# Make sure row names are unique
rownames(count_data_ashr_raw) <- make.unique(as.character(count_data_ashr_raw[, 1]))

# Remove the column that was originally used for row names (if necessary)
count_data_ashr <- count_data_ashr_raw[, -1]  # Adjust if the first column is not gene names

# Convert to a proper format (numeric matrix)
count_data_ashr <- as.matrix(count_data_ashr)

count_data_normal <- read.csv("~/Downloads/mut_vs_ctrl_normal_filtered.csv", row.names = 1)
count_data_apelgm <- read.csv("~/Downloads/mut_vs_ctrl_apeglm_filtered.csv", row.names = 1)

# Create metadata for conditions (adjust if necessary)
col_data <- data.frame(condition = factor(c("mut", "ctrl", "mut", "ctrl")))

# Ensure no duplicate row names
count_data_ashr <- make_unique_rows(count_data_ashr)
count_data_normal <- make_unique_rows(count_data_normal)
count_data_apelgm <- make_unique_rows(count_data_apelgm)

# Run DESeq2 for ASHR, Normal, and Apeglm datasets
dds_ashr <- DESeqDataSetFromMatrix(countData = count_data_ashr, colData = col_data, design = ~ condition)
dds_ashr <- DESeq(dds_ashr)
res_ashr <- results(dds_ashr)

dds_normal <- DESeqDataSetFromMatrix(countData = count_data_normal, colData = col_data, design = ~ condition)
dds_normal <- DESeq(dds_normal)
res_normal <- results(dds_normal)

dds_apelgm <- DESeqDataSetFromMatrix(countData = count_data_apelgm, colData = col_data, design = ~ condition)
dds_apelgm <- DESeq(dds_apelgm)
res_apelgm <- results(dds_apelgm)

# Convert gene IDs to Entrez IDs using clusterProfiler
convert_to_entrez <- function(gene_vector) {
  bitr(gene_vector,
       fromType = "ENSEMBL",
       toType = c("ENTREZID", "SYMBOL", "ENSEMBL"),
       OrgDb = org.Hs.eg.db)
}

# Extract upregulated and downregulated genes for ASHR, Normal, and Apeglm
upreg_genes_ashr <- convert_to_entrez(rownames(res_ashr[res_ashr$padj < 0.05 & res_ashr$log2FoldChange > 1,]))$ENTREZID
downreg_genes_ashr <- convert_to_entrez(rownames(res_ashr[res_ashr$padj < 0.05 & res_ashr$log2FoldChange < -1,]))$ENTREZID

upreg_genes_normal <- convert_to_entrez(rownames(res_normal[res_normal$padj < 0.05 & res_normal$log2FoldChange > 1,]))$ENTREZID
downreg_genes_normal <- convert_to_entrez(rownames(res_normal[res_normal$padj < 0.05 & res_normal$log2FoldChange < -1,]))$ENTREZID

upreg_genes_apelgm <- convert_to_entrez(rownames(res_apelgm[res_apelgm$padj < 0.05 & res_apelgm$log2FoldChange > 1,]))$ENTREZID
downreg_genes_apelgm <- convert_to_entrez(rownames(res_apelgm[res_apelgm$padj < 0.05 & res_apelgm$log2FoldChange < -1,]))$ENTREZID

# Enrichment analysis function
run_enrichment <- function(entrez_ids, organism = "hsa") {
  # KEGG enrichment
  enrich_kegg <- enrichKEGG(
    gene = entrez_ids,
    organism = organism,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  # GO Biological Process enrichment
  enrich_go_bp <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db, 
    ont = "BP",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  # GO Cellular Component enrichment
  enrich_go_cc <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    ont = "CC",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  # GO Molecular Function enrichment
  enrich_go_mf <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    ont = "MF",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  # Return a list of enrichment results
  return(list(
    KEGG = enrich_kegg,
    GO_BP = enrich_go_bp,
    GO_CC = enrich_go_cc,
    GO_MF = enrich_go_mf
  ))
}

# Run enrichment for each set of upregulated and downregulated genes
enrich_ashr_up <- run_enrichment(upreg_genes_ashr)
enrich_ashr_down <- run_enrichment(downreg_genes_ashr)

enrich_normal_up <- run_enrichment(upreg_genes_normal)
enrich_normal_down <- run_enrichment(downreg_genes_normal)

enrich_apelgm_up <- run_enrichment(upreg_genes_apelgm)
enrich_apelgm_down <- run_enrichment(downreg_genes_apelgm)

# Create dotplot function for enrichment results
create_dotplot <- function(enrich_result, condition, type, fill_color) {
  # Check if the enrichment result has any rows
  if (nrow(as.data.frame(enrich_result)) == 0) {
    warning("No enrichment data for ", condition, " ", type)
    return(NULL)
  }
  
  dotplot(enrich_result, showCategory = 20) +
    scale_fill_gradient(low = fill_color[1], high = fill_color[2], name = "p-value") +
    scale_size_continuous(name = "Enrichment") +
    ggtitle(paste0(condition, " ", type, " Enrichment")) +
    theme(plot.title = element_text(face = "bold", size = 14))
}

# Define color gradients for upregulated and downregulated genes
up_color_gradient <- c("#FF0000", "#7B68EE")   # Red to SlateBlue for upregulated
down_color_gradient <- c("#FFA500", "#00BFFF")   # Orange to DeepSkyBlue for downregulated

# Generate enrichment plots for each method and condition (upregulated and downregulated)
create_dotplot(enrich_ashr_up$KEGG, "ASHR", "KEGG", c("#FF0000", "#7B68EE"))
create_dotplot(enrich_ashr_down$KEGG, "ASHR", "KEGG", c("#FFA500", "#00BFFF"))

create_dotplot(enrich_normal_up$KEGG, "Normal", "KEGG", c("#FF0000", "#7B68EE"))
create_dotplot(enrich_normal_down$KEGG, "Normal", "KEGG", c("#FFA500", "#00BFFF"))

create_dotplot(enrich_apelgm_up$KEGG, "Apeglm", "KEGG", c("#FF0000", "#7B68EE"))
create_dotplot(enrich_apelgm_down$KEGG, "Apeglm", "KEGG", c("#FFA500", "#00BFFF"))
