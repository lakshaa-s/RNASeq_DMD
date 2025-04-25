########################################
########################################

setwd ("/Users/lakshaasrishankar/Downloads")

library(GenomicRanges)
install.packages("BiocManager") 

# genomation package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("genomation", force = TRUE)
library(genomation)

# rtracklayer package
if(!require("rtracklayer", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("rtracklayer")
}
library(rtracklayer)

# GenomicRanges package
if(!require("GenomicRanges", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("GenomicRanges")
}
library(GenomicRanges)

# DESeq2 package
if(!require("DESeq2", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("DESeq2")
}
library(DESeq2)

if(!require("EnhancedVolcano", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("EnhancedVolcano")
}
library(EnhancedVolcano)

########################################
# DESeq2 analysis
########################################


# Adjusted sample files and conditions
sampleFiles <- c("SRR16969752.txt", "SRR16969757.txt", "SRR16969762.txt", "SRR16969767.txt", "SRR12798466.txt", "SRR12798472.txt", "SRR12798478.txt",
                 "SRR16969732.txt", "SRR16969737.txt", "SRR16969742.txt", "SRR16969747.txt", "SRR12798469.txt", "SRR12798475.txt", "SRR12798481.txt")



# Assigning conditions based on colors in the table
sampleCondition <- c("ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl",
                     "mut", "mut", "mut", "mut", "mut", "mut", "mut")


# Directory where your htseq files are stored
directory <- "~/Downloads"


sampleTable <- data.frame(
  sampleName = gsub(".txt", "", sampleFiles),
  fileName = sampleFiles,
  condition = sampleCondition
)


dds <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable,
  directory = directory,
  design = ~ condition
)

# View dataset to confirm
dds
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Perform DESeq2 analysis
dds <- DESeq(dds, sfType="ratio")

# Contrasts for each comparison
contrast_mut_vs_ctrl <- c("condition", "mut", "ctrl")
res_mut_vs_ctrl <- results(dds, contrast = contrast_mut_vs_ctrl, alpha = 0.05)

# Shrink results for better visualization
res_mut_vs_ctrl_ashr <- lfcShrink(dds, contrast = contrast_mut_vs_ctrl, res = res_mut_vs_ctrl, type = "ashr")
res_mut_vs_ctrl_apeglm <- lfcShrink(dds, coef = "condition_mut_vs_ctrl", res = res_mut_vs_ctrl, type = "apeglm")
res_mut_vs_ctrl_normal <- lfcShrink(dds, coef = "condition_mut_vs_ctrl", res = res_mut_vs_ctrl, type = "normal")


# Save results
save(dds, res_mut_vs_ctrl, res_mut_vs_ctrl_ashr, res_mut_vs_ctrl_apeglm, res_mut_vs_ctrl_normal, file="DESeq2_results.RData")

# Filter for genes that are both statistically significant and have a strong fold change
deg_filtered_ashr <- res_mut_vs_ctrl_ashr[!is.na(res_mut_vs_ctrl_ashr$padj)
                                             & abs(res_mut_vs_ctrl_ashr$log2FoldChange) >= 2, ]

deg_filtered_apeglm <- res_mut_vs_ctrl_apeglm[!is.na(res_mut_vs_ctrl_apeglm$padj)
                                               & abs(res_mut_vs_ctrl_apeglm$log2FoldChange) >= 2, ]

deg_filtered_normal <- res_mut_vs_ctrl_normal[!is.na(res_mut_vs_ctrl_normal$padj)
                                               & abs(res_mut_vs_ctrl_normal$log2FoldChange) >= 2, ]

# Save filtered results
write.table(deg_filtered_ashr, file = "mut_vs_ctrl_ashr_filtered.csv", sep = ",", row.names = FALSE)
write.table(deg_filtered_apeglm, file = "mut_vs_ctrl_apeglm_filtered.csv", sep = ",", row.names = FALSE)
write.table(deg_filtered_normal, file = "mut_vs_ctrl_normal_filtered.csv", sep = ",", row.names = FALSE)

################################################################################
# plots
################################################################################

pdf("RNAseqVolcanoPlot_mut_vs_ctrl_ashr.pdf", width = 10, height = 8, pointsize = 10)
par(cex = 1.1)
par(las = 1)
par(mar = c(4, 4.0, 3, 1) + 0.1)
EnhancedVolcano(res_mut_vs_ctrl_ashr,
                lab = rep("", nrow(res_mut_vs_ctrl_ashr)),
                x ='log2FoldChange',
                y = 'padj',
                title = "Volcano Plot: Ashr",
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'adjusted p-value'),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                colAlpha = 1,
                xlim = c(-10,10),
                ylim = c(0,25),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = FALSE,
                widthConnectors = 0.75)
dev.off()

pdf("plotCounts_ashr.pdf", width = 8, height = 6)
plotCounts(dds, gene = which.min(res_mut_vs_ctrl_normal$padj), intgroup = "condition")
title(main = "Ashr", adj=0.73)
dev.off()

pdf("plotCounts_normal_ashr.pdf", width = 8, height = 6)
plotCounts(dds, gene="ENSG00000198947", intgroup="condition")
title(main = "Ashr", adj=0.73)
dev.off()

pdf("RNAseqVolcanoPlot_mut_vs_ctrl_apelgm.pdf", width = 10, height = 8, pointsize = 10)
par(cex = 1.1)
par(las = 1)
par(mar = c(4, 4.0, 3, 1) + 0.1)
EnhancedVolcano(res_mut_vs_ctrl_apeglm,
                lab = rep("", nrow(res_mut_vs_ctrl_apeglm)),
                x ='log2FoldChange',
                y = 'padj',
                title = "Volcano Plot: Apeglm",
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'adjusted p-value'),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                colAlpha = 1,
                xlim = c(-10,10),
                ylim = c(0,25),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = FALSE,
                widthConnectors = 0.75)
dev.off()

pdf("plotCounts_apelgm.pdf", width = 8, height = 6)
plotCounts(dds, gene = which.min(res_mut_vs_ctrl_normal$padj), intgroup = "condition")
title(main = "Apelgm", adj=0.73)
dev.off()

pdf("plotCounts_apelgm_DMD.pdf", width = 8, height = 6)
plotCounts(dds, gene="ENSG00000198947", intgroup="condition")
title(main = "Apelgm", adj=0.73)
dev.off()

pdf("RNAseqVolcanoPlot_mut_vs_ctrl_normal.pdf", width = 10, height = 8, pointsize = 10)
par(cex = 1.1)
par(las = 1)
par(mar = c(4, 4.0, 3, 1) + 0.1)
EnhancedVolcano(res_mut_vs_ctrl_normal,
                lab = rep("", nrow(res_mut_vs_ctrl_normal)),
                x ='log2FoldChange',
                y = 'padj',
                title = "Volcano Plot: Normal",
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'adjusted p-value'),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                colAlpha = 1,
                xlim = c(-10,10),
                ylim = c(0,25),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = FALSE,
                widthConnectors = 0.75)
dev.off()

pdf("plotCounts_normal.pdf", width = 8, height = 6)
plotCounts(dds, gene = which.min(res_mut_vs_ctrl_normal$padj), intgroup = "condition")
title(main = "Normal", adj=0.73)
dev.off()

pdf("plotCounts_normal_DMD.pdf", width = 8, height = 6)
plotCounts(dds, gene="ENSG00000198947", intgroup="condition")
title(main = "Normal", adj=0.73)
dev.off()
################################################################################
# DEG for ashr
################################################################################
library(rtracklayer)


genes <- import("~/Downloads/genome.gtf")
genes <- genes[genes$type == "gene"]



################################################################################
# ashr
################################################################################
positions <- match(rownames(res_mut_vs_ctrl_ashr), genes$gene_id)

# Removing NAs
valid <- !is.na(positions)

genes_ashr <- genes

genes_ashr$baseMeanCtrl <- 0
genes_ashr$baseMeanCtrl[positions[valid]] <- res_mut_vs_ctrl_ashr$baseMean[valid]

genes_ashr$log2FoldChange_mutCtrl <- 0
genes_ashr$log2FoldChange_mutCtrl[positions[valid]] <- res_mut_vs_ctrl_ashr$log2FoldChange[valid]

genes_ashr$padj_mutCtrl <- 1
genes_ashr$padj_mutCtrl[positions[valid]] <- res_mut_vs_ctrl_ashr$padj[valid]

# Filtering DEG
genes_ashr_DEG <- genes_ashr[!is.na(genes_ashr$padj_mutCtrl)]

genes_ashr_DEG <- genes_ashr_DEG[genes_ashr_DEG$padj_mutCtrl < 0.05
                       & abs(genes_ashr_DEG$log2FoldChange_mutCtrl) > 1]

genes_ashr_DEG_UP <- genes_ashr_DEG[genes_ashr_DEG$padj_mutCtrl < 0.05
                          & genes_ashr_DEG$log2FoldChange_mutCtrl > 1]
genes_ashr_DEG_DOWN <- genes_ashr_DEG[genes_ashr_DEG$padj_mutCtrl < 0.05
                            & genes_ashr_DEG$log2FoldChange_mutCtrl < -1]


write.table(genes_ashr_DEG, file = "mut_vs_ctrl_ashr_DEG.csv", sep = ",", row.names = FALSE)


write.table(genes_ashr_DEG_UP, file = "mut_vs_ctrl_ashr_UP.csv", sep = ",", row.names = FALSE)


write.table(genes_ashr_DEG_DOWN, file = "mut_vs_ctrl_ashr_DOWN.csv", sep = ",", row.names = FALSE)

write.table(genes_ashr, file = "mut_vs_ctrl_ashr_filtered.csv", sep = ",", row.names = FALSE)

################################################################################
# DEG for apeglm
################################################################################
library(rtracklayer)


genes <- import("~/Downloads/genome.gtf")
genes <- genes[genes$type == "gene"]



# For apelgm
positions_apelgm <- match(rownames(res_mut_vs_ctrl_apeglm), genes$gene_id)

# Removing NAs
valid_apelgm <- !is.na(positions_apelgm)

genes_apelgm <- genes
genes_apelgm$baseMeanCtrl <- 0
genes_apelgm$baseMeanCtrl[positions_apelgm[valid_apelgm]] <- res_mut_vs_ctrl_apeglm$baseMean[valid_apelgm]

genes_apelgm$log2FoldChange_mutCtrl <- 0
genes_apelgm$log2FoldChange_mutCtrl[positions_apelgm[valid_apelgm]] <- res_mut_vs_ctrl_apeglm$log2FoldChange[valid_apelgm]

genes_apelgm$padj_mutCtrl <- 1
genes_apelgm$padj_mutCtrl[positions_apelgm[valid_apelgm]] <- res_mut_vs_ctrl_apeglm$padj[valid_apelgm]

# Filtering DEG
genes_apelgm_DEG <- genes_apelgm[!is.na(genes_apelgm$padj_mutCtrl)]
genes_apelgm_DEG <- genes_apelgm_DEG[genes_apelgm_DEG$padj_mutCtrl < 0.05
                                     & abs(genes_apelgm_DEG$log2FoldChange_mutCtrl) > 1]

genes_apelgm_DEG_UP <- genes_apelgm_DEG[genes_apelgm_DEG$padj_mutCtrl < 0.05
                                        & genes_apelgm_DEG$log2FoldChange_mutCtrl > 1]
genes_apelgm_DEG_DOWN <- genes_apelgm_DEG[genes_apelgm_DEG$padj_mutCtrl < 0.05
                                          & genes_apelgm_DEG$log2FoldChange_mutCtrl < -1]

# Save results for apelgm
write.table(genes_apelgm_DEG, file = "mut_vs_ctrl_apeglm_DEG.csv", sep = ",", row.names = FALSE)
write.table(genes_apelgm_DEG_UP, file = "mut_vs_ctrl_apeglm_UP.csv", sep = ",", row.names = FALSE)
write.table(genes_apelgm_DEG_DOWN, file = "mut_vs_ctrl_apeglm_DOWN.csv", sep = ",", row.names = FALSE)
write.table(genes_apelgm, file = "mut_vs_ctrl_apeglm_filtered.csv", sep = ",", row.names = FALSE)


################################################################################
# DEG for normal
################################################################################
library(rtracklayer)


genes <- import("~/Downloads/genome.gtf")
genes <- genes[genes$type == "gene"]



# Positions matching gene IDs in the normal results
positions_normal <- match(rownames(res_mut_vs_ctrl_normal), genes$gene_id)

# Removing NAs
valid_normal <- !is.na(positions_normal)

genes_normal <- genes
genes_normal$baseMeanCtrl <- 0
genes_normal$baseMeanCtrl[positions_normal[valid_normal]] <- res_mut_vs_ctrl_normal$baseMean[valid_normal]

genes_normal$log2FoldChange_mutCtrl <- 0
genes_normal$log2FoldChange_mutCtrl[positions_normal[valid_normal]] <- res_mut_vs_ctrl_normal$log2FoldChange[valid_normal]

genes_normal$padj_mutCtrl <- 1
genes_normal$padj_mutCtrl[positions_normal[valid_normal]] <- res_mut_vs_ctrl_normal$padj[valid_normal]

# Filtering DEG for normal
genes_normal_DEG <- genes_normal[!is.na(genes_normal$padj_mutCtrl)]

genes_normal_DEG <- genes_normal_DEG[genes_normal_DEG$padj_mutCtrl < 0.05
                                     & abs(genes_normal_DEG$log2FoldChange_mutCtrl) > 1]

genes_normal_DEG_UP <- genes_normal_DEG[genes_normal_DEG$padj_mutCtrl < 0.05
                                        & genes_normal_DEG$log2FoldChange_mutCtrl > 1]
genes_normal_DEG_DOWN <- genes_normal_DEG[genes_normal_DEG$padj_mutCtrl < 0.05
                                          & genes_normal_DEG$log2FoldChange_mutCtrl < -1]

# Save results for normal
write.table(genes_normal_DEG, file = "mut_vs_ctrl_normal_DEG.csv", sep = ",", row.names = FALSE)
write.table(genes_normal_DEG_UP, file = "mut_vs_ctrl_normal_UP.csv", sep = ",", row.names = FALSE)
write.table(genes_normal_DEG_DOWN, file = "mut_vs_ctrl_normal_DOWN.csv", sep = ",", row.names = FALSE)
write.table(genes_normal, file = "mut_vs_ctrl_normal_filtered.csv", sep = ",", row.names = FALSE)
