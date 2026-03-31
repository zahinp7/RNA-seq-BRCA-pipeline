# Stage 6: Differential Expression Analysis with DESeq2
# Dataset: GSE183947 - BRCA tumor vs normal (6 samples)

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(EnhancedVolcano)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

#  1. Load count matrix 
counts_raw <- read.table("data/counts/counts.txt",
                         header = TRUE, skip = 1, row.names = 1)

counts <- counts_raw[, 6:ncol(counts_raw)]

colnames(counts) <- c("SRR15852399", "SRR15852407", "SRR15852416",
                      "SRR15852426", "SRR15852432", "SRR15852438")

#  2. Sample metadata 
coldata <- data.frame(
    sample    = colnames(counts),
    condition = c("tumor", "tumor", "tumor", "normal", "normal", "normal"),
    row.names = colnames(counts)
)
coldata$condition <- factor(coldata$condition, levels = c("normal", "tumor"))

#  3. DESeq2 
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = coldata,
                              design    = ~ condition)

dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "tumor", "normal"))
res <- lfcShrink(dds, coef = "condition_tumor_vs_normal", type = "apeglm")

dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

#  4. Add gene symbols 
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys      = rownames(res),
                       column    = "SYMBOL",
                       keytype   = "ENSEMBL",
                       multiVals = "first")

res_df <- as.data.frame(res)
res_df$gene   <- rownames(res_df)
res_df$symbol <- gene_symbols[rownames(res_df)]

# Save full annotated results
write.csv(res_df, "results/deseq2/deseq2_results_annotated.csv")

# Significant genes
sig <- res_df %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
    arrange(padj)
write.csv(sig, "results/deseq2/significant_genes.csv")
cat("Significant DEGs (padj<0.05, |LFC|>1):", nrow(sig), "\n")

#  5. PCA 
vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
    geom_point(size = 4) +
    geom_text_repel(size = 3) +
    scale_color_manual(values = c("normal" = "#2166AC", "tumor" = "#D6604D")) +
    xlab(paste0("PC1: ", pct_var[1], "% variance")) +
    ylab(paste0("PC2: ", pct_var[2], "% variance")) +
    ggtitle("PCA - BRCA Tumor vs Normal") +
    theme_classic()
ggsave("results/figures/pca.png", p_pca, width = 6, height = 5, dpi = 300)

#  6. Volcano plot 
png("results/figures/volcano.png", width = 2400, height = 2000, res = 300)
EnhancedVolcano(res_df,
    lab            = ifelse(is.na(res_df$symbol), res_df$gene, res_df$symbol),
    x              = "log2FoldChange",
    y              = "padj",
    pCutoff        = 0.05,
    FCcutoff       = 1,
    pointSize      = 2,
    labSize        = 3,
    title          = "BRCA Tumor vs Normal",
    subtitle       = "DESeq2 | padj < 0.05 | |LFC| > 1",
    col            = c("grey70", "grey70", "#2166AC", "#D6604D"),
    colAlpha       = 0.7,
    legendPosition = "right")
dev.off()

# 7. Heatmap 
top50_symbols <- ifelse(is.na(sig$symbol[1:50]),
                        sig$gene[1:50],
                        sig$symbol[1:50])

mat <- assay(vsd)[sig$gene[1:50], ]
mat <- mat - rowMeans(mat)
rownames(mat) <- top50_symbols

ann  <- data.frame(condition = coldata$condition, row.names = colnames(mat))
cols <- list(condition = c(normal = "#2166AC", tumor = "#D6604D"))

png("results/figures/heatmap_top50.png", width = 2000, height = 2800, res = 300)
pheatmap(mat,
         annotation_col    = ann,
         annotation_colors = cols,
         color             = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         show_rownames     = TRUE,
         show_colnames     = TRUE,
         fontsize_row      = 7,
         main              = "Top 50 DEGs - BRCA Tumor vs Normal")
dev.off()

cat("DESeq2 analysis complete. Figures saved to results/figures/\n")