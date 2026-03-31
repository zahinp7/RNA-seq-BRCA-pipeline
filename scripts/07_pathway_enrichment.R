# Stage 7: Pathway Enrichment Analysis
# GO enrichment + KEGG + GSEA

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(dplyr)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/enrichment", recursive = TRUE, showWarnings = FALSE)

# 1. Load DESeq2 results 
res_df <- read.csv("results/deseq2/deseq2_results_annotated.csv", row.names = 1)
sig    <- read.csv("results/deseq2/significant_genes.csv", row.names = 1)

# Convert Ensembl IDs to Entrez IDs
entrez_all <- mapIds(org.Hs.eg.db,
                     keys    = res_df$gene,
                     column  = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")

entrez_sig <- mapIds(org.Hs.eg.db,
                     keys    = sig$gene,
                     column  = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")

# 2. GO Enrichment (over-representation) 
go_bp <- enrichGO(gene          = na.omit(entrez_sig),
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)

write.csv(as.data.frame(go_bp), "results/enrichment/go_bp_results.csv")

# GO dotplot
p_go <- dotplot(go_bp, showCategory = 20) +
    ggtitle("GO Biological Process - BRCA Tumor vs Normal") +
    theme(axis.text.y = element_text(size = 8))
ggsave("results/figures/go_dotplot.png", p_go, width = 10, height = 8, dpi = 300)

# 3. KEGG Enrichment 
kegg <- enrichKEGG(gene         = na.omit(entrez_sig),
                   organism     = "hsa",
                   pvalueCutoff = 0.05)

write.csv(as.data.frame(kegg), "results/enrichment/kegg_results.csv")

p_kegg <- dotplot(kegg, showCategory = 20) +
    ggtitle("KEGG Pathways - BRCA Tumor vs Normal") +
    theme(axis.text.y = element_text(size = 8))
ggsave("results/figures/kegg_dotplot.png", p_kegg, width = 10, height = 8, dpi = 300)

# 4. GSEA 
# Ranked gene list by log2FoldChange
res_ranked <- res_df %>%
    filter(!is.na(log2FoldChange)) %>%
    mutate(entrez = entrez_all[gene]) %>%
    filter(!is.na(entrez)) %>%
    arrange(desc(log2FoldChange))

# Remove duplicate entrez IDs - keep highest absolute LFC
res_ranked <- res_ranked %>%
    group_by(entrez) %>%
    slice_max(abs(log2FoldChange), n = 1) %>%
    ungroup() %>%
    arrange(desc(log2FoldChange))

gene_list <- res_ranked$log2FoldChange
names(gene_list) <- res_ranked$entrez

gsea_go <- gseGO(geneList     = gene_list,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 minGSSize    = 10,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)

write.csv(as.data.frame(gsea_go), "results/enrichment/gsea_go_results.csv")

p_gsea <- dotplot(gsea_go, showCategory = 20, split = ".sign") +
    facet_grid(. ~ .sign) +
    ggtitle("GSEA GO BP - BRCA Tumor vs Normal") +
    theme(axis.text.y = element_text(size = 7))
ggsave("results/figures/gsea_dotplot.png", p_gsea, width = 14, height = 8, dpi = 300)

cat("Pathway enrichment complete\n")
cat("GO BP terms:", nrow(as.data.frame(go_bp)), "\n")
cat("KEGG pathways:", nrow(as.data.frame(kegg)), "\n")
cat("GSEA GO terms:", nrow(as.data.frame(gsea_go)), "\n")