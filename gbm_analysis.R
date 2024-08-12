# GBM RNA-Seq Analysis Script

# Setup -------------------------------------------------------------------
setwd("/Users/Derleen/Desktop/Hackbio_courses/Projects/GBM_Project")
.libPaths("/Users/Derleen/Desktop/Hackbio_courses/Hackbio_test")

# Install and load packages -----------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

packages <- c("TCGAbiolinks", "SummarizedExperiment", "edgeR", "limma", 
              "ggplot2", "pheatmap", "org.Hs.eg.db", "clusterProfiler")

BiocManager::install(packages)
sapply(packages, library, character.only = TRUE)

# Data Acquisition --------------------------------------------------------
query <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query)
data <- GDCprepare(query)

# Data Preprocessing ------------------------------------------------------
counts <- assay(data)
sample_info <- colData(data)

dge <- DGEList(counts = counts, genes = rowData(data))
keep <- filterByExpr(dge, group = sample_info$sample_type)
dge <- dge[keep, ]
dge <- calcNormFactors(dge)

# Differential Expression Analysis ----------------------------------------
design <- model.matrix(~ sample_info$sample_type)
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit)

results <- topTags(lrt, n = Inf)
results_df <- as.data.frame(results)
sig_genes <- results_df[results_df$FDR < 0.05 & abs(results_df$logFC) > 1, ]

# Pathway Enrichment Analysis ---------------------------------------------
gene_list <- sig_genes$logFC
names(gene_list) <- rownames(sig_genes)
gene_list <- sort(gene_list, decreasing = TRUE)

remove_version <- function(ids) {
  gsub("\\.[0-9]+$", "", ids)
}
names(gene_list) <- remove_version(names(gene_list))

gene_df <- bitr(names(gene_list), 
                fromType = "ENSEMBL", 
                toType = c("ENTREZID", "SYMBOL"), 
                OrgDb = org.Hs.eg.db)
gene_list <- gene_list[gene_df$ENSEMBL]

go_results <- enrichGO(gene = gene_df$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

kegg_results <- enrichKEGG(gene = gene_df$ENTREZID,
                           organism = 'hsa',
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2)

# Visualization -----------------------------------------------------------
top_gene_ids <- rownames(sig_genes)[1:20]
valid_top_genes <- intersect(top_gene_ids, rownames(counts))

if (length(valid_top_genes) == 0) {
  stop("No valid top genes found in the expression data.")
}

heatmap_data <- counts[valid_top_genes, ]
heatmap_data[is.na(heatmap_data) | is.nan(heatmap_data) | is.infinite(heatmap_data)] <- 0

if (ncol(heatmap_data) != length(sample_info$sample_type)) {
  stop("Mismatch between the number of samples in heatmap_data and sample_info$sample_type.")
}

heatmap_data <- t(scale(t(heatmap_data)))
annotation_col <- data.frame(sample_type = sample_info$sample_type)
rownames(annotation_col) <- colnames(heatmap_data)

pheatmap(heatmap_data,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         filename = "heatmap_top_genes.png")

pdf("GO_dotplot.pdf", width = 12, height = 8)
dotplot(go_results, showCategory = 20)
dev.off()

pdf("KEGG_dotplot.pdf", width = 12, height = 8)
dotplot(kegg_results, showCategory = 20)
dev.off()

# Save Results ------------------------------------------------------------
write.csv(results_df, "differential_expression_results.csv")
write.csv(as.data.frame(go_results), "go_enrichment_results.csv")
write.csv(as.data.frame(kegg_results), "kegg_pathway_enrichment_results.csv")

print("Analysis complete. Results and visualizations have been saved.")