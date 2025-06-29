# Load required libraries
library(readr)
library(clusterProfiler)
library(org.Rn.eg.db)
library(enrichplot)
library(ggplot2)
library(DOSE)
library(plotly)

# Read in the significant differential expression results
results <- read_csv("limma_ldopa_significant_genes.csv")

# Rename the gene symbol column if it's under ...1
results$Gene.symbol <- results$`...1`

# Extract significant genes (adjusted p-value < 0.05)
sig_genes <- subset(results, adj.P.Val < 0.05)$Gene.symbol

# Convert gene symbols to Entrez IDs
gene_ids <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Rn.eg.db)

# Perform GO enrichment analysis
go_enrichBP <- enrichGO(gene          = gene_ids$ENTREZID,
                        OrgDb         = org.Rn.eg.db,
                        ont           = "BP", 
                        readable      = TRUE,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05)

go_enrichCC <- enrichGO(gene          = gene_ids$ENTREZID,
                        OrgDb         = org.Rn.eg.db,
                        ont           = "CC", 
                        readable      = TRUE,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05)

go_enrichMF <- enrichGO(gene          = gene_ids$ENTREZID,
                        OrgDb         = org.Rn.eg.db,
                        ont           = "MF", 
                        readable      = TRUE,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05)

# Perform KEGG pathway enrichment analysis
kegg_enrich <- enrichKEGG(gene         = gene_ids$ENTREZID,
                          organism     = "rno", 
                          pvalueCutoff = 0.05)

# Ensure pairwise similarity matrices for visualization
go_enrichBP <- pairwise_termsim(go_enrichBP)
kegg_enrich <- pairwise_termsim(kegg_enrich)

################# PLOTS #################

# Dot Plot for GO
go_dotBP <- dotplot(go_enrichBP, showCategory = 10) + ggtitle("GO Enrichment (BP) - Dot Plot")
ggsave("go_dotBP.png", plot = go_dotBP, width = 8, height = 6, dpi = 300)
ggplotly(go_dotBP)

# Dot Plot for KEGG
kegg_dot <- dotplot(kegg_enrich, showCategory = 10) + ggtitle("KEGG Enrichment - Dot Plot")
ggsave("kegg_dot.png", plot = kegg_dot, width = 8, height = 6, dpi = 300)
ggplotly(kegg_dot)

# Heatmap for GO
go_heatBP <- heatplot(go_enrichBP, showCategory = 10) + ggtitle("GO Enrichment (BP) - Heatmap")
ggplotly(go_heatBP)

# Heatmap for KEGG
kegg_heat <- heatplot(kegg_enrich, showCategory = 10) + ggtitle("KEGG Enrichment - Heatmap")
ggplotly(kegg_heat)

# Static Network-Based Plots
pdf("/home/nora/Semester6/BSP6/static_network_plots.pdf", width = 8, height = 6)

# Cnet Plot for GO
cnetplot(go_enrichBP, showCategory = 5) + ggtitle("GO Enrichment (BP) Network Plot")

# Cnet Plot for KEGG
cnetplot(kegg_enrich, showCategory = 5) + ggtitle("KEGG Enrichment Network Plot")

# Tree Plot for GO
treeplot(go_enrichBP) + ggtitle("GO Enrichment (BP) Tree Plot")

# Tree Plot for KEGG
treeplot(kegg_enrich) + ggtitle("KEGG Enrichment Tree Plot")

# Enrichment Map for GO
emapplot(go_enrichBP, showCategory = 10) + ggtitle("GO Enrichment (BP) Map")

# Enrichment Map for KEGG
emapplot(kegg_enrich, showCategory = 10) + ggtitle("KEGG Enrichment Map")

# Close PDF device
dev.off()

################# SORTING #################

# Sort and filter GO results
go_results <- go_enrichBP@result
sorted_go_results <- go_results[order(go_results$p.adjust), ]
significant_go_results <- subset(sorted_go_results, p.adjust < 0.05)
significant_go_results$core_genes <- sapply(significant_go_results$geneID, function(g) paste(strsplit(g, "/")[[1]], collapse = ", "))
head(significant_go_results[, c("Description", "p.adjust", "core_genes")], 10)

# Sort and filter KEGG results
kegg_results <- kegg_enrich@result
sorted_kegg_results <- kegg_results[order(kegg_results$p.adjust), ]
significant_kegg_results <- subset(sorted_kegg_results, p.adjust < 0.05)
significant_kegg_results$core_genes <- sapply(significant_kegg_results$geneID, function(g) paste(strsplit(g, "/")[[1]], collapse = ", "))
head(significant_kegg_results[, c("Description", "p.adjust", "core_genes")], 10)


# Extract significant GO pathways (BP)
sig_go_pathways <- subset(go_enrichBP@result, p.adjust < 0.05)
sig_go_pathway_names <- sig_go_pathways$Description  # or ID if you prefer
sig_go_pathway_genes <- sapply(sig_go_pathways$geneID, function(g) unlist(strsplit(g, "/")))

# Extract significant KEGG pathways
sig_kegg_pathways <- subset(kegg_enrich@result, p.adjust < 0.05)
sig_kegg_pathway_names <- sig_kegg_pathways$Description
sig_kegg_pathway_genes <- sapply(sig_kegg_pathways$geneID, function(g) unlist(strsplit(g, "/")))

# Save these for later (optional)
save(sig_go_pathways, sig_kegg_pathways, file = "dataset1_pathway_results.RData")

################# FILES FOR NETWORK ANALYSIS #################

# 1. Gene symbols only
gene_symbols <- results["Gene.symbol"]
write.table(gene_symbols, "/home/nora/Semester6/BSP6/gene_symbols.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

gene_symbols_adj_P_Val <- results[, c("Gene.symbol", "adj.P.Val")]
write.table(gene_symbols_adj_P_Val, "/home/nora/Semester6/BSP6/gene_symbols_adj_P_Val.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


# 2. Gene symbols with logFC
gene_symbols_logfc <- results[, c("Gene.symbol", "logFC")]
write.table(gene_symbols_logfc, "/home/nora/Semester6/BSP6/gene_symbols_logfc.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# 3. Gene symbols with negative logFC
gene_symbols_neg_logfc <- gene_symbols_logfc
gene_symbols_neg_logfc$logFC <- -gene_symbols_neg_logfc$logFC
write.table(gene_symbols_neg_logfc, "/home/nora/Semester6/BSP6/gene_symbols_neg_logfc.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# 4. Gene symbols with binary regulation (1 = up, 0 = down)
gene_symbols_up_down <- gene_symbols_logfc
gene_symbols_up_down$Regulation <- ifelse(gene_symbols_up_down$logFC > 0, 1, 0)
gene_symbols_up_down <- gene_symbols_up_down[, c("Gene.symbol", "Regulation")]
write.table(gene_symbols_up_down, "/home/nora/Semester6/BSP6/gene_symbols_up_down.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
