library(GEOquery)
library(dplyr)

# Load datasets
gset1 <- getGEO("GSE139438", GSEMatrix = TRUE, AnnotGPL = TRUE)
if (length(gset1) > 1) idx1 <- grep("GPL6101", attr(gset1, "names")) else idx1 <- 1
gset1 <- gset1[[idx1]]
feature_data1 <- fData(gset1)
full_genes1 <- unique(na.omit(as.character(feature_data1$`Gene symbol`)))

gset2 <- getGEO("GSE55096", GSEMatrix = TRUE, AnnotGPL = TRUE)
if (length(gset2) > 1) idx2 <- grep("GPL1261", attr(gset2, "names")) else idx2 <- 1
gset2 <- gset2[[idx2]]
feature_data2 <- fData(gset2)
full_genes2 <- unique(na.omit(as.character(feature_data2$`Gene symbol`)))

# Load significant genes with statistics
sig_genes1_df <- read.delim("/home/nora/Semester5/BSP5/FINAL/Code/differential_expression_results.tsv", 
                            stringsAsFactors = FALSE, header = TRUE)
sig_genes2_df <- read.csv("/home/nora/Semester6/BSP6/CODE/limma_ldopa_significant_genes.csv", 
                          stringsAsFactors = FALSE, header = TRUE)

# Filter to background genes
sig_genes1_df <- sig_genes1_df[sig_genes1_df$Gene.symbol %in% full_genes1, ]
sig_genes2_df <- sig_genes2_df[sig_genes2_df$Gene.title %in% full_genes2, ]

# Define significant gene lists
sig_genes1 <- unique(na.omit(as.character(sig_genes1_df$Gene.symbol)))
sig_genes2 <- unique(na.omit(as.character(sig_genes2_df$Gene.title)))

# Define the universe for Fisher's test (union of both dataset backgrounds)
universe <- unique(c(full_genes1, full_genes2))

# Calculate overlap and exclusive sets for contingency table
overlap <- length(intersect(sig_genes1, sig_genes2))                 # Sig in both
only1   <- length(setdiff(sig_genes1, sig_genes2))                   # Sig only in file1
only2   <- length(setdiff(sig_genes2, sig_genes1))                   # Sig only in file2
rest    <- length(setdiff(universe, union(sig_genes1, sig_genes2)))  # Not sig in either

# Construct contingency table
contingency_table <- matrix(c(overlap, only1, only2, rest),
                            nrow = 2,
                            byrow = TRUE)
rownames(contingency_table) <- c("File1_Sig", "File1_Not_Sig")
colnames(contingency_table) <- c("File2_Sig", "File2_Not_Sig")

# Perform Fisher's exact test
fisher_res <- fisher.test(contingency_table)

# Print results
print("Contingency Table:")
print(contingency_table)
print("Fisher's Exact Test Result:")
print(fisher_res)

# Extract overlapping genes with their stats
overlap_genes <- intersect(sig_genes1, sig_genes2)

overlap_df1 <- sig_genes1_df %>%
  filter(Gene.symbol %in% overlap_genes) %>%
  select(Gene.symbol, adj.P.Val, logFC) %>%
  rename(AdjPVal_File1 = adj.P.Val, LogFC_File1 = logFC)

overlap_df2 <- sig_genes2_df %>%
  filter(Gene.title %in% overlap_genes) %>%
  select(Gene.title, adj.P.Val, logFC) %>%
  rename(Gene.symbol = Gene.title, AdjPVal_File2 = adj.P.Val, LogFC_File2 = logFC)

# Merge stats from both files
merged_overlap <- merge(overlap_df1, overlap_df2, by = "Gene.symbol")

# Save merged overlap with Fisher test summary
write.table(merged_overlap, "overlapping_genes_with_stats.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
