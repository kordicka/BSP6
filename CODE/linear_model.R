# Load required libraries
library(GEOquery)
library(limma)
library(edgeR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# -------------------------------------------------------------------------------
# 1. Load dataset
# -------------------------------------------------------------------------------

# Download and load the dataset from GEO
gset <- getGEO("GSE55096", GSEMatrix = TRUE, AnnotGPL = TRUE)
if (length(gset) > 1) idx <- grep("GPL1261", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# Ensure proper column names
fvarLabels(gset) <- make.names(fvarLabels(gset))

# Extract expression data and feature (gene) data
exprs_data <- exprs(gset)
feature_data <- fData(gset)

# -------------------------------------------------------------------------------
# 2. Prepare sample metadata
# -------------------------------------------------------------------------------

sample_info <- pData(gset)

# Exclude the ascorbate lesion type
sample_info <- sample_info[sample_info$`lesion:ch1` != "ascorbate", ]
exprs_data <- exprs_data[, rownames(sample_info)]  # keep only matching columns

# Define treatment groups
sample_info$group <- factor(sample_info$`drug treatment:ch1`,
                            levels = c("Chronic saline", "Chronic low levodopa", "Chronic high levodopa"),
                            labels = c("PD model", "low L-DOPA", "high L-DOPA"))

# Assign numeric L-DOPA concentrations(based on the doses given to the mice)
sample_info$ldopa_conc <- as.numeric(recode(sample_info$group,
                                            "PD model" = 0,
                                            "low L-DOPA" = 2,
                                            "high L-DOPA" = 6))

# -------------------------------------------------------------------------------
# 3. Select transcript with highest average expression per gene
# -------------------------------------------------------------------------------

gene_column <- "Gene.symbol"

# Create a dataframe with expression and gene symbols
exprs_df <- as.data.frame(exprs_data)
exprs_df$Gene.symbol <- feature_data[rownames(exprs_data), gene_column]

# Remove rows with missing or empty gene symbols
exprs_df <- exprs_df[!is.na(exprs_df$Gene.symbol) & exprs_df$Gene.symbol != "", ]

# Add rowMeans to identify the most expressed probe per gene
exprs_df$avg_expr <- rowMeans(exprs_df[, 1:(ncol(exprs_df) - 2)])  # exclude Gene.symbol & avg_expr columns

# Arrange by descending average expression and keep the top probe per gene
exprs_top_probe <- exprs_df %>%
  arrange(desc(avg_expr)) %>%
  distinct(Gene.symbol, .keep_all = TRUE)

# Set gene symbols as rownames and keep only expression columns
rownames(exprs_top_probe) <- exprs_top_probe$Gene.symbol
exprs_data <- as.matrix(exprs_top_probe[, 1:(ncol(exprs_top_probe) - 2)])  # drop Gene.symbol and avg_expr

# -------------------------------------------------------------------------------
# 4. Create DGEList object and filter low-expressed genes
# -------------------------------------------------------------------------------

dge <- DGEList(counts = exprs_data)

# Filter genes expressed in at least 6 samples (cpm > 1)
keep <- rowSums(cpm(dge) > 1) >= 6
dge <- dge[keep, ]

# Normalize using TMM
dge <- calcNormFactors(dge, method = "TMM")

# -------------------------------------------------------------------------------
# 5. Design matrix using L-DOPA concentration
# -------------------------------------------------------------------------------

design <- model.matrix(~ ldopa_conc, data = sample_info)
colnames(design) <- c("Intercept", "L_DOPA_Linear")

# -------------------------------------------------------------------------------
# 6. Fit linear model using limma-voom
# -------------------------------------------------------------------------------

v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

# -------------------------------------------------------------------------------
# 7. Extract results and identify differentially expressed genes
# -------------------------------------------------------------------------------

results <- topTable(fit, coef = "L_DOPA_Linear", n = Inf, sort.by = "P")
results$direction <- ifelse(results$adj.P.Val >= 0.05, "Not Significant",
                            ifelse(results$logFC > 0, "Upregulated", "Downregulated"))

results$Gene.title <- rownames(results)
sig_genes <- results[results$adj.P.Val < 0.05, ]

cat("Total genes tested:", nrow(results), "\n")
cat("Significantly differentially expressed:", nrow(sig_genes), "\n")
cat("Upregulated:", sum(sig_genes$logFC > 0), " | Downregulated:", sum(sig_genes$logFC < 0), "\n")

# -------------------------------------------------------------------------------
# 8. Box plot for initial distribution of gene expression levels
# -------------------------------------------------------------------------------

exprs_df <- data.frame(exprs_data) 
exprs_df$Gene.symbol <- rownames(exprs_df) 

long_exprs_df <- exprs_df %>%
  gather(key = "Sample", value = "Expression", -Gene.symbol) %>%
  mutate(Group = rep(sample_info$group, each = nrow(exprs_df)))

box_plot <- ggplot(long_exprs_df, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot() +
  labs(
    title = "Initial Gene Expression Distribution by Treatment Group",
    x = "Treatment Group",
    y = "Gene Expression (Log2 Counts)"
  ) +
  theme_bw() +
  theme(legend.position = "none")

print(box_plot)

# -------------------------------------------------------------------------------
# 9. Boxplot visualization for top genes across groups
# -------------------------------------------------------------------------------

top_genes <- rownames(sig_genes)[1:9]

boxplot_data <- v$E[top_genes, ] %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Expression")

boxplot_data$Group <- sample_info[boxplot_data$Sample, "group"]

ggplot(boxplot_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot() +
  facet_wrap(~Gene, scales = "free_y") +
  labs(
    title = "Expression of Top DE Genes Across Treatment Groups",
    x = "Group",
    y = "Normalized Expression (log-CPM)"
  ) +
  theme_bw() +
  theme(legend.position = "none")

# -------------------------------------------------------------------------------
# 10. Volcano Plot
# -------------------------------------------------------------------------------

volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = direction)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c("Downregulated" = "blue", "Upregulated" = "red", "Not Significant" = "gray")) +
  labs(
    title = "Volcano Plot of Differential Expression",
    x = "Log2 Fold Change",
    y = "-log10(adj. P-Value)",
    color = "Direction"
  ) +
  theme_minimal()

print(volcano_plot)

# -------------------------------------------------------------------------------
# 11. Save results
# -------------------------------------------------------------------------------

write.csv(results, "limma_ldopa_concentration_results.csv", row.names = TRUE)
write.csv(sig_genes, "limma_ldopa_significant_genes.csv", row.names = TRUE)

# -------------------------------------------------------------------------------
# 12. Session Info
# -------------------------------------------------------------------------------

sessionInfo()
# -------------------------------------------------------------------------------
# CHECKING FOR BATCH EFFECTS OR CONFOUNDERS
# -------------------------------------------------------------------------------

# Load additional library
library(pheatmap)

# Check which potential batch/confounder columns exist in sample_info
colnames(sample_info)

# Example: plot PCA by treatment group
v <- voom(dge, design, plot = FALSE)

# Calculate PCA
pca <- prcomp(t(v$E), scale. = TRUE)

# Combine PCA with sample info
pca_df <- data.frame(PC1 = pca$x[,1],
                     PC2 = pca$x[,2],
                     group = sample_info$group)

# Plot PCA colored by treatment group
ggplot(pca_df, aes(x=PC1, y=PC2, color=group)) +
  geom_point(size=3) +
  labs(title="PCA plot colored by Treatment Group") +
  theme_bw()

# -------------------------------------------------------------------------------
# Example: plot PCA by a potential batch variable (change 'batch_column' below)
# -------------------------------------------------------------------------------

# Replace 'batch_column' with actual column names in your sample_info,
# e.g. 'batch' or 'date:ch1' if available.

batch_column <- "growth_protocol_ch1"  # CHANGE THIS to your actual batch/confounder column name

if (batch_column %in% colnames(sample_info)) {
  pca_df$batch <- sample_info[[batch_column]]
  
  ggplot(pca_df, aes(x=PC1, y=PC2, color=batch)) +
    geom_point(size=3) +
    labs(title=paste("PCA plot colored by", batch_column)) +
    theme_bw()
}

# -------------------------------------------------------------------------------
# Sample-to-sample hierarchical clustering heatmap
# -------------------------------------------------------------------------------

sample_dists <- dist(t(v$E))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- colnames(v$E)
colnames(sample_dist_matrix) <- colnames(v$E)

pheatmap(sample_dist_matrix,
         clustering_distance_rows=sample_dists,
         clustering_distance_cols=sample_dists,
         main="Sample-to-sample distance heatmap")

# -------------------------------------------------------------------------------
# Interpretation:
# - If samples cluster by treatment group: desired
# - If samples cluster by an unintended variable: batch effect suspected
# -------------------------------------------------------------------------------