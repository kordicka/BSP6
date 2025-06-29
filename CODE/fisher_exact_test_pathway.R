# Load dataset 1
load("/home/nora/Semester6/BSP6/CODE/dataset1_pathway_results.RData")  # loads sig_go_pathways and sig_kegg_pathways
sig_go_pathways1 <- sig_go_pathways
sig_kegg_pathways1 <- sig_kegg_pathways

# Load dataset 2
load("/home/nora/Semester6/BSP6/CODE/dataset2_pathway_results.RData")  # loads sig_go_pathways and sig_kegg_pathways
sig_go_pathways2 <- sig_go_pathways
sig_kegg_pathways2 <- sig_kegg_pathways

##### Function to extract overlaps #####
extract_overlaps <- function(df1, df2, analysis_name) {
  # Find overlapping pathway IDs
  overlap_ids <- intersect(df1$ID, df2$ID)
  
  # Extract relevant columns from both datasets
  df1_overlap <- df1[df1$ID %in% overlap_ids, c("ID", "Description", "GeneRatio", "p.adjust")]
  df2_overlap <- df2[df2$ID %in% overlap_ids, c("ID", "Description", "GeneRatio", "p.adjust")]
  
  # Merge the two dataframes by ID
  merged <- merge(
    df1_overlap,
    df2_overlap,
    by = "ID",
    suffixes = c("_Analysis1", "_Analysis2")
  )
  
  # Add Type column
  merged$Type <- analysis_name
  
  # Reorder columns
  merged <- merged[, c(
    "Type",
    "ID",
    "Description_Analysis1",
    "Description_Analysis2",
    "GeneRatio_Analysis1",
    "p.adjust_Analysis1",
    "GeneRatio_Analysis2",
    "p.adjust_Analysis2"
  )]
  
  return(merged)
}

# Extract overlapping GO pathways
overlap_go <- extract_overlaps(sig_go_pathways1, sig_go_pathways2, "GO")

# Extract overlapping KEGG pathways
overlap_kegg <- extract_overlaps(sig_kegg_pathways1, sig_kegg_pathways2, "KEGG")

# Combine both into one dataframe
overlap_all <- rbind(overlap_go, overlap_kegg)

# Print the overlapping results for inspection
print(overlap_all)

# Write to CSV for your thesis tables
write.csv(overlap_all, "overlapping_GO_KEGG_pathways_comparison.csv", row.names = FALSE)
