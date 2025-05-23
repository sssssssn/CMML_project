library(SCINA)

# Read marker gene data
pan_reference <- read.csv("/Users/shennuo/Desktop/pan reference.csv")

# Convert marker genes to a vector
marker_genes_vector <- as.vector(t(pan_reference))

# Get the normalized matrix
normalized_matrix <- as.matrix(GetAssayData(panc8, layer = "data"))

# Extract the subset matrix containing only marker genes
subset_matrix <- normalized_matrix[rownames(normalized_matrix) %in% marker_genes_vector, ]

# Write the subset matrix to a CSV file
write.csv(subset_matrix, "/Users/shennuo/Desktop/scina.csv")

# Load the SCINA analysis results
analysis_results <- load("/Users/shennuo/Downloads/results.RData")

# Run the SCINA analysis
scina_analysis <- SCINA(
  subset_matrix,
  pan_reference,
  max_iter = 100,
  convergence_n = 10,
  convergence_rate = 0.99,
  sensitivity_cutoff = 1,
  rm_overlap = 1,
  allow_unknown = 1,
  log_file = "SCINA.log"
)

# Add the analysis results to the metadata
panc8@meta.data[["scina_analysis"]] <- scina_analysis[["cell_labels"]]