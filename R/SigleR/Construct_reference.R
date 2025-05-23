ref <- readRDS("./singleRpan.rds")

# Data normalization steps for the reference
ref <- NormalizeData(ref)
ref <- FindVariableFeatures(ref)
ref <- ScaleData(ref)

# Add index to reference metadata for merging
num_rows <- nrow(ref@meta.data)
new_values <- seq(0, num_rows - 1)
ref@meta.data[["Index"]] <- new_values

# Extract count matrix and metadata for reference
ref_count <- ref[["RNA"]]@counts
pdata <- ref@meta.data[, c("Index", "cell_type")]
rownames(pdata) <- pdata$Index
pdata$Index <- NULL
colnames(pdata) <- "ref_label"

# Create SummarizedExperiment object
ref_SE <- SummarizedExperiment(
  assays = list(counts = ref_count),
  colData = pdata
)
ref_SE <- logNormCounts(ref_SE)

# Prepare test dataset
panc8 <- readRDS("./zheng68k.rds")
pancounts <- panc8[["RNA"]]$counts

# Find common genes between test and reference
common_genes <- intersect(rownames(panc8), rownames(ref_SE))
pancounts <- pancounts[common_genes, ]
ref_SE <- ref_SE[common_genes, ]

# Create SummarizedExperiment for test data
pan_se <- SummarizedExperiment(
  assays = list(counts = pancounts)
)
pan_se <- logNormCounts(pan_se)

# Run SingleR with custom reference
singleR_ref <- SingleR(
  test = pan_se, 
  ref = ref_SE, 
  labels = ref_SE$ref_label
)

# Prepare annotation results for merging
anno_df <- as.data.frame(singleR_ref$labels)
anno_df$Index <- rownames(singleR_ref)
colnames(anno_df)[1] <- 'ref_label_from_pan8'

# Add annotations to Seurat object metadata
panc8@meta.data <- merge(
  panc8@meta.data,
  anno_df,
  by = "Index",
  all.x = TRUE
)

# Set row names for metadata
rownames(panc8@meta.data) <- panc8@meta.data$Index

# Visualize annotations
DimPlot(panc8, reduction = "umap", group.by = "ref_label_from_pan8")