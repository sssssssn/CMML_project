library(SingleR)
library(celldex)
library(Seurat)

# Reference datasets
hpca.se <- HumanPrimaryCellAtlasData() # Load Human Primary Cell Atlas dataset
BP <- BlueprintEncodeData()            # Load Blueprint Encode dataset

# Get normalized data matrix from Seurat object
pan_normalized_matrix <- GetAssayData(panc8, layer = "data")

# Run SingleR using HPCA as reference
pan.hesc <- SingleR(
  test = pan_normalized_matrix, 
  ref = hpca.se, 
  labels = hpca.se$label.main
)

# Add SingleR labels to Seurat object metadata
panc8@meta.data$labels <- pan.hesc$labels

# Plot UMAP with annotations
print(DimPlot(panc8, group.by = "labels", reduction = "umap"))

# Table of predicted labels vs original cell types
meta <- panc8@meta.data
table(pan.hesc$labels, meta$celltype)

