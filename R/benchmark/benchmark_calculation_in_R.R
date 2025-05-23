library(dplyr)
library(purrr)
library(stringr)

data = read.csv("./pancreatic_obs.csv")

# Create a mapping list to map the predicted values to the real labels
#example
mapping <- list(
  "Gamma_Cell" = "gamma",
  "Alpha_Cell" = "alpha",
  "Ductal_Cell" = "ductal",
  "Beta_Cell" = "beta",
  "Acinar_Cell" = "acinar",
  "Delta_Cell" = "delta",
  "Endothelial_Cell" = "endothelial",
  "Macrophage" = "macrophage",
  "Epsilon_Cell" = "epsilon"
)


# Map the predicted labels back to the corresponding true labels
data$predicted_mapped <- sapply(data$celltype_hint, function(x) {
  mapped_value <- mapping[[x]]
  if (is.null(mapped_value)) {
    return(NA)  
  }
  return(mapped_value)
})

# Obtain all the actual labels
actual_labels <- unique(data$combined_label)

# Initialize the list to store the metrics of each category
metrics <- data.frame(
  CellType = character(),
  Precision = numeric(),
  Recall = numeric(),
  F1 = numeric(),
  stringsAsFactors = FALSE
)

# Calculate the precision rate, recall rate and F1 value of each category
for (i in seq_along(actual_labels)) {
  actual_label <- actual_labels[i]
  
  # TP：the actual label is actual_label and the predicted label is actual_label
  TP <- sum(data$combined_label == actual_label & data$predicted_mapped == actual_label, na.rm = TRUE)
  
  # FP：the predicted label is actual_label, but the actual label is not actual_label
  FP <- sum(data$predicted_mapped == actual_label & data$combined_label != actual_label, na.rm = TRUE)
  
  # FN：the actual label is actual_label, but the predicted label is not actual_label
  FN <- sum(data$combined_label == actual_label & data$predicted_mapped != actual_label, na.rm = TRUE)
  
  
  # Calculate the precision rate and recall rate
  precision <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
  recall <- ifelse(TP + FN == 0, 0, TP / (TP + FN))
  
  # Calculate the F1-score
  f1 <- ifelse(precision + recall == 0, 0, 2 * (precision * recall) / (precision + recall))
  
  
  # add results
  metrics <- rbind(metrics, data.frame(
    CellType = actual_label,
    Precision = precision,
    Recall = recall,
    F1 = f1,
  ))
}

# Calculate the overall average metrics
overall_metrics <- data.frame(
  CellType = "Overall",
  Precision = mean(metrics$Precision, na.rm = TRUE),
  Recall = mean(metrics$Recall, na.rm = TRUE),
  F1 = mean(metrics$F1, na.rm = TRUE)
)

metrics <- rbind(metrics, overall_metrics)

