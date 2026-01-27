# Global variable declarations to silence R CMD check NOTES for non-standard evaluation
# Used across dplyr/tidy pipelines and internal helpers
utils::globalVariables(c(
  ".data", "score.x", "score.y", "id.x", "id.y", "score", "score0", "true",
  "celltype", "group", "dissimilarity_method", "consistency.method",
  "consistency.metric", "data.type", "clustering_method.", "product", "parent",
  "optimal", "optimal0", "sil_width", "mean_sil_width", "measure", "value", "k",
  "rate", "pval", "score_celltype", "score_mean", "label", "pruned.labels", "x",
  "y", "FDR", "summary.logFC"
))
