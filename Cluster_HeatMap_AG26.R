library (pheatmap)
library(RColorBrewer)

# - 1. Load Data --------------------------------------------------------------
path <- "/Users/aramintagray/Documents/DMG_BulkRNA/Gene_Lists/"
setwd(path)
df <- read.csv(paste0(path, "mCherry_vs_DoubleNeg_DEG_Filtered_Senescence_Scott.csv"))

# remove fully duplicated rows (same gene AND same cell type)
df <- df[!duplicated(df[, c("gene", "cell_type")]), ]

# Create unique row identifiers: gene name + cell type (handles duplicates)
df$row_id <- ifelse(
  duplicated(df$gene) | duplicated(df$gene, fromLast = TRUE),
  paste0(df$gene, "_", df$cell_type),
  df$gene
)

# - 2. DataFrame1: Gene -> Cell Identity --------------------------------------
gene_df <- data.frame(Cell_Identity = df$cell_type, row.names = df$row_id)

# - 3. DataFrame2: Samples -> Group + Normalised Values -----------------------
#extract normalised column counts
norm_cols <- c("norm.2", "norm.5", "norm.8", "norm.11",
               "norm.3", "norm.6", "norm.9", "norm.12")

#build matrix: genes x samples
mat <- as.matrix(df[, norm_cols])
rownames(mat) <- df$row_id

#rename columns to meaningful mouse IDs
colnames(mat) <- c("mCherry_M1", "mCherry_M2", "mCherry_M3", "mCherry_M4",
                  "DoubleNeg_M1", "DoubleNeg_M2", "DoubleNeg_M3", "DoubleNeg_M4")

#log-transform for visualisation (log2+1 to handle zeroes)
mat_log <- log2(mat+1)

#sample annotation dataframe
ann_col <- data.frame(
  Group = c("mCherry", "mCherry", "mCherry", "mCherry",
            "DoubleNeg", "DoubleNeg", "DoubleNeg", "DoubleNeg")
)
rownames(ann_col) <- colnames(mat_log)

# - 4. Annotation colours -----------------------------------------------------
ann_colours <-list(
  Cell_Identity = c(
    "Negative Regulation of the Cell Cycle" = "#4575B4",
    "Positive Regulation of the Cell Cycle" = "#74ADD1",
    "Lysosomal Compartment Expansion" = "#ABD9E9",
    "Positive Regulation of Cellular Senescence" = "#708238",
    "Inflammatory Response" = "#D73027",
    "Stress Response" = "#F46D43"
  ),
  Group = c(
    "mCherry" = "#C51B7D",
    "DoubleNeg" = "#4D4D4D"
  )
)

# - 5. Plot heatmap -----------------------------------------------------------
# Sort rows by cell identity so genes are grouped together
row_order <- order(gene_df$Cell_Identity)
mat_log <- mat_log[row_order, ]
gene_df <- gene_df[row_order, , drop = FALSE]

heat_plot <- pheatmap(
  mat_log,
  col = colorRampPalette(rev(brewer.pal(10, "RdYlBu"))) (100),
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  clustering_distance_cols = "euclidean",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D",
  annotation_row = gene_df, #gene -> cell identity
  annotation_col = ann_col, #sample -> group
  annotation_colors = ann_colours,
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  fontsize_row = 6,
  fontsize_col = 6, 
  angle_col = 45,
  show_colnames = TRUE,
  show_rownames = TRUE,
  scale = "row", #z-score per gene to highlight relative differences
  legend_breaks = c(-2, 0, 2),
  legend_labels = c("Low", "Medium", "High"),
  fontsize_legend = 6,
  main = "Normalised Expression per Mouse: mCherry vs. Double Neg"
)