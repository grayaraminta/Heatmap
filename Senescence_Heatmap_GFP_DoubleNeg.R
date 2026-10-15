library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)

# - 1. Load Data --------------------------------------------------------------
path <- "/Users/aramintagray/Documents/DMG_BulkRNA/Gene_Lists/"
setwd(path)

df <- read.csv(
  paste0(path, "GFP_vs_DoubleNeg_DEG_Filtered_DMG_Immune.csv")
)

df <- df[!duplicated(df[, c("gene", "cell_type")]), ]

# - SAFE UNIQUE ROW IDS -------------------------------------------------------
df$row_id <- make.unique(
  ifelse(
    duplicated(df$gene) | duplicated(df$gene, fromLast = TRUE),
    paste0(df$gene, "_", df$cell_type),
    df$gene
  )
)

# - Gene categories -----------------------------------------------------------
category_order <- c(
  "B Cell",
  "Dendritic Cell",
  "M1 Macrophage",
  "M2 Macrophage",
  "Microglia",
  "Monocyte",
  "Myeloid (General)"
)

category_labels <- c(
  "B Cell",
  "Dendritic Cell",
  "M1\nMacrophage",
  "M2\nMacrophage",
  "Microglia",
  "Monocyte",
  "Myeloid\n(General)"
)

gene_df <- data.frame(
  Cell_Identity = factor(df$cell_type, levels = category_order),
  row.names = df$row_id
)

gene_df <- gene_df[order(gene_df$Cell_Identity), , drop = FALSE]

# - Expression matrix ---------------------------------------------------------
norm_cols <- c(
  "norm.1", "norm.4", "norm.10",
  "norm.3", "norm.6", "norm.9", "norm.12"
)

mat <- as.matrix(df[, norm_cols])
rownames(mat) <- df$row_id

colnames(mat) <- c(
  "GFP_M1", "GFP_M2", "GFP_M4",
  "DoubleNeg_M1", "DoubleNeg_M2", "DoubleNeg_M3", "DoubleNeg_M4"
)

# reorder genes
mat <- mat[rownames(gene_df), ]

# log transform
mat_log <- log2(mat + 1)

# sample order
sample_order <- c(
  "DoubleNeg_M1", "DoubleNeg_M2", "DoubleNeg_M3", "DoubleNeg_M4",
  "GFP_M1", "GFP_M2", "GFP_M4"
)

mat_log <- mat_log[, sample_order]

# clean sample labels
colnames(mat_log) <- c(
  "M1", "M2", "M3", "M4",
  "M1", "M2", "M4"
)

group <- factor(
  c(rep("DoubleNeg", 4), rep("GFP", 3)),
  levels = c("DoubleNeg", "GFP")
)

# - TRANSPOSE + Z-SCORE -------------------------------------------------------
mat_t <- t(mat_log)
mat_z <- as.matrix(scale(mat_t))

# remove zero-variance genes
keep <- colSums(is.na(mat_z)) == 0
mat_z <- mat_z[, keep, drop = FALSE]

gene_df <- gene_df[colnames(mat_z), , drop = FALSE]

# -------------------- Chil3 CLEAN LABEL FIX ---------------------------------
clean_names <- colnames(mat_z)

clean_names <- ifelse(
  grepl("^Itgam", clean_names),
  "Itgam",
  clean_names
)

clean_names <- ifelse(
  grepl("^Aif1", clean_names),
  "Aif1",
  clean_names
)

colnames(mat_z) <- clean_names

# - COLORS --------------------------------------------------------------------
col_fun <- colorRamp2(
  seq(-2, 2, length.out = 11),
  rev(brewer.pal(11, "RdYlBu"))
)

# - PLOT SIZE (SQUARE CELL CONTROL) -------------------------------------------
n_genes <- ncol(mat_z)

cell_size <- 0.25  # tweak this (0.2–0.35)

png(
  paste0(path, "heatmap_transposed.png"),
  width = n_genes * cell_size + 3,
  height = 3.5,
  units = "in",
  res = 600
)

ht <- Heatmap(
  mat_z,
  name = "Z-score",
  col = col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  row_split = group,
  
  row_title = c(
    "Double Negative",
    expression(GFP^"+")
  ),
  
  row_title_gp = gpar(
    fontsize = 10,
    fontface = "bold",
    col = c("black", "#2E6F40")
  ),
  
  column_split = gene_df$Cell_Identity,
  
  column_title = category_labels,
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  
  row_gap = unit(1.5, "mm"),
  column_gap = unit(1.5, "mm"),
  
  border = TRUE,
  
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 90,
  
  heatmap_legend_param = list(
    title = "Z-score",
    legend_height = unit(3, "cm")
  )
)

draw(ht)
dev.off()