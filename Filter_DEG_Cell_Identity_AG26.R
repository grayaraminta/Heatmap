#==============================================================================
# filter DEGs by Cell Identity Marker Genes
#==============================================================================

library(readr)
library(dplyr)
library(tidyr)

# set input paths -------------------------------------------------------------
path <- "/Users/aramintagray/Documents/DMG_BulkRNA/Raw_Data/"
path2 <- "/Users/aramintagray/Documents/DMG_BulkRNA/Gene_Lists/"

# 1. load DEG results =========================================================
deg_raw <- read_csv(paste0(path, "mCherry_vs_GFP_Complete.csv"))

deg <- deg_raw %>%
  mutate(gene = sub("^[^_]+_", "", Id))

cat("DEG file loaded:", nrow(deg), "genes\n")

# 2. load cell identity marker genes ==========================================
cell_id_raw <- read_csv(
  paste0(path2, "DMG_Immune.csv"),
  skip = 1,
  col_names = TRUE
)

# cell types present
cell_types <- names(cell_id_raw)[!is.na(names(cell_id_raw))]

cat("Cell types found:", paste(cell_types, collapse = ", "), "\n")

# reshape to long format (THIS PRESERVES EXCEL ORDER)
marker_long <- cell_id_raw %>%
  select(all_of(cell_types)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "cell_type",
    values_to = "marker_gene"
  ) %>%
  filter(!is.na(marker_gene))

# 🔥 CREATE EXCEL ORDER INDEX (CRITICAL FIX)
marker_long <- marker_long %>%
  mutate(order_id = row_number())

cat("Total marker gene entries:", nrow(marker_long), "\n")
cat("Unique marker genes:", n_distinct(marker_long$marker_gene), "\n\n")

# 3. match genes (case-insensitive) ===========================================
deg <- deg %>%
  mutate(gene_lower = tolower(gene))

marker_long <- marker_long %>%
  mutate(gene_lower = tolower(marker_gene))

# 4. filter DEGs to marker genes ==============================================
filtered_deg <- deg %>%
  inner_join(
    marker_long %>% select(gene_lower, cell_type, marker_gene, order_id),
    by = "gene_lower",
    relationship = "many-to-many"
  ) %>%
  select(Id, gene, marker_gene, cell_type, everything(), -gene_lower) %>%
  
  # 🔥 FORCE Excel order AFTER join
  arrange(order_id) %>%
  
  # optional cleanup
  select(-order_id)

# 5. summaries ================================================================
cat("---Results---\n")
cat("DEG genes matched to at least one cell type:", n_distinct(filtered_deg$gene), "\n")
cat("Total rows (gene x cell_type pairs):        ", nrow(filtered_deg), "\n\n")

cat("Genes matched per cell type:\n")
filtered_deg %>%
  count(cell_type, name = "n_genes") %>%
  print()

# 6. save output ==============================================================
write_csv(filtered_deg, "mCherry_vs_GFP_DEG_Filtered_DMG_Immune.csv")

cat("\nSaved: mCherry_vs_GFP_DEG_Filtered_DMG_Immune.csv\n")