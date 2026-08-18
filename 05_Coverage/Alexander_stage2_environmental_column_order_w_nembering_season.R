###############################################################################
# ALEXANDER ESTUARY — STAGE 2
# VISUALIZATION WITH MANUAL COLUMN ORDER (NO CLUSTERING)
###############################################################################

library(dplyr)
library(readr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

###############################################################################
# 1. WORKING DIRECTORY
###############################################################################

setwd("C:/Users/user/Desktop/Research Project Oded/metadata/Alexander/New")

###############################################################################
# 2. LOAD EXPORTED OBJECTS FROM STAGE 1
###############################################################################

coverage_mat <- read.delim(
  "alexander_coverage_filt_for_clustering.tsv",
  check.names = FALSE,
  row.names = 1
)

coverage_mat <- as.matrix(coverage_mat)

# TRUE row order
row_order_list <- readRDS("alexander_row_order_by_cluster_slice.rds")

# TRUE original cluster object  ← HERE
plasmid_clusters <- readRDS("alexander_plasmid_clusters_original.rds")

###############################################################################
# 3. RECONSTRUCT TRUE ROW ORDER (CORRECT WAY)
###############################################################################

# row_order_list contains rownames per slice
row_order_names <- unlist(
  lapply(row_order_list, function(x) {
    rownames(coverage_mat)[x]
  })
)

# Reorder matrix by rownames
coverage_mat <- coverage_mat[row_order_names, , drop = FALSE]

# Reorder cluster vector by rownames
plasmid_clusters <- plasmid_clusters[row_order_names]
plasmid_clusters <- factor(
  plasmid_clusters,
  levels = unique(plasmid_clusters)
)
###############################################################################
# 4. LOAD METADATA
###############################################################################

alex_meta <- read.csv(
  "alex_clean_metadata.csv",
  stringsAsFactors = FALSE
)

###############################################################################
# 5. DEFINE NEW COLUMN ORDER: SEASON → REGION
###############################################################################

region_levels <- c(
  "Upstream",
  "Midstream",
  "Urban_Tributary",
  "Lower_River",
  "Estuary"
)

season_levels <- c(
  "Winter",
  "Spring",
  "LateSummer"
)

alex_meta <- alex_meta %>%
  filter(sample_id %in% colnames(coverage_mat)) %>%
  mutate(
    season = factor(season, levels = season_levels),
    region = factor(region, levels = region_levels),
    Temperature = factor(Temperature, levels = c(19, 30))
  ) %>%
  arrange(season, region)

# Apply manual column order (ROWS NOT TOUCHED)
coverage_mat <- coverage_mat[, alex_meta$sample_id, drop = FALSE]

###############################################################################
# 6. REBUILD ANNOTATIONS
###############################################################################

n_clusters <- nlevels(plasmid_clusters)

cluster_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(n_clusters),
  levels(plasmid_clusters)
)

row_ha <- rowAnnotation(
  Cluster = anno_block(
    gp = gpar(fill = cluster_colors[levels(plasmid_clusters)]),
    labels = levels(plasmid_clusters),
    labels_gp = gpar(col = "black", fontsize = 10)
  )
)

temp_cols <- setNames(
  brewer.pal(max(3, nlevels(alex_meta$Temperature)), "Set1")[1:nlevels(alex_meta$Temperature)],
  levels(alex_meta$Temperature)
)

region_cols <- setNames(
  brewer.pal(max(3, nlevels(alex_meta$region)), "Set2")[1:nlevels(alex_meta$region)],
  levels(alex_meta$region)
)

season_cols <- setNames(
  brewer.pal(max(3, nlevels(alex_meta$season)), "Dark2")[1:nlevels(alex_meta$season)],
  levels(alex_meta$season)
)

top_ha <- HeatmapAnnotation(
  Temperature = alex_meta$Temperature,
  Region      = alex_meta$region,
  Season      = alex_meta$season,
  col = list(
    Temperature = temp_cols,
    Region      = region_cols,
    Season      = season_cols
  ),
  annotation_name_side = "left"
)

cov_col <- colorRamp2(
  c(0, 20, 50, 70, 100),
  c("white", "#deebf7", "#9ecae1", "#3182bd", "#08519c")
)

###############################################################################
# 7. DRAW FINAL HEATMAP (NO CLUSTERING)
###############################################################################

png(
  "Alexander_environmental_column_order_w_numbering_season.png",
  width = 4200,
  height = 5200,
  res = 300
)

Heatmap(
  coverage_mat,
  name = "Coverage (%)",
  col = cov_col,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  row_split = plasmid_clusters,
  row_title = NULL,   # ← THIS removes duplicated numbers
  
  show_row_names = FALSE,
  show_column_names = TRUE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  
  left_annotation = row_ha,
  top_annotation = top_ha,
  
  use_raster = TRUE,
  raster_quality = 2,
  column_names_rot = 45
)

dev.off()