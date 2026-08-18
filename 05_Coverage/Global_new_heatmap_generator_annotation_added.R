###############################################################################
# GLOBAL PLASMIDOME — COVERAGE-BASED PLASMID HEATMAP
# WITH STATION + SEASON ANNOTATIONS
###############################################################################

library(dplyr)
library(readr)
library(ComplexHeatmap)
library(circlize)
library(dynamicTreeCut)
library(RColorBrewer)

###############################################################################
# 1. SET WORKING DIRECTORY
###############################################################################

setwd("C:/Users/user/Desktop/Research Project Oded/metadata/Global")

###############################################################################
# 2. LOAD COVERAGE MATRIX
###############################################################################

coverage <- read.delim(
  "coverage_non_singletons.tsv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rownames(coverage) <- coverage[, 1]
coverage <- coverage[, -1]

###############################################################################
# 3. LOAD & MERGE METADATA (ALEXANDER + BOGOTA)
###############################################################################

alex_meta <- read.csv(
  "C:/Users/user/Desktop/Research Project Oded/metadata/Alexander/New/alex_clean_metadata.csv",
  stringsAsFactors = FALSE
)

bogo_meta <- read.csv(
  "C:/Users/user/Desktop/Research Project Oded/metadata/Bogota/bogo_clean_metadata.csv",
  stringsAsFactors = FALSE
)

# Keep only required columns
alex_meta <- alex_meta %>%
  select(sample_id, station_name, season)

bogo_meta <- bogo_meta %>%
  select(sample_id, station_name, season)

global_meta <- bind_rows(alex_meta, bogo_meta)

###############################################################################
# 4. GLOBAL DATASET — NO FILTERING
###############################################################################

coverage_filt <- coverage

###############################################################################
# 5. ALIGN METADATA TO COVERAGE MATRIX
###############################################################################

global_meta <- global_meta %>%
  filter(sample_id %in% colnames(coverage_filt)) %>%
  arrange(match(sample_id, colnames(coverage_filt)))

###############################################################################
# 6. COLUMN CLUSTERING (SAMPLES) — PEARSON / UPGMA
###############################################################################

sample_dist <- as.dist(
  1 - cor(
    coverage_filt,
    method = "pearson",
    use = "pairwise.complete.obs"
  )
)

sample_hc <- hclust(sample_dist, method = "average")

###############################################################################
# 7. ROW CLUSTERING (PLASMIDS) — PEARSON / UPGMA
###############################################################################

plasmid_dist <- as.dist(
  1 - cor(
    t(coverage_filt),
    method = "pearson",
    use = "pairwise.complete.obs"
  )
)

plasmid_hc <- hclust(plasmid_dist, method = "average")

###############################################################################
# 8. AUTOMATIC PLASMID CLUSTERING
###############################################################################

plasmid_clusters <- cutreeDynamic(
  dendro = plasmid_hc,
  distM = as.matrix(plasmid_dist),
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = 50
)

plasmid_clusters <- factor(plasmid_clusters)

###############################################################################
# 9. EXPORT PLASMID → CLUSTER TABLE
###############################################################################

write.csv(
  data.frame(
    Plasmid = rownames(coverage_filt),
    Cluster = plasmid_clusters
  ),
  "global_plasmid_clusters.csv",
  row.names = FALSE
)

###############################################################################
# 10. ROW ANNOTATION — CLUSTER BLOCKS
###############################################################################

n_clusters <- nlevels(plasmid_clusters)

cluster_cols <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(n_clusters),
  levels(plasmid_clusters)
)

row_ha <- rowAnnotation(
  Cluster = plasmid_clusters,
  col = list(Cluster = cluster_cols),
  show_annotation_name = TRUE
)

###############################################################################
# 11. DEFINE COLOR MAPPINGS (STATION + SEASON)
###############################################################################

station_levels <- unique(global_meta$station_name)
season_levels  <- unique(global_meta$season)

station_cols <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(station_levels)),
  station_levels
)

season_cols <- setNames(
  colorRampPalette(brewer.pal(8, "Dark2"))(length(season_levels)),
  season_levels
)

###############################################################################
# 12. TOP ANNOTATION — STATION + SEASON
###############################################################################

top_ha <- HeatmapAnnotation(
  Station = global_meta$station_name,
  Season  = global_meta$season,
  col = list(
    Station = station_cols,
    Season  = season_cols
  ),
  annotation_name_side = "left"
)

###############################################################################
# 13. COVERAGE COLOR SCALE (NO NORMALIZATION)
###############################################################################

cov_col <- colorRamp2(
  c(0, 20, 50, 70, 100),
  c("white", "#deebf7", "#9ecae1", "#3182bd", "#08519c")
)

###############################################################################
# 14. DRAW HEATMAP (PNG, THESIS SCALE)
###############################################################################

png(
  "Global_coverage_heatmap_UPGMA_station_season.png",
  width = 4200,
  height = 5200,
  res = 300
)

Heatmap(
  as.matrix(coverage_filt),
  name = "Coverage (%)",
  col = cov_col,
  cluster_columns = sample_hc,
  show_row_names = FALSE,
  show_column_names = TRUE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_split = plasmid_clusters,
  left_annotation = row_ha,
  top_annotation = top_ha,
  use_raster = TRUE,
  raster_quality = 2,
  column_names_rot = 45
)

dev.off()
