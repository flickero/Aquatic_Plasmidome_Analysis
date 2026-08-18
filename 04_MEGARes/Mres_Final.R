###############################################################################
# SAMPLE × RESISTANCE TARGET (PLASMID-NORMALIZED) HEATMAP
# Bray–Curtis ecological clustering
###############################################################################

library(dplyr)
library(readr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(vegan)

###############################################################################
# 1. WORKING DIRECTORY
###############################################################################

setwd("C:/Users/user/Desktop/Research Project Oded/MEGAres")

###############################################################################
# 2. LOAD SAMPLE × TARGET RICHNESS MATRIX
###############################################################################

amr_mat <- read.delim(
  "sample_AMR_target_richness.tsv",
  check.names = FALSE
)

rownames(amr_mat) <- amr_mat[, 1]
amr_mat <- amr_mat[, -1]

# Transpose → rows = resistance targets, columns = samples
amr_mat <- t(as.matrix(amr_mat))

###############################################################################
# 3. DERIVE PLASMID COUNTS PER SAMPLE (≥70% rule)
###############################################################################

merged <- read.delim("merged.tsv", check.names = FALSE)

rownames(merged) <- merged$Contig
merged <- merged[, -1]

plasmid_presence <- merged >= 70
plasmid_counts <- colSums(plasmid_presence)

# Align order with heatmap
plasmid_counts <- plasmid_counts[colnames(amr_mat)]

###############################################################################
# 3A. LOAD SAMPLE METADATA
###############################################################################

alex_meta <- read.csv(
  "alex_clean_metadata.csv",
  stringsAsFactors = FALSE
)

bogo_meta <- read.csv(
  "bogo_clean_metadata.csv",
  stringsAsFactors = FALSE
)

###############################################################################
# Build publication sample labels
###############################################################################

alex_labels <- alex_meta %>%
  transmute(
    sample_id,
    sample_label = paste0(
      "S",
      station_numeric,
      "-",
      recode(
        season,
        Winter = "W",
        Spring = "Sp",
        LateSummer = "LS"
      )
    ),
    sample_order = paste0(
      "A_",
      sprintf("%02d", station_numeric),
      "_",
      match(
        season,
        c("Winter", "Spring", "LateSummer")
      )
    )
  )

bogo_labels <- bogo_meta %>%
  transmute(
    sample_id,
    sample_label = station_name,
    sample_order = paste0(
      "B_",
      match(
        region,
        c(
          "River Source",
          "River Middle",
          "River Lower",
          "Hospital A",
          "Hospital B",
          "Hospital C"
        )
      ),
      "_",
      sampling_year
    )
  )

sample_metadata <- bind_rows(
  alex_labels,
  bogo_labels
)

###############################################################################
# 4. NORMALIZE RICHNESS BY PLASMID COUNT
###############################################################################

amr_mat_norm <- sweep(amr_mat, 2, plasmid_counts, FUN = "/")

amr_mat_norm[is.na(amr_mat_norm)] <- 0

###############################################################################
# 4A. REMOVE VERY RARE RESISTANCE TARGETS
###############################################################################

# Minimum total richness across all samples
MIN_TOTAL_RICHNESS <- 10

# Calculate total richness per resistance target
target_totals <- rowSums(amr_mat)

# Store omitted targets for reporting
removed_targets <- names(target_totals[target_totals < MIN_TOTAL_RICHNESS])

# Filter matrices
amr_mat <- amr_mat[
  target_totals >= MIN_TOTAL_RICHNESS,
]

amr_mat_norm <- amr_mat_norm[
  target_totals >= MIN_TOTAL_RICHNESS,
]

cat(
  length(removed_targets),
  "rare resistance targets removed from heatmap.\n"
)

print(removed_targets)

###############################################################################
# 4A. REPLACE ACCESSIONS WITH PUBLICATION SAMPLE LABELS
###############################################################################

# Keep only metadata for samples actually present
sample_metadata <- sample_metadata %>%
  filter(sample_id %in% colnames(amr_mat_norm)) %>%
  arrange(match(sample_id, colnames(amr_mat_norm)))

# Safety check
if(any(is.na(sample_metadata$sample_label))){
  stop("Missing publication labels for one or more samples.")
}

# Rename matrix columns
colnames(amr_mat_norm) <- sample_metadata$sample_label

# Rename plasmid count vector
names(plasmid_counts) <- sample_metadata$sample_label

###############################################################################
# 5. BUILD COLUMN ANNOTATION (ORFs + Reads)
###############################################################################

# ORF counts
orf_counts <- read.delim(
  "ORF_count_per_sample_filtered.tsv",
  stringsAsFactors = FALSE
)

orf_counts <- orf_counts %>%
  left_join(
    sample_metadata,
    by = c("Sample" = "sample_id")
  ) %>%
  filter(!is.na(sample_label)) %>%
  arrange(match(sample_label, colnames(amr_mat_norm)))

# Read counts
read_counts <- read.csv(
  "read_per_sample.csv",
  stringsAsFactors = FALSE
)

read_counts <- read_counts %>%
  left_join(
    sample_metadata,
    by = c("run_accession" = "sample_id")
  ) %>%
  filter(!is.na(sample_label)) %>%
  arrange(match(sample_label, colnames(amr_mat_norm)))

col_annot <- HeatmapAnnotation(
  ORFs = anno_barplot(
    orf_counts$ORF_count,
    gp = gpar(fill = "grey40"),
    border = FALSE,
    height = unit(2, "cm")
  ),
  Reads = anno_barplot(
    read_counts$read_count,
    gp = gpar(fill = "grey70"),
    border = FALSE,
    height = unit(2, "cm")
  ),
  annotation_name_side = "left"
)

###############################################################################
# 6. DEFINE RESISTANCE DOMAIN (ROW ANNOTATION)
###############################################################################

resistance_domain <- ifelse(
  grepl("metal|Zinc|Copper|Mercury|Cadmium|Nickel|Lead|Arsenic|Gold|Iron|Tellurium|Tungsten",
        rownames(amr_mat_norm), ignore.case = TRUE),
  "Metals",
  ifelse(
    grepl("biocide|QAC|Peroxide|Phenolic|Biguanide|Cationic",
          rownames(amr_mat_norm), ignore.case = TRUE),
    "Biocides",
    ifelse(
      grepl("multi", rownames(amr_mat_norm), ignore.case = TRUE),
      "Multi-compound",
      "Drugs"
    )
  )
)

resistance_domain <- factor(
  resistance_domain,
  levels = c("Drugs", "Metals", "Biocides", "Multi-compound")
)

domain_colors <- c(
  "Drugs" = "#377eb8",
  "Metals" = "#4daf4a",
  "Biocides" = "#984ea3",
  "Multi-compound" = "#e41a1c"
)

row_annot <- rowAnnotation(
  Domain = resistance_domain,
  col = list(Domain = domain_colors)
)

###############################################################################
# 7. COLOR FUNCTION
###############################################################################

max_val <- max(amr_mat_norm)

col_fun <- colorRamp2(
  c(0, max_val/5, max_val/2, max_val),
  c("#f7fbff", "#9ecae1", "#3182bd", "#08306b")
)

###############################################################################
# 8. BRAY–CURTIS CLUSTERING
###############################################################################

column_dist <- vegdist(t(amr_mat_norm), method = "bray")
row_dist <- vegdist(amr_mat_norm, method = "bray")

###############################################################################
# 9. BUILD HEATMAP
###############################################################################

ht <- Heatmap(
  amr_mat_norm,
  name = "Target proportion",
  col = col_fun,
  
  cluster_rows = as.dendrogram(hclust(row_dist, method = "average")),
  cluster_columns = as.dendrogram(hclust(column_dist, method = "average")),
  
  left_annotation = row_annot,
  top_annotation = col_annot,
  
  show_row_names = TRUE,
  show_column_names = TRUE,
  
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 60,
  
  row_title = "Resistance targets",
  column_title = "Samples",
  
  heatmap_legend_param = list(
    title = "Plasmid-normalized\nresistance richness",
    legend_height = unit(4, "cm")
  )
)

###############################################################################
# 10. EXPORT PNG
###############################################################################

png(
  "Mres_Final.png",
  width = 2400,
  height = 1600,
  res = 200
)

draw(ht,
     heatmap_legend_side = "right",
     annotation_legend_side = "right")

dev.off()