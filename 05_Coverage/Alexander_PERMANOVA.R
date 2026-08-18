###############################################################################

# ALEXANDER RIVER PLASMIDOME

# PERMANOVA + DISPERSION DIAGNOSTICS + PCoA

#

# Purpose:

# Evaluate whether plasmid composition varies:

# 1. Along the river continuum (Region)

# 2. Between seasons (Season)

# 3. Under a combined Region + Season model

#

# Biological definition:

# Plasmid present = coverage >= 70%

#

# Distance metric:

# Jaccard (binary)

#

# Statistical framework:

# - PERMANOVA (adonis2)

# - betadisper

# - permutation-based dispersion test

# - PCoA visualization

#

# Outputs:

# Alexander_PERMANOVA_results.txt

# Alexander_PCoA_coordinates.csv

# Alexander_PCoA_Jaccard.png

###############################################################################

###############################################################################

# 1. LIBRARIES

###############################################################################

library(dplyr)
library(readr)
library(vegan)
library(ggplot2)

###############################################################################

# 2. WORKING DIRECTORY

###############################################################################

setwd("C:/Users/user/Desktop/Research Project Oded/metadata/Alexander/New")

###############################################################################

# 3. INPUT FILES

###############################################################################

coverage_mat <- read.delim(
  "alexander_coverage_filt_for_clustering.tsv",
  check.names = FALSE,
  row.names = 1
)

coverage_mat <- as.matrix(coverage_mat)

alex_meta <- read.csv(
  "alex_clean_metadata.csv",
  stringsAsFactors = FALSE
)

###############################################################################

# 4. ALIGN METADATA TO MATRIX

###############################################################################

alex_meta <- alex_meta %>%
  filter(sample_id %in% colnames(coverage_mat)) %>%
  arrange(match(sample_id, colnames(coverage_mat)))

coverage_mat <- coverage_mat[, alex_meta$sample_id]

###############################################################################

# 5. CONVERT TO PRESENCE / ABSENCE

###############################################################################

coverage_pa <- coverage_mat >= 70

###############################################################################

# 6. FACTOR LEVELS

###############################################################################

alex_meta$region <- factor(
  alex_meta$region,
  levels = c(
    "Upstream",
    "Midstream",
    "Urban_Tributary",
    "Lower_River",
    "Estuary"
  )
)

alex_meta$season <- factor(
  alex_meta$season,
  levels = c(
    "Winter",
    "Spring",
    "LateSummer"
  )
)

###############################################################################

# 7. JACCARD DISTANCE MATRIX

###############################################################################

dist_mat <- vegdist(
  t(coverage_pa),
  method = "jaccard",
  binary = TRUE
)

###############################################################################
# 8. PERMANOVA
###############################################################################

set.seed(42)

#############################
# REGION ONLY
#############################

permanova_region <- adonis2(
  dist_mat ~ region,
  data = alex_meta,
  permutations = 9999
)

#############################
# SEASON ONLY
#############################

permanova_season <- adonis2(
  dist_mat ~ season,
  data = alex_meta,
  permutations = 9999
)

#############################
# REGION + SEASON
#############################

permanova_combined <- adonis2(
  dist_mat ~ region + season,
  data = alex_meta,
  permutations = 9999
)

#############################
# MARGINAL EFFECTS
#############################

permanova_marginal <- adonis2(
  dist_mat ~ region + season,
  data = alex_meta,
  permutations = 9999,
  by = "margin"
)

###############################################################################

# 9. DISPERSION DIAGNOSTICS

###############################################################################

#############################

# REGION

#############################

bd_region <- betadisper(
  dist_mat,
  alex_meta$region
)

anova_region <- anova(bd_region)

perm_region <- permutest(
  bd_region,
  permutations = 9999
)

#############################

# SEASON

#############################

bd_season <- betadisper(
  dist_mat,
  alex_meta$season
)

anova_season <- anova(bd_season)

perm_season <- permutest(
  bd_season,
  permutations = 9999
)

###############################################################################

# 10. PCoA

###############################################################################

pcoa <- cmdscale(
  dist_mat,
  eig = TRUE,
  k = 2
)

pcoa_df <- data.frame(
  sample_id = rownames(pcoa$points),
  PCoA1 = pcoa$points[, 1],
  PCoA2 = pcoa$points[, 2]
)

pcoa_df <- left_join(
  pcoa_df,
  alex_meta,
  by = "sample_id"
)

###############################################################################

# 11. EXPORT PCoA COORDINATES

###############################################################################

write.csv(
  pcoa_df,
  "Alexander_PCoA_coordinates.csv",
  row.names = FALSE
)

###############################################################################

# 12. PCoA FIGURE

###############################################################################

variance_explained <- round(
  100 * pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]),
  2
)

p <- ggplot(
  pcoa_df,
  aes(
    x = PCoA1,
    y = PCoA2,
    color = region,
    shape = season
  )
) +
  geom_point(size = 4) +
  theme_bw() +
  labs(
    title = "Alexander River Plasmid Composition",
    subtitle = "Jaccard distance (presence/absence, ≥70% coverage)",
    x = paste0(
      "PCoA1 (",
      variance_explained[1],
      "%)"
    ),
    y = paste0(
      "PCoA2 (",
      variance_explained[2],
      "%)"
    )
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  "Alexander_PCoA_Jaccard.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

###############################################################################

# 13. EXPORT REPORT

###############################################################################

sink("Alexander_PERMANOVA_results.txt")

cat(
  "========================================================\n",
  "ALEXANDER RIVER PLASMIDOME\n",
  "PERMANOVA ANALYSIS\n",
  "========================================================\n\n",
  
  "Distance metric: Jaccard\n",
  "Matrix type: Presence / Absence\n",
  "Occurrence threshold: 70% coverage\n",
  "Permutations: 9999\n\n",
  
  "========================================================\n",
  "MODEL 1: REGION ONLY\n",
  "========================================================\n\n"
)

print(permanova_region)

cat(
  "\n\n========================================================\n",
  "MODEL 2: SEASON ONLY\n",
  "========================================================\n\n"
)

print(permanova_season)

cat(
  "\n\n========================================================\n",
  "MODEL 3: REGION + SEASON\n",
  "========================================================\n\n"
)

print(permanova_combined)

cat(
  "\n\n========================================================\n",
  "MARGINAL EFFECTS MODEL\n",
  "========================================================\n\n"
)

print(permanova_marginal)

cat(
  "\n\n========================================================\n",
  "DISPERSION DIAGNOSTICS: REGION\n",
  "========================================================\n\n",
  
  "ANOVA:\n\n"
)

print(anova_region)

cat(
  "\n\nPERMUTATION TEST:\n\n"
)

print(perm_region)

cat(
  "\n\n========================================================\n",
  "DISPERSION DIAGNOSTICS: SEASON\n",
  "========================================================\n\n",
  
  "ANOVA:\n\n"
)

print(anova_season)

cat(
  "\n\nPERMUTATION TEST:\n\n"
)

print(perm_season)

sink()

###############################################################################

# 14. CONSOLE SUMMARY

###############################################################################

cat("\n")
cat("========================================\n")
cat("Analysis completed successfully\n")
cat("========================================\n")
cat("Generated files:\n")
cat(" - Alexander_PERMANOVA_results.txt\n")
cat(" - Alexander_PCoA_coordinates.csv\n")
cat(" - Alexander_PCoA_Jaccard.png\n")
cat("========================================\n")
cat("\n")
###############################################################################

# END

###############################################################################
