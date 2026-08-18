###############################################################################
#
# BOGOTA RIVER PLASMIDOME
#
# COMMUNITY ECOLOGY ANALYSIS
#
# Purpose
#
# Evaluate whether plasmid community composition differs among
# ecological compartments of the Bogotá River system and determine
# whether measured environmental variables are associated with
# these community patterns.
#
# Biological definition
#
# Plasmid present = coverage ≥70%
#
# Community metric
#
# Jaccard distance (presence / absence)
#
# Statistical framework
#
# • PERMANOVA
# • Dispersion diagnostics
# • Principal Coordinates Analysis (PCoA)
# • Environmental fitting (envfit)
#
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

setwd(
  "C:/Users/user/Desktop/Research Project Oded/metadata/Bogota"
)

###############################################################################
# 3. LOAD COVERAGE MATRIX
###############################################################################

coverage <- read.delim(
  
  "coverage_non_singletons.tsv",
  
  check.names = FALSE,
  
  stringsAsFactors = FALSE
  
)

rownames(coverage) <- coverage[,1]

coverage <- coverage[,-1]

###############################################################################
# 4. LOAD METADATA
###############################################################################

bogo_meta <- read.csv(
  
  "bogo_clean_metadata.csv",
  
  stringsAsFactors = FALSE
  
)

###############################################################################
# 5. SUBSET TO BOGOTA SAMPLES
###############################################################################

coverage_bogo <-
  
  coverage[,
           
           colnames(coverage) %in%
             
             bogo_meta$sample_id
           
  ]

###############################################################################
# 6. IDENTICAL PLASMID FILTERING USED FOR HEATMAP
###############################################################################

is_ERR <- grepl(
  
  "^ERR",
  
  rownames(coverage_bogo)
  
)

is_SRR <- grepl(
  
  "^SRR",
  
  rownames(coverage_bogo)
  
)

SRR_pass <- rownames(coverage_bogo)[
  
  is_SRR &
    
    apply(
      
      coverage_bogo >= 70,
      
      1,
      
      sum
      
    ) >= 2
  
]

keep_plasmids <- rownames(coverage_bogo)[
  
  is_ERR |
    
    rownames(coverage_bogo) %in%
    
    SRR_pass
  
]

coverage_filt <-
  
  coverage_bogo[
    
    keep_plasmids,
    
  ]

###############################################################################
# 7. ALIGN METADATA
###############################################################################

bogo_meta <-
  
  bogo_meta %>%
  
  filter(
    
    sample_id %in%
      
      colnames(coverage_filt)
    
  ) %>%
  
  arrange(
    
    match(
      
      sample_id,
      
      colnames(coverage_filt)
      
    )
    
  )

coverage_filt <-
  
  coverage_filt[,
                
                bogo_meta$sample_id
                
  ]

###############################################################################
# 8. DEFINE ECOLOGICAL GROUPS
###############################################################################

bogo_meta$ecological_group <-
  
  case_when(
    
    bogo_meta$region %in%
      
      c(
        
        "River Source",
        
        "River Middle"
        
      )
    
    ~
      
      "Upstream",
    
    bogo_meta$region ==
      
      "River Lower"
    
    ~
      
      "Lower River",
    
    TRUE
    
    ~
      
      "Hospital"
    
  )

bogo_meta$ecological_group <-
  
  factor(
    
    bogo_meta$ecological_group,
    
    levels = c(
      
      "Upstream",
      
      "Lower River",
      
      "Hospital"
      
    )
    
  )

###############################################################################
# 9. CREATE PRESENCE / ABSENCE MATRIX
###############################################################################

coverage_pa <-
  
  coverage_filt >= 70

###############################################################################
# 10. BUILD JACCARD DISTANCE MATRIX
###############################################################################

dist_mat <- vegdist(
  
  t(coverage_pa),
  
  method = "jaccard",
  
  binary = TRUE
  
)

###############################################################################
# 11. INITIAL SUMMARY
###############################################################################

cat("\n")

cat("========================================\n")
cat("BOGOTA COMMUNITY ANALYSIS\n")
cat("========================================\n\n")

cat(
  
  "Samples:",
  
  ncol(coverage_pa),
  
  "\n"
  
)

cat(
  
  "Plasmids:",
  
  nrow(coverage_pa),
  
  "\n"
  
)

cat(
  
  "Ecological groups:\n"
  
)

print(
  
  table(
    
    bogo_meta$ecological_group
    
  )
  
)

cat("\n")

cat(
  
  "Distance metric: Jaccard (presence/absence)\n"
  
)

cat(
  
  "Coverage threshold: 70%\n"
  
)

cat("\n")

###############################################################################
# END OF CHUNK 1
###############################################################################
###############################################################################
# 12. PERMANOVA
###############################################################################

set.seed(42)

permanova_ecology <- adonis2(
  
  dist_mat ~ ecological_group,
  
  data = bogo_meta,
  
  permutations = 9999
  
)

###############################################################################
# 13. DISPERSION DIAGNOSTICS
###############################################################################

bd_ecology <- betadisper(
  
  dist_mat,
  
  bogo_meta$ecological_group
  
)

anova_ecology <- anova(
  
  bd_ecology
  
)

perm_ecology <- permutest(
  
  bd_ecology,
  
  permutations = 9999
  
)

###############################################################################
# 14. PRINCIPAL COORDINATES ANALYSIS
###############################################################################

pcoa <- cmdscale(
  
  dist_mat,
  
  eig = TRUE,
  
  k = 2
  
)

###############################################################################
# 15. CREATE PCOA DATA FRAME
###############################################################################

pcoa_df <- data.frame(
  
  sample_id = rownames(pcoa$points),
  
  PCoA1 = pcoa$points[,1],
  
  PCoA2 = pcoa$points[,2]
  
)

pcoa_df <- left_join(
  
  pcoa_df,
  
  bogo_meta,
  
  by = "sample_id"
  
)

###############################################################################
# 16. VARIANCE EXPLAINED
###############################################################################

positive_axes <- pcoa$eig[
  pcoa$eig > 0
]

variance_explained <- round(
  
  100 *
    
    positive_axes /
    
    sum(positive_axes),
  
  2
  
)

###############################################################################
# 17. EXPORT PCOA COORDINATES
###############################################################################

write.csv(
  
  pcoa_df,
  
  "Bogota_PCoA_coordinates.csv",
  
  row.names = FALSE
  
)

###############################################################################
# 18. CONSOLE SUMMARY
###############################################################################

cat("\n")
cat("=========================================\n")
cat("COMMUNITY ANALYSIS\n")
cat("=========================================\n\n")

cat("PERMANOVA\n\n")

print(permanova_ecology)

cat("\n")

cat("Observed R² =",
    
    round(
      
      permanova_ecology$R2[1],
      
      3
      
    ),
    
    "\n\n"
    
)

cat("DISPERSION ANOVA\n\n")

print(anova_ecology)

cat("\n")

cat("PERMUTATION TEST\n\n")

print(perm_ecology)

cat("\n")

cat("Variance explained\n")

cat(
  
  "PCoA1:",
  
  variance_explained[1],
  
  "%\n"
  
)

cat(
  
  "PCoA2:",
  
  variance_explained[2],
  
  "%\n"
  
)

cat("\n")

###############################################################################
# END OF CHUNK 2
###############################################################################
###############################################################################
# 19. SELECT SAMPLES WITH ENVIRONMENTAL DATA
###############################################################################

env_meta <-
  
  bogo_meta %>%
  
  filter(
    
    include_env == TRUE
    
  )

###############################################################################
# Keep the corresponding PCoA coordinates
###############################################################################

env_pcoa <-
  
  pcoa_df %>%
  
  filter(
    
    sample_id %in%
      
      env_meta$sample_id
    
  ) %>%
  
  arrange(
    
    match(
      
      sample_id,
      
      env_meta$sample_id
      
    )
    
  )

###############################################################################
# Verify ordering
###############################################################################

stopifnot(
  
  identical(
    
    env_pcoa$sample_id,
    
    env_meta$sample_id
    
  )
  
)

###############################################################################
# 20. CREATE ENVIRONMENTAL MATRIX
###############################################################################

abx_cols <- c(
  
  "MNZ",
  
  "SMX",
  
  "TMP",
  
  "NOR",
  
  "CIP",
  
  "CLI",
  
  "CLR",
  
  "ERY",
  
  "AZM"
  
)

env_meta$API <-
  
  rowMeans(
    
    env_meta[,abx_cols],
    
    na.rm = TRUE
    
  )

env_data <-
  
  env_meta %>%
  
  dplyr::select(
    
    BOD,
    
    COD,
    
    TSS,
    
    TN,
    
    API
    
  )

###############################################################################
# 21. ENVIRONMENTAL SUMMARY
###############################################################################

env_summary <- data.frame(
  
  Variable = names(env_data),
  
  Mean = sapply(env_data, mean),
  
  SD = sapply(env_data, sd),
  
  Minimum = sapply(env_data, min),
  
  Maximum = sapply(env_data, max)
  
)

write.csv(
  
  env_summary,
  
  "Bogota_environment_summary.csv",
  
  row.names = FALSE
  
)

###############################################################################
# 22. ENVIRONMENTAL CORRELATION
###############################################################################

env_cor <- cor(
  
  env_data,
  
  method = "pearson"
  
)

write.csv(
  
  round(env_cor,3),
  
  "Bogota_environment_correlation.csv"
  
)

###############################################################################
# 23. ENVFIT
###############################################################################

set.seed(42)

env_fit <- envfit(
  
  env_pcoa[,c("PCoA1","PCoA2")],
  
  env_data,
  
  permutations = 9999
  
)
###############################################################################
# 23B. PERMUTATION RESOLUTION
###############################################################################

n_env_samples <- nrow(env_data)

total_permutations <- factorial(n_env_samples)

minimum_p <- 1 / total_permutations

cat("\n")
cat("=========================================\n")
cat("PERMUTATION INFORMATION\n")
cat("=========================================\n\n")

cat("Environmental samples:", n_env_samples, "\n")
cat("Unique permutations:", total_permutations, "\n")
cat("Minimum attainable p-value:",
    signif(minimum_p, 3), "\n\n")
###############################################################################
# 24. EXTRACT RESULTS
###############################################################################

env_vectors <- data.frame(
  
  Variable = rownames(
    
    env_fit$vectors$arrows
    
  ),
  
  Axis1 = env_fit$vectors$arrows[,1],
  
  Axis2 = env_fit$vectors$arrows[,2],
  
  R2 = env_fit$vectors$r,
  
  Pvalue = env_fit$vectors$pvals
  
)

env_vectors <-
  
  env_vectors %>%
  
  arrange(
    
    Pvalue,
    
    desc(R2)
    
  )

write.csv(
  
  env_vectors,
  
  "Bogota_envfit_vectors.csv",
  
  row.names = FALSE
  
)

###############################################################################
# Significant variables
###############################################################################

env_vectors_sig <-
  
  env_vectors %>%
  
  filter(
    
    Pvalue <= 0.05
    
  )

###############################################################################
# 25. CONSOLE SUMMARY
###############################################################################

cat("\n")
cat("=========================================\n")
cat("ENVIRONMENTAL FITTING\n")
cat("=========================================\n\n")

cat(
  
  "Samples with environmental data:",
  
  nrow(env_meta),
  
  "\n\n"
  
)

cat(
  
  "Environmental variables:\n\n"
  
)

print(
  
  env_summary
  
)

cat("\n")

cat(
  
  "Correlation matrix\n\n"
  
)

print(
  
  round(env_cor,3)
  
)

cat("\n")

cat(
  
  "ENVFIT RESULTS\n\n"
  
)

print(
  
  env_vectors
  
)

cat("\n")

cat(
  
  "Significant variables:",
  
  nrow(env_vectors_sig),
  
  "\n"
  
)

###############################################################################
# END OF CHUNK 3
###############################################################################
###############################################################################
# 26. PUBLICATION QUALITY PCoA
###############################################################################

library(ggrepel)

vector_scale <- 0.45

plot_vectors <- env_vectors

plot_vectors$Axis1 <-
  plot_vectors$Axis1 * vector_scale

plot_vectors$Axis2 <-
  plot_vectors$Axis2 * vector_scale

p <- ggplot(
  
  pcoa_df,
  
  aes(
    
    PCoA1,
    
    PCoA2,
    
    colour = ecological_group
    
  )
  
) +
  
  geom_point(
    
    size = 4
    
  ) +
  
  geom_text_repel(
    
    aes(label = station_name),
    
    size = 3.8,
    
    max.overlaps = Inf
    
  ) +
  
  geom_segment(
    
    data = plot_vectors,
    
    aes(
      
      x = 0,
      
      y = 0,
      
      xend = Axis1,
      
      yend = Axis2
      
    ),
    
    inherit.aes = FALSE,
    
    arrow = arrow(length = unit(0.25,"cm")),
    
    linewidth = 0.7,
    
    colour = "black"
    
  ) +
  
  geom_text_repel(
    
    data = plot_vectors,
    
    aes(
      
      Axis1,
      
      Axis2,
      
      label = Variable
      
    ),
    
    inherit.aes = FALSE,
    
    colour = "black",
    
    size = 4
    
  ) +
  
  theme_bw(base_size = 13) +
  
  labs(
    
    title = "Bogotá plasmid community structure",
    
    subtitle = "Jaccard distance (presence / absence ≥70%)",
    
    x = paste0(
      
      "PCoA1 (",
      
      variance_explained[1],
      
      "%)"
      
    ),
    
    y = paste0(
      
      "PCoA2 (",
      
      variance_explained[2],
      
      "%)"
      
    ),
    
    colour = "Ecological group"
    
  )

ggsave(
  
  "Bogota_PCoA_envfit.png",
  
  p,
  
  width = 8,
  
  height = 6,
  
  dpi = 300
  
)

###############################################################################
# 27. ENVIRONMENTAL PCA EXPORTS
###############################################################################

env_loading <- data.frame(
  
  Variable = rownames(
    
    env_pca$rotation
    
  ),
  
  PC1 = env_pca$rotation[,1],
  
  PC2 = env_pca$rotation[,2]
  
)

write.csv(
  
  env_loading,
  
  "Bogota_environmental_PCA_loadings.csv",
  
  row.names = FALSE
  
)

write.csv(
  
  env_scores,
  
  "Bogota_environmental_PC1_scores.csv",
  
  row.names = FALSE
  
)

###############################################################################
# 28. GRADIENT CORRELATION
###############################################################################

gradient_cor <- cor.test(
  
  comparison$Pollution_PC1,
  
  comparison$PCoA1,
  
  method = "pearson"
  
)

capture.output(
  
  gradient_cor,
  
  file = "Bogota_gradient_correlation.txt"
  
)

###############################################################################
# 29. CENTROIDS
###############################################################################

group_centroids <- aggregate(
  
  cbind(
    
    PCoA1,
    
    PCoA2
    
  )
  
  ~
    
    ecological_group,
  
  data = pcoa_df,
  
  mean
  
)

write.csv(
  
  group_centroids,
  
  "Bogota_group_centroids.csv",
  
  row.names = FALSE
  
)

###############################################################################
# 30. EXPORT COMMUNITY SUMMARY
###############################################################################

community_summary <- data.frame(
  
  Metric = c(
    
    "PERMANOVA_R2",
    
    "PERMANOVA_P",
    
    "Dispersion_P",
    
    "PCoA1",
    
    "PCoA2"
    
  ),
  
  Value = c(
    
    permanova_ecology$R2[1],
    
    permanova_ecology$`Pr(>F)`[1],
    
    perm_ecology$tab[1,"Pr(>F)"],
    
    variance_explained[1],
    
    variance_explained[2]
    
  )
  
)

write.csv(
  
  community_summary,
  
  "Bogota_community_summary.csv",
  
  row.names = FALSE
  
)

###############################################################################
# 31. FINAL CONSOLE OUTPUT
###############################################################################

cat("\n")
cat("=================================================\n")
cat("BOGOTA ANALYSIS COMPLETE\n")
cat("=================================================\n\n")

cat("Generated outputs\n\n")

cat("Community\n")
cat("---------\n")
cat("Bogota_PCoA_coordinates.csv\n")
cat("Bogota_group_centroids.csv\n")
cat("Bogota_community_summary.csv\n\n")

cat("Environment\n")
cat("-----------\n")
cat("Bogota_environment_summary.csv\n")
cat("Bogota_environment_correlation.csv\n")
cat("Bogota_envfit_vectors.csv\n")
cat("Bogota_environmental_PCA_loadings.csv\n")
cat("Bogota_environmental_PC1_scores.csv\n")
cat("Bogota_gradient_correlation.txt\n\n")

cat("Figure\n")
cat("------\n")
cat("Bogota_PCoA_envfit.png\n\n")

cat("=================================================\n")