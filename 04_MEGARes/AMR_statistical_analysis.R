###############################################################################
#
# PLASMID-ASSOCIATED RESISTANCE PROFILE ANALYSIS
#
# Purpose:
#
# Characterize ecological variation in plasmid-associated resistance target
# composition across the Alexander River and Bogotá River datasets.
#
# Response matrix:
#
# Plasmid-normalized resistance target richness
#
# Community distance:
#
# Bray–Curtis
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
  "C:/Users/user/Desktop/Research Project Oded/MEGAres"
)

###############################################################################
# 3. LOAD SAMPLE × TARGET RICHNESS MATRIX
###############################################################################

amr_mat <- read.delim(
  "sample_AMR_target_richness.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

rownames(amr_mat) <- amr_mat[,1]

amr_mat <- amr_mat[,-1]

###############################################################################
# Rows = resistance targets
# Columns = samples
###############################################################################

amr_mat <- t(as.matrix(amr_mat))

###############################################################################
# 4. LOAD COVERAGE MATRIX
#
# Used to calculate plasmid counts per sample
# (presence defined as ≥70% coverage)
###############################################################################

merged <- read.delim(
  "merged.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

rownames(merged) <- merged[,1]

merged <- merged[,-1]

###############################################################################
# 5. DERIVE PLASMID COUNTS
###############################################################################

plasmid_presence <-
  merged >= 70

plasmid_counts <-
  colSums(plasmid_presence)

###############################################################################
# Match sample order
###############################################################################

plasmid_counts <-
  plasmid_counts[
    colnames(amr_mat)
  ]

###############################################################################
# 6. NORMALIZE AMR RICHNESS
###############################################################################

amr_mat_norm <-
  sweep(
    amr_mat,
    2,
    plasmid_counts,
    FUN = "/"
  )

amr_mat_norm[
  is.na(amr_mat_norm)
] <- 0

###############################################################################
# 7. LOAD METADATA
###############################################################################

alex_meta <- read.csv(
  "C:/Users/user/Desktop/Research Project Oded/metadata/Alexander/New/alex_clean_metadata.csv",
  stringsAsFactors = FALSE
)

bogo_meta <- read.csv(
  "C:/Users/user/Desktop/Research Project Oded/metadata/Bogota/bogo_clean_metadata.csv",
  stringsAsFactors = FALSE
)

###############################################################################
# 8. DERIVE UNIFIED METADATA
###############################################################################

#############################
# Alexander
#############################

alex_meta <-
  alex_meta %>%
  mutate(
    
    Study = "Alexander",
    
    Sample_type = "River",
    
    Ecological_group = region
    
  )

#############################
# Bogotá
#############################

bogo_meta <-
  bogo_meta %>%
  mutate(
    
    Study = "Bogota",
    
    Sample_type =
      ifelse(
        grepl("^Hospital", region),
        "Hospital",
        "River"
      ),
    
    Ecological_group =
      case_when(
        
        region %in%
          c(
            "River Source",
            "River Middle"
          ) ~ "Upstream",
        
        region ==
          "River Lower" ~
          "Lower_River",
        
        grepl(
          "Hospital",
          region
        ) ~
          "Hospital"
        
      )
    
  )

###############################################################################
# 9. KEEP COMMON VARIABLES
###############################################################################

common_cols <- c(
  
  "sample_id",
  
  "Study",
  
  "Sample_type",
  
  "Ecological_group",
  
  "station_name",
  
  "region",
  
  "season"
  
)

alex_meta <-
  alex_meta[
    ,
    common_cols
  ]

bogo_meta <-
  bogo_meta[
    ,
    common_cols
  ]

###############################################################################
# 10. MERGE METADATA
###############################################################################

global_meta <-
  bind_rows(
    
    alex_meta,
    
    bogo_meta
    
  )

###############################################################################
# 11. ALIGN TO AMR MATRIX
###############################################################################

global_meta <-
  global_meta %>%
  filter(
    
    sample_id %in%
      colnames(amr_mat_norm)
    
  ) %>%
  arrange(
    
    match(
      
      sample_id,
      
      colnames(amr_mat_norm)
      
    )
    
  )

amr_mat_norm <-
  amr_mat_norm[
    ,
    global_meta$sample_id
  ]

###############################################################################
# 12. DATA SUMMARY
###############################################################################

cat("\n")

cat("====================================================\n")
cat("PLASMID-ASSOCIATED RESISTANCE ANALYSIS\n")
cat("====================================================\n\n")

cat(
  "Samples:",
  ncol(amr_mat_norm),
  "\n"
)

cat(
  "Resistance targets:",
  nrow(amr_mat_norm),
  "\n\n"
)

cat("Studies\n")

print(
  table(global_meta$Study)
)

cat("\n")

cat("Sample type\n")

print(
  table(global_meta$Sample_type)
)

cat("\n")

cat("Ecological groups\n")

print(
  table(global_meta$Ecological_group)
)

cat("\n")

cat(
  "Response matrix:",
  "Plasmid-normalized resistance target richness\n"
)

cat(
  "Community distance:",
  "Bray-Curtis\n"
)

###############################################################################
# END OF CHUNK 1
###############################################################################
###############################################################################
# 13. BRAY–CURTIS DISTANCE MATRIX
###############################################################################

dist_mat <-
  vegdist(
    
    t(amr_mat_norm),
    
    method = "bray"
    
  )

###############################################################################
# 14. RECONSTRUCT SAMPLE DENDROGRAM
#
# Same clustering used in the heatmap:
#
# Bray–Curtis
# ↓
# UPGMA (average linkage)
###############################################################################

sample_hc <-
  hclust(
    
    dist_mat,
    
    method = "average"
    
  )

###############################################################################
# 15. DENDROGRAM VALIDATION
#
# Cophenetic correlation
###############################################################################

cophenetic_r <-
  cor(
    
    dist_mat,
    
    cophenetic(sample_hc)
    
  )

###############################################################################
# 16. PRINCIPAL COORDINATE ANALYSIS (PCoA)
###############################################################################

pcoa <-
  cmdscale(
    
    dist_mat,
    
    eig = TRUE,
    
    k = 2
    
  )

###############################################################################
# Variance explained
###############################################################################

variance_explained <-
  round(
    
    100 *
      pcoa$eig /
      sum(pcoa$eig[pcoa$eig > 0]),
    
    2
    
  )

###############################################################################
# 17. BUILD PCoA DATA FRAME
###############################################################################

pcoa_df <-
  data.frame(
    
    sample_id = rownames(pcoa$points),
    
    PCoA1 = pcoa$points[,1],
    
    PCoA2 = pcoa$points[,2]
    
  )

###############################################################################
# Add metadata
###############################################################################

pcoa_df <-
  left_join(
    
    pcoa_df,
    
    global_meta,
    
    by = "sample_id"
    
  )

###############################################################################
# Export coordinates
###############################################################################

write.csv(
  
  pcoa_df,
  
  "AMR_PCoA_coordinates.csv",
  
  row.names = FALSE
  
)

###############################################################################
# 18. PCoA FIGURE
#
# Validation figure
# (not necessarily intended for the thesis)
###############################################################################

p <-
  
  ggplot(
    
    pcoa_df,
    
    aes(
      
      x = PCoA1,
      
      y = PCoA2,
      
      colour = Study,
      
      shape = Sample_type
      
    )
    
  ) +
  
  geom_point(
    
    size = 4
    
  ) +
  
  theme_bw(
    
    base_size = 16
    
  ) +
  
  labs(
    
    title =
      "Global plasmid-associated resistance structure",
    
    subtitle =
      "Bray–Curtis dissimilarity represented by PCoA",
    
    x =
      paste0(
        "PCoA1 (",
        variance_explained[1],
        "%)"
      ),
    
    y =
      paste0(
        "PCoA2 (",
        variance_explained[2],
        "%)"
      )
    
  )

ggsave(
  
  "AMR_Global_PCoA.png",
  
  p,
  
  width = 9,
  
  height = 7,
  
  dpi = 300
  
)

###############################################################################
# 19. SAMPLE DENDROGRAM
###############################################################################

png(
  
  "AMR_Global_Dendrogram.png",
  
  width = 2200,
  
  height = 1200,
  
  res = 300
  
)

plot(
  
  sample_hc,
  
  main = "Global AMR Sample Dendrogram",
  
  xlab = "",
  
  sub = "",
  
  cex = 0.8
  
)

dev.off()

###############################################################################
# 20. EXPORT DENDROGRAM ORDER
###############################################################################

dendrogram_order <-
  data.frame(
    
    Position = seq_along(sample_hc$order),
    
    Sample = sample_hc$labels[sample_hc$order]
    
  )

write.csv(
  
  dendrogram_order,
  
  "AMR_Global_Dendrogram_Order.csv",
  
  row.names = FALSE
  
)

###############################################################################
# 21. CONSOLE SUMMARY
###############################################################################

cat("\n")

cat("====================================================\n")
cat("BRAY–CURTIS COMMUNITY STRUCTURE\n")
cat("====================================================\n\n")

cat(
  
  "Cophenetic correlation coefficient:\n"
  
)

print(
  
  round(cophenetic_r,3)
  
)

cat("\n")

cat(
  
  "Variance explained by first two PCoA axes:\n\n"
  
)

print(
  
  variance_explained[1:2]
  
)

cat("\n")

cat(
  
  "PCoA coordinates exported:\n",
  
  "AMR_PCoA_coordinates.csv\n"
  
)

cat("\n")

cat(
  
  "Figures created:\n",
  
  "- AMR_Global_PCoA.png\n",
  
  "- AMR_Global_Dendrogram.png\n"
  
)

###############################################################################
# END OF CHUNK 2
###############################################################################
###############################################################################
# 22. PERMANOVA
###############################################################################

set.seed(42)

###############################################################################
# MODEL 1 — STUDY
###############################################################################

permanova_study <-
  adonis2(
    
    dist_mat ~ Study,
    
    data = global_meta,
    
    permutations = 9999
    
  )

###############################################################################
# MODEL 2 — SAMPLE TYPE
###############################################################################

permanova_sampletype <-
  adonis2(
    
    dist_mat ~ Sample_type,
    
    data = global_meta,
    
    permutations = 9999
    
  )

###############################################################################
# MODEL 3 — ECOLOGICAL GROUP
###############################################################################

permanova_ecology <-
  adonis2(
    
    dist_mat ~ Ecological_group,
    
    data = global_meta,
    
    permutations = 9999
    
  )

###############################################################################
# MODEL 4 — SEASON (ALEXANDER ONLY)
###############################################################################

alex_samples <-
  global_meta %>%
  filter(
    Study == "Alexander"
  )

alex_mat <-
  amr_mat_norm[
    ,
    alex_samples$sample_id
  ]

alex_dist <-
  vegdist(
    
    t(alex_mat),
    
    method = "bray"
    
  )

permanova_season <-
  adonis2(
    
    alex_dist ~ season,
    
    data = alex_samples,
    
    permutations = 9999
    
  )

###############################################################################
# 23. HOMOGENEITY OF DISPERSION
###############################################################################

#############################
# STUDY
#############################

bd_study <-
  betadisper(
    
    dist_mat,
    
    global_meta$Study
    
  )

perm_study <-
  permutest(
    
    bd_study,
    
    permutations = 9999
    
  )

#############################
# SAMPLE TYPE
#############################

bd_sampletype <-
  betadisper(
    
    dist_mat,
    
    global_meta$Sample_type
    
  )

perm_sampletype <-
  permutest(
    
    bd_sampletype,
    
    permutations = 9999
    
  )

#############################
# ECOLOGICAL GROUP
#############################

bd_ecology <-
  betadisper(
    
    dist_mat,
    
    global_meta$Ecological_group
    
  )

perm_ecology <-
  permutest(
    
    bd_ecology,
    
    permutations = 9999
    
  )

#############################
# SEASON
#############################

bd_season <-
  betadisper(
    
    alex_dist,
    
    alex_samples$season
    
  )

perm_season <-
  permutest(
    
    bd_season,
    
    permutations = 9999
    
  )

###############################################################################
# 24. EXPORT REPORT
###############################################################################

sink(
  "AMR_PERMANOVA_results.txt"
)

cat(
  "========================================================\n",
  "PLASMID-ASSOCIATED RESISTANCE PROFILE ANALYSIS\n",
  "========================================================\n\n"
)

cat(
  "Response matrix: plasmid-normalized resistance target richness\n",
  "Distance metric: Bray–Curtis\n",
  "Clustering: UPGMA\n",
  "Ordination: PCoA\n\n"
)

###############################################################################
cat(
  "========================================================\n",
  "MODEL 1: STUDY\n",
  "========================================================\n\n"
)

print(permanova_study)

cat("\nDispersion test\n\n")

print(perm_study)

###############################################################################
cat(
  "\n========================================================\n",
  "MODEL 2: SAMPLE TYPE\n",
  "========================================================\n\n"
)

print(permanova_sampletype)

cat("\nDispersion test\n\n")

print(perm_sampletype)

###############################################################################
cat(
  "\n========================================================\n",
  "MODEL 3: ECOLOGICAL GROUP\n",
  "========================================================\n\n"
)

print(permanova_ecology)

cat("\nDispersion test\n\n")

print(perm_ecology)

###############################################################################
cat(
  "\n========================================================\n",
  "MODEL 4: ALEXANDER SEASON\n",
  "========================================================\n\n"
)

print(permanova_season)

cat("\nDispersion test\n\n")

print(perm_season)

###############################################################################
cat(
  "\n========================================================\n",
  "DENDROGRAM VALIDATION\n",
  "========================================================\n\n"
)

cat("Cophenetic correlation:\n")

print(cophenetic_r)

cat("\n")

cat("Variance explained by first two PCoA axes:\n")

print(
  variance_explained[1:2]
)

sink()

###############################################################################
# 25. CONSOLE SUMMARY
###############################################################################

cat("\n")

cat("====================================================\n")
cat("PERMANOVA ANALYSIS COMPLETE\n")
cat("====================================================\n\n")

summary_table <-
  data.frame(
    
    Model = c(
      "Study",
      "Sample type",
      "Ecological group",
      "Alexander season"
    ),
    
    R2 = c(
      permanova_study$R2[1],
      permanova_sampletype$R2[1],
      permanova_ecology$R2[1],
      permanova_season$R2[1]
    ),
    
    Pvalue = c(
      permanova_study$`Pr(>F)`[1],
      permanova_sampletype$`Pr(>F)`[1],
      permanova_ecology$`Pr(>F)`[1],
      permanova_season$`Pr(>F)`[1]
    ),
    
    Dispersion_P = c(
      perm_study$tab[1,"Pr(>F)"],
      perm_sampletype$tab[1,"Pr(>F)"],
      perm_ecology$tab[1,"Pr(>F)"],
      perm_season$tab[1,"Pr(>F)"]
    )
    
  )

print(
  summary_table,
  row.names = FALSE
)

cat("\n")

cat(
  "Detailed report written to:\n",
  "AMR_PERMANOVA_results.txt\n"
)

###############################################################################
# END OF CHUNK 3
###############################################################################
###############################################################################
# 26. PAIRWISE PERMANOVA
###############################################################################

group_levels <- levels(
  factor(global_meta$Ecological_group)
)

pairwise_results <- data.frame()

###############################################################################
# Compare every ecological-group pair
###############################################################################

for(i in 1:(length(group_levels)-1)){
  
  for(j in (i+1):length(group_levels)){
    
    g1 <- group_levels[i]
    g2 <- group_levels[j]
    
    ###########################################################################
    # Metadata subset
    ###########################################################################
    
    meta_sub <-
      global_meta %>%
      filter(
        Ecological_group %in%
          c(g1, g2)
      )
    
    ###########################################################################
    # Community matrix
    ###########################################################################
    
    comm_sub <-
      amr_mat_norm[
        ,
        meta_sub$sample_id
      ]
    
    ###########################################################################
    # Bray–Curtis distance
    ###########################################################################
    
    dist_sub <-
      vegdist(
        t(comm_sub),
        method = "bray"
      )
    
    ###########################################################################
    # PERMANOVA
    ###########################################################################
    
    model <-
      adonis2(
        
        dist_sub ~ Ecological_group,
        
        data = meta_sub,
        
        permutations = 9999
        
      )
    
    ###########################################################################
    # Save results
    ###########################################################################
    
    pairwise_results <-
      
      rbind(
        
        pairwise_results,
        
        data.frame(
          
          Group1 = g1,
          
          Group2 = g2,
          
          F = model$F[1],
          
          R2 = model$R2[1],
          
          P = model$`Pr(>F)`[1]
          
        )
        
      )
    
  }
  
}

###############################################################################
# Multiple testing correction
###############################################################################

pairwise_results$P_adjusted <-
  
  p.adjust(
    
    pairwise_results$P,
    
    method = "BH"
    
  )

###############################################################################
# Order by significance
###############################################################################

pairwise_results <-
  
  pairwise_results %>%
  
  arrange(
    
    P_adjusted,
    
    desc(R2)
    
  )

###############################################################################
# Export
###############################################################################

write.csv(
  
  pairwise_results,
  
  "AMR_pairwise_PERMANOVA.csv",
  
  row.names = FALSE
  
)

capture.output(
  
  pairwise_results,
  
  file = "AMR_pairwise_PERMANOVA.txt"
  
)

###############################################################################
# Console summary
###############################################################################

cat("\n")

cat("=====================================================\n")
cat("PAIRWISE PERMANOVA\n")
cat("=====================================================\n\n")

print(pairwise_results)

cat("\n")

cat(
  
  "Significant comparisons (BH-adjusted p ≤ 0.05):",
  
  sum(pairwise_results$P_adjusted <= 0.05),
  
  "\n"
  
)

###############################################################################
# END OF CHUNK 4
###############################################################################
###############################################################################
# 27. ADJUSTED PERMANOVA MODELS
#
# Evaluate whether Sample type and Ecological group explain
# additional variation after accounting for Study.
###############################################################################

###############################################################################
# MODEL A
# Sequential test
#
# Study -> Sample type
###############################################################################

model_study_sample_seq <-
  adonis2(
    
    dist_mat ~ Study + Sample_type,
    
    data = global_meta,
    
    permutations = 9999
    
  )

###############################################################################
# MODEL B
# Sequential test
#
# Study -> Ecological group
###############################################################################

model_study_ecology_seq <-
  adonis2(
    
    dist_mat ~ Study + Ecological_group,
    
    data = global_meta,
    
    permutations = 9999
    
  )

###############################################################################
# MODEL C
# Marginal effects
#
# Study + Sample type
###############################################################################

model_study_sample_margin <-
  adonis2(
    
    dist_mat ~ Study + Sample_type,
    
    data = global_meta,
    
    permutations = 9999,
    
    by = "margin"
    
  )

###############################################################################
# MODEL D
# Marginal effects
#
# Study + Ecological group
###############################################################################

model_study_ecology_margin <-
  adonis2(
    
    dist_mat ~ Study + Ecological_group,
    
    data = global_meta,
    
    permutations = 9999,
    
    by = "margin"
    
  )

###############################################################################
# EXPORT
###############################################################################

sink(
  "AMR_Adjusted_PERMANOVA.txt"
)

cat(
  "=============================================================\n",
  "ADJUSTED PERMANOVA MODELS\n",
  "=============================================================\n\n"
)

###############################################################################
cat(
  "MODEL A\n",
  "Sequential\n",
  "Study + Sample type\n\n"
)

print(model_study_sample_seq)

###############################################################################
cat(
  "\n=============================================================\n",
  "MODEL B\n",
  "Sequential\n",
  "Study + Ecological group\n\n"
)

print(model_study_ecology_seq)

###############################################################################
cat(
  "\n=============================================================\n",
  "MODEL C\n",
  "Marginal effects\n",
  "Study + Sample type\n\n"
)

print(model_study_sample_margin)

###############################################################################
cat(
  "\n=============================================================\n",
  "MODEL D\n",
  "Marginal effects\n",
  "Study + Ecological group\n\n"
)

print(model_study_ecology_margin)

sink()

###############################################################################
# CONSOLE SUMMARY
###############################################################################

cat("\n")

cat("====================================================\n")
cat("ADJUSTED PERMANOVA COMPLETE\n")
cat("====================================================\n\n")

cat("Output file:\n")
cat("AMR_Adjusted_PERMANOVA.txt\n\n")

cat("Sequential models:\n")
cat("- Study -> Sample type\n")
cat("- Study -> Ecological group\n\n")

cat("Marginal models:\n")
cat("- Study + Sample type\n")
cat("- Study + Ecological group\n")

###############################################################################
# END OF CHUNK 5
###############################################################################
###############################################################################
# 28. FINAL STATISTICAL REPORT
###############################################################################

sink("AMR_Statistical_Report.txt")

###############################################################################
# HEADER
###############################################################################

cat(
  "=======================================================================\n",
  "PLASMID-ASSOCIATED RESISTANCE PROFILE ANALYSIS\n",
  "FINAL STATISTICAL REPORT\n",
  "=======================================================================\n\n"
)

cat(
  "Response matrix : Plasmid-normalized resistance target richness\n",
  "Distance metric : Bray-Curtis\n",
  "Clustering      : UPGMA\n",
  "Ordination      : PCoA\n\n"
)

###############################################################################
# VALIDATION
###############################################################################

cat(
  "=======================================================================\n",
  "1. DATA VALIDATION\n",
  "=======================================================================\n\n"
)

cat(
  "Samples analysed: ",
  ncol(amr_mat_norm),
  "\n"
)

cat(
  "Resistance targets: ",
  nrow(amr_mat_norm),
  "\n\n"
)

cat(
  "Cophenetic correlation:\n"
)

print(
  round(cophenetic_r,3)
)

cat("\n")

cat(
  "Variance explained by first two PCoA axes:\n"
)

print(
  variance_explained[1:2]
)

cat("\n")

###############################################################################
# PRIMARY PERMANOVA
###############################################################################

cat(
  "=======================================================================\n",
  "2. PRIMARY PERMANOVA\n",
  "=======================================================================\n\n"
)

summary_table <- data.frame(
  
  Model = c(
    
    "Study",
    
    "Sample type",
    
    "Ecological group",
    
    "Alexander season"
    
  ),
  
  R2 = c(
    
    permanova_study$R2[1],
    
    permanova_sampletype$R2[1],
    
    permanova_ecology$R2[1],
    
    permanova_season$R2[1]
    
  ),
  
  Pvalue = c(
    
    permanova_study$`Pr(>F)`[1],
    
    permanova_sampletype$`Pr(>F)`[1],
    
    permanova_ecology$`Pr(>F)`[1],
    
    permanova_season$`Pr(>F)`[1]
    
  ),
  
  Dispersion_P = c(
    
    perm_study$tab[1,"Pr(>F)"],
    
    perm_sampletype$tab[1,"Pr(>F)"],
    
    perm_ecology$tab[1,"Pr(>F)"],
    
    perm_season$tab[1,"Pr(>F)"]
    
  )
  
)

summary_table$Interpretation <- c(
  
  "Robust",
  
  "Robust",
  
  "Interpret with caution (heterogeneous dispersion)",
  
  "Not significant"
  
)

print(
  summary_table,
  row.names=FALSE
)

cat("\n")

###############################################################################
# MARGINAL PERMANOVA
###############################################################################

cat(
  "=======================================================================\n",
  "3. MARGINAL PERMANOVA\n",
  "=======================================================================\n\n"
)

cat(
  "Study + Sample type\n\n"
)

print(
  model_study_sample_margin
)

cat("\n")

cat(
  "Study + Ecological group\n\n"
)

print(
  model_study_ecology_margin
)

cat("\n")

###############################################################################
# PAIRWISE PERMANOVA
###############################################################################

cat(
  "=======================================================================\n",
  "4. PAIRWISE PERMANOVA\n",
  "=======================================================================\n\n"
)

print(
  pairwise_results
)

cat("\n")

cat(
  "Significant comparisons after BH correction:\n\n"
)

print(
  
  pairwise_results %>%
    
    filter(
      P_adjusted <= 0.05
    )
  
)

cat("\n")

###############################################################################
# BIOLOGICAL SUMMARY
###############################################################################

cat(
  "=======================================================================\n",
  "5. BIOLOGICAL SUMMARY\n",
  "=======================================================================\n\n"
)

cat(
  
  "- Resistance composition differed significantly between the Alexander and Bogotá datasets.\n",
  
  "- River and Hospital samples retained significant differences after accounting for Study.\n",
  
  "- Ecological group explained additional variation after controlling for Study,\n",
  
  "  although interpretation should consider heterogeneous multivariate dispersion.\n",
  
  "- No significant seasonal effect was detected within the Alexander River dataset.\n",
  
  "- Following multiple-testing correction, Hospital vs Midstream was the only\n",
  
  "  significant pairwise ecological comparison.\n"
  
)

cat("\n")

###############################################################################
# OUTPUT FILES
###############################################################################

cat(
  "=======================================================================\n",
  "GENERATED OUTPUTS\n",
  "=======================================================================\n\n"
)

cat(
  
  "AMR_Global_PCoA.png\n",
  
  "AMR_Global_Dendrogram.png\n",
  
  "AMR_PCoA_coordinates.csv\n",
  
  "AMR_PERMANOVA_results.txt\n",
  
  "AMR_Adjusted_PERMANOVA.txt\n",
  
  "AMR_pairwise_PERMANOVA.csv\n",
  
  "AMR_pairwise_PERMANOVA.txt\n"
  
)

sink()

###############################################################################
# CONSOLE
###############################################################################

cat("\n")

cat("=====================================================\n")
cat("FINAL REPORT CREATED\n")
cat("=====================================================\n")
cat("Output:\n")
cat("AMR_Statistical_Report.txt\n")
cat("=====================================================\n")

###############################################################################
# END OF ANALYSIS
###############################################################################
###############################################################################
# 29. SUPPLEMENTARY TABLE
#
# Table Sx
#
# Summary of statistical analyses for plasmid-associated
# resistance profile composition.
###############################################################################

supp_table <- data.frame(
  
  Analysis = c(
    
    "PERMANOVA",
    "PERMANOVA",
    "PERMANOVA",
    "PERMANOVA",
    
    "Marginal PERMANOVA",
    "Marginal PERMANOVA",
    "Marginal PERMANOVA",
    
    "Pairwise PERMANOVA"
    
  ),
  
  Factor = c(
    
    "Study",
    
    "Sample type",
    
    "Ecological group",
    
    "Alexander season",
    
    "Study",
    
    "Sample type",
    
    "Ecological group",
    
    "Hospital vs Midstream"
    
  ),
  
  R2 = c(
    
    permanova_study$R2[1],
    
    permanova_sampletype$R2[1],
    
    permanova_ecology$R2[1],
    
    permanova_season$R2[1],
    
    model_study_sample_margin$R2[1],
    
    model_study_sample_margin$R2[2],
    
    model_study_ecology_margin$R2[2],
    
    pairwise_results$R2[
      pairwise_results$Group1=="Hospital" &
        pairwise_results$Group2=="Midstream"
    ]
    
  ),
  
  P_value = c(
    
    permanova_study$`Pr(>F)`[1],
    
    permanova_sampletype$`Pr(>F)`[1],
    
    permanova_ecology$`Pr(>F)`[1],
    
    permanova_season$`Pr(>F)`[1],
    
    model_study_sample_margin$`Pr(>F)`[1],
    
    model_study_sample_margin$`Pr(>F)`[2],
    
    model_study_ecology_margin$`Pr(>F)`[2],
    
    pairwise_results$P_adjusted[
      pairwise_results$Group1=="Hospital" &
        pairwise_results$Group2=="Midstream"
    ]
    
  ),
  
  Dispersion_P = c(
    
    perm_study$tab[1,"Pr(>F)"],
    
    perm_sampletype$tab[1,"Pr(>F)"],
    
    perm_ecology$tab[1,"Pr(>F)"],
    
    perm_season$tab[1,"Pr(>F)"],
    
    NA,
    
    NA,
    
    NA,
    
    NA
    
  ),
  
  Interpretation = c(
    
    "Robust",
    
    "Robust",
    
    "Interpret with caution (heterogeneous dispersion)",
    
    "Not significant",
    
    "Unique study effect",
    
    "Unique sample-type effect",
    
    "Unique ecological effect",
    
    "Significant after BH correction"
    
  ),
  
  stringsAsFactors = FALSE
  
)

###############################################################################
# Round numerical columns
###############################################################################

supp_table$R2 <-
  round(
    supp_table$R2,
    3
  )

supp_table$P_value <-
  signif(
    supp_table$P_value,
    3
  )

supp_table$Dispersion_P <-
  signif(
    supp_table$Dispersion_P,
    3
  )

###############################################################################
# Export
###############################################################################

write.csv(
  
  supp_table,
  
  "Supplementary_Table_AMR_Statistics.csv",
  
  row.names = FALSE
  
)

###############################################################################
# Console
###############################################################################

cat("\n")

cat("=====================================================\n")
cat("SUPPLEMENTARY TABLE CREATED\n")
cat("=====================================================\n\n")

print(
  supp_table,
  row.names = FALSE
)

cat("\n")

cat(
  "Output:\n",
  "Supplementary_Table_AMR_Statistics.csv\n"
)

###############################################################################
# END
###############################################################################