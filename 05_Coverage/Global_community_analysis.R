###############################################################################
#
# GLOBAL PLASMIDOME
#
# COMMUNITY ANALYSIS
#
# Purpose:
#
# Characterize the global organization of plasmid community composition
# across Alexander River and Bogotá River samples using an unsupervised
# community ecology framework.
#
# Community definition:
#
# Plasmid present = coverage ≥ 70%
#
# Community distance:
#
# Jaccard (binary)
#
###############################################################################

###############################################################################
# 1. LIBRARIES
###############################################################################

library(dplyr)
library(readr)
library(vegan)
library(ggplot2)
library(dynamicTreeCut)

###############################################################################
# 2. WORKING DIRECTORY
###############################################################################

setwd(
  "C:/Users/user/Desktop/Research Project Oded/metadata/Global"
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

alex_meta <- read.csv(
  
  "C:/Users/user/Desktop/Research Project Oded/metadata/Alexander/New/alex_clean_metadata.csv",
  
  stringsAsFactors = FALSE
  
)

bogo_meta <- read.csv(
  
  "C:/Users/user/Desktop/Research Project Oded/metadata/Bogota/bogo_clean_metadata.csv",
  
  stringsAsFactors = FALSE
  
)

###############################################################################
# 5. DERIVE UNIFIED METADATA
###############################################################################

#############################
# Alexander
#############################

alex_meta <- alex_meta %>%
  
  mutate(
    
    Study = "Alexander",
    
    Sample_type = "River",
    
    Ecological_group = region
    
  )

#############################
# Bogotá
#############################

bogo_meta <- bogo_meta %>%
  
  mutate(
    
    Study = "Bogota",
    
    Sample_type = ifelse(
      
      grepl("^Hospital", region),
      
      "Hospital",
      
      "River"
      
    ),
    
    Ecological_group = case_when(
      
      region %in% c(
        
        "River Source",
        
        "River Middle"
        
      ) ~
        
        "Upstream",
      
      region == "River Lower" ~
        
        "Lower_River",
      
      grepl(
        
        "Hospital",
        
        region
        
      ) ~
        
        "Hospital"
      
    )
    
  )

###############################################################################
# 6. KEEP COMMON VARIABLES
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

alex_meta <- alex_meta[,common_cols]

bogo_meta <- bogo_meta[,common_cols]

###############################################################################
# 7. MERGE METADATA
###############################################################################

global_meta <- bind_rows(
  
  alex_meta,
  
  bogo_meta
  
)

###############################################################################
# 8. ALIGN TO COVERAGE MATRIX
###############################################################################

global_meta <- global_meta %>%
  
  filter(
    
    sample_id %in%
      
      colnames(coverage)
    
  ) %>%
  
  arrange(
    
    match(
      
      sample_id,
      
      colnames(coverage)
      
    )
    
  )

coverage <- coverage[
  
  ,
  
  global_meta$sample_id
  
]

###############################################################################
# 9. PRESENCE / ABSENCE MATRIX
###############################################################################

coverage_pa <-
  
  coverage >= 70

###############################################################################
# 10. BUILD JACCARD DISTANCE MATRIX
###############################################################################

dist_mat <- vegdist(
  
  t(coverage_pa),
  
  method = "jaccard",
  
  binary = TRUE
  
)

###############################################################################
# 11. DATA SUMMARY
###############################################################################

cat("\n")

cat("====================================================\n")
cat("GLOBAL COMMUNITY ANALYSIS\n")
cat("====================================================\n\n")

cat(
  
  "Samples:",
  
  ncol(coverage_pa),
  
  "\n"
  
)

cat(
  
  "Plasmids:",
  
  nrow(coverage_pa),
  
  "\n\n"
  
)

cat(
  
  "Studies\n"
  
)

print(
  
  table(global_meta$Study)
  
)

cat("\n")

cat(
  
  "Sample type\n"
  
)

print(
  
  table(global_meta$Sample_type)
  
)

cat("\n")

cat(
  
  "Ecological groups\n"
  
)

print(
  
  table(global_meta$Ecological_group)
  
)

cat("\n")

cat(
  
  "Regions\n"
  
)

print(
  
  table(global_meta$region)
  
)

cat("\n")

cat(
  
  "Distance metric: Jaccard (binary)\n"
  
)

cat(
  
  "Coverage threshold: 70%\n"
  
)

###############################################################################
# END OF CHUNK 1
###############################################################################
###############################################################################
# 12. RECONSTRUCT SAMPLE DENDROGRAM
#
# Same methodology used in the global heatmap:
#
# Pearson correlation
# ↓
# Distance = 1 - Pearson
# ↓
# UPGMA clustering
###############################################################################

sample_dist <- as.dist(
  
  1 - cor(
    
    coverage,
    
    method = "pearson",
    
    use = "pairwise.complete.obs"
    
  )
  
)

sample_hc <- hclust(
  
  sample_dist,
  
  method = "average"
  
)

###############################################################################
# 13. COPHENETIC CORRELATION
#
# Quantifies how faithfully the dendrogram preserves
# the original pairwise sample distances.
###############################################################################

cophenetic_r <- cor(
  
  sample_dist,
  
  cophenetic(sample_hc)
  
)

###############################################################################
# 14. EXPORT DENDROGRAM ORDER
###############################################################################

dendrogram_order <- data.frame(
  
  Position = seq_along(sample_hc$order),
  
  Sample = sample_hc$labels[sample_hc$order]
  
)

write.csv(
  
  dendrogram_order,
  
  "Global_sample_dendrogram_order.csv",
  
  row.names = FALSE
  
)

###############################################################################
# 15. SAVE DENDROGRAM
###############################################################################

png(
  
  "Global_sample_dendrogram.png",
  
  width = 2200,
  
  height = 1200,
  
  res = 300
  
)

plot(
  
  sample_hc,
  
  main = "Global Sample Dendrogram",
  
  xlab = "",
  
  sub = "",
  
  cex = 0.8
  
)

dev.off()

###############################################################################
# 16. PRINCIPAL COORDINATE ANALYSIS
###############################################################################

pcoa <- cmdscale(
  
  sample_dist,
  
  eig = TRUE,
  
  k = 2
  
)

###############################################################################
# 17. PCoA COORDINATES
###############################################################################

pcoa_df <- data.frame(
  
  sample_id = rownames(pcoa$points),
  
  PCoA1 = pcoa$points[,1],
  
  PCoA2 = pcoa$points[,2]
  
)

pcoa_df <- left_join(
  
  pcoa_df,
  
  global_meta,
  
  by = "sample_id"
  
)

###############################################################################
# 18. VARIANCE EXPLAINED
###############################################################################

variance_explained <- round(
  
  100 *
    
    pcoa$eig /
    
    sum(
      
      pcoa$eig[pcoa$eig > 0]
      
    ),
  
  2
  
)

###############################################################################
# 19. EXPORT PCoA COORDINATES
###############################################################################

write.csv(
  
  pcoa_df,
  
  "Global_PCoA_coordinates.csv",
  
  row.names = FALSE
  
)

###############################################################################
# 20. PCoA FIGURE
###############################################################################

p <- ggplot(
  
  pcoa_df,
  
  aes(
    
    x = PCoA1,
    
    y = PCoA2,
    
    color = Study,
    
    shape = Sample_type
    
  )
  
) +
  
  geom_point(
    
    size = 4
    
  ) +
  
  theme_bw() +
  
  labs(
    
    title = "Global plasmid community structure",
    
    subtitle = "Pearson similarity represented by PCoA",
    
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
  
  "Global_PCoA.png",
  
  p,
  
  width = 8,
  
  height = 6,
  
  dpi = 300
  
)

###############################################################################
# 21. CONSOLE SUMMARY
###############################################################################

cat("\n")

cat("=====================================================\n")
cat("GLOBAL COMMUNITY STRUCTURE\n")
cat("=====================================================\n\n")

cat(
  
  "Cophenetic correlation coefficient:\n\n"
  
)

print(
  
  round(cophenetic_r,3)
  
)

cat("\n")

cat(
  
  "Variance explained\n\n"
  
)

print(
  
  variance_explained[1:2]
  
)

cat("\n")

cat(
  
  "PCoA coordinate summary\n\n"
  
)

print(
  
  summary(
    
    pcoa_df[,c("PCoA1","PCoA2")]
    
  )
  
)

###############################################################################
# END OF CHUNK 2
###############################################################################
###############################################################################
# 22. PERMANOVA
###############################################################################

set.seed(42)

#############################
# STUDY
#############################

permanova_study <- adonis2(
  
  dist_mat ~ Study,
  
  data = global_meta,
  
  permutations = 9999
  
)

#############################
# SAMPLE TYPE
#############################

permanova_sampletype <- adonis2(
  
  dist_mat ~ Sample_type,
  
  data = global_meta,
  
  permutations = 9999
  
)

#############################
# ECOLOGICAL GROUP
#############################

permanova_ecology <- adonis2(
  
  dist_mat ~ Ecological_group,
  
  data = global_meta,
  
  permutations = 9999
  
)

###############################################################################
# 23. DISPERSION DIAGNOSTICS
###############################################################################

#############################
# STUDY
#############################

bd_study <- betadisper(
  
  dist_mat,
  
  global_meta$Study
  
)

anova_study <- anova(bd_study)

perm_study <- permutest(
  
  bd_study,
  
  permutations = 9999
  
)

#############################
# SAMPLE TYPE
#############################

bd_sampletype <- betadisper(
  
  dist_mat,
  
  global_meta$Sample_type
  
)

anova_sampletype <- anova(bd_sampletype)

perm_sampletype <- permutest(
  
  bd_sampletype,
  
  permutations = 9999
  
)

#############################
# ECOLOGICAL GROUP
#############################

bd_ecology <- betadisper(
  
  dist_mat,
  
  global_meta$Ecological_group
  
)

anova_ecology <- anova(bd_ecology)

perm_ecology <- permutest(
  
  bd_ecology,
  
  permutations = 9999
  
)

###############################################################################
# 24. EXPORT REPORT
###############################################################################

sink("Global_PERMANOVA_results.txt")

cat(
  "========================================================\n",
  "GLOBAL PLASMID COMMUNITY ANALYSIS\n",
  "========================================================\n\n"
)

cat(
  "Community definition: Presence/absence (coverage ≥70%)\n",
  "Community distance: Jaccard\n",
  "Clustering: Pearson + UPGMA (heatmap)\n",
  "Ordination: PCoA\n\n"
)

####################################################

cat(
  "========================================================\n",
  "MODEL 1: STUDY\n",
  "========================================================\n\n"
)

print(permanova_study)

cat("\n\nDispersion ANOVA\n\n")

print(anova_study)

cat("\n\nDispersion permutation test\n\n")

print(perm_study)

####################################################

cat(
  "\n\n========================================================\n",
  "MODEL 2: SAMPLE TYPE\n",
  "========================================================\n\n"
)

print(permanova_sampletype)

cat("\n\nDispersion ANOVA\n\n")

print(anova_sampletype)

cat("\n\nDispersion permutation test\n\n")

print(perm_sampletype)

####################################################

cat(
  "\n\n========================================================\n",
  "MODEL 3: ECOLOGICAL GROUP\n",
  "========================================================\n\n"
)

print(permanova_ecology)

cat("\n\nDispersion ANOVA\n\n")

print(anova_ecology)

cat("\n\nDispersion permutation test\n\n")

print(perm_ecology)

####################################################

cat(
  "\n\n========================================================\n",
  "DENDROGRAM VALIDATION\n",
  "========================================================\n\n"
)

cat(
  
  "Cophenetic correlation coefficient:\n"
  
)

print(cophenetic_r)

cat("\n")

cat(
  
  "Variance explained by PCoA axes:\n"
  
)

print(
  
  variance_explained[1:2]
  
)

sink()

###############################################################################
# 25. CONSOLE SUMMARY
###############################################################################

cat("\n")
cat("=========================================\n")
cat("GLOBAL ANALYSIS COMPLETE\n")
cat("=========================================\n")
cat("Generated files:\n")
cat(" - Global_PERMANOVA_results.txt\n")
cat(" - Global_PCoA.png\n")
cat(" - Global_sample_dendrogram.png\n")
cat("=========================================\n")
###############################################################################
# END OF CHUNK 3
###############################################################################
###############################################################################
# 26. PAIRWISE PERMANOVA
#
# Compare all ecological groups pairwise
# using adonis2()
###############################################################################

group_levels <- levels(
  factor(global_meta$Ecological_group)
)

pairwise_results <- data.frame()

###############################################################################
# Loop over every pair of ecological groups
###############################################################################

for(i in 1:(length(group_levels)-1)){
  
  for(j in (i+1):length(group_levels)){
    
    g1 <- group_levels[i]
    g2 <- group_levels[j]
    
    ###########################################################################
    # Subset metadata
    ###########################################################################
    
    meta_sub <-
      
      global_meta %>%
      
      filter(
        
        Ecological_group %in%
          
          c(g1, g2)
        
      )
    
    ###########################################################################
    # Subset community matrix
    ###########################################################################
    
    comm_sub <-
      
      coverage_pa[,
                  
                  meta_sub$sample_id
                  
      ]
    
    ###########################################################################
    # Distance matrix
    ###########################################################################
    
    dist_sub <-
      
      vegdist(
        
        t(comm_sub),
        
        method = "jaccard",
        
        binary = TRUE
        
      )
    
    ###########################################################################
    # PERMANOVA
    ###########################################################################
    
    model <-
      
      adonis2(
        
        dist_sub ~
          
          Ecological_group,
        
        data = meta_sub,
        
        permutations = 9999
        
      )
    
    ###########################################################################
    # Store results
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
# Multiple-testing correction
###############################################################################

pairwise_results$P_adjusted <-
  
  p.adjust(
    
    pairwise_results$P,
    
    method = "BH"
    
  )

###############################################################################
# Order by adjusted p-value
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
  
  "Global_pairwise_PERMANOVA.csv",
  
  row.names = FALSE
  
)

capture.output(
  
  pairwise_results,
  
  file = "Global_pairwise_PERMANOVA.txt"
  
)

###############################################################################
# Console summary
###############################################################################

cat("\n")
cat("=========================================\n")
cat("PAIRWISE PERMANOVA\n")
cat("=========================================\n\n")

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
# 27. PAIRWISE DISPERSION DIAGNOSTICS
###############################################################################

###############################################################################
# Significant comparisons from pairwise PERMANOVA
###############################################################################

sig_pairs <- pairwise_results %>%
  
  filter(
    
    P_adjusted <= 0.05
    
  )

###############################################################################
# Empty results table
###############################################################################

pairwise_dispersion <- data.frame()

###############################################################################
# Loop over significant comparisons
###############################################################################

for(i in seq_len(nrow(sig_pairs))){
  
  g1 <- sig_pairs$Group1[i]
  
  g2 <- sig_pairs$Group2[i]
  
  ###########################################################################
  # Subset metadata
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
    
    coverage_pa[,
                
                meta_sub$sample_id
                
    ]
  
  ###########################################################################
  # Distance matrix
  ###########################################################################
  
  dist_sub <-
    
    vegdist(
      
      t(comm_sub),
      
      method = "jaccard",
      
      binary = TRUE
      
    )
  
  ###########################################################################
  # betadisper
  ###########################################################################
  
  bd <-
    
    betadisper(
      
      dist_sub,
      
      meta_sub$Ecological_group
      
    )
  
  perm <-
    
    permutest(
      
      bd,
      
      permutations = 9999
      
    )
  
  ###########################################################################
  # Store results
  ###########################################################################
  
  pairwise_dispersion <-
    
    rbind(
      
      pairwise_dispersion,
      
      data.frame(
        
        Group1 = g1,
        
        Group2 = g2,
        
        F = perm$tab[1,"F"],
        
        P = perm$tab[1,"Pr(>F)"]
        
      )
      
    )
  
}

###############################################################################
# Export
###############################################################################

write.csv(
  
  pairwise_dispersion,
  
  "Global_pairwise_dispersion.csv",
  
  row.names = FALSE
  
)

###############################################################################
# Console
###############################################################################

cat("\n")
cat("=========================================\n")
cat("PAIRWISE DISPERSION TESTS\n")
cat("=========================================\n\n")

print(pairwise_dispersion)

cat("\n")

cat(
  
  "Comparisons violating homogeneity (p < 0.05):",
  
  sum(pairwise_dispersion$P < 0.05),
  
  "\n"
  
)

###############################################################################
# END OF CHUNK 5
###############################################################################
###############################################################################
# 28. FINAL STATISTICAL REPORT
###############################################################################

sink("Global_Statistical_Report.txt")

###############################################################################
# HEADER
###############################################################################

cat(
  "=======================================================================\n",
  "GLOBAL PLASMIDOME COMMUNITY ANALYSIS\n",
  "=======================================================================\n\n"
)

cat(
  "Purpose:\n",
  "Characterize the global organization of plasmid communities across\n",
  "Alexander River and Bogotá River samples and evaluate whether\n",
  "community structure is associated with ecological context.\n\n"
)

###############################################################################
# DATASET SUMMARY
###############################################################################

cat(
  "=======================================================================\n",
  "DATASET SUMMARY\n",
  "=======================================================================\n\n"
)

cat(
  "Samples: ", ncol(coverage_pa), "\n",
  "Plasmids: ", nrow(coverage_pa), "\n",
  "Coverage threshold: 70%\n",
  "Community matrix: Presence / Absence\n",
  "Community distance: Jaccard\n",
  "Heatmap clustering: Pearson similarity + UPGMA\n",
  "Ordination: Principal Coordinate Analysis (PCoA)\n\n",
  sep=""
)

cat("Studies\n")
print(table(global_meta$Study))

cat("\nSample type\n")
print(table(global_meta$Sample_type))

cat("\nEcological groups\n")
print(table(global_meta$Ecological_group))

###############################################################################
# DENDROGRAM VALIDATION
###############################################################################

cat(
  "\n=======================================================================\n",
  "DENDROGRAM VALIDATION\n",
  "=======================================================================\n\n"
)

cat(
  "Cophenetic correlation coefficient:\n\n"
)

print(round(cophenetic_r,3))

cat("\n")

if(cophenetic_r >= 0.90){
  
  cat(
    "Interpretation:\n",
    "Excellent agreement between the dendrogram and\n",
    "the original Pearson distance matrix.\n\n"
  )
  
}else if(cophenetic_r >= 0.80){
  
  cat(
    "Interpretation:\n",
    "Very good agreement between dendrogram and\n",
    "original distance matrix.\n\n"
  )
  
}else{
  
  cat(
    "Interpretation:\n",
    "Moderate agreement.\n\n"
  )
  
}

###############################################################################
# PCoA
###############################################################################

cat(
  "=======================================================================\n",
  "PCoA SUMMARY\n",
  "=======================================================================\n\n"
)

cat(
  "Variance explained\n\n"
)

print(round(variance_explained[1:2],2))

cat("\nTotal variance explained:\n")

print(round(sum(variance_explained[1:2]),2))

###############################################################################
# PERMANOVA
###############################################################################

cat(
  "\n=======================================================================\n",
  "PERMANOVA\n",
  "=======================================================================\n\n"
)

cat("--------------------------\n")
cat("MODEL 1 : STUDY\n")
cat("--------------------------\n\n")

print(permanova_study)

cat("\nDispersion\n\n")

print(perm_study)

cat("\n")

cat("--------------------------\n")
cat("MODEL 2 : SAMPLE TYPE\n")
cat("--------------------------\n\n")

print(permanova_sampletype)

cat("\nDispersion\n\n")

print(perm_sampletype)

cat("\n")

cat("--------------------------\n")
cat("MODEL 3 : ECOLOGICAL GROUP\n")
cat("--------------------------\n\n")

print(permanova_ecology)

cat("\nDispersion\n\n")

print(perm_ecology)

###############################################################################
# PAIRWISE PERMANOVA
###############################################################################

cat(
  "\n=======================================================================\n",
  "PAIRWISE PERMANOVA\n",
  "=======================================================================\n\n"
)

print(pairwise_results)

###############################################################################
# PAIRWISE DISPERSION
###############################################################################

cat(
  "\n=======================================================================\n",
  "PAIRWISE DISPERSION TESTS\n",
  "=======================================================================\n\n"
)

print(pairwise_dispersion)

###############################################################################
# FINAL STATISTICAL SUMMARY
###############################################################################

cat(
  "\n=======================================================================\n",
  "FINAL STATISTICAL SUMMARY\n",
  "=======================================================================\n\n"
)

summary_table <- data.frame(
  
  Analysis=c(
    
    "Study",
    "Sample_type",
    "Ecological_group"
    
  ),
  
  PERMANOVA_P=c(
    
    permanova_study$`Pr(>F)`[1],
    permanova_sampletype$`Pr(>F)`[1],
    permanova_ecology$`Pr(>F)`[1]
    
  ),
  
  R2=c(
    
    permanova_study$R2[1],
    permanova_sampletype$R2[1],
    permanova_ecology$R2[1]
    
  ),
  
  Dispersion_P=c(
    
    perm_study$tab[1,"Pr(>F)"],
    perm_sampletype$tab[1,"Pr(>F)"],
    perm_ecology$tab[1,"Pr(>F)"]
    
  )
  
)

summary_table$Interpretation <- ifelse(
  
  summary_table$PERMANOVA_P <= 0.05 &
    summary_table$Dispersion_P > 0.05,
  
  "Robust ecological signal",
  
  ifelse(
    
    summary_table$PERMANOVA_P <= 0.05 &
      summary_table$Dispersion_P <= 0.05,
    
    "Interpret with caution (dispersion differs)",
    
    "Not significant"
    
  )
  
)

print(summary_table,row.names=FALSE)

###############################################################################
# WORKFLOW
###############################################################################

cat(
  "\n=======================================================================\n",
  "ANALYSIS WORKFLOW\n",
  "=======================================================================\n\n"
)

cat(
  "Coverage matrix\n",
  "      ↓\n",
  "Presence / Absence conversion (≥70%)\n",
  "      ↓\n",
  "Jaccard community distance\n",
  "      ↓\n",
  "Pearson sample clustering (UPGMA)\n",
  "      ↓\n",
  "Cophenetic validation\n",
  "      ↓\n",
  "PCoA\n",
  "      ↓\n",
  "PERMANOVA\n",
  "      ↓\n",
  "Dispersion diagnostics\n",
  "      ↓\n",
  "Pairwise PERMANOVA\n",
  "      ↓\n",
  "Pairwise dispersion diagnostics\n"
)

sink()

###############################################################################
# CONSOLE
###############################################################################

cat("\n")
cat("=============================================================\n")
cat("GLOBAL STATISTICAL REPORT CREATED\n")
cat("=============================================================\n")
cat("Output file:\n")
cat(" - Global_Statistical_Report.txt\n")
cat("=============================================================\n")