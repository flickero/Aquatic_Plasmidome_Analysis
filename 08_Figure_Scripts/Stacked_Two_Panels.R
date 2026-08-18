###############################################################################
# Figure 4.6
# Sample-level COG functional composition
#
# Generates a two-panel publication figure:
# (A) Including Function Unknown (S)
# (B) Excluding Function Unknown (S)
###############################################################################

# ============================================================================
# 1. Libraries
# ============================================================================

library(tidyverse)
library(scales)
library(patchwork)

# ============================================================================
# 2. Global parameters
# ============================================================================

WORKDIR <- "C:/Users/user/Desktop/Enog"

EGGNOG_FILE <- "all_ORFs.emapper.annotations.filtered"
COVERAGE_FILE <- "merged.tsv"

ALEX_METADATA <- "alex_clean_metadata.csv"
BOGO_METADATA <- "bogo_clean_metadata.csv"

COVERAGE_THRESHOLD <- 70
RARE_CATEGORY_THRESHOLD <- 0.02

setwd(WORKDIR)

# ============================================================================
# 3. Metadata
# ============================================================================

alex_meta <- read.csv(
  ALEX_METADATA,
  stringsAsFactors = FALSE
)

bogo_meta <- read.csv(
  BOGO_METADATA,
  stringsAsFactors = FALSE
)

# ============================================================================
# 4. eggNOG annotation import
# ============================================================================

lines <- readLines(EGGNOG_FILE)

header_line <- grep("^#query", lines)[1]

header <- strsplit(lines[header_line], "\t")[[1]]
header[1] <- sub("^#", "", header[1])

annotations <- read.delim(
  EGGNOG_FILE,
  skip = header_line,
  header = FALSE,
  sep = "\t",
  fill = TRUE,
  quote = "",
  stringsAsFactors = FALSE
)

colnames(annotations) <- header

annotations <- annotations %>%
  filter(!grepl("^#", query)) %>%
  filter(!query %in% c("query", "#query"))

# ============================================================================
# 5. Coverage matrix
# ============================================================================

coverage <- read.delim(
  COVERAGE_FILE,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

colnames(coverage)[1] <- "plasmid"

coverage <- coverage %>%
  filter(!grepl("^#", plasmid)) %>%
  filter(!plasmid %in% c("query", "#query"))
# ============================================================================
# 6. Publication sample metadata
# ============================================================================

# ---------- Alexander River ----------

alex_labels <- alex_meta %>%
  transmute(
    sample_id,
    dataset = "Alexander River",
    
    # Publication label
    sample_label = paste0(
      "S",
      station_numeric,
      "\n",
      season
    ),
    
    # Biological ordering
    sample_order = paste0(
      "A_",
      sprintf("%02d", station_numeric),
      "_",
      match(
        season,
        c(
          "Winter",
          "Spring",
          "LateSummer"
        )
      )
    )
  )

# ---------- Bogotá River ----------

bogo_labels <- bogo_meta %>%
  transmute(
    sample_id,
    dataset = "Bogotá River",
    
    sample_label = paste(
      station_name,
      sampling_year,
      sep = "\n"
    ),
    
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

# ---------- Combined lookup table ----------

sample_metadata <- bind_rows(
  alex_labels,
  bogo_labels
)

# Verify uniqueness

stopifnot(
  !anyDuplicated(sample_metadata$sample_id)
)

# ============================================================================
# 7. Coverage matrix (long format)
# ============================================================================

coverage_long <- coverage %>%
  pivot_longer(
    cols = -plasmid,
    names_to = "sample_cov",
    values_to = "coverage"
  ) %>%
  filter(
    coverage >= COVERAGE_THRESHOLD
  )

# ============================================================================
# 8. Parse eggNOG annotations
# ============================================================================

cog_long <- annotations %>%
  mutate(
    
    sample_origin = sub(
      "_RNODE.*",
      "",
      query
    ),
    
    plasmid = sub(
      "_[0-9]+$",
      "",
      query
    )
    
  ) %>%
  filter(
    !is.na(COG_category),
    COG_category != "",
    COG_category != "-"
  ) %>%
  mutate(
    COG_category = toupper(COG_category),
    COG_category = strsplit(
      COG_category,
      ""
    )
  ) %>%
  unnest(COG_category)
# ============================================================================
# COG category names
# ============================================================================

COG_NAMES <- c(
  "A" = "RNA processing and modification",
  "B" = "Chromatin structure and dynamics",
  "C" = "Energy production and conversion",
  "D" = "Cell cycle control, cell division, chromosome partitioning",
  "E" = "Amino acid transport and metabolism",
  "F" = "Nucleotide transport and metabolism",
  "G" = "Carbohydrate transport and metabolism",
  "H" = "Coenzyme transport and metabolism",
  "I" = "Lipid transport and metabolism",
  "J" = "Translation, ribosomal structure and biogenesis",
  "K" = "Transcription",
  "L" = "Replication, recombination and repair",
  "M" = "Cell wall/membrane/envelope biogenesis",
  "N" = "Cell motility",
  "O" = "Posttranslational modification, protein turnover and chaperones",
  "P" = "Inorganic ion transport and metabolism",
  "Q" = "Secondary metabolites biosynthesis, transport and catabolism",
  "R" = "General function prediction only",
  "S" = "Function unknown",
  "T" = "Signal transduction mechanisms",
  "U" = "Intracellular trafficking, secretion and vesicular transport",
  "V" = "Defense mechanisms",
  "W" = "Extracellular structures",
  "Y" = "Nuclear structure",
  "Z" = "Cytoskeleton",
  "Other" = "Other (<2% each category)"
)
# ============================================================================
# 9. Function: prepare plotting table
# ============================================================================

prepare_plot_data <- function(cog_data) {
  
  # --------------------------------------------------------------------------
  # Assign plasmids to samples
  # --------------------------------------------------------------------------
  
  plot_data <- cog_data %>%
    left_join(
      coverage_long,
      by = "plasmid",
      relationship = "many-to-many"
    ) %>%
    mutate(
      final_sample = if_else(
        is.na(sample_cov),
        sample_origin,
        sample_cov
      )
    ) %>%
    filter(
      !grepl("^#", final_sample),
      final_sample != "query"
    )
  
  # --------------------------------------------------------------------------
  # Count COG categories per sample
  # --------------------------------------------------------------------------
  
  plot_data <- plot_data %>%
    count(
      final_sample,
      COG_category,
      name = "count"
    ) %>%
    group_by(final_sample) %>%
    mutate(
      prop = count / sum(count)
    ) %>%
    ungroup()
  
  # --------------------------------------------------------------------------
  # Merge rare categories
  # --------------------------------------------------------------------------
  
  plot_data <- plot_data %>%
    mutate(
      Category = if_else(
        prop < RARE_CATEGORY_THRESHOLD,
        "Other (<2% each category)",
        COG_category
      )
    ) %>%
    group_by(
      final_sample,
      Category
    ) %>%
    summarise(
      prop = sum(prop),
      .groups = "drop"
    )
  plot_data <- plot_data %>%
    mutate(
      Category = recode(
        Category,
        !!!COG_NAMES
      )
    )
      plot_data$Category <- factor(
        plot_data$Category,
        levels = COG_ORDER
      )     
    
  
  # --------------------------------------------------------------------------
  # Attach publication metadata
  # --------------------------------------------------------------------------
  
  plot_data <- plot_data %>%
    left_join(
      sample_metadata,
      by = c("final_sample" = "sample_id")
    )
  
  # --------------------------------------------------------------------------
  # Sanity check
  # --------------------------------------------------------------------------
  
  if (any(is.na(plot_data$sample_label))) {
    
    stop(
      "One or more samples could not be matched to the metadata."
    )
    
  }
  
  # --------------------------------------------------------------------------
  # Biological ordering
  # --------------------------------------------------------------------------
  
  plot_data <- plot_data %>%
    arrange(
      dataset,
      sample_order
    )
  
  plot_data$sample_label <- factor(
    plot_data$sample_label,
    levels = unique(plot_data$sample_label)
  )
  
  plot_data$dataset <- factor(
    plot_data$dataset,
    levels = c(
      "Alexander River",
      "Bogotá River"
    )
  )
  
  plot_data
  
}
# ============================================================================
# 10. Fixed COG color palette
# ============================================================================

# Use the EXACT palette from Figure 4.5
# (Replace this placeholder with your finalized named vector.)

COG_COLORS <- c(
  
  
  "Function unknown" = "#696166",
  
  "Replication, recombination and repair" = "#D6D2D3",
  
  "Cell wall/membrane/envelope biogenesis" = "#FF1D2D",
  
  "Transcription" = "#F000E8",
  
  "Posttranslational modification, protein turnover and chaperones" = "#19FF2A",
  
  "General function prediction only" = "#3F7DE3",
  
  "Inorganic ion transport and metabolism" = "#FFB020",
  
  "Defense mechanisms" = "#BC0979",
  
  "Energy production and conversion" = "#2DE4C4",
  
  "Amino acid transport and metabolism" = "#97B516",
  
  "Signal transduction mechanisms" = "#3BBCE2",
  
  "Lipid transport and metabolism" = "#C38DEB",
  
  "Translation, ribosomal structure and biogenesis" = "#9D15F5",
  
  "Carbohydrate transport and metabolism" = "#F59B9B",
  
  "Coenzyme transport and metabolism" = "#3F64A8",
  
  "Cell cycle control, cell division, chromosome partitioning" = "#CB4D1C",
  
  "Intracellular trafficking, secretion and vesicular transport" = "#27895B",
  
  "Secondary metabolites biosynthesis, transport and catabolism" = "#9B7A10",
  
  "Nucleotide transport and metabolism" = "#0B8C8C",
  
  "RNA processing and modification" = "#8C4F0B",
  
  "Chromatin structure and dynamics" = "#7F7F7F",
  
  "Cell motility" = "#00A6A6",
  
  "Cytoskeleton" = "#8E44AD",
  
  "Extracellular structures" = "#B8860B",
  
  "Nuclear structure" = "#2C3E50",
  
  "Other (<2% each category)" = "#444444"
  
)
COG_ORDER <- c(
  
  "Function unknown",
  
  "Replication, recombination and repair",
  
  "Cell wall/membrane/envelope biogenesis",
  
  "Transcription",
  
  "Posttranslational modification, protein turnover and chaperones",
  
  "General function prediction only",
  
  "Inorganic ion transport and metabolism",
  
  "Defense mechanisms",
  
  "Energy production and conversion",
  
  "Amino acid transport and metabolism",
  
  "Signal transduction mechanisms",
  
  "Lipid transport and metabolism",
  
  "Translation, ribosomal structure and biogenesis",
  
  "Carbohydrate transport and metabolism",
  
  "Coenzyme transport and metabolism",
  
  "Cell cycle control, cell division, chromosome partitioning",
  
  "Intracellular trafficking, secretion and vesicular transport",
  
  "Secondary metabolites biosynthesis, transport and catabolism",
  
  "Nucleotide transport and metabolism",
  
  "RNA processing and modification",
  
  "Chromatin structure and dynamics",
  
  "Cell motility",
  
  "Cytoskeleton",
  
  "Extracellular structures",
  
  "Nuclear structure",
  
  "Other (<2% each category)"
  
)

# ============================================================================
# 11. Function: create publication plot
# ============================================================================

make_plot <- function(plot_data) {
  
  ggplot(
    plot_data,
    aes(
      x = sample_label,
      y = prop,
      fill = Category
    )
  ) +
    
    geom_col(
      width = 0.88
    ) +
    
    facet_grid(
      ~dataset,
      scales = "free_x",
      space = "free_x"
    ) +
    
    scale_fill_manual(
      values = COG_COLORS,
      drop = TRUE
    ) +
    
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    
    labs(
      x = NULL,
      y = "Relative abundance"
    ) +
    
    theme_bw(base_size = 13) +
    
    theme(
      
      plot.title = element_text(
        face = "bold",
        size = 14,
        hjust = 0
      ),
      
      strip.background = element_blank(),
      
      strip.text = element_text(
        face = "bold",
        size = 12
      ),
      
      axis.text.x = element_text(
        angle = 60,
        hjust = 1,
        size = 7.5
      ),
      
      axis.title.y = element_text(
        face = "bold"
      ),
      
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      
      legend.title = element_text(
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 8
      ),
      
      legend.key.height = unit(
        0.45,
        "cm"
      )
    )
  
}
plot_data_with_unknown <-
  prepare_plot_data(cog_long)

plot_data_without_unknown <-
  prepare_plot_data(
    filter(
      cog_long,
      COG_category != "S"
    )
  )
plot_with_unknown <-
  make_plot(plot_data_with_unknown) +
  ggtitle("(A)")

plot_without_unknown <-
  make_plot(plot_data_without_unknown) +
  ggtitle("(B)") +
  theme(
    legend.position = "none"
  )
final_figure <-
  
  plot_with_unknown /
  
  plot_without_unknown

print(final_figure)
ggsave(
  filename = "Figure_4_6_COG_per_sample.png",
  plot = final_figure,
  width = 15,
  height = 12,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "Figure_4_6_COG_per_sample.pdf",
  plot = final_figure,
  width = 15,
  height = 12,
  bg = "white"
)