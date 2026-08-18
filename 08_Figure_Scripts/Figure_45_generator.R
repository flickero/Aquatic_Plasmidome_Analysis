###############################################################################
# Figure 4.5
# Global Functional Composition of the Plasmidome
# Panel A - Including Function Unknown
# Panel B - Excluding Function Unknown
###############################################################################

library(tidyverse)
library(scales)
library(patchwork)
library(grid)

###############################################################################
# User settings
###############################################################################

# Categories below this fraction will not receive a percentage label
LABEL_THRESHOLD <- 0.02

###############################################################################
# Working directory
###############################################################################

setwd("C:/Users/user/Desktop/Enog")

###############################################################################
# Read eggNOG annotations
###############################################################################

lines <- readLines("all_ORFs.emapper.annotations.filtered")

header_line <- grep("^#query", lines)[1]

header <- strsplit(lines[header_line], "\t")[[1]]
header[1] <- sub("^#", "", header[1])

annotations <- read.delim(
  "all_ORFs.emapper.annotations.filtered",
  skip = header_line,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  fill = TRUE,
  quote = ""
)

colnames(annotations) <- header

annotations <- annotations %>%
  filter(!grepl("^#", query)) %>%
  filter(query != "query")

###############################################################################
# Expand multi-letter COG assignments
###############################################################################

cog_long <- annotations %>%
  filter(
    !is.na(COG_category),
    COG_category != "",
    COG_category != "-"
  ) %>%
  mutate(
    COG_category = toupper(COG_category),
    COG_category = strsplit(COG_category, "")
  ) %>%
  unnest(COG_category)

###############################################################################
# COG dictionary
###############################################################################

cog_dict <- tribble(
  ~COG_category, ~Description,
  "J","Translation, ribosomal structure and biogenesis",
  "A","RNA processing and modification",
  "K","Transcription",
  "L","Replication, recombination and repair",
  "B","Chromatin structure and dynamics",
  "D","Cell cycle control, cell division, chromosome partitioning",
  "Y","Nuclear structure",
  "V","Defense mechanisms",
  "T","Signal transduction mechanisms",
  "M","Cell wall/membrane/envelope biogenesis",
  "N","Cell motility",
  "Z","Cytoskeleton",
  "W","Extracellular structures",
  "U","Intracellular trafficking, secretion and vesicular transport",
  "O","Posttranslational modification, protein turnover and chaperones",
  "C","Energy production and conversion",
  "G","Carbohydrate transport and metabolism",
  "E","Amino acid transport and metabolism",
  "F","Nucleotide transport and metabolism",
  "H","Coenzyme transport and metabolism",
  "I","Lipid transport and metabolism",
  "P","Inorganic ion transport and metabolism",
  "Q","Secondary metabolites biosynthesis, transport and catabolism",
  "R","General function prediction only",
  "S","Function unknown"
)

###############################################################################
# Fixed color palette
# (Colors remain identical between panels)
###############################################################################

fixed_colors <- c(
  
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
###############################################################################
# Function producing one donut
###############################################################################

make_donut <- function(remove_unknown = FALSE,
                       panel_label = "(A)") {
  
  dat <- cog_long
  
  if (remove_unknown) {
    dat <- dat %>%
      filter(COG_category != "S")
  }
  
  ###########################################################################
  # Count categories
  ###########################################################################
  
  counts <- dat %>%
    count(COG_category) %>%
    inner_join(cog_dict, by = "COG_category") %>%
    mutate(prop = n / sum(n))
  
  ###########################################################################
  # Merge small categories
  ###########################################################################
  
  counts <- counts %>%
    mutate(
      Description = ifelse(
        prop < 0.02,
        "Other (<2% each category)",
        Description
      )
    ) %>%
    group_by(Description) %>%
    summarise(
      n = sum(n),
      prop = sum(prop),
      .groups = "drop"
    )
  
  ###########################################################################
  # Preserve identical color mapping between panels
  ###########################################################################
  
  counts <- counts %>%
    arrange(desc(prop))
  counts$Description <- factor(
    counts$Description,
    levels = counts$Description
  )
  
  ###########################################################################
  # Plot
  ###########################################################################
  
  ggplot(
    counts,
    aes(
      x = 1.7,
      y = prop,
      fill = Description
    )
  ) +
    
    geom_col(
      colour = "white",
      linewidth = 0.35,
      width = 1
    ) +
    
    coord_polar(theta = "y") +
    
    xlim(0.55, 2.45) +
    
    #########################################################################
  # Percentage labels
  #########################################################################
  
  geom_text(
    aes(
      label = ifelse(
        prop >= LABEL_THRESHOLD,
        percent(prop, accuracy = 0.1),
        ""
      )
    ),
    position = position_stack(vjust = 0.5),
    size = 3.3,
    fontface = "bold",
    colour = "black"
  ) +
    
    #########################################################################
  # Fixed colors
  #########################################################################
  
  scale_fill_manual(
    values = fixed_colors,
    drop = TRUE,
    name = "COG functional category"
  ) +
    
    #########################################################################
  # Theme
  #########################################################################
  
  theme_void(base_size = 13) +
    
    theme(
      
      plot.title = element_text(
        face = "bold",
        size = 16,
        hjust = 0.02
      ),
      
      legend.position = c(1.12, 0.8),
      legend.justification = c("left", "top"),
      
      legend.title = element_text(
        face = "bold",
        size = 11
      ),
      
      legend.text = element_text(
        size = 9
      ),
      
      legend.key.size = unit(0.45, "cm"),
      
      plot.margin = margin(8,8,8,8)
    ) +
    
    labs(
      title = panel_label
    )
  
}
###############################################################################
# Create panels
###############################################################################

panelA <- make_donut(
  remove_unknown = FALSE,
  panel_label = "(A)"
)

panelB <- make_donut(
  remove_unknown = TRUE,
  panel_label = "(B)"
) +
  theme(legend.position = "none")

###############################################################################
# Collect legend and combine
###############################################################################

combined <-
  
  (panelA | panelB) +
  
  plot_layout(
    widths = c(1,1)
  )

###############################################################################
# Export
###############################################################################

ggsave(
  
  filename = "Figure4.5_Global_COG_Functions_Tall.tiff",
  
  plot = combined,
  
  dpi = 300,
  
  compression = "lzw",
  
  width = 26,
  
  height = 9,

  bg = "white"
  
)

ggsave(
  
  filename = "Figure4.5_Global_COG_Functions_Tall.pdf",
  
  plot = combined,
  
  device = cairo_pdf,
  
  width = 26,
  
  height = 9,
  
  bg = "white"
  
)

###############################################################################
# Preview
###############################################################################

print(combined)