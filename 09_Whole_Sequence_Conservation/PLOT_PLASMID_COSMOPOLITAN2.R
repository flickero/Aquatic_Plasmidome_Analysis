##############################################################################
#
# plot_cosmopolitan_plasmids.R
#
# Cosmopolitan plasmid visualization
#
# Inputs
#   feature_summary.tsv
#   master_annotations.tsv
#   feature_labels.tsv
#   sample_metadata.tsv
#
# Output
#   One publication-quality PNG per plasmid
#
###############################################################################

rm(list = ls())

###############################################################################
# Libraries
###############################################################################

library(ggplot2)
library(gggenes)
library(dplyr)
library(stringr)
library(grid)

###############################################################################
# USER SETTINGS
###############################################################################

WORKDIR <- "C:/Users/user/Downloads"

FEATURE_SUMMARY <-
  file.path(
    WORKDIR,
    "feature_summary.tsv"
  )

MASTER_ANNOTATIONS <-
  file.path(
    WORKDIR,
    "master_annotations.tsv"
  )

FEATURE_LABELS <-
  file.path(
    WORKDIR,
    "feature_labels_curated.tsv"
  )

SAMPLE_METADATA <-
  file.path(
    WORKDIR,
    "sample_metadata.tsv"
  )

OUTPUT_DIR <- WORKDIR

###############################################################################
# Read input tables
###############################################################################

feature_summary <-
  
  read.delim(
    
    FEATURE_SUMMARY,
    
    sep="\t",
    
    stringsAsFactors=FALSE,
    
    check.names=FALSE
    
  )

master_annotations <-
  
  read.delim(
    
    MASTER_ANNOTATIONS,
    
    sep="\t",
    
    stringsAsFactors=FALSE,
    
    check.names=FALSE
    
  )

feature_labels <-
  
  read.delim(
    
    FEATURE_LABELS,
    
    sep="\t",
    
    stringsAsFactors=FALSE,
    
    check.names=FALSE
    
  )

sample_metadata <-
  
  read.delim(
    
    SAMPLE_METADATA,
    
    sep="\t",
    
    stringsAsFactors=FALSE,
    
    check.names=FALSE
    
  )

###############################################################################
# Merge curated labels into master annotations
###############################################################################

master_annotations <-
  
  left_join(
    
    master_annotations,
    
    feature_labels,
    
    by="Feature"
    
  )

###############################################################################
# Merge metadata into feature summary
###############################################################################

feature_summary <-
  
  left_join(
    
    feature_summary,
    
    sample_metadata,
    
    by="Sample"
    
  )

###############################################################################
# Biological colour palette
###############################################################################

category_colours <- c(
  
  Origin        = "#228B22",
  
  Replication   = "#2F6BFF",
  
  Mobility      = "#F39C12",
  
  Maintenance   = "#16A085",
  
  Recombination = "#8E44AD",
  
  Addiction     = "#C2185B",
  
  Resistance    = "#C0392B",
  
  Other         = "#A9A9A9"
  
)

###############################################################################
# Determine plasmids
###############################################################################

plasmids <-
  
  sort(
    
    unique(
      
      feature_summary$Plasmid
      
    )
    
  )

cat("\n")

cat("Detected plasmids:\n\n")

print(plasmids)

cat("\n")

cat(length(plasmids)," plasmids will be plotted.\n")

###############################################################################
# Main plotting loop begins here
###############################################################################

current_plasmid <-
  "SRR23967301_RNODE_242_length_5723_cov_76.02694"
  ###############################################################################
  # Current plasmid
  ###############################################################################
  
  current_features <-
    
    feature_summary %>%
    
    filter(
      
      Plasmid == current_plasmid
      
    )
  
  ###############################################################################
  # Determine plasmid length
  ###############################################################################
  
plasmid_length <-
  
  as.numeric(
    
    sub(
      
      ".*_length_([0-9]+)_cov_.*",
      
      "\\1",
      
      current_plasmid
      
    )
    
  )
  
  ###############################################################################
  # Retrieve reference annotation
  ###############################################################################
  
  reference <-
    
    master_annotations %>%
    
    filter(
      
      Plasmid == current_plasmid
      
    )
  
  ###############################################################################
  # Keep only features represented in feature_summary
  ###############################################################################
  
  reference <-
    
    reference %>%
    
    semi_join(
      
      current_features %>%
        
        distinct(Feature),
      
      by="Feature"
      
    )
  
  ###############################################################################
  # Retrieve display names from feature_summary
  ###############################################################################
  
  display_lookup <-
    
    current_features %>%
    
    distinct(
      
      Feature,
      
      Display_Name,
      
      FeatureType,
      
      Source
      
    )
  
  reference <-
    
    left_join(
      
      reference,
      
      display_lookup,
      
      by="Feature"
      
    )
  
  ###############################################################################
  # Fill missing labels
  ###############################################################################
  
  reference$Display_Name[
    
    is.na(reference$Display_Name)
    
  ] <-
    
    reference$Feature
  
  reference$Category[
    is.na(reference$Category)
  ] <- "Other"
  
  reference$Category[
    !is.na(reference$AMR_gene) &
      reference$AMR_gene != ""
  ] <- "Resistance"
  
  ###############################################################################
  # Plot labels
  ###############################################################################
  
  # Keep curated labels from feature_labels.tsv
  # Fall back to Display_Name only if no curated label exists
  
  reference$PlotLabel <-
    ifelse(
      is.na(reference$PlotLabel) | reference$PlotLabel == "",
      reference$Display_Name,
      reference$PlotLabel
    )
  ###############################################################################
  # Remove labels that should never be displayed
  ###############################################################################
  
  reference$PlotLabel[
    grepl(
      "^SRR.*_length_|^ERR.*_length_|^Coil$|^consensus disorder prediction$|^Hypothetical ORF$",
      reference$PlotLabel,
      ignore.case = TRUE
    )
  ] <- ""
  
  ###############################################################################
  # Plot colours
  ###############################################################################
  
  # Start with biological colours
  reference$Colour <- category_colours[reference$Category]
  
  print(
    reference %>%
      select(
        PlotLabel,
        Category,
        Colour
      )
  )
  
  # Unknown proteins:
  # ORF predicted but no informative annotation
  unknown_idx <-
    reference$PlotLabel == ""
  
  reference$Colour[unknown_idx] <-
    "#9E9E9E"
  
  # Structural-only predictions
  structural_idx <-
    grepl(
      "coil|disorder",
      reference$Display_Name,
      ignore.case = TRUE
    )
  
  reference$Colour[structural_idx] <-
    "#E8E8E8"
  
  ###############################################################################
  # Sort left to right
  ###############################################################################
  
  reference <-
    
    reference %>%
    
    arrange(
      
      Start
      
    )
  
  ###############################################################################
  # gggenes requires these columns
  ###############################################################################
  
  reference <-
    
    reference %>%
    
    mutate(
      
      molecule="Reference",
      
      xmin=Start,
      
      xmax=End,
      
      forward=Strand %in% c(
        
        1,
        
        "+",
        
        "+1"
        
      )
      
    )
  
  ###############################################################################
  # Quick diagnostics
  ###############################################################################
  
  cat("\n")
  
  cat("Reference features:\n")
  
  print(
    
    reference %>%
      
      select(
        
        Feature,
        
        PlotLabel,
        
        Category,
        
        Start,
        
        End,
        
        Strand
        
      )
    
  )
  
  cat("\n")
  ###############################################################################
  # Build plotting dataframe
  ###############################################################################
  
  plot_reference <-
    
    reference %>%
    
    mutate(
      
      ymin = -0.18,
      
      ymax = 0.18
      
    )
  
  ###############################################################################
  # OriV separately
  ###############################################################################
  
  ori_df <-
    
    plot_reference %>%
    
    filter(
      
      FeatureType == "OriV"
      
    )
  wrap_ori <- ori_df %>%
    filter(Start > End)
  
  normal_ori <- ori_df %>%
    filter(Start <= End)
  
  if (nrow(wrap_ori) > 0) {
    
    wrap_part1 <- wrap_ori
    wrap_part1$xmin <- wrap_ori$Start
    wrap_part1$xmax <- plasmid_length
    
    wrap_part2 <- wrap_ori
    wrap_part2$xmin <- 1
    wrap_part2$xmax <- wrap_ori$End
    
    normal_ori$xmin <- normal_ori$Start
    normal_ori$xmax <- normal_ori$End
    
    ori_df <- bind_rows(
      normal_ori,
      wrap_part1,
      wrap_part2
    )
    
  } else {
    
    ori_df$xmin <- ori_df$Start
    ori_df$xmax <- ori_df$End
    
  }
  
  ###############################################################################
  # Genes only
  ###############################################################################
  
  gene_df <-
    
    plot_reference %>%
    
    filter(
      
      FeatureType != "OriV"
      
    )
  hypo <- which(
    
    grepl(
      
      "hypothetical",
      
      gene_df$PlotLabel,
      
      ignore.case = TRUE
      
    )
    
  )
  
  if(length(hypo) > 1){
    
    gene_df$PlotLabel[hypo[-1]] <- ""
    
  }
  
  cat("\nGene labels after suppression:\n")
  print(gene_df[, c("Start", "End", "PlotLabel")])
  
  ###############################################################################
  # Label positions
  ###############################################################################
  
  gene_df$label_y <- 0.82
  
  gene_df$label_y[
    
    seq(2, nrow(gene_df), 2)
    
  ] <- 1.05
  
  ###############################################################################
  # Start plot
  ###############################################################################
  
  p <-
    
    ggplot()
  
  ###############################################################################
  # Backbone
  ###############################################################################
  
  p <-
    
    p +
    
    geom_segment(
      
      aes(
        
        x=1,
        
        xend=plasmid_length,
        
        y=0,
        
        yend=0
        
      ),
      
      linewidth=1.3,
      
      colour="black"
      
    )
  
  ###############################################################################
  # OriV
  ###############################################################################
  
  if(nrow(ori_df)>0){
    
    p <-
      
      p +
      
      geom_rect(
        
        data=ori_df,
        
        aes(
          
          xmin=xmin,
          
          xmax=xmax,
          
          ymin=-0.15,
          
          ymax=0.15
          
        ),
        
        fill=category_colours["Origin"],
        
        colour="black",
        
        linewidth=0.4
        
      )
    
  }
  ###############################################################################
  # Display current plot
  ###############################################################################
  
  print(p)
  ###############################################################################
  # Draw gene arrows
  ###############################################################################
  
  p <-
    
    p +
    
    geom_gene_arrow(
      
      data = gene_df,
      
      aes(
        
        xmin = xmin,
        xmax = xmax,
        y = 0,
        forward = forward,
        fill = Colour
        
      ),
      
      colour = "black",
      linewidth = 0.35,
      arrowhead_height = unit(4,"mm"),
      arrow_body_height = unit(4,"mm"),
      show.legend = FALSE
      
    )
  
  ###############################################################################
  # Gene labels
  ###############################################################################
  
  p <-
    
    p +
    
    geom_segment(
      
      data = subset(gene_df, PlotLabel != ""),
      
      aes(
        
        x = (Start + End)/2,
        
        xend = (Start + End)/2,
        
        y = 0.20,
        
        yend = 0.74
        
      ),
      
      colour = "grey40",
      
      linewidth = 0.25,
      
      linetype = "dotted"
      
    ) +
    
    geom_text(
      
      data = subset(gene_df, PlotLabel != ""),
      
      aes(
        
        x = (Start + End)/2,
        
        y = label_y,
        
        label = PlotLabel
        
      ),
      
      size = 3.3,
      
      fontface = "bold"
      
    )
  
  ###############################################################################
  # Colour palette
  ###############################################################################
  
  p <-
    
    p +
    
    scale_fill_identity()

  
  ###############################################################################
  # Remove unnecessary elements
  ###############################################################################
  
  p <-
    
    p +
    
    theme_classic(
      
      base_size = 13
      
    )
  
  ###############################################################################
  # Final formatting
  ###############################################################################
  
  p <-
    
    p +
    
    theme(
      
      axis.title.x = element_text(
        size = 12,
        face = "bold"
      ),
      
      axis.title.y = element_blank(),
      
      axis.text.x = element_text(
        size = 10
      ),
      
      axis.text.y = element_blank(),
      
      axis.ticks.x = element_line(),
      
      axis.ticks.y = element_blank(),
      
      panel.grid = element_blank(),
      
      legend.position = "none",
      
      plot.title = element_text(
        
        face = "bold",
        
        hjust = 0.5,
        
        size = 16
        
      )
      
    )
  
  ###############################################################################
  # Title
  ###############################################################################
  
  short_name <-
    
    sub(
      
      "_length_.*",
      
      "",
      
      current_plasmid
      
    )
  
  p <-
    
    p +
    
    ggtitle(
      
      paste0(
        
        "Plasmid ",
        
        short_name,
        
        "   ",
        
        plasmid_length,
        
        " bp"
        
      )
      
    )
  
  ###############################################################################
  # Display
  ###############################################################################
  
  print(p)
  ###############################################################################
  # Build sample coverage tracks
  ###############################################################################
  
  track_df <-
    
    current_features %>%
    
    filter(
      
      CoveredBP > 0
      
    )
  
  ###############################################################################
  # Order samples
  ###############################################################################
  
  sample_order <-
    
    track_df %>%
    
    distinct(
      
      Sample,
      
      CoveragePercent
      
    ) %>%
    
    arrange(
      
      desc(CoveragePercent),
      
      Sample
      
    )
  
  ###############################################################################
  # Assign y positions
  ###############################################################################
  
  sample_order$track_y <-
    
    seq(
      
      from = -1.2,
      
      by = -0.75,
      
      length.out = nrow(sample_order)
      
    )
  
  ###############################################################################
  # Add y positions
  ###############################################################################
  
  track_df <-
    
    left_join(
      
      track_df,
      
      sample_order,
      
      by = c(
        
        "Sample",
        
        "CoveragePercent"
        
      )
      
    )
  
  ###############################################################################
  # Expand CoveredSegments
  ###############################################################################
  
  coverage_boxes <-
    
    list()
  
  ###############################################################################
  # Parse every feature
  ###############################################################################
  
  for(i in seq_len(nrow(track_df))){
    
    seg_string <-
      
      track_df$CoveredSegments[i]
    
    if(
      
      is.na(seg_string) ||
      
      seg_string == ""
      
    ){
      
      next
      
    }
    
    segments <-
      
      unlist(
        
        strsplit(
          
          seg_string,
          
          ";"
          
        )
        
      )
    
    for(seg in segments){
      
      coords <-
        
        as.numeric(
          
          unlist(
            
            strsplit(
              
              seg,
              
              "-"
              
            )
            
          )
          
        )
      
      if(length(coords)!=2){
        
        next
        
      }
      
      coverage_boxes[[length(coverage_boxes)+1]] <-
        
        data.frame(
          
          Sample = track_df$Sample[i],
          
          CoveragePercent = track_df$CoveragePercent[i],
          
          xmin = coords[1],
          
          xmax = coords[2],
          
          ymin = track_df$track_y[i]-0.12,
          
          ymax = track_df$track_y[i]+0.12
          
        )
      
    }
    
  }
  
  ###############################################################################
  # Convert to dataframe
  ###############################################################################
  
  coverage_boxes <-
    
    bind_rows(
      
      coverage_boxes
      
    )
  ###############################################################################
  # Draw coverage rectangles
  ###############################################################################
  
  ###############################################################################
  # Sample backbones
  ###############################################################################
  
  p <-
    
    p +
    
    geom_segment(
      
      data = sample_order,
      
      aes(
        
        x = 1,
        
        xend = plasmid_length,
        
        y = track_y,
        
        yend = track_y
        
      ),
      
      colour = "grey80",
      
      linewidth = 0.6
      
    ) +
    
    geom_rect(
      
      data = coverage_boxes,
      
      aes(
        
        xmin = xmin,
        
        xmax = xmax,
        
        ymin = ymin,
        
        ymax = ymax
        
      ),
      
      fill = "#5B8FD1",
      
      colour = "#5B8FD1"
      
    )
  
  ###############################################################################
  # Sample names
  ###############################################################################
  
  p <-
    
    p +
    
    geom_text(
      
      data = sample_order,
      
      aes(
        
        x = -380,
        
        y = track_y,
        
        label = Sample
        
      ),
      
      hjust = 1,
      
      fontface = "bold",
      
      size = 3.4
      
    )
  
  ###############################################################################
  # Coverage percentage
  ###############################################################################
  
  p <-
    
    p +
    
    geom_text(
      
      data = sample_order,
      
      aes(
        
        x = plasmid_length + 250,
        
        y = track_y,
        
        label = paste0(
          
          round(CoveragePercent,1),
          
          "%"
          
        )
        
      ),
      
      hjust = 0,
      
      size = 3.3
      
    )
  ###############################################################################
  # Legend
  ###############################################################################
  
  legend_x <- plasmid_length + 420
  
  ###############################################################################
  # Determine legend contents
  ###############################################################################
  
  module_order <- c(
    
    "Origin",
    "Replication",
    "Maintenance",
    "Mobility",
    "Recombination",
    "Addiction",
    "Resistance"
    
  )
  
  biological_categories <-
    unique(reference$Category)
  
  modules_present <-
    
    module_order[
      module_order %in% biological_categories
    ]
  
  annotation_labels <- c(
    "Unknown function",
    "Structural prediction"
  )
  
  legend_labels <-
    c(
      modules_present,
      annotation_labels
    )
  
  legend_fills <-
    c(
      category_colours[modules_present],
      "#9E9E9E",
      "#E8E8E8"
    )
  
  legend_y <-
    seq(
      from = 1.55,
      by = -0.15,
      length.out = length(modules_present)
    )
  
  annotation_y <-
    c(
      min(legend_y) - 0.30,
      min(legend_y) - 0.45
    )
  
  legend_df <-
    data.frame(
      
      x = legend_x,
      
      y = c(legend_y, annotation_y),
      
      label = legend_labels,
      
      fill = legend_fills,
      
      stringsAsFactors = FALSE
      
    )
  p <-
    
    p +
    
    geom_tile(
      
      data = legend_df,
      
      aes(
        x = x,
        y = y,
        fill = fill
      ),
      
      width = plasmid_length * 0.015,
      height = 0.05,
      
      colour = "black",
      
      inherit.aes = FALSE
      
    )
  p <-
    
    p +
    
    geom_text(
      
      data = legend_df,
      
      aes(
        
        x = x + plasmid_length * 0.03,
        
        y = y,
        
        label = label
        
      ),
      
      hjust = 0,
      
      size = 3.0,
      
      inherit.aes = FALSE
      
    )
  p <-
    
    p +
    
    annotate(
      
      "text",
      
      x = legend_x,
      
      functional_title_y <- max(legend_df$y) + 0.14,
      
      y = functional_title_y,
      
      label = "Functional modules",
      
      hjust = 0,
      
      fontface = "bold",
      
      size = 3.4
      
    ) +
    
    annotate(
      
      "text",
      
      x = legend_x,
      
      annotation_title_y <-
        
        annotation_y[1] + 0.18,
      
      y = annotation_title_y,
      
      label = "Annotation status",
      
      hjust = 0,
      
      fontface = "bold",
      
      size = 3.4
      
    )
  
  ###############################################################################
  # Expand plotting region
  ###############################################################################
  
  p <-
    
    p +
    
    scale_x_continuous(
      
      name = "Position (bp)",
      
      breaks = seq(
        
        0,
        
        ceiling(plasmid_length / 1000) * 1000,
        
        by = 1000
        
      ),
      
      expand = c(0,0)
      
    ) +
    
    coord_cartesian(
      
      xlim = c(
        
        -1100,
        
        plasmid_length + 1700
        
      ),
      
      ylim = c(
        
        min(sample_order$track_y)-0.6,
        
        1.95
        
      ),
      
      clip = "off"
      
    )
  
  ###############################################################################
  # Export figure
  ###############################################################################
  
  output_file <-
    
    file.path(
      
      OUTPUT_DIR,
      
      paste0(
        
        gsub(
          
          "_length_.*",
          
          "",
          
          current_plasmid
          
        ),
        
        ".png"
        
      )
      
    )
  
  ggsave(
    
    filename = output_file,
    
    plot = p,
    
    width = 12,
    
    height = 7,
    
    dpi = 600,
    
    bg = "white"
    
  )
  
  ###############################################################################
  # Display
  ###############################################################################
  
  print(p)