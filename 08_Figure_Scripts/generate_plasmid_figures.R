###############################################################################
# Reproducible Thesis Figures – Plasmidome Project
# Author: Oded Flicker
# Description: Regenerates 4 thesis-quality figures programmatically
###############################################################################

# ==============================
# 0. Load libraries
# ==============================

library(tidyverse)
library(readxl)
library(scales)

# ==============================
# 1. Set Working Directory
# ==============================

setwd("C:/Users/user/Desktop/Research Project Oded/Figures and their sources")

# ==============================
# 2. Global Theme (Thesis Style)
# ==============================

theme_thesis <- theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank()
  )

###############################################################################
# FIGURE 1 – Plasmid Length Histogram
###############################################################################

plasmid_lengths <- read_excel(
  "input_for_plasmid_length_histo.xlsx",
  col_names = FALSE
)

colnames(plasmid_lengths) <- c("Plasmid", "Length_bp", "Log10_Length")

fig1 <- ggplot(plasmid_lengths, aes(x = Log10_Length)) +
  geom_histogram(
    binwidth = 0.1,
    boundary = 3.0,
    fill = "#1f4e79",
    color = "white"
  ) +
  labs(
    title = "Length Histogram",
    x = "Length(logarithmic scale)",
    y = "Count"
  ) +
  theme_thesis

ggsave("Figure1_Length_Histogram.tiff", fig1, dpi = 300, width = 8, height = 6)
ggsave("Figure1_Length_Histogram.pdf", fig1, width = 8, height = 6)

###############################################################################
# FIGURE 2 – ORF Count per Plasmid Candidate (Truncated Display)
###############################################################################

orf_counts <- read_excel(
  "predictions_fixed.xlsx",
  sheet = "ORF_Counts"
)

# Count frequency of each ORF count value
orf_counts_summary <- orf_counts %>%
  count(ORF_Count)

# Define truncation limit (adjust if needed)
x_limit <- 120

fig2 <- ggplot(orf_counts_summary, aes(x = ORF_Count, y = n)) +
  geom_col(fill = "#3b82c4") +
  scale_y_log10() +
  coord_cartesian(xlim = c(0, x_limit)) +  # <<< Truncation happens HERE
  labs(
    title = "ORF Count per Plasmid Candidate",
    x = "Number of ORFs",
    y = "Number of candidates(logarithmic scale)"
  ) +
  theme_thesis

ggsave("Figure2_ORF_Count_per_Plasmid_truncated.tiff",
       fig2, dpi = 300, width = 8, height = 6)

ggsave("Figure2_ORF_Count_per_Plasmid_truncated.pdf",
       fig2, width = 8, height = 6)

###############################################################################
# FIGURE 3 – ORF Length Distribution (Truncated at 16,000 bp)
###############################################################################

orf_lengths <- read_excel(
  "predictions_fixed.xlsx",
  sheet = "ORF_Lengths"
)

orf_lengths_truncated <- orf_lengths %>%
  filter(ORF_Length <= 16000)

bins <- seq(0, 16000, by = 1000)

fig3 <- ggplot(orf_lengths_truncated, aes(x = ORF_Length)) +
  geom_histogram(
    breaks = bins,
    fill = "#9db8d6",
    color = "black"
  ) +
  scale_y_log10() +
  labs(
    title = "Distribution of ORF Lengths in Plasmid Candidates",
    x = "ORF LENGTH BINS (BP)",
    y = "Number of ORFs(logarithmic scale)"
  ) +
  theme_thesis

ggsave("Figure3_ORF_Length_Distribution.tiff", fig3, dpi = 300, width = 9, height = 6)
ggsave("Figure3_ORF_Length_Distribution.pdf", fig3, width = 9, height = 6)

###############################################################################
# FIGURE 4 – Total ORF Length / Plasmid Length Ratio
###############################################################################

ratio_data <- read_excel(
  "predictions_fixed.xlsx",
  sheet = "ORF_to_Plasmid_Ratio"
)

fig4 <- ggplot(ratio_data, aes(x = ORF_to_Plasmid_Ratio)) +
  geom_histogram(
    bins = 14,
    fill = "#4f81bd",
    color = "white"
  ) +
  scale_y_log10() +
  labs(
    title = "Total ORF Length / Plasmid Length Ratio",
    x = "ORF-to-Plasmid Length Ratio",
    y = "Number of Plasmid Candidates(logarithmic scale)"
  ) +
  theme_thesis

ggsave("Figure4_ORF_to_Plasmid_Ratio.tiff", fig4, dpi = 300, width = 8, height = 6)
ggsave("Figure4_ORF_to_Plasmid_Ratio.pdf", fig4, width = 8, height = 6)

###############################################################################
# STATISTICAL SUMMARIES – Clean, Separated by Figure
###############################################################################

# Helper function to save stats cleanly
save_stats <- function(data, filename) {
  write.csv(data, filename, row.names = FALSE)
  print(data)
  cat(paste0("\nSaved: ", filename, "\n\n"))
}

# -----------------------------
# FIGURE 1 – Plasmid Length Stats
# -----------------------------

stats_length <- tibble(
  Metric = c("N",
             "Min_bp",
             "Max_bp",
             "Mean_bp",
             "Median_bp",
             "SD_bp",
             "IQR_bp",
             "P5_bp",
             "P95_bp"),
  Value = c(
    nrow(plasmid_lengths),
    min(plasmid_lengths$Length_bp),
    max(plasmid_lengths$Length_bp),
    mean(plasmid_lengths$Length_bp),
    median(plasmid_lengths$Length_bp),
    sd(plasmid_lengths$Length_bp),
    IQR(plasmid_lengths$Length_bp),
    quantile(plasmid_lengths$Length_bp, 0.05),
    quantile(plasmid_lengths$Length_bp, 0.95)
  )
)

save_stats(stats_length, "Figure1_Length_Statistics.csv")


# -----------------------------
# FIGURE 2 – ORF Count Stats
# -----------------------------

stats_orf_count <- tibble(
  Metric = c("N",
             "Min_ORF",
             "Max_ORF",
             "Mean_ORF",
             "Median_ORF",
             "SD_ORF",
             "P95_ORF"),
  Value = c(
    nrow(orf_counts),
    min(orf_counts$ORF_Count),
    max(orf_counts$ORF_Count),
    mean(orf_counts$ORF_Count),
    median(orf_counts$ORF_Count),
    sd(orf_counts$ORF_Count),
    quantile(orf_counts$ORF_Count, 0.95)
  )
)

save_stats(stats_orf_count, "Figure2_ORF_Count_Statistics.csv")


# -----------------------------
# FIGURE 3 – ORF Length Stats
# -----------------------------

stats_orf_length <- tibble(
  Metric = c("Total_ORFs",
             "Min_ORF_bp",
             "Max_ORF_bp",
             "Mean_ORF_bp",
             "Median_ORF_bp",
             "SD_ORF_bp",
             "P5_ORF_bp",
             "P95_ORF_bp"),
  Value = c(
    nrow(orf_lengths),
    min(orf_lengths$ORF_Length),
    max(orf_lengths$ORF_Length),
    mean(orf_lengths$ORF_Length),
    median(orf_lengths$ORF_Length),
    sd(orf_lengths$ORF_Length),
    quantile(orf_lengths$ORF_Length, 0.05),
    quantile(orf_lengths$ORF_Length, 0.95)
  )
)

save_stats(stats_orf_length, "Figure3_ORF_Length_Statistics.csv")


# -----------------------------
# FIGURE 4 – Ratio Stats
# -----------------------------

stats_ratio <- tibble(
  Metric = c("N",
             "Min_ratio",
             "Max_ratio",
             "Mean_ratio",
             "Median_ratio",
             "SD_ratio",
             "P95_ratio",
             "Percent_over_1"),
  Value = c(
    nrow(ratio_data),
    min(ratio_data$ORF_to_Plasmid_Ratio),
    max(ratio_data$ORF_to_Plasmid_Ratio),
    mean(ratio_data$ORF_to_Plasmid_Ratio),
    median(ratio_data$ORF_to_Plasmid_Ratio),
    sd(ratio_data$ORF_to_Plasmid_Ratio),
    quantile(ratio_data$ORF_to_Plasmid_Ratio, 0.95),
    mean(ratio_data$ORF_to_Plasmid_Ratio > 1) * 100
  )
)

save_stats(stats_ratio, "Figure4_Ratio_Statistics.csv")

###############################################################################
# END
###############################################################################

cat("All figures generated successfully.\n")
