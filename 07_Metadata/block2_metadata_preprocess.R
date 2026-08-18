setwd("C:/Users/user/Desktop/Research Project Oded/metadata")
install.packages(c("readxl", "dplyr", "tidyr", "lubridate"))

library(readxl)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)

###############################################################################
# ============================  PART A — ALEXANDER  ============================
###############################################################################

alex_raw <- read_excel("Alexander_est_samples_metadata.xlsx")

# ---------------------------------------------------------------------------
# A.1 — Parse basic fields
# Extract sample_id as run_accession
alex <- alex_raw %>%
  rename(sample_id = run_accession) %>%
  mutate(
    # Ensure date is parsed (e.g., "Feb-22")
    date = my(date),
    
    # Extract numeric station (1–9)
    station_numeric = as.integer(str_extract(Station, "(?<=station\\s)[0-9]+")),
    
    # Raw label after number
    raw_label = str_trim(str_replace(Station, "^station\\s*[0-9]+\\s*", "")),
    
    # Fine-grained station name (Option A — Final)
    station_name = case_when(
      station_numeric == 1 ~ "Tsur Natan",
      station_numeric == 2 ~ "Tira",
      station_numeric == 3 ~ "Avrech",
      station_numeric == 4 ~ "Qalansuwa",
      station_numeric == 5 ~ "Viker",
      station_numeric == 6 ~ "Nablus River",
      station_numeric == 7 ~ "Avihayl",
      station_numeric == 8 ~ "Estuary_50m",
      station_numeric == 9 ~ "Estuary_Shoreline",
      TRUE ~ raw_label
    ),
    
    # Broad ecological region (used for grouping)
    region = case_when(
      station_numeric %in% c(1, 2) ~ "Upstream",
      station_numeric %in% c(3, 4, 5) ~ "Midstream",
      station_numeric == 6 ~ "Urban_Tributary",
      station_numeric == 7 ~ "Lower_River",
      station_numeric %in% c(8, 9) ~ "Estuary",
      TRUE ~ "Unknown"
    )
  )

# ---------------------------------------------------------------------------
# A.2 — Standardize seasons for the three sampling campaigns
alex <- alex %>%
  mutate(
    season = case_when(
      month(date) %in% c(2) ~ "Winter",
      month(date) %in% c(4) ~ "Spring",
      month(date) %in% c(9) ~ "LateSummer",
      TRUE ~ "Unknown"
    )
  )

# ---------------------------------------------------------------------------
# A.3 — Final clean Alexander metadata
alex_clean <- alex %>%
  select(
    sample_id, sample_title, station_numeric, station_name, region,
    date, season, read_count, Temperature
  )

# =============================================================================
# A.4 — Export Alexander metadata
# =============================================================================

# CSV version
write.csv(alex_clean, "alex_clean_metadata.csv", row.names = FALSE)

# Excel version
library(writexl)
write_xlsx(alex_clean, "alex_clean_metadata.xlsx")

cat("\nAlexander metadata exported as:\n - alex_clean_metadata.csv\n - alex_clean_metadata.xlsx\n")


###############################################################################
# ==============================  PART B — BOGOTA  =============================
###############################################################################

library(readxl)
library(dplyr)
library(stringr)

# ---------------------  B.0 — Read Bogotá Excel Sheets  ----------------------

bogo_raw <- read_excel("Bogoata_river_samples_metadata.xlsx")
antib_raw <- read_excel("Bogoata_river_antibiotics.xlsx")
phys_raw  <- read_excel("Bogoata_river_physicochemical.xlsx")


# =============================================================================
# B.1 — CLEAN SAMPLE METADATA TABLE
# =============================================================================

bogo <- bogo_raw %>%
  rename(sample_id = run_accession) %>%
  mutate(
    # Convert date ranges to year:
    date = case_when(
      grepl("Aug", date)  ~ "2015",
      grepl("Nov", date)  ~ "2015",
      grepl("May", date)  ~ "2016",
      grepl("June", date) ~ "2016",
      TRUE ~ NA_character_
    ),
    
    sampling_year = as.integer(date),
    
    # Station code stays as-is (RS15 / RM15 / RL15 / HA15 / HB15 / HC15)
    station_name = Station,
    
    # Assign ecological region
    region = case_when(
      grepl("^RS", Station) ~ "River Source",
      grepl("^RM", Station) ~ "River Middle",
      grepl("^RL", Station) ~ "River Lower",
      grepl("^HA", Station) ~ "Hospital A",
      grepl("^HB", Station) ~ "Hospital B",
      grepl("^HC", Station) ~ "Hospital C",
      TRUE ~ "Unknown"
    ),
    
    # Only samples from 2015 have antibiotics & physicochemical metadata
    include_env = (sampling_year == 2015)
  )


# =============================================================================
# B.2 — PROCESS ANTIBIOTIC DATA (CORRECTED IMPORT FOR RS15)
# =============================================================================

# Load WITH skipping first header row
antib_raw2 <- read_excel(
  "Bogoata_river_antibiotics.xlsx",
  skip = 1   # <--- THIS FIXES RS15
)

# Now rename properly based on your exact structure:
antib <- antib_raw2 %>%
  rename(
    Station = 1,
    MNZ = 2,
    SMX = 3,
    TMP = 4,
    NOR = 5,
    CIP = 6,
    CLI = 7,
    CLR = 8,
    ERY = 9,
    AZM = 10
  ) %>%
  filter(!is.na(Station)) %>%   # RS15 now passes
  mutate(
    across(MNZ:AZM, ~ ifelse(. == "ND", 0, as.numeric(.)))
  ) %>%
  mutate(
    station_group = str_extract(Station, "^[A-Za-z]+[0-9]+")
  ) %>%
  group_by(station_group) %>%
  summarise(across(MNZ:AZM, ~ mean(.x, na.rm = TRUE))) %>%
  ungroup()


# =============================================================================
# B.3 — PROCESS PHYSICOCHEMICAL DATA (AVERAGE DUPLICATES)
# =============================================================================

phys_clean <- phys_raw %>%
  rename(
    Station = `Triplicate samples used for physicochemical analyses`,
    BOD     = `BOD (mgO2/L)`,
    COD     = `COD (mgO2/L)`,
    TSS     = `TSS (mg/L)`,
    TN      = `TN (%)`
  ) %>%
  filter(grepl("^(RS|RM|RL|HA|HB|HC)", Station)) %>%   # keep only real samples
  mutate(
    across(c(BOD, COD, TSS, TN), as.numeric)
  )

phys <- phys_clean %>%
  mutate(
    station_group = str_extract(Station, "^[A-Za-z]+[0-9]+")   # HB15-1, HB15-2 → HB15
  ) %>%
  group_by(station_group) %>%
  summarise(across(c(BOD, COD, TSS, TN), ~ mean(.x, na.rm = TRUE))) %>%
  ungroup()


# =============================================================================
# B.4 — MERGE EVERYTHING INTO FINAL BOGOTA METADATA TABLE
# =============================================================================

bogo_clean <- bogo %>%
  mutate(station_group = str_extract(Station, "^[A-Za-z]+[0-9]+")) %>%
  left_join(antib, by = "station_group") %>%
  left_join(phys, by = "station_group")


# =============================================================================
# B.5 — Inspect final output
# =============================================================================

cat("\n==== Bogotá Metadata (bogo_clean) Preview ====\n")
print(head(bogo_clean, 10))

cat("\nColumns available in bogo_clean:\n")
print(names(bogo_clean))

# =============================================================================
# B.6 — Export final Bogotá metadata tables
# =============================================================================

# CSV version
write.csv(bogo_clean, "bogo_clean_metadata.csv", row.names = FALSE)

# Excel version (requires writexl)
install.packages("writexl")
library(writexl)

write_xlsx(bogo_clean, "bogo_clean_metadata.xlsx")

cat("\nBogotá metadata exported as:\n - bogo_clean_metadata.csv\n - bogo_clean_metadata.xlsx\n")

