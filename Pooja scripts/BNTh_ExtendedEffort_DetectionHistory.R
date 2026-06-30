# ====================================================================
# BNTh EXTENDED SAMPLING EFFORT: DATA PREPARATION
# ====================================================================
library(tidyverse)
library(lubridate)

# ---- 1. File Path Configurations ----
dir_extended    <- file.path("Data", "Raw_Data", "BN_Data", "Extended_Rec_SelectionFiles_All") 
master_sp_file  <- file.path("Data", "Raw_Data", "Species_Metadata", "species_metadata_masterfile.csv")
output_ext_csv  <- file.path("Data", "Clean_Data", "BNTh_extendedSampling.csv")

# List all hourly BirdNET selection files across all subfolders recursively
extBN_files <- list.files(path = dir_extended, pattern = "*BirdNET.selection.table.txt", full.names = TRUE, recursive = TRUE)

# ---- 2. Parse Raw Hourly BirdNET Selection Tables ----
parse_extended_file <- function(file_path) {
  temp_df <- read.table(file_path, header = TRUE, sep = "\t", check.names = FALSE)
  if (nrow(temp_df) == 0) return(NULL)
  
  # Example 1: ACAD3005_20220531_050000.BirdNET.selection.table.txt
  # Example 2: HBEFV308_20220531_050000.BirdNET.selection.table.txt
  base_name <- tools::file_path_sans_ext(basename(file_path))
  clean_name <- sub(".BirdNET.selection.table", "", base_name)
  
  # Extract string elements safely
  parts <- str_split(clean_name, "_")[[1]]
  raw_site <- parts[1]                           # e.g., "HBEFV007" or "ACAD3005"
  raw_date <- parts[2]                           # e.g., "20220531"
  raw_time <- parts[3]                           # e.g., "050000"
  
  # Strip only the specific project location prefixes
  site_num <- str_remove(raw_site, "^(ACAD|HBEF|MABI|KAWW)") # "HBEFV007" -> "V007"
  
  year_str <- substr(raw_date, 1, 4)             # Extract "2022" or "2023"
  site_year_id <- paste0(site_num, "_", year_str) # Results in "V007_2022" (or "3005_2022")
  
  temp_df %>%
    mutate(
      fileName     = clean_name,
      site_year_id = site_year_id,
      EventDate    = raw_date,
      Hour         = raw_time,
      Confidence   = as.numeric(Confidence)
    ) %>%
    select(fileName, site_year_id, EventDate, Hour, `Species Code`, Confidence)
}

# Combine files dynamically
if(length(extBN_files) > 0) {
  raw_extBN_combined <- map_df(extBN_files, parse_extended_file, .progress = TRUE)
} else {
  stop("Critical Error: No raw selection text files located in your target directory.")
}

# ---- 3. Establish Core Sampling Calendar Framework ----
# Extract unique dates and assign visit numbers sequentially
target_dates <- sort(unique(raw_extBN_combined$EventDate))

# Create a robust mapping dataframe that handles multiple years seamlessly
date_visit_map <- raw_extBN_combined %>%
  distinct(EventDate) %>%
  mutate(Year = substr(EventDate, 1, 4)) %>%
  group_by(Year) %>%
  arrange(EventDate) %>%
  mutate(visit_no = paste0("V", row_number())) %>%
  ungroup() %>%
  select(EventDate, visit_no)

# ---- 4. Apply Logit Conversion and Threshold Filtering ----
master_spp_meta <- read_csv(master_sp_file, col_types = cols(.default = "c")) %>%
  mutate(mean_lTh90_Uniform200 = as.numeric(mean_lTh90_Uniform200))

extBN_filtered <- raw_extBN_combined %>%
  filter(EventDate %in% target_dates) %>%
  left_join(date_visit_map, by = "EventDate") %>%
  left_join(
    select(master_spp_meta, Spp6, Spp4, mean_lTh90_Uniform200, ForOccModel),
    by = c("Species Code" = "Spp6")
  ) %>%
  # Drop data if species lacks optimization configurations or is rejected
  filter(!is.na(mean_lTh90_Uniform200), str_to_lower(ForOccModel) == "yes") %>%
  # Apply logit probability transformation formula
  mutate(
    Confidence_adj = case_when(
      Confidence >= 1.0 ~ 1.0 - 1e-5,
      Confidence <= 0.0 ~ 0.0 + 1e-5,
      TRUE ~ Confidence
    ),
    logit_score = log(Confidence_adj / (1 - Confidence_adj))
  ) %>%
  filter(logit_score >= mean_lTh90_Uniform200)

# ---- 5. Merge Hourly Data into Day-Level Detections ----
BNTh_extendedSampling <- extBN_filtered %>%
  group_by(site_year_id, EventDate, visit_no, Spp4) %>%
  summarize(
    voc_count = n(), 
    method    = "BNTh_extendedsamplingeffort",
    .groups   = 'drop'
  ) %>%
  mutate(
    Survey_ID = paste0(site_year_id, "_", EventDate, "_", visit_no)
  ) %>%
  rename(AOU_Code = Spp4, count = voc_count)

# Save clean detection history to your database
write_csv(BNTh_extendedSampling, output_ext_csv)

# ---- 6. Generate an Explicit Effort Tracking Matrix ----
effort_grid <- raw_extBN_combined %>%
  filter(EventDate %in% target_dates) %>%
  distinct(site_year_id, EventDate) %>%
  left_join(date_visit_map, by = "EventDate")

save(effort_grid, date_visit_map, file = file.path("Data","Clean_Data","bnth_extended_effort_grid.RData"))
