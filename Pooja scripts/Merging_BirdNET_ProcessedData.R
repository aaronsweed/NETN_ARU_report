####loading required libraries
library("librarian")
#Packages name
librarian::shelf(ggplot2, tidyverse, lubridate, tibble, hms, stringr, data.table, purrr, tools)

##-------------------------------------------------------------------------------------------------------------------
## Method 4 v1 : Get detection history for BN_FP model
##-------------------------------------------------------------------------------------------------------------------
manifest_files <- file.path("Data", "Raw_Data", "BN_Data", 
                            c("soundManifest_NPS_2022_09Apr2025.csv", 
                              "soundManifest_NPS_2023_08Apr2025.csv"))

all_rec_data <- rbindlist(lapply(manifest_files, fread))

# Pre-filter and process effort in one go
all_rec_data_clean <- all_rec_data[, `:=`(
  date_time = mdy_hm(startTime.Excel),
  site = as.character(plot)
)][, `:=`(
  year = year(date_time),
  month = month(date_time),
  week_num = week(date_time),
  hour = hour(date_time),
  date = as.Date(date_time)
)]

june_weeks <- unique(all_rec_data_clean[month == 6, week_num])

daily_effort <- all_rec_data_clean[
  week_num %in% june_weeks & hour >= 5 & hour < 8 & area != "SAGA",
  .(total_minutes = sum(fileLength.min, na.rm = TRUE)),
  by = .(area, site, year, date, week_num)
][, day_of_year := yday(date)]

# --- 2. Load BirdNET Runs (The Bottleneck) ---
# We use fread and filter IMMEDIATELY during load to save memory
#spp_Th90_valid <- fread(file.path("Data", "Raw_Data", "BN_Data", "BN_Threshold_90.csv"))[Include == "1"]
#species_thresholds <- setNames(spp_Th90_valid$mean_Th, spp_Th90_valid$Spp6)
#species4BNOcc <- names(species_thresholds)

all_files <- file.path("Data","Raw_Data", "BN_Data", "All_BirdNET_Run", "Merged", 
                       c("ACAD_2022_BN_All.csv", "HBEF_2022_BN_All.csv", "KAWW_2022_BN_All.csv", "MABI_2022_BN_All.csv",
                         "ACAD_2023_BN_All.csv", "HBEF_2023_BN_All.csv", "KAWW_2023_BN_All.csv", "MABI_2023_BN_All.csv"))

combined_final_data <- rbindlist(lapply(all_files, function(f) {
  fread(f)#[Species.Code %in% species4BNOcc] # Filter species WHILE reading
}))

# Process DateTime after combining
combined_final_data[, `:=`(
  DateTime = as.POSIXct(DateTime, format = "%Y-%m-%d %H:%M:%S"),
  date = as.Date(DateTime),
  hour = hour(DateTime),
  week_num = week(DateTime)
)]

# --- 3. Optimized Detection Summary ---
# Convert to data.table for high-speed grouping
setDT(combined_final_data)
BNFP_detection_summary <- combined_final_data[
  week_num %in% june_weeks & hour >= 5 & hour < 8,
  .(
    n_detections_all = .N,
    # Summing logicals is faster than sum(Confidence >= th)
    #n_detections_above_th = sum(Confidence >= species_thresholds[Species.Code], na.rm = TRUE),
    max_confidence = as.numeric(suppressWarnings(
      max(Confidence, na.rm = TRUE) #[Confidence >= species_thresholds[Species.Code]]
    ))
  ),
  by = .(SiteName, date, Species.Code)
]
# Clean up -Inf from max()
BNFP_detection_summary[is.infinite(max_confidence), max_confidence := 0]

BNFP_Spp <- unique(BNFP_detection_summary$Species.Code)

# --- 4. The Optimized Grid (Sparse Join) ---
# Instead of crossing(), we only keep records where effort exists
setnames(BNFP_detection_summary, c("SiteName", "Species.Code"), c("site", "species"))

BNFP_df <- as.data.table(crossing(daily_effort, species = BNFP_Spp))

# Perform a rolling/keyed join (very fast)
setkey(BNFP_df, site, date, species)
setkey(BNFP_detection_summary, site, date, species)

BNFP_df <- BNFP_detection_summary[BNFP_df]

# Fill NAs in bulk
cols_to_fix <- c("n_detections_all", "max_confidence")
for (j in cols_to_fix) set(BNFP_df, i = which(is.na(BNFP_df[[j]])), j = j, v = 0)
BNFP_df[, site_year_id := paste0(site, "_", year)]

# --- Save Results ---
fwrite(BNFP_df, file.path("Data", "Clean_Data", "BNFP_Combined_maxcf_Allspp.csv"))