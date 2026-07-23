# =============================================================================
# Script name: clean_input_csv.R
# Author: Estephe Kana
# Date created: 2026-07-01
# Purpose: Drop unused columns from the two raw input CSV files.
#          Only columns referenced in make_db_object.R, codePredomics_prev.R,
#          scalenet_network_inference.R, save_table_data.R, and figure scripts
#          are retained.
# Inputs:  data/eDNA_SEAMOUNTS_REEF3.0_merged_Environmental_Variables_raw.csv
#          data/eDNA_Data_SEAMOUNTS_REEF3.0_merged.csv
# Outputs: data/eDNA_SEAMOUNTS_REEF3.0_merged_Environmental_Variables_final.csv
#          data/eDNA_Data_SEAMOUNTS_REEF3.0_merged_final.csv
# Run:     Rscript analyses/scripts/clean_input_csv.R
# =============================================================================

library(readr)

# Derive paths from script location
script_path <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
script_dir  <- dirname(script_path)
repo_root   <- dirname(dirname(script_dir))
data_dir    <- file.path(repo_root, "data")

# Habitat value renaming
habitat_recode <- c(
  "Soft_back_reef"   = "BackReef",
  "Reef_outer_slope" = "OuterSlope",
  "Summit50"         = "Seamount50",
  "DeepSlope"        = "DeepSlope150",
  "Summit250"        = "Seamount250",
  "Summit500"        = "Seamount500"
)
recode_habitat <- function(x) {
  ifelse(x %in% names(habitat_recode), unname(habitat_recode[x]), x)
}

# ================================================================================
# 1. eDNA_SEAMOUNTS_REEF3.0_merged_Environmental_Variables_raw.csv pre-processing
# ================================================================================

env_cols_used <- c(
  "Station",                
  "Latitude",                
  "Longitude",               
  "Site",                    
  "Habitat",                 
  "Depth",                  
  "Salinity",         
  "SSTmean",
  "seafloorTemp",
  "Chla",
  "TravelTime",
  "ReefMinDist.m"
)

env_path <- file.path(data_dir, "eDNA_SEAMOUNTS_REEF3.0_merged_Environmental_Variables_raw.csv")
env_raw  <- as.data.frame(suppressMessages(read_csv(env_path)))

env_missing <- setdiff(env_cols_used, colnames(env_raw))
if (length(env_missing) > 0)
  warning("Columns listed as used but not found in env file: ",
          paste(env_missing, collapse = ", "))

env_dropped <- setdiff(colnames(env_raw), env_cols_used)
env_clean   <- env_raw[, intersect(env_cols_used, colnames(env_raw))]
env_clean$Habitat <- recode_habitat(env_clean$Habitat)

message("\n--- Environmental variables file ---")
message("Original : ", ncol(env_raw),   " columns, ", nrow(env_raw),   " rows")
message("Cleaned  : ", ncol(env_clean), " columns, ", nrow(env_clean), " rows")
message("Dropped  : ", paste(env_dropped, collapse = ", "))

out_env <- file.path(data_dir, "eDNA_SEAMOUNTS_REEF3.0_merged_Environmental_Variables_final_V2.csv")
write_csv(env_clean, out_env)
message("Saved    : ", out_env)

# =============================================================================
# 2. eDNA_Data_SEAMOUNTS_REEF3.0_merged.csv pre-processing
# =============================================================================

edna_cols_used <- c(
  "Station",
  "Site",
  "Habitat",
  "Depth",
  "Latitude",
  "Longitude",
  "sequence",
  "new_scientific_name_ncbi", 
  "mean_pcr_count_reads",
  "SpeciesFB"
)

edna_path <- file.path(data_dir, "eDNA_Data_SEAMOUNTS_REEF3.0_merged.csv")
edna_raw  <- as.data.frame(suppressMessages(read_csv(edna_path)))

edna_missing <- setdiff(edna_cols_used, colnames(edna_raw))
if (length(edna_missing) > 0)
  warning("Columns listed as used but not found in eDNA file: ",
          paste(edna_missing, collapse = ", "))

edna_dropped <- setdiff(colnames(edna_raw), edna_cols_used)
edna_clean   <- edna_raw[, intersect(edna_cols_used, colnames(edna_raw))]
edna_clean$Habitat <- recode_habitat(edna_clean$Habitat)
edna_clean$new_scientific_name_ncbi <- gsub(" ", ".", edna_clean$new_scientific_name_ncbi)

# Merge with clean env to keep station order
cols_to_drop <- intersect(c("Site", "Habitat", "Depth", "Latitude", "Longitude"),
                          colnames(env_clean))
edna_clean <- merge(x = edna_clean,
                    y = env_clean[, !colnames(env_clean) %in% cols_to_drop],
                    by = "Station", all.x = TRUE)

# Derive best_taxonomic_assignment  from (sequence, new_scientific_name_ncbi, SpeciesFB), assign a sequential MOTU label
edna_annot <- unique(edna_clean[, c("sequence", "new_scientific_name_ncbi", "SpeciesFB")])

valid <- !is.na(edna_annot$new_scientific_name_ncbi)
edna_annot <- edna_annot[valid,]
edna_annot$best_taxonomic_assignment <- paste0(
  "MOTU", sprintf("%03d", 1:nrow(edna_annot)), "_",
  edna_annot$new_scientific_name_ncbi
)
edna_clean <- merge(edna_clean, edna_annot[, c("sequence", "best_taxonomic_assignment")],
                    by = "sequence", all.x = TRUE)

# Drop columns used only temporarily
env_cols_added <- setdiff(colnames(env_clean),
                          c("Station", "Site", "Habitat", "Depth", "Latitude", "Longitude"))
cols_to_remove <- c("sequence", "new_scientific_name_ncbi", "SpeciesFB", env_cols_added)
edna_clean <- edna_clean[, !colnames(edna_clean) %in% cols_to_remove]

message("\n--- eDNA data file ---")
message("Original : ", ncol(edna_raw),   " columns, ", nrow(edna_raw),   " rows")
message("Cleaned  : ", ncol(edna_clean), " columns (incl. best_taxonomic_assignment), ",
        nrow(edna_clean), " rows")
message("Dropped  : ", paste(edna_dropped, collapse = ", "))

out_edna <- file.path(data_dir, "eDNA_Data_SEAMOUNTS_REEF3.0_merged_final_V2.csv")
write_csv(edna_clean, out_edna)
message("Saved    : ", out_edna)
