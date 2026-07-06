# =============================================================================
# Script name: clean_input_csv.R
# Author: Estephe Kana
# Date created: 2026-07-01
# Purpose: Drop unused columns from the two raw input CSV files.
#          Only columns referenced in make_db_object.R, codePredomics_prev.R,
#          scalenet_network_inference.R, save_table_data.R, and figure scripts
#          are retained. Cleaned files are saved with a _clean suffix.
#          Order mirrors make_db_object.R: env is cleaned first so it is ready
#          for the merge that precedes OTU annotation in the main pipeline.
# Inputs:  data/eDNA_SEAMOUNTS_REEF3.0_merged_Environmental_Variables_raw.csv
#          data/eDNA_Data_SEAMOUNTS_REEF3.0_merged.csv
# Outputs: data/eDNA_SEAMOUNTS_REEF3.0_merged_Environmental_Variables_clean.csv
#          data/eDNA_Data_SEAMOUNTS_REEF3.0_merged_clean.csv
# Run:     Rscript analyses/scripts/clean_input_csv.R
# =============================================================================

library(readr)

# Habitat value renaming applied to both cleaned files (Bay and Lagoon are
# left unchanged).
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

# Derive paths from script location
script_path <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
script_dir  <- dirname(script_path)
repo_root   <- dirname(dirname(script_dir))
data_dir    <- file.path(repo_root, "data")

# =============================================================================
# 1. eDNA_SEAMOUNTS_REEF3.0_merged_Environmental_Variables_raw.csv
#    Cleaned first so it is available for the merge in section 2, mirroring
#    the order in make_db_object.R (env loaded before db merge).
# =============================================================================

env_cols_used <- c(
  "Station",                   # merge key
  "Latitude",                  # geographic coordinate (dropped at merge, kept for reference)
  "Longitude",                 # geographic coordinate (dropped at merge, kept for reference)
  "Site",                      # site name (dropped at merge, kept for reference)
  "Habitat",                   # habitat (dropped at merge, kept for reference)
  "Depth",                     # depth (dropped at merge, kept for reference)
  "EastwardVelocity",          # current velocity
  "NorthwardVelocity",         # current velocity
  "Salinity",                  # salinity (PCoA biplot)
  "SuspendedParticulateMatter",
  # "SSTmax",                  # sea-surface temperature max
  "SSTmean",                   # sea-surface temperature mean (PCoA biplot)
  # "SSTmin",                  # sea-surface temperature min
  # "SSTsd",                   # sea-surface temperature sd
  "seafloorTemp",              # seafloor temperature (PCoA biplot)
  "Chla",                      # chlorophyll a (PCoA biplot)
  "TravelTime",                # travel time (PCoA biplot)
  "ReefMinDist.m"              # reef minimum distance (PCoA biplot)
  # LandMinDist.m intentionally excluded: present in raw file but never used
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

out_env <- file.path(data_dir, "eDNA_SEAMOUNTS_REEF3.0_merged_Environmental_Variables_clean.csv")
write_csv(env_clean, out_env)
message("Saved    : ", out_env)

# =============================================================================
# 2. eDNA_Data_SEAMOUNTS_REEF3.0_merged.csv
#    Loaded after env so the merge with env_clean can happen before the
#    best_taxonomic_assignment annotation, exactly as in make_db_object.R.
# =============================================================================

edna_cols_used <- c(
  "Spygen",                   # sample identifier
  "Site",                     # site name
  "Habitat",                  # habitat category
  # "Substrate",              # substrate type
  # "Stratum",                # depth stratum
  "Station",                  # station ID (merge key with env file)
  # "Date",                   # collection date
  # "Time",                   # collection time
  "Depth",                    # depth
  # "Sampling_depth",         # actual sampling depth
  "Latitude",                 # geographic coordinate
  "Longitude",                # geographic coordinate
  "sequence",                 # DNA sequence (used to build OTU IDs)
  "new_scientific_name_ncbi", # taxonomic name (used to build OTU labels)
  # "Project",                # source project
  # "Method",                 # sampling method
  "mean_pcr_count_reads",     # abundance value pivoted into the OTU matrix
  "SpeciesFB"                 # FishBase species name
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

# Merge with clean env to bring in environmental variables, dropping columns
# already present in edna_clean (mirrors make_db_object.R line 41).
cols_to_drop <- intersect(c("Site", "Habitat", "Depth", "Latitude", "Longitude"),
                          colnames(env_clean))
edna_clean <- merge(x = edna_clean,
                    y = env_clean[, !colnames(env_clean) %in% cols_to_drop],
                    by = "Station", all.x = TRUE)

# Derive best_taxonomic_assignment using the same logic as the id column in
# make_db_object.R: build a unique annotation table from (sequence,
# new_scientific_name_ncbi, SpeciesFB), assign a sequential OTU label, then
# merge back by sequence so every row gets its label.
edna_annot <- unique(edna_clean[, c("sequence", "new_scientific_name_ncbi", "SpeciesFB")])
# Number only rows with a valid species name, mirroring make_db_object.R
# which filters !is.na(sequence) before building the annotation table.
# NA-named rows get NA so no slot is consumed, and the OTU numeric prefix
# matches the original pipeline exactly.
valid <- !is.na(edna_annot$new_scientific_name_ncbi)
edna_annot$best_taxonomic_assignment <- NA_character_
edna_annot$best_taxonomic_assignment[valid] <- paste0(
  "MOTU", sprintf("%03d", seq_len(sum(valid))), "_",
  edna_annot$new_scientific_name_ncbi[valid]
)
edna_clean <- merge(edna_clean, edna_annot[, c("sequence", "best_taxonomic_assignment")],
                    by = "sequence", all.x = TRUE)

# Drop columns used only temporarily: source annotation columns and env
# variables brought in solely to ensure species-name consistency during the
# merge (the env file is saved separately as its own clean output).
env_cols_added <- setdiff(colnames(env_clean),
                          c("Station", "Site", "Habitat", "Depth", "Latitude", "Longitude"))
cols_to_remove <- c("sequence", "new_scientific_name_ncbi", "SpeciesFB", env_cols_added)
edna_clean <- edna_clean[, !colnames(edna_clean) %in% cols_to_remove]

message("\n--- eDNA data file ---")
message("Original : ", ncol(edna_raw),   " columns, ", nrow(edna_raw),   " rows")
message("Cleaned  : ", ncol(edna_clean), " columns (incl. best_taxonomic_assignment), ",
        nrow(edna_clean), " rows")
message("Dropped  : ", paste(edna_dropped, collapse = ", "))

out_edna <- file.path(data_dir, "eDNA_Data_SEAMOUNTS_REEF3.0_merged_clean.csv")
write_csv(edna_clean, out_edna)
message("Saved    : ", out_edna)
