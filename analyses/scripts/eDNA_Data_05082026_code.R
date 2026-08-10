# =============================================================================
# Script name: eDNA_Data_05082026_code.R
# Author: Estephe Kana & Eugeni Belda
# Date created: 2026-08-05
# Purpose: Drop unused columns from the two raw input CSV files.
#          Only columns referenced in make_db_object_clean.R, codePredomics_prev.R,
#          scalenet_network_inference.R, save_table_data.R, and figure scripts
#          are retained.
# Inputs:  data/eDNA_SEAMOUNTS_REEF3.0_merged_Environmental_Variables_raw.csv
#          data/eDNA_Data_SEAMOUNTS_REEF3.0_merged.csv
# Outputs: data/Environmental_Variables_Data.csv
#          data/eDNA_Data.csv
# Run:     Rscript analyses/scripts/eDNA_Data_05082026_code.R
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
## check habitat nenaming
env_clean$Habitat_new <- recode_habitat(env_clean$Habitat)
table(env_clean$Habitat, env_clean$Habitat_new)
## rename the original habitat for save
env_clean$Habitat <- recode_habitat(env_clean$Habitat)
env_clean <- env_clean[,-ncol(env_clean)]

message("\n--- Environmental variables file ---")
message("Original : ", ncol(env_raw),   " columns, ", nrow(env_raw),   " rows")
message("Cleaned  : ", ncol(env_clean), " columns, ", nrow(env_clean), " rows")
message("Dropped  : ", paste(env_dropped, collapse = ", "))

out_env <- file.path(data_dir, paste0("Environmental_Variables_Data_", Sys.Date(), ".csv"))
write_csv(env_clean, out_env)
write.table(env_clean, file = gsub(".csv$",".tsv", out_env), sep = "\t", quote = FALSE, append = FALSE, row.names = FALSE, col.names = TRUE)
save(env_clean, file = gsub(".csv$",".Rda", out_env))
message("Saved    : ", out_env)

# =============================================================================
# 2. eDNA_Data_SEAMOUNTS_REEF3.0_merged.csv pre-processing
# =============================================================================

edna_cols_used <- c(
  "Station",
  "Spygen",
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
head(edna_raw)


edna_missing <- setdiff(edna_cols_used, colnames(edna_raw))
if (length(edna_missing) > 0)
  warning("Columns listed as used but not found in eDNA file: ",
          paste(edna_missing, collapse = ", "))

edna_dropped <- setdiff(colnames(edna_raw), edna_cols_used)
edna_clean   <- edna_raw[, intersect(edna_cols_used, colnames(edna_raw))]
##check habitat renaming
edna_clean$Habitat_new <- recode_habitat(edna_clean$Habitat)
table(edna_clean$Habitat, edna_clean$Habitat_new)
##ik, rename the original variable
edna_clean$Habitat <- recode_habitat(edna_clean$Habitat)
edna_clean <- edna_clean[,-ncol(edna_clean)]
edna_clean[str_detect(string = edna_clean$new_scientific_name_ncbi, pattern = "  ") & !is.na(edna_clean$new_scientific_name_ncbi),]
##recoding needed for scalenet format
edna_clean$new_scientific_name_ncbi <- gsub(" ", ".", edna_clean$new_scientific_name_ncbi)

# Merge with clean env to keep station order
cols_to_drop <- intersect(c("Site", "Habitat", "Depth", "Latitude", "Longitude"),
                          colnames(env_clean))
dim(edna_clean)
edna_clean <- merge(x = edna_clean,
                    y = env_clean[, !colnames(env_clean) %in% cols_to_drop],
                    by = "Station", all.x = TRUE)
dim(edna_clean)
# Derive best_taxonomic_assignment  from (sequence, new_scientific_name_ncbi, SpeciesFB), assign a sequential MOTU label
edna_annot <- unique(edna_clean[, c("sequence", "new_scientific_name_ncbi", "SpeciesFB")])
dim(edna_annot)
length(which(is.na(edna_annot$sequence)))
edna_annot[is.na(edna_annot$sequence),] ##1NA

dim(edna_clean[is.na(edna_clean$sequence),]) ##11 records=11 stations
unique(edna_clean[is.na(edna_clean$sequence),"Station"]) ##11 records=11 stations
unique(edna_clean[is.na(edna_clean$sequence),"mean_pcr_count_reads"]) ##11 records=11 stations; all zero mean_pcr_count_reads
##check the 11 stations
tstations <- unique(edna_clean[is.na(edna_clean$sequence),"Station"])
edna_clean[edna_clean$Station %in% tstations,] ##all zero in mean_pcr_count_reads, all NA in sequences

##check dims after excluding the 11 stations
length(unique(edna_clean[!edna_clean$Station %in% tstations,"Station"])) ##217 stations
length(unique(edna_clean[!edna_clean$Station %in% tstations,"sequence"])) ##967 sequences


valid <- !is.na(edna_annot$new_scientific_name_ncbi)
edna_annot$best_taxonomic_assignment <- NA_character_
edna_annot$best_taxonomic_assignment[valid] <- paste0(
  "MOTU", sprintf("%03d", seq_len(sum(valid))), "_",
  edna_annot$new_scientific_name_ncbi[valid]
)
edna_annot[is.na(edna_annot$best_taxonomic_assignment),] ##1 line (the NA sequence above)
edna_annot[str_detect(string = edna_annot$best_taxonomic_assignment, pattern = "MOTU037")  & !is.na(edna_annot$best_taxonomic_assignment),] ##1 single line for MOTU037

dim(edna_clean)
edna_clean <- merge(edna_clean, edna_annot[, c("sequence", "best_taxonomic_assignment")],
                    by = "sequence", all.x = TRUE)
dim(edna_clean)
head(edna_clean)


# Drop columns used only temporarily
env_cols_added <- setdiff(colnames(env_clean),
                          c("Station", "Spygen", "Site", "Habitat", "Depth", "Latitude", "Longitude"))
cols_to_remove <- c("sequence", "new_scientific_name_ncbi", "SpeciesFB", env_cols_added)
edna_clean <- edna_clean[, !colnames(edna_clean) %in% cols_to_remove]

message("\n--- eDNA data file ---")
message("Original : ", ncol(edna_raw),   " columns, ", nrow(edna_raw),   " rows")
message("Cleaned  : ", ncol(edna_clean), " columns (incl. best_taxonomic_assignment), ",
        nrow(edna_clean), " rows")
message("Dropped  : ", paste(edna_dropped, collapse = ", "))

head(edna_clean)
length(unique(edna_clean$Station))
length(unique(edna_clean$best_taxonomic_assignment))
length(which(is.na(edna_clean$best_taxonomic_assignment)))
tstat <- unique(edna_clean[is.na(edna_clean$best_taxonomic_assignment),"Station"])

tstat.df <- edna_clean[edna_clean$Station %in% tstat,]
dim(tstat.df)

length(unique(edna_clean$best_taxonomic_assignment)) ##968 unique sequences
length(unique(edna_clean$best_taxonomic_assignment)[!is.na(unique(edna_clean$best_taxonomic_assignment))]) ##967 unique non-na sequences

unique(edna_clean[str_detect(string = edna_clean$best_taxonomic_assignment, pattern = "MOTU037") & !is.na(edna_clean$best_taxonomic_assignment),"best_taxonomic_assignment"])
unique(edna_clean[str_detect(string = edna_clean$best_taxonomic_assignment, pattern = "\\.\\.") & !is.na(edna_clean$best_taxonomic_assignment),"best_taxonomic_assignment"])

##last modification (05/08/2026)=> change Spygen colname by Sample + replace SPY by SMP in same column
colnames(edna_clean) <- gsub("Spygen","Sample", colnames(edna_clean))
edna_clean$Sample <- gsub("^SPY","SMP", edna_clean$Sample)

out_edna <- file.path(data_dir, paste0("eDNA_Data_", Sys.Date(), ".csv"))
write_csv(edna_clean, out_edna)
write.table(edna_clean, file = gsub(".csv$",".tsv", out_edna), sep = "\t", quote = FALSE, append = FALSE, row.names = FALSE, col.names = TRUE)
save(edna_clean, file = gsub(".csv$",".Rda", out_edna))
message("Saved    : ", out_edna)
