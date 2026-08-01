# ====================================================================================================
# Script name: save_table_data.R
# Author: Estephe Kana & Edi Prifti & Eugeni Belda
# Date created: 2025-12-10
# Purpose: Filter eDNA presence/absence data by prevalence and save
#          presence/absence table(s) for downstream analyses.
# Inputs: data/smX_pres_abs_matrix.rda
# Outputs:
#   (default) analyses/files/txt/presanceAbsence_table_prev_<N>.txt                  — full dataset, used by ScaleNet
#   (optional) analyses/files/txt/presanceAbsence_table_prev_<N>_<stratum>.txt         — per stratum
#              analyses/files/txt/presanceAbsence_table_prev_<N>_<stratum1>_<stratum2>.txt — pairwise
# ====================================================================================================

# ----------------------------------------------------------------------------------------------------
# Commands to run the script (from repo root):
#   export REPO_ROOT="/data/projects/aime/analyses/seamount/metabarcoding/"
#   export DATA_DIR="/data/projects/aime/data/seamount/metabarcoding"
#
#   From REPO_ROOT, run the following commands:
#
#   Rscript analyses/scripts/save_table_data.R 3          # ScaleNet table only (default)
#   Rscript analyses/scripts/save_table_data.R 3 TRUE     # also save stratified versions
# ----------------------------------------------------------------------------------------------------

# define arguments for the script
args <- commandArgs(trailingOnly = TRUE)

# Prevalence rate used to label the output filenames (see comment below)
species_prev_rate <- as.numeric(args[1])

# Whether to also save per-stratum and pairwise tables (default: FALSE)
save_stratified <- if (length(args) >= 2) as.logical(args[2]) else FALSE

# Derive repo root from script location
script_path  <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
script_dir   <- dirname(script_path)         # <repo>/analyses/scripts
repo_root    <- dirname(dirname(script_dir)) # <repo>/
data_dir     <- file.path(repo_root, "data")

# Allow override via environment variable (set on the remote server)
repo_root <- Sys.getenv("REPO_ROOT", unset = repo_root)
data_dir  <- Sys.getenv("DATA_DIR",  unset = file.path(repo_root, "data"))

analyses_dir <- file.path(repo_root, "analyses")
source(file.path(script_dir, "utils.R"))

# load dataset: samples x MOTU presence/absence, with Habitat/Zone/hab_inoff
load(file.path(data_dir, "smX_pres_abs_matrix.rda"))
meta_cols <- c("Station", "Habitat", "Zone", "hab_inoff")
edna_presenceAbsence <- as.matrix(smX_pres_abs_matrix[, !colnames(smX_pres_abs_matrix) %in% meta_cols])
rownames(edna_presenceAbsence) <- smX_pres_abs_matrix$Station

# Already prevalence-filtered upstream (make_db_object_clean.R); species_prev_rate only labels output filenames.
# filtered_edna_presenceAbsence <- get_sample_by_prevalence(edna_presenceAbsence, species_prev_rate)
filtered_edna_presenceAbsence <- edna_presenceAbsence

# Save full table for ScaleNet (default output)
write.table(
  filtered_edna_presenceAbsence,
  file  = file.path(analyses_dir, "files", "txt", paste0("presanceAbsence_table_prev_", species_prev_rate, ".txt")),
  quote = FALSE,
  sep   = "\t",
  row.names = TRUE,
  col.names = TRUE
)
message("Saved: presanceAbsence_table_prev_", species_prev_rate, ".txt (",
        nrow(filtered_edna_presenceAbsence), " samples x ",
        ncol(filtered_edna_presenceAbsence), " species)")

# Optionally save per-stratum and pairwise tables
if (save_stratified) {
  sample.info <- smX_pres_abs_matrix[, c("Station", "Zone")]

  presabs.df          <- as.data.frame(filtered_edna_presenceAbsence)
  presabs.df$Station   <- rownames(presabs.df)
  presabs.info.df     <- merge(presabs.df, sample.info[, c("Station", "Zone")], by = "Station", all.x = TRUE)

  save_table <- function(df, path) {
    df <- df[rowSums(df) != 0, ]
    df <- df[, colSums(df) != 0, drop = FALSE]
    write.table(df, file = path, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
    message("Saved: ", basename(path), " (", nrow(df), " samples × ", ncol(df), " species)")
  }

  # Per-stratum tables
  for (zone in unique(na.omit(presabs.info.df$Zone))) {
    tbl <- presabs.info.df[presabs.info.df$Zone %in% zone, ]
    rownames(tbl) <- tbl$Station
    tbl <- tbl[, -match(c("Station", "Zone"), colnames(tbl))]
    save_table(tbl, file.path(analyses_dir, "files", "txt",
                              paste0("presanceAbsence_table_prev_", species_prev_rate, "_", zone, ".txt")))
  }

  # Pairwise comparison tables
  comp <- combn(x = unique(na.omit(presabs.info.df$Zone)), m = 2, simplify = FALSE)
  for (pair in comp) {
    tbl <- presabs.info.df[presabs.info.df$Zone %in% pair, ]
    rownames(tbl) <- tbl$Station
    tbl <- tbl[, -match(c("Station", "Zone"), colnames(tbl))]
    save_table(tbl, file.path(analyses_dir, "files", "txt",
                              paste0("presanceAbsence_table_prev_", species_prev_rate,
                                     "_", pair[[1]], "_", pair[[2]], ".txt")))
  }
}

