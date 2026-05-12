# ====================================================================================================
# Script name: save_table_data.R
# Author: Estephe Kana & Edi Prifti & Eugeni Belda
# Date created: 2025-12-10
# Purpose: Filter eDNA abundance data by prevalence and save presence/absence
#          table(s) for downstream analyses.
# Inputs: abundance_data_matrix
# Outputs:
#   (default) presanceAbsence_table.txt                              — full dataset, used by ScaleNet
#   (optional) presanceAbsence_table_prev_<N>_<stratum>.txt         — per stratum
#              presanceAbsence_table_prev_<N>_<stratum1>_<stratum2>.txt — pairwise
# ====================================================================================================

# ----------------------------------------------------------------------------------------------------
# Commands to run the script (from repo root):
#   Rscript analyses/scripts/save_table_data.R 3          # ScaleNet table only (default)
#   Rscript analyses/scripts/save_table_data.R 3 TRUE     # also save stratified versions
# ----------------------------------------------------------------------------------------------------

# Check for required packages
required_pkgs <- c("vegan")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))

library(vegan)

# define arguments for the script
args <- commandArgs(trailingOnly = TRUE)

# Define threshold variables to select edges
species_prev_rate <- as.numeric(args[1])

# Whether to also save per-stratum and pairwise tables (default: FALSE)
save_stratified <- if (length(args) >= 2) as.logical(args[2]) else FALSE

# Derive repo root from script location
script_path  <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
script_dir   <- dirname(script_path)         # <repo>/analyses/scripts
repo_root    <- dirname(dirname(script_dir)) # <repo>/
data_dir     <- file.path(repo_root, "data")
analyses_dir <- file.path(repo_root, "analyses")
source(file.path(script_dir, "utils.R"))

# load dataset
load(file.path(data_dir, "seamount_integrated_dataset.rda"))

# get data abundance table
edna_abundance <- t(sm$X)

# # get samples filtered at species_prev_rate % of prevalence of the total of samples
filtered_edna_abundance <- get_sample_by_prevalence(edna_abundance, species_prev_rate)

# presence/absence table
filtered_edna_presenceAbsence <- decostand(filtered_edna_abundance, method = "pa")

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
  sample.info <- sm$sample_info
  sample.info$Zone <- ifelse(sample.info$Habitat %in% c("Bay","Lagoon","Reef_outer_slope", "Soft_back_reef"), "Shallow",
                             ifelse(sample.info$Habitat %in% c("Summit50", "DeepSlope"), "Middle",
                                    ifelse(sample.info$Habitat %in% c("Summit250", "Summit500"), "Deep",
                                           NA)))

  presabs.df          <- as.data.frame(filtered_edna_presenceAbsence)
  presabs.df$Spygen   <- rownames(presabs.df)
  presabs.info.df     <- merge(presabs.df, sample.info[, c("Spygen", "Zone")], by = "Spygen", all.x = TRUE)

  save_table <- function(df, path) {
    df <- df[rowSums(df) != 0, ]
    df <- df[, colSums(df) != 0, drop = FALSE]
    write.table(df, file = path, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
    message("Saved: ", basename(path), " (", nrow(df), " samples × ", ncol(df), " species)")
  }

  # Per-stratum tables
  for (zone in unique(na.omit(presabs.info.df$Zone))) {
    tbl <- presabs.info.df[presabs.info.df$Zone %in% zone, ]
    rownames(tbl) <- tbl$Spygen
    tbl <- tbl[, -match(c("Spygen", "Zone"), colnames(tbl))]
    save_table(tbl, file.path(analyses_dir, "files", "txt",
                              paste0("presanceAbsence_table_prev_", species_prev_rate, "_", zone, ".txt")))
  }

  # Pairwise comparison tables
  comp <- combn(x = unique(na.omit(presabs.info.df$Zone)), m = 2, simplify = FALSE)
  for (pair in comp) {
    tbl <- presabs.info.df[presabs.info.df$Zone %in% pair, ]
    rownames(tbl) <- tbl$Spygen
    tbl <- tbl[, -match(c("Spygen", "Zone"), colnames(tbl))]
    save_table(tbl, file.path(analyses_dir, "files", "txt",
                              paste0("presanceAbsence_table_prev_", species_prev_rate,
                                     "_", pair[[1]], "_", pair[[2]], ".txt")))
  }
}

