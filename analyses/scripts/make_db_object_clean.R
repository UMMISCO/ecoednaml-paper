# =============================================================================
# Script name: make_db_object_clean.R
# Author: Estephe Kana
# Date created: 2026-07-01
# Purpose: Build the integrated eDNA dataset object (sm) from the cleaned CSV
#          inputs produced by clean_input_csv.R. The clean eDNA file contains
#          sample/abundance/OTU columns only; environmental variables are
#          loaded from the separate clean env file and merged by Station,
#          mirroring the logic of make_db_object.R. sm$X stores MOTU
#          presence/absence (0/1), not normalized abundance.
# Inputs:  data/eDNA_Data_SEAMOUNTS_REEF3.0_merged_clean.csv
#          data/eDNA_SEAMOUNTS_REEF3.0_merged_Environmental_Variables_clean.csv
#          data/DBtaxonomy.xlsx
# Outputs: data/seamount_integrated_dataset_clean.rda   (sm object, sm$X is presence/absence)
#          data/SmX_pres_abs_matrix.rda                 (samples x MOTU presence/absence,
#                                                         with Habitat/Zone/hab_inoff columns
#                                                         used by Predomics, ScaleNet, GSEA)
# Run:     Rscript analyses/scripts/make_db_object_clean.R
# =============================================================================

required_pkgs <- c("readr", "readxl", "data.table")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))

library(readr)
library(readxl)
library(data.table)

# Derive repo root from script location
script_path <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
script_dir  <- dirname(script_path)
repo_root   <- dirname(dirname(script_dir))
data_dir    <- file.path(repo_root, "data")

# -----------------------------------------------------------------------------
# Load clean inputs
# -----------------------------------------------------------------------------
db  <- as.data.frame(suppressMessages(
  read_csv(file.path(data_dir, "eDNA_Data_SEAMOUNTS_REEF3.0_merged_clean.csv"))
))
env <- as.data.frame(suppressMessages(
  read_csv(file.path(data_dir, "eDNA_SEAMOUNTS_REEF3.0_merged_Environmental_Variables_clean.csv"))
))

# Drop rows with no taxonomic assignment, mirroring the !is.na(sequence)
# filter in make_db_object.R line 38.
db <- db[!is.na(db$best_taxonomic_assignment), ]

# -----------------------------------------------------------------------------
# Merge env into db by Station, dropping columns already present in db
# (mirrors make_db_object.R line 41)
# -----------------------------------------------------------------------------
cols_to_drop <- intersect(c("Site", "Habitat", "Depth", "Latitude", "Longitude"),
                          colnames(env))
dba <- merge(x = db,
             y = env[, !colnames(env) %in% cols_to_drop],
             by = "Station", all.x = TRUE)

# -----------------------------------------------------------------------------
# Pivot long -> wide (samples x MOTUs)
# -----------------------------------------------------------------------------
meta_cols <- c("Spygen", "Site", "Habitat", "Station", "Depth",
               "Latitude", "Longitude", "EastwardVelocity", "NorthwardVelocity",
               "Salinity", "SuspendedParticulateMatter", "SSTmean",
               "seafloorTemp", "Chla", "TravelTime", "ReefMinDist.m")

setDT(dba)
db.wide <- dcast(dba,
  Spygen + Site + Habitat + Station + Depth + Latitude + Longitude +
  EastwardVelocity + NorthwardVelocity + Salinity + SuspendedParticulateMatter +
  SSTmean + seafloorTemp + Chla + TravelTime + ReefMinDist.m ~ best_taxonomic_assignment,
  value.var = "mean_pcr_count_reads", na.rm = FALSE)
db.wide <- as.data.frame(db.wide)

# Extract sample metadata and build OTU presence/absence matrix
sample.info <- db.wide[, meta_cols]

X <- db.wide[, !colnames(db.wide) %in% meta_cols]
rownames(X) <- db.wide$Spygen
X <- t(X)
X[is.na(X)] <- 0
X <- (X > 0) * 1

# -----------------------------------------------------------------------------
# Build taxonomy lookup from DBtaxonomy.xlsx
# Scientific names are extracted from best_taxonomic_assignment by stripping
# the "MOTU001_" prefix.
# -----------------------------------------------------------------------------
reftax     <- readxl::read_xlsx(path = file.path(data_dir, "DBtaxonomy.xlsx"))
motu_ids   <- rownames(X)
motu_names <- sub("^MOTU[0-9]+_", "", motu_ids)
motu_df    <- data.frame(id = motu_ids, name = motu_names, stringsAsFactors = FALSE)

reftax.list <- vector("list", nrow(motu_df))
for (i in seq_len(nrow(motu_df))) {
  i.id  <- motu_df$id[i]
  iname <- motu_df$name[i]

  idf <- data.frame(id = i.id, tax_name = NA_character_, tax_id = NA_character_,
                    genus = NA_character_, family = NA_character_,
                    order = NA_character_, class = NA_character_,
                    stringsAsFactors = FALSE)

  for (lvl in c("tax_name", "genus", "family", "order", "class")) {
    col <- if (lvl == "tax_name") "tax_name" else lvl
    idx <- match(iname, reftax[[col]])
    if (!is.na(idx)) {
      idf$tax_id <- as.character(reftax$tax_id[idx])
      for (fill in c("tax_name", "genus", "family", "order", "class"))
        if (fill %in% colnames(reftax)) idf[[fill]] <- as.character(reftax[[fill]][idx])
      break
    }
  }
  reftax.list[[i]] <- idf
}

reftax.list <- do.call("rbind", reftax.list)
rownames(reftax.list) <- reftax.list$id
reftax.list <- reftax.list[, -1]
reftax.list <- reftax.list[, rev(colnames(reftax.list))]

# -----------------------------------------------------------------------------
# Enrich sample metadata
# -----------------------------------------------------------------------------
richness    <- data.frame(sample = colnames(X), otu_richness = colSums(X))
sample.info <- merge(sample.info, richness, by.x = "Spygen", by.y = "sample")

ind7 <- nchar(sample.info$Station) == 7
ind9 <- nchar(sample.info$Station) == 9
sample.info$Station_id <- NA_character_
sample.info$Station_id[ind7] <- substr(sample.info$Station[ind7], 1, 5)
sample.info$Station_id[ind9] <- substr(sample.info$Station[ind9], 1, 7)
sample.info$Filter_id <- NA_character_
sample.info$Filter_id[ind7] <- substr(sample.info$Station[ind7], 6, 7)
sample.info$Filter_id[ind9] <- substr(sample.info$Station[ind9], 8, 9)

# -----------------------------------------------------------------------------
# Output 1: sm object (.rda)
# -----------------------------------------------------------------------------
sm <- list(db_long = dba, sample_info = sample.info, X = X, taxonomy = reftax.list)
out_rda <- file.path(data_dir, "seamount_integrated_dataset_clean.rda")
save(sm, file = out_rda)
message("Saved: seamount_integrated_dataset_clean.rda (",
        nrow(sm$X), " species x ", ncol(sm$X), " samples, presence/absence)")

# -----------------------------------------------------------------------------
# Output 2: SmX_pres_abs_matrix (samples x MOTU presence/absence), with
# Habitat/Zone/hab_inoff columns so downstream analyses (Predomics, ScaleNet,
# GSEA) don't each need to re-derive these groupings from Habitat.
# -----------------------------------------------------------------------------
zone_lookup <- c(
  "Bay" = "Shallow", "Lagoon" = "Shallow", "OuterSlope" = "Shallow", "BackReef" = "Shallow",
  "Seamount50" = "Middle", "DeepSlope150" = "Middle",
  "Seamount250" = "Deep", "Seamount500" = "Deep"
)
hab_inoff_lookup <- c(
  "Bay" = "INSHORE", "Lagoon" = "INSHORE", "OuterSlope" = "INSHORE", "BackReef" = "INSHORE",
  "Seamount50" = "OFFSHORE", "DeepSlope150" = "OFFSHORE",
  "Seamount250" = "OFFSHORE", "Seamount500" = "OFFSHORE"
)

SmX_pres_abs_matrix <- merge(sample.info[, c("Spygen", "Habitat")],
                             as.data.frame(t(X)),
                             by.x = "Spygen", by.y = "row.names")
SmX_pres_abs_matrix$Zone      <- unname(zone_lookup[SmX_pres_abs_matrix$Habitat])
SmX_pres_abs_matrix$hab_inoff <- unname(hab_inoff_lookup[SmX_pres_abs_matrix$Habitat])
rownames(SmX_pres_abs_matrix) <- SmX_pres_abs_matrix$Spygen

motu_cols <- setdiff(colnames(SmX_pres_abs_matrix),
                     c("Spygen", "Habitat", "Zone", "hab_inoff"))
SmX_pres_abs_matrix <- SmX_pres_abs_matrix[, c("Spygen", "Habitat", "Zone", "hab_inoff", motu_cols)]

out_rda_X <- file.path(data_dir, "SmX_pres_abs_matrix.rda")
save(SmX_pres_abs_matrix, file = out_rda_X)
message("Saved: smX_pres_abs_matrix.rda (",
        nrow(SmX_pres_abs_matrix), " samples x ", length(motu_cols), " MOTUs)")
