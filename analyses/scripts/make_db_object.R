# =============================================================================
# Script name: make_db_object.R
# Author: Estephe Kana & Eugeni Belda & Edi Prifti
# Date created: 2025-12-10
# Purpose: Build the integrated eDNA dataset object (sm) from raw CSV inputs
# Inputs: eDNA_Data_SEAMOUNTS_REEF3.0_merged.csv,
#         eDNA_SEAMOUNTS_REEF3.0_merged_Environmental_Variables_raw.csv,
#         DBtaxonomy.xlsx
# Outputs: data/seamount_integrated_dataset.rda
# =============================================================================

# -----------------------------------------------------------------------------
# Commands to run the script (from repo root):
# Rscript analyses/scripts/make_db_object.R
# -----------------------------------------------------------------------------

# Check for required packages
# Note: momr — remotes::install_github("eprifti/momr")
required_pkgs <- c("readr", "readxl", "data.table", "momr")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))

library(readr)
library(readxl)
library(data.table)
library(momr)

# Derive repo root from script location
script_path  <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
script_dir   <- dirname(script_path)         # <repo>/analyses/scripts
repo_root    <- dirname(dirname(script_dir)) # <repo>/
data_dir     <- file.path(repo_root, "data")

db  <- as.data.frame(suppressMessages(read_csv(file.path(data_dir, "eDNA_Data_SEAMOUNTS_REEF3.0_merged.csv"))))
env <- as.data.frame(suppressMessages(read_csv(file.path(data_dir, "eDNA_SEAMOUNTS_REEF3.0_merged_Environmental_Variables_raw.csv"))))

db <- db[!is.na(db$sequence), ]

# Add Station annotation to the main dataset
dba <- merge(x = db, y = env[, -match(c("Site", "Habitat", "Depth", "Latitude", "Longitude"), colnames(env))], by = "Station")

# Rename sequences with OTU ids
dba_annot <- unique(dba[, c("sequence", "new_scientific_name_ncbi", "SpeciesFB")])
dba_annot$id <- paste0("OTU", sprintf("%03d", seq_len(nrow(dba_annot))), "_", dba_annot$new_scientific_name_ncbi)

dba <- merge(dba, dba_annot[, c("sequence", "id")], by = "sequence", all.x = TRUE)

# Pivot to wide samples × species table
setDT(dba)
db.wide <- dcast(dba,
  Spygen + Site + Habitat + Substrate + Stratum + Station +
  Date + Time + Depth + Sampling_depth + Latitude + Longitude +
  Project + Method + EastwardVelocity + NorthwardVelocity +
  Salinity + SuspendedParticulateMatter + SSTmax + SSTmean +
  SSTmin + SSTsd + seafloorTemp + Chla + TravelTime + ReefMinDist.m ~ id,
  value.var = "mean_pcr_count_reads", na.rm = FALSE)
db.wide <- as.data.frame(db.wide)

# Extract sample metadata
meta_cols <- c("Spygen", "Site", "Habitat", "Substrate", "Stratum",
               "Station", "Date", "Time", "Depth", "Sampling_depth",
               "Latitude", "Longitude", "Project", "Method", "EastwardVelocity",
               "NorthwardVelocity", "Salinity", "SuspendedParticulateMatter",
               "SSTmax", "SSTmean", "SSTmin", "SSTsd", "seafloorTemp", "Chla",
               "TravelTime", "ReefMinDist.m")
sample.info <- db.wide[, meta_cols]

# Extract OTU abundance matrix (species × samples)
X <- db.wide[, -match(meta_cols, colnames(db.wide))]
rownames(X) <- db.wide$Spygen
X <- t(X)
X[is.na(X)] <- 0

# Build taxonomy lookup from DBtaxonomy.xlsx
reftax <- readxl::read_xlsx(path = file.path(data_dir, "DBtaxonomy.xlsx"))
reftax.list <- vector("list", nrow(dba_annot))

for (i in seq_len(nrow(dba_annot))) {
  i.id   <- as.character(dba_annot$id[i])
  iname  <- dba_annot$new_scientific_name_ncbi[i]
  levels <- c("tax_name", "genus", "family", "order", "class")

  idf <- data.frame(id = i.id, tax_name = NA_character_, tax_id = NA_character_,
                    genus = NA_character_, family = NA_character_,
                    order = NA_character_, class = NA_character_,
                    stringsAsFactors = FALSE)

  for (lvl in levels) {
    col <- if (lvl == "tax_name") "tax_name" else lvl
    idx <- match(iname, reftax[[col]])
    if (!is.na(idx)) {
      if (lvl == "tax_name") idf$tax_id <- as.character(reftax$tax_id[idx])
      for (fill in c("tax_name", "genus", "family", "order", "class")) {
        if (fill %in% colnames(reftax)) idf[[fill]] <- as.character(reftax[[fill]][idx])
      }
      break
    }
  }
  reftax.list[[i]] <- idf
}

reftax.list <- do.call("rbind", reftax.list)
rownames(reftax.list) <- reftax.list$id
reftax.list <- reftax.list[, -1]
reftax.list <- reftax.list[, rev(colnames(reftax.list))]

# Normalise abundance data
X <- normFreqTC(as.matrix(X))

richness    <- data.frame(sample = colnames(X), otu_richness = colSums(X > 0))
sample.info <- merge(sample.info, richness, by.x = "Spygen", by.y = "sample")

ind7 <- nchar(sample.info$Station) == 7
ind9 <- nchar(sample.info$Station) == 9
sample.info$Station_id <- NA_character_
sample.info$Station_id[ind7] <- substr(sample.info$Station[ind7], 1, 5)
sample.info$Station_id[ind9] <- substr(sample.info$Station[ind9], 1, 7)
sample.info$Filter_id <- NA_character_
sample.info$Filter_id[ind7] <- substr(sample.info$Station[ind7], 6, 7)
sample.info$Filter_id[ind9] <- substr(sample.info$Station[ind9], 8, 9)

# Assemble and save the unified data object
sm <- list(db_long = dba, sample_info = sample.info, X = X, taxonomy = reftax.list)

save(sm, file = file.path(data_dir, "seamount_integrated_dataset.rda"))
message("Saved: seamount_integrated_dataset.rda (",
        nrow(sm$X), " species x ", ncol(sm$X), " samples)")
