# =============================================================================
# Script name: make_db_object_clean.R
# Author: Estephe Kana
# Date created: 2026-07-01
# Purpose: Build the integrated eDNA dataset object (sm) from the raw data
# Inputs:  data/eDNA_Data_SEAMOUNTS_REEF3.0_merged_final.csv
#          data/eDNA_SEAMOUNTS_REEF3.0_merged_Environmental_Variables_final.csv
#          data/highlighted_taxa_list.csv (built by extract_highlighted_taxa.R)
#          NCBI Taxonomy, queried live via taxize (no local reference file)
# Outputs: data/seamount_integrated_dataset_final.rda   (sm object, sm$X is presence/absence)
#          data/smX_pres_abs_matrix.rda                 (samples x MOTU presence/absence,
#                                                         with Habitat/Zone/hab_inoff columns
#                                                         used by Predomics, ScaleNet, GSEA)
#          data/sm_taxonomy.csv                         (sm$taxonomy, one row per MOTU)
#          data/sm_sample_info.csv                      (sm$sample_info, one row per sample)
# Run:     Rscript analyses/scripts/make_db_object_clean.R
# =============================================================================

required_pkgs <- c("readr", "data.table", "taxize", "reshape2")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))

library(readr)
library(data.table)
library(taxize)

# Derive repo root from script location
script_path <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
script_dir  <- dirname(script_path)
repo_root   <- dirname(dirname(script_dir))
data_dir    <- file.path(repo_root, "data")

# -----------------------------------------------------------------------------
# Load clean inputs
# -----------------------------------------------------------------------------
db  <- as.data.frame(suppressMessages(
  read_csv(file.path(data_dir, "eDNA_Data_SEAMOUNTS_REEF3.0_merged_final.csv"))
))
env <- as.data.frame(suppressMessages(
  read_csv(file.path(data_dir, "eDNA_SEAMOUNTS_REEF3.0_merged_Environmental_Variables_final.csv"))
))

# Drop rows with no taxonomic assignment
table(is.na(db$best_taxonomic_assignment))
db <- db[!is.na(db$best_taxonomic_assignment), ]

# Keep MOTU names only for taxa in the manuscript
highlighted_tax_list <- read.csv(file.path(data_dir, "highlighted_taxa_list.csv"),
                                 stringsAsFactors = FALSE)
motu_codes_to_keep <- highlighted_tax_list$feature

db$MOTU_code <- ifelse(
  db$best_taxonomic_assignment %in% motu_codes_to_keep,
  db$best_taxonomic_assignment,          # keep as-is
  sub("_.*", "", db$best_taxonomic_assignment)  # else strip to MOTU code
)

# -----------------------------------------------------------------------------
# Merge env into db by Station
# -----------------------------------------------------------------------------
cols_to_drop <- intersect(c("Site", "Habitat", "Depth", "Latitude", "Longitude"),
                          colnames(env))
dba <- merge(x = db,
             y = env[, !colnames(env) %in% cols_to_drop],
             by = "Station", all.x = TRUE)

# -----------------------------------------------------------------------------
# Pivot long -> wide (samples x MOTUs)
# -----------------------------------------------------------------------------
meta_cols <- c("Station", "Site", "Habitat", "Depth",
               "Latitude", "Longitude","Salinity", "SSTmean",
               "seafloorTemp", "Chla", "TravelTime", "ReefMinDist.m")

setDT(dba)
db.wide <- dcast(dba,
  Station + Site + Habitat + Depth + Latitude + Longitude +
  + Salinity  + SSTmean + seafloorTemp + Chla + TravelTime + ReefMinDist.m ~ MOTU_code,
  value.var = "mean_pcr_count_reads", na.rm = FALSE)

db.wide <- as.data.frame(db.wide)
# order by Station
db.wide <- db.wide[order(db.wide$Station), ]

# Keep species at 3% prevalence, matching scalenet_network_inference.R /
# codePredomics_prev.R / save_table_data.R
source(file.path(script_dir, "utils.R"))
motu_mat <- as.matrix(db.wide[, !colnames(db.wide) %in% meta_cols])
motu_mat[is.na(motu_mat)] <- 0
rownames(motu_mat) <- db.wide$Station
motu_mat <- get_sample_by_prevalence(motu_mat, 3)
db.wide <- db.wide[db.wide$Station %in% rownames(motu_mat), c(meta_cols, colnames(motu_mat))]

# Extract sample metadata and build OTU presence/absence matrix
sample.info <- db.wide[, meta_cols]

X <- db.wide[, !colnames(db.wide) %in% meta_cols]
rownames(X) <- db.wide$Station
X <- t(X)
X[is.na(X)] <- 0
X <- (X > 0) * 1

# -----------------------------------------------------------------------------
# Build taxonomy lookup by querying NCBI Taxonomy directly (via taxize)
# -----------------------------------------------------------------------------
motu_ids     <- rownames(X)
motu_lookup  <- unique(dba[, c("MOTU_code", "best_taxonomic_assignment")])
motu_names   <- motu_lookup$best_taxonomic_assignment[match(motu_ids, motu_lookup$MOTU_code)]
motu_names   <- sub("^MOTU[0-9]+_", "", motu_names)
motu_df      <- data.frame(id = motu_ids, name = gsub("\\.", " ", motu_names),
                           stringsAsFactors = FALSE)

taxa_names <- unique(motu_df$name)

# Retrieve NCBI taxonomy IDs, then the full lineage for each taxon
uids <- get_uid(taxa_names, ask = FALSE, messages = TRUE)
tax_classification <- classification(uids, db = "ncbi")

# Tag each lineage with the NCBI taxid it was retrieved for
for (i in names(tax_classification))
  if (!is.na(i)) tax_classification[[i]]$taxid_source <- i

tax_classification.df <- do.call("rbind", tax_classification)
tax_classification.df <- tax_classification.df[
  tax_classification.df$rank %in% c("class", "order",
                                    "family", "genus", "species"), ]

tax_classification.df_wide <- reshape2::dcast(
  data = tax_classification.df, formula = taxid_source ~ rank, value.var = "name")
tax_classification.df_wide <- tax_classification.df_wide[
  , c("taxid_source", "class", "order", "family", "genus", "species")]

# Link each resolved taxid back to the name it was queried with
uids.df <- as.data.frame(uids)
uids.df$sourceTax <- taxa_names
tax_classification.df_wide <- merge(tax_classification.df_wide, uids.df[, c("ids", "sourceTax")],
                                    by.x = "taxid_source", by.y = "ids", all.x = TRUE)

reftax.list <- merge(motu_df, tax_classification.df_wide,
                     by.x = "name", by.y = "sourceTax", all.x = TRUE)
reftax.list <- reftax.list[match(motu_df$id, reftax.list$id), ]
rownames(reftax.list) <- reftax.list$id

reftax.list$tax_name <- reftax.list$species
reftax.list$tax_id   <- reftax.list$taxid_source
reftax.list$raw_tax_name <- motu_lookup$best_taxonomic_assignment[match(reftax.list$id, motu_lookup$MOTU_code)]
reftax.list$MOTU_name   <- reftax.list$name
# Start from the original name
reftax.list$MOTU_code <- reftax.list$name

# Flag all occurrences of names that are duplicated
dup <- reftax.list$MOTU_code %in% reftax.list$MOTU_code[duplicated(reftax.list$MOTU_code)]

# Per-group sequential counter
suffix <- ave(seq_along(reftax.list$MOTU_code), reftax.list$MOTU_code,
              FUN = function(i) sprintf("MOTU%02d", seq_along(i)))

# Append suffix only to duplicated names
reftax.list$MOTU_code[dup] <- paste0(reftax.list$MOTU_code[dup], "_", suffix[dup])

# Replace spaces with "."
reftax.list$MOTU_code <- gsub(" ", ".", reftax.list$MOTU_code)
reftax.list <- reftax.list[, c("MOTU_name", "MOTU_code", "class", "order", "family",
                              "genus", "tax_id", "tax_name", "raw_tax_name")]

# -----------------------------------------------------------------------------
# Enrich sample metadata
# -----------------------------------------------------------------------------
richness    <- data.frame(sample = colnames(X), MOTU_richness = colSums(X))
sample.info <- merge(sample.info, richness, by.x = "Station", by.y = "sample")

# -----------------------------------------------------------------------------
# Output 1: sm object (.rda)
# -----------------------------------------------------------------------------
sm <- list(db_long = dba, sample_info = sample.info, X = X, taxonomy = reftax.list)
out_rda <- file.path(data_dir, "seamount_integrated_dataset_clean.rda")
save(sm, file = out_rda)
message("Saved: seamount_integrated_dataset_final.rda (",
        nrow(sm$X), " species x ", ncol(sm$X), " samples, presence/absence)")

# -----------------------------------------------------------------------------
# Output 2: smX_pres_abs_matrix (samples x MOTU presence/absence), with
# Habitat/Zone/hab_inoff columns so downstream analyses (Predomics, ScaleNet,
# GSEA)
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

smX_pres_abs_matrix <- merge(sample.info[, c("Station", "Habitat")],
                             as.data.frame(t(X)),
                             by.x = "Station", by.y = "row.names")
smX_pres_abs_matrix$Zone      <- unname(zone_lookup[smX_pres_abs_matrix$Habitat])
smX_pres_abs_matrix$hab_inoff <- unname(hab_inoff_lookup[smX_pres_abs_matrix$Habitat])
rownames(smX_pres_abs_matrix) <- smX_pres_abs_matrix$Station

motu_cols <- setdiff(colnames(smX_pres_abs_matrix),
                     c("Station", "Habitat", "Zone", "hab_inoff"))
smX_pres_abs_matrix <- smX_pres_abs_matrix[, c("Station", "Habitat", "Zone", "hab_inoff", motu_cols)]

out_rda_X <- file.path(data_dir, "smX_pres_abs_matrix_final.rda")
save(smX_pres_abs_matrix, file = out_rda_X)
message("Saved: smX_pres_abs_matrix_final.rda (",
        nrow(smX_pres_abs_matrix), " samples x ", length(motu_cols), " MOTUs)")

# -----------------------------------------------------------------------------
# Output 3: sm$taxonomy as CSV, so downstream figure scripts (e.g. Figure6,
# FigureS2_S3)
# -----------------------------------------------------------------------------
taxo_out <- sm$taxonomy
taxo_out$feature <- rownames(taxo_out)
taxo_out <- taxo_out[, c("feature", setdiff(colnames(taxo_out), "feature"))]

out_csv_taxo <- file.path(data_dir, "sm_taxonomy_final.csv")
write.csv(taxo_out, out_csv_taxo, row.names = FALSE)
message("Saved: sm_taxonomy.csv (", nrow(taxo_out), " MOTUs)")

# -----------------------------------------------------------------------------
# Output 4: sm$sample_info as CSV, so downstream figure scripts (e.g. Figure1)
# -----------------------------------------------------------------------------
out_csv_sample_info <- file.path(data_dir, "sm_sample_info_final.csv")
write.csv(sample.info, out_csv_sample_info, row.names = FALSE)
message("Saved: sm_sample_info_final.csv (", nrow(sample.info), " samples)")

