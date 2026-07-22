# =============================================================================
# Script name: extract_highlighted_taxa.R
# Authors: Estephe Kana
# Purpose: Programmatically identify every MOTU that is specifically labelled
#          or singled out in the manuscript figures, then filter the raw eDNA
#          detection table down to just those MOTUs.
#
#          Sources of "highlighted" status, each reproducing the exact
#          selection rule used by its figure script:
#            - Figure 3    : the 55 indicator MOTUs identified by Predomics
#                            (bininter/terinter, IsIndSp == 1) across the
#                            three pairwise zone contrasts.
#            - Figure 4/5  : the subset of those 55 that survive in the
#                            |ecorr| > 0.5 co-occurrence network.
#            - Figure S2-D : the 3D-scatter outliers, using the same
#                            (betw_cent > 0.9 | log10_mfi > log10(3.62+1))
#                            & degr_cent > 7 condition as FigureS2_S3_code.R.
#            - Figure S3   : the union of the top-20 MOTUs by degree
#                            centrality and by betweenness centrality.
#
#          Outputs (written to data/):
#            - highlighted_taxa_list.csv                      : one row per
#              highlighted MOTU with its source-figure flags.
#            - eDNA_Data_SEAMOUNTS_REEF3.0_highlighted_taxa.csv : the rows of
#              eDNA_Data_SEAMOUNTS_REEF3.0_merged_clean.csv whose
#              best_taxonomic_assignment matches one of these MOTUs.
#
#          Note: matching against eDNA_Data_SEAMOUNTS_REEF3.0_merged_clean.csv
#          is done on the MOTU### code prefix only (not the full feature
#          string), since best_taxonomic_assignment separates multi-word
#          taxon names with a space while the figures/network data use a dot.
# =============================================================================

# -----------------------------------------------------------------------------
# Commands to run the script (from repo root):
#   export REPO_ROOT="/data/projects/aime/analyses/seamount/metabarcoding/"
#   export DATA_DIR="/data/projects/aime/data/seamount/metabarcoding"
#
#   From REPO_ROOT, run the following commands:
#
#   Rscript analyses/scripts/extract_highlighted_taxa.R
# -----------------------------------------------------------------------------

required_pkgs <- c("dplyr")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))

library(dplyr)

# Derive repo root from script location
script_path  <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
script_dir   <- dirname(script_path)         # <repo>/analyses/scripts
repo_root    <- dirname(dirname(script_dir)) # <repo>/

# Allow override via environment variable (set on the remote server)
repo_root <- Sys.getenv("REPO_ROOT", unset = repo_root)
data_dir  <- Sys.getenv("DATA_DIR",  unset = file.path(repo_root, "data"))

analyses_dir <- file.path(repo_root, "analyses")

# =============================================================================
# LOAD CO-OCCURRENCE NETWORK (Figures 4, 5, S2, S3)
# =============================================================================

rda_files       <- list.files(file.path(analyses_dir, "files", "rdata", "graph_data"),
                               pattern = "^graph_data_ecorr50_all_strat_.*\\.rda$", full.names = TRUE)
graph_data_path <- rda_files[which.max(file.mtime(rda_files))]
load(graph_data_path)   # -> network, lay, nodes, nodes.annot, sp.chisq_posthoc, sp.chisq_posthoc.habitat

# =============================================================================
# LOAD PREDOMICS INDICATOR-TAXA RESULTS (Figure 3)
# =============================================================================

load(file.path(analyses_dir, "analysis_outputs", "bininter_output_data",
                "bininter_Predomics_all_analyses_overall_data_strat_group_prev_3.Rda"))
bin.df <- predout.bin.sub
rm(adonis_pred.bin, predout.bin, predout.bin.sub)

load(file.path(analyses_dir, "analysis_outputs", "terinter_output_data",
                "terinter_Predomics_all_analyses_overall_data_strat_group_prev_3.Rda"))
ter.df <- predout.bin.sub
rm(adonis_pred.bin, predout.bin, predout.bin.sub)

mergedf.all <- rbind(bin.df, ter.df)

# 55 unique indicator MOTUs across the 3 pairwise zone contrasts (Figure 3)
species_predomics <- unique(as.character(
  mergedf.all$feature[mergedf.all$IsIndSp == 1 & mergedf.all$data == "pres/abs"]
))

# The subset of those 55 that survive in the |ecorr| > 0.5 network (Figures 4-5)
species_network <- intersect(species_predomics, nodes.annot$name)

# =============================================================================
# CENTRALITY / FEATURE-IMPORTANCE DERIVED SETS (Figures S2, S3)
# =============================================================================

feat_summary <- mergedf.all %>%
  dplyr::group_by(feature) %>%
  dplyr::summarise(mean_featureImportance = mean(featureImportance, na.rm = TRUE),
                    .groups = "drop")

nodes_merged <- merge(nodes.annot, feat_summary, by.x = "name", by.y = "feature", all.x = TRUE)
nodes_merged$mean_featureImportance[is.na(nodes_merged$mean_featureImportance)] <- 0
nodes_merged$log10_mfi <- log10(nodes_merged$mean_featureImportance + 1)

# Top-20 by degree / betweenness centrality (Figure S3)
top20_degree <- nodes_merged$name[order(-nodes_merged$degr_cent)][1:20]
top20_betw   <- nodes_merged$name[order(-nodes_merged$betw_cent)][1:20]

# 3D-plot outliers (Figure S2-D) — same selection rule as FigureS2_S3_code.R
outlier_condition <- (nodes_merged$betw_cent > 0.9 |
                         nodes_merged$log10_mfi > log10(3.62 + 1)) &
  nodes_merged$degr_cent > 7
outliers <- nodes_merged$name[outlier_condition]

# =============================================================================
# COMBINE INTO A SINGLE "HIGHLIGHTED TAXA" LOOKUP TABLE
# =============================================================================

highlighted <- unique(c(species_predomics, species_network, top20_degree, top20_betw, outliers))

taxa_lookup <- data.frame(
  motu_code    = sub("_.*", "", highlighted),
  feature      = highlighted,
  in_Figure3   = highlighted %in% species_predomics,
  in_Figure4_5 = highlighted %in% species_network,
  in_FigureS2  = highlighted %in% outliers,
  in_FigureS3  = highlighted %in% union(top20_degree, top20_betw),
  stringsAsFactors = FALSE
)
taxa_lookup <- taxa_lookup[order(as.numeric(sub("MOTU", "", taxa_lookup$motu_code))), ]

lookup_out <- file.path(data_dir, "highlighted_taxa_list.csv")
write.csv(taxa_lookup, lookup_out, row.names = FALSE)
message("Saved taxa lookup table (", nrow(taxa_lookup), " MOTUs) to ", lookup_out)