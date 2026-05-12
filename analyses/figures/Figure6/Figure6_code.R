# =============================================================================
# Script name: Figure6_code.R
# Authors: Estephe Kana & Edi Prifti & Eugeni Belda
# Purpose: Build an iTOL-style circular phylogenetic tree of marine fish
#   species from the network inferred by ScaleNet, annotated with:
#     - Concentric rings for Zone, and Habitat associations
#     - Bar plots for mean feature importance, betweenness centrality, and degree centrality
#     - Indicator taxa highlighted in red

# =============================================================================

# -----------------------------------------------------------------------------
# Commands to run the script (from repo root):
# Rscript analyses/figures/Figure6/Figure6_code.R
# -----------------------------------------------------------------------------

# Check for required packages
# Note: ggtree, ggtreeExtra — BiocManager::install(c("ggtree", "ggtreeExtra"))
required_pkgs <- c("ape", "ggtree", "ggtreeExtra", "ggnewscale",
                   "ggplot2", "dplyr", "tidyr", "igraph")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))

library(ape)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(ggplot2)
library(dplyr)
library(tidyr)
library(igraph)

# Derive repo root from script location
script_path  <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
script_dir   <- dirname(script_path)                      # <repo>/analyses/figures/Figure6
repo_root    <- dirname(dirname(dirname(script_dir)))     # <repo>/
data_dir     <- file.path(repo_root, "data")
analyses_dir <- file.path(repo_root, "analyses")
source(file.path(repo_root, "analyses", "scripts", "utils.R"))

objlist <- list()

# load bininter results
x <- load(file.path(analyses_dir, "analysis_outputs", "bininter_output_data", "bininter_Predomics_all_analyses_overall_data_strat_group_prev_3.Rda"))
for(i in x)
{
  objlist[["bininter"]][[i]] <- get(i)
}
rm(list = x)

# load terinter results
x <- load(file.path(analyses_dir, "analysis_outputs", "terinter_output_data", "terinter_Predomics_all_analyses_overall_data_strat_group_prev_3.Rda"))
for(i in x)
{
  objlist[["terinter"]][[i]] <- get(i)
}
rm(list = x)

mergedf.all <- rbind(objlist$bininter$predout.bin.sub,
                     objlist$terinter$predout.bin.sub)
mergedf.all.bplot <- mergedf.all[mergedf.all$IsIndSp == 1,]
mergedf.all.bplot <- data.frame(table(mergedf.all.bplot$source, mergedf.all.bplot$data, mergedf.all.bplot$comparison))

# =============================================================================
# PATHS
# =============================================================================

rda_files       <- list.files(file.path(analyses_dir, "files", "rdata", "graph_data"),
                               pattern = "^graph_data_ecorr50_all_strat_", full.names = TRUE)
graph_data_path <- rda_files[which.max(file.mtime(rda_files))]
dataset_path    <- file.path(data_dir, "seamount_integrated_dataset.rda")

species_prev_rate <- 3

# out_pdf <- file.path(script_dir, paste0("Figure6_", Sys.Date(), ".pdf"))

# =============================================================================
# LOAD DATA
# =============================================================================

load(graph_data_path)   # -> network, lay, nodes.annot

# -- eDNA abundance table -----------------------------------------------------
load(dataset_path)   # -> sm

edna_abundance <- t(sm$X)

filtered_edna_abundance       <- get_sample_by_prevalence(edna_abundance, species_prev_rate)
filtered_edna_presenceAbsence <- vegan::decostand(filtered_edna_abundance, method = "pa")

# -- Sample metadata ----------------------------------------------------------
sample.info <- sm$sample_info
colnames(sample.info)[colnames(sample.info) == "DeepSlope"] <- "DeepSlope150"
sample.info$Zone <- ifelse(
  sample.info$Habitat %in% c("Bay", "Lagoon", "Reef_outer_slope", "Soft_back_reef"), "Shallow",
  ifelse(sample.info$Habitat %in% c("Summit50", "DeepSlope150"), "Middle",
         ifelse(sample.info$Habitat %in% c("Summit250", "Summit500"), "Deep", NA))
)

# =============================================================================
# MODULE DETECTION — fast greedy on undirected graph
# =============================================================================

network_undirected  <- as_undirected(network, mode = "collapse")

fast_greedy_cluster <- cluster_fast_greedy(network_undirected)

modularity.df <- data.frame(
  species     = V(network_undirected)$label,
  name        = V(network_undirected)$name,
  fast_greedy = fast_greedy_cluster$membership,
  stringsAsFactors = FALSE
)

fast_greedy_modules <- split(modularity.df$name, modularity.df$fast_greedy)

# -- Cluster colours (one per module) -----------------------------------------
cluster_ids    <- sort(unique(modularity.df$fast_greedy))
cluster_labels <- paste0("Cluster_", sprintf("%02d", cluster_ids))

distinct_13 <- c(
  "#E6194B", "#3CB44B", "#4363D8", "#F58231",
  "#911EB4", "#42D4F4", "#F032E6", "#BFEF45",
  "#FABED4", "#469990", "#DCBEFF", "#9A6324",
  "#800000"
)

cluster_palette <- setNames(
  distinct_13[seq_along(cluster_ids)],
  cluster_labels
)

species_to_cluster <- setNames(
  paste0("Cluster_", sprintf("%02d", modularity.df$fast_greedy)),
  modularity.df$name
)

node_cluster_colors <- sapply(V(network)$name, function(node) {
  if (node %in% names(species_to_cluster)) cluster_palette[species_to_cluster[node]]
  else "gray80"
})

# get results of chi-square post Hoc
sp.chisq_zone <- sp.chisq_posthoc
sp.chisq_habitat <- sp.chisq_posthoc.habitat

# Step 1 — fix column names first
cols_to_rename <- c("pval_chisq", "padj_chisq", "mean_presence_count",
                    "enriched_zones", "depleted_zones", "assigned_class")
colnames(sp.chisq_habitat)[colnames(sp.chisq_habitat) %in% paste0(cols_to_rename, ".habitat")] <- cols_to_rename

# Step 2 — recode DeepSlope AFTER assigned_class column exists under its correct name
sp.chisq_habitat$assigned_class[sp.chisq_habitat$assigned_class == "DeepSlope"] <- "DeepSlope150"

# Step 3 — also recode any DeepSlope column name suffixes (padj_PH_, residual_)
colnames(sp.chisq_habitat) <- gsub("DeepSlope", "DeepSlope150", colnames(sp.chisq_habitat))

nodes_tree.df <- modularity.df[, c("name", "fast_greedy")] %>%
  # Join habitat-level chi-sq assignment
  dplyr::left_join(
    sp.chisq_habitat %>% dplyr::select(feature, assigned_habitat = assigned_class),
    by = c("name" = "feature")
  ) %>%
  # Join zone-level chi-sq assignment
  dplyr::left_join(
    sp.chisq_zone %>% dplyr::select(feature, assigned_zone = assigned_class),
    by = c("name" = "feature")
  ) %>%
  dplyr::mutate(
    Habitat = dplyr::case_when(
      !is.na(assigned_habitat) & assigned_habitat != "NS" ~ assigned_habitat,
      TRUE                                                ~ "NS.habitat"
    ),
    Zone = dplyr::case_when(
      !is.na(assigned_zone) & assigned_zone != "NS" ~ assigned_zone,
      TRUE                                          ~ "NS.zone"
    ),
    Module = paste0("Cluster_", sprintf("%02d", fast_greedy))
  )

nodes_tree.df <- nodes_tree.df[, c("name", "Habitat", "Zone", "Module")]

colnames(nodes_tree.df)[colnames(nodes_tree.df)=="name"] <- "feature"

taxonomy.df <- sm$taxonomy
rownames(taxonomy.df) <- gsub(" ", ".", rownames(taxonomy.df))

# merge with taxonomy
nodes_tree_taxo.df <- merge(nodes_tree.df, taxonomy.df, by.x = "feature", by.y = 0, all.x = TRUE)

# merge with centrality metrics
nodes_tree_taxo.df <- merge(
  nodes_tree_taxo.df,
  nodes.annot[, c("name", "betw_cent", "degr_cent")],
  by.x = "feature", by.y = "name", all.x = TRUE
)

# mean feature importance + indicator-species flag across all predomics models
feat_summary <- mergedf.all %>%
  dplyr::group_by(feature) %>%
  dplyr::summarise(
    mean_featureImportance = mean(featureImportance, na.rm = TRUE),
    IsIndSp                = any(IsIndSp == 1, na.rm = TRUE),
    .groups = "drop"
  )

nodes_tree_taxo.df <- nodes_tree_taxo.df %>%
  dplyr::left_join(feat_summary, by = "feature")

nodes_tree_taxo.df$mean_featureImportance[is.na(nodes_tree_taxo.df$mean_featureImportance)] <- 0
nodes_tree_taxo.df$IsIndSp[is.na(nodes_tree_taxo.df$IsIndSp)]                               <- FALSE

# ---------------------------------------------------------------------------
# 1. Build phylo object (unchanged)
# ---------------------------------------------------------------------------

tree_df <- nodes_tree_taxo.df %>%
  dplyr::select(feature, class, order, family, genus, Module, Zone, Habitat,
                betw_cent, degr_cent, mean_featureImportance, IsIndSp) %>%
  dplyr::mutate(
    class  = ifelse(is.na(class),  "Unk_class",  class),
    order  = ifelse(is.na(order),  paste0(class,  "_unk_order"),  order),
    family = ifelse(is.na(family), paste0(order,  "_unk_family"), family),
    genus  = ifelse(is.na(genus),  paste0(family, "_unk_genus"),  genus)
  )

edges_df <- bind_rows(
  tree_df %>% distinct(class)          %>% transmute(parent = "root", child = class),
  tree_df %>% distinct(class, order)   %>% transmute(parent = class,  child = order),
  tree_df %>% distinct(order, family)  %>% transmute(parent = order,  child = family),
  tree_df %>% distinct(family, genus)  %>% transmute(parent = family, child = genus),
  tree_df %>% distinct(genus, feature) %>% transmute(parent = genus,  child = feature)
) %>%
  distinct() %>%
  filter(parent != child)

tips      <- unique(tree_df$feature)
internal  <- setdiff(unique(c(edges_df$parent, edges_df$child)), tips)
all_nodes <- c(tips, internal)
node_map  <- setNames(seq_along(all_nodes), all_nodes)

edge_mat  <- matrix(
  c(node_map[edges_df$parent],
    node_map[edges_df$child]),
  ncol = 2
)

phylo_obj <- structure(
  list(
    edge        = edge_mat,
    tip.label   = tips,
    node.label  = internal,
    Nnode       = length(internal),
    edge.length = rep(0.3, nrow(edge_mat))
  ),
  class = "phylo"
)
phylo_obj <- ape::reorder.phylo(phylo_obj)

# ---------------------------------------------------------------------------
# 2. Annotation data frame (geom_fruit needs feature as a column)
# ---------------------------------------------------------------------------

annot <- tree_df %>%
  dplyr::select(feature, Module, Zone, Habitat,
                mean_featureImportance, betw_cent, degr_cent)

# ---------------------------------------------------------------------------
# 3. Palettes
# ---------------------------------------------------------------------------

zone_pal <- c(
  "Shallow" = "#5ae6ab",
  "Middle"  = "#ffe699",
  "Deep"    = "#25456B",
  "NS.zone" = "gray"
)

habitat_pal <- c(
  "Bay"              = "#5ae6ab",
  "Lagoon"           = "#88d941",
  "Soft_back_reef"   = "#2e7d00",
  "Reef_outer_slope" = "#4e8273",
  "Summit50"         = "#ffe699",
  "DeepSlope150"     = "#d79c3b",
  "Summit250"        = "#4a6a94",
  "Summit500"        = "#1a3250",
  "NS.habitat"       = "gray"
)

# ---------------------------------------------------------------------------
# 4. Base circular tree with family highlights
# ---------------------------------------------------------------------------

p_base <- suppressWarnings(
  ggtree(phylo_obj, layout = "circular", linewidth = 0.25, color = "grey50")
) %<+% (tree_df %>% dplyr::select(feature, IsIndSp) %>% dplyr::rename(label = feature))

tree_data    <- p_base$data
tree_radius  <- max(tree_data$x, na.rm = TRUE)

# Extract internal nodes for a given taxonomic rank and attach hjust.
.rank_nodes <- function(labels_vec, tree_data) {
  tree_data %>%
    dplyr::filter(!isTip, label %in% labels_vec) %>%
    dplyr::select(node, label) %>%
    dplyr::left_join(tree_data %>% dplyr::select(node, angle), by = "node") %>%
    dplyr::mutate(
      hjust = ifelse(angle > 90 & angle < 260, 1, 0),
      label = gsub("_unk_family", "\nunk_family", label) 
    )
}

class_nodes  <- .rank_nodes(unique(tree_df$class),  tree_data)
order_nodes  <- .rank_nodes(unique(tree_df$order),  tree_data)
family_nodes <- .rank_nodes(unique(tree_df$family), tree_data)
genus_nodes  <- .rank_nodes(unique(tree_df$genus),  tree_data)

# ---------------------------------------------------------------------------
# 5. Add concentric annotation rings and tip labels (with indicator species in red)
# ---------------------------------------------------------------------------

p3 <- p_base +
  geom_tiplab(
    aes(
      fontface = ifelse(IsIndSp, "bold.italic", "italic"),
      color    = IsIndSp  # map to logical directly
    ),
    size   = 2.4,
    offset = tree_radius * 1.15,
    align  = FALSE,
    family = "serif"
  ) +
  scale_color_manual(
    values = c("TRUE" = "firebrick", "FALSE" = "black"),
    guide  = "none"   # hides the legend
  ) +
  # Bar 1 — Mean feature importance (innermost bar ring)
  ggtreeExtra::geom_fruit(
    data    = annot %>% dplyr::filter(!is.na(mean_featureImportance)),
    geom    = geom_bar,
    mapping = aes(y = feature, x = mean_featureImportance,
                  fill = "Feature Importance"),
    stat        = "identity",
    orientation = "y",
    # colour      = NA,
    offset      = 0.05,
    axis.params = list(axis = "x", hjust = 1, vjust=1.2, nbreak=3, text.size  = 1.5),
    grid.params = list(linetype = 2, colour = "grey80")
  ) +
  # Bar 2 — Betweenness centrality
  ggtreeExtra::geom_fruit(
    data    = annot %>% dplyr::filter(!is.na(betw_cent)),
    geom    = geom_bar,
    mapping = aes(y = feature, x = betw_cent,
                  fill = "Betweenness centrality"),
    stat        = "identity",
    orientation = "y",
    # colour      = NA,
    offset      = 0.1,
    axis.params = list(axis = "x", hjust = 1, vjust=1.2, nbreak=3, text.size  = 1.5),
    grid.params = list(linetype = 2, colour = "grey80")
  ) +
  # Bar 3 — Degree centrality
  ggtreeExtra::geom_fruit(
    data    = annot %>% dplyr::filter(!is.na(degr_cent)),
    geom    = geom_bar,
    mapping = aes(y = feature, x = degr_cent,
                  fill = "Degree centrality"),
    stat        = "identity",
    orientation = "y",
    # colour      = NA,
    offset      = 0.1,
    axis.params = list(axis = "x", hjust = 1, vjust=1.2, nbreak=3, text.size  = 1.5),
    grid.params = list(linetype = 2, colour = "grey80")
  ) +
  scale_fill_manual(
    name   = "Metrics",
    values = c("Feature Importance"     = "#E07B54",
               "Betweenness centrality" = "#4A90C4",
               "Degree centrality"      = "#7B4F9E"),
    guide  = guide_legend(order = 10)
  ) +
  new_scale_fill() +
  # Ring 1 — Zone (innermost tile ring)
  ggtreeExtra::geom_fruit(
    data    = annot,
    geom    = geom_tile,
    mapping = aes(y = feature, fill = Zone),
    width   = 0.12,
    offset  = 0.1,
    color   = NA
  ) +
  scale_fill_manual(name = "Zone", values = zone_pal,
                    breaks = names(zone_pal), na.value = "grey85") +
  new_scale_fill() +
  # Ring 2 — Habitat (outermost)
  ggtreeExtra::geom_fruit(
    data    = annot,
    geom    = geom_tile,
    mapping = aes(y = feature, fill = Habitat),
    width   = 0.12,
    offset  = 0.1,
    color   = NA
  ) +
  scale_fill_manual(name = "Habitat", values = habitat_pal, breaks = names(habitat_pal), na.value = "grey85")

# ---------------------------------------------------------------------------
# 6. Final tree (tip labels and points already added in p3)
# ---------------------------------------------------------------------------

p_itol <- p3 +
  theme(
    legend.position  = c(1.0,0.88),
    legend.text      = element_text(size = 8),
    legend.title     = element_text(size = 10, face = "bold"),
    legend.key.size  = unit(0.4,  "cm"),
    legend.spacing.y = unit(0.12, "cm"),
    plot.margin      = margin(30, 30, 30, 30)
  )

# ---------------------------------------------------------------------------
# 7. Save
# ---------------------------------------------------------------------------

out_dir <- script_dir

ggsave(
  filename = file.path(out_dir, "Figure6.pdf"),
  plot     = p_itol,
  width    = 16, height = 15, units = "in",
  device   = cairo_pdf
)


message("iTOL-style circular tree saved.")