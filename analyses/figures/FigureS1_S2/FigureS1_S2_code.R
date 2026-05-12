# =================================================================================
# Script name: FigureS1_S2.R
# Authors: Estephe Kana & Edi Prifti & Eugeni Belda
# Purpose: Build two figures
#            - Figure S1: 4-panel figure for pairwise correlation between degree 
#                         centrality, betweeness centrality and feature importance; 
#                         and a scatterplot3D showing all these metrics
#            - Figure S2: Get the top 20 MOTUs based on betweeness centrality and 
#                         degree centraily
# =================================================================================

# -----------------------------------------------------------------------------
# Commands to run the script (from repo root):
# Rscript analyses/figures/FigureS1_S2/FigureS1_S2_code.R
# -----------------------------------------------------------------------------

# Check for required packages
# Note: scatterplot3d, ggplotify — install.packages(c("scatterplot3d", "ggplotify"))
required_pkgs <- c("ggplot2", "dplyr", "tidyr", "patchwork", "igraph",
                   "scatterplot3d", "ggplotify")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(igraph)


# Derive repo root from script location
script_path  <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
script_dir   <- dirname(script_path)                      # <repo>/analyses/figures/FigureS1_S2
repo_root    <- dirname(dirname(dirname(script_dir)))     # <repo>/
data_dir     <- file.path(repo_root, "data")
analyses_dir <- file.path(repo_root, "analyses")
source(file.path(repo_root, "analyses", "scripts", "utils.R"))

out_dir <- script_dir

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
# MODULARITY DETECTION — fast greedy on undirected graph
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
  # Join habitat-level chi-sq assignment — rename inside select to avoid clash
  dplyr::left_join(
    sp.chisq_habitat %>% dplyr::select(feature, assigned_habitat = assigned_class),
    by = c("name" = "feature")
  ) %>%
  # Join zone-level chi-sq assignment — same pattern
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
  dplyr::select(feature, family, genus, Module, Zone, Habitat,
                betw_cent, degr_cent, mean_featureImportance, IsIndSp) %>%
  dplyr::mutate(
    family = ifelse(is.na(family), "Unk_family", family),
    genus  = ifelse(is.na(genus),  paste0(family, "_unk_genus"), genus)
  )

# ---------------------------------------------------------------------------
# 4-panel correlation figure
# ---------------------------------------------------------------------------
# install.packages(c("scatterplot3d", "ggplotify"))
library(scatterplot3d)
library(ggplotify)

scatter_df <- tree_df %>%
  dplyr::filter(!is.na(degr_cent), !is.na(betw_cent)) %>%
  dplyr::mutate(
    IsIndSp_lbl  = factor(ifelse(IsIndSp, "Indicator", "Non-indicator"),
                          levels = c("Indicator", "Non-indicator")),
    Habitat      = factor(Habitat, levels = names(habitat_pal)),
    log10_mfi    = log10(mean_featureImportance + 1)
  )

indsp_shapes <- c("Indicator" = 17, "Non-indicator" = 16)

# Spearman correlation label helper
.cor_label <- function(x, y) {
  ct <- cor.test(x, y, method = "spearman", exact = FALSE)
  sprintf("rho = %.2f,  p = %.2e", ct$estimate, ct$p.value)
}

corr_theme <- theme_bw(base_size = 12) +
  theme(legend.position  = "right",
        legend.title     = element_text(face = "bold"),
        plot.title       = element_text(face = "bold", size = 11))

# Panel A — degree vs mean feature importance
pA <- ggplot(scatter_df, aes(x = log10_mfi, y = degr_cent)) +
  geom_point(aes(colour = Habitat, shape = IsIndSp_lbl), alpha = 0.75, size = 2) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.7) +
  annotate("text", x = Inf, y = Inf,
           label = .cor_label(scatter_df$log10_mfi, scatter_df$degr_cent),
           hjust = 1.05, vjust = 1.6, size = 3.5, fontface = "italic") +
  scale_colour_manual("Habitat", values = habitat_pal) +
  scale_shape_manual("Species type", values = indsp_shapes) +
  labs(title = "A", x = "Log10(Mean Feature Importance)", y = "Degree Centrality") +
  corr_theme

# Panel B — betweenness vs mean feature importance
pB <- ggplot(scatter_df, aes(x = log10_mfi, y = betw_cent)) +
  geom_point(aes(colour = Habitat, shape = IsIndSp_lbl), alpha = 0.75, size = 2) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.7) +
  annotate("text", x = Inf, y = Inf,
           label = .cor_label(scatter_df$log10_mfi, scatter_df$betw_cent),
           hjust = 1.05, vjust = 1.6, size = 3.5, fontface = "italic") +
  scale_colour_manual("Habitat", values = habitat_pal) +
  scale_shape_manual("Species type", values = indsp_shapes) +
  labs(title = "B", x = "Log10(Mean Feature Importance)", y = "Betweenness Centrality") +
  corr_theme

# Panel C — degree vs betweenness
pC <- ggplot(scatter_df, aes(x = degr_cent, y = betw_cent)) +
  geom_point(aes(colour = Habitat, shape = IsIndSp_lbl), alpha = 0.75, size = 2) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.7) +
  annotate("text", x = Inf, y = Inf,
           label = .cor_label(scatter_df$degr_cent, scatter_df$betw_cent),
           hjust = 1.05, vjust = 1.6, size = 3.5, fontface = "italic") +
  scale_colour_manual("Habitat", values = habitat_pal) +
  scale_shape_manual("Species type", values = indsp_shapes) +
  labs(title = "C", x = "Degree Centrality", y = "Betweenness Centrality") +
  corr_theme

# Panel D — 3D scatter (converted to ggplot via ggplotify)
pD <- as.ggplot(function() {
  plot3d_df   <- scatter_df
  point_colors <- habitat_pal[as.character(plot3d_df$Habitat)]
  point_pch    <- indsp_shapes[as.character(plot3d_df$IsIndSp_lbl)]
  
  s3d <- scatterplot3d(
    x     = plot3d_df$log10_mfi,
    y     = plot3d_df$degr_cent,
    z     = plot3d_df$betw_cent,
    color = point_colors,
    pch   = point_pch,
    xlim  = range(plot3d_df$log10_mfi, na.rm = TRUE),
    ylim  = range(plot3d_df$degr_cent, na.rm = TRUE),
    zlim  = range(plot3d_df$betw_cent, na.rm = TRUE),
    xlab  = "Log10(Mean Feature Importance)",
    ylab  = "Degree Centrality",
    zlab  = "Betweenness Centrality",
    main  = "",
    grid  = TRUE, box = FALSE, cex.main = 1.1, font.main = 2
  )
  title(main = "D", adj = 0, font.main = 2, cex.main = 1, line = 1.8)
  # Label outlier points
  condition <- (plot3d_df$betw_cent > 0.9 | 
                  plot3d_df$log10_mfi > log10(3.62 + 1)) & 
    plot3d_df$degr_cent > 7
  
  outliers       <- plot3d_df[condition, ]
  outlier_coords <- s3d$xyz.convert(
    outliers$log10_mfi,
    outliers$degr_cent,
    outliers$betw_cent
  )
  text(outlier_coords$x, outlier_coords$y,
       labels = outliers$feature,
       pos    = 4,
       cex    = 0.6,
       font   = 2,
       col    = "firebrick")
})

four_panel <- (pA | pB) / (pC | pD) +
  patchwork::plot_layout(guides = "collect",
                         widths  = c(0.6, 1.4),
                         heights = c(1, 1)) &
  theme(legend.position = "right")

ggsave(
  filename = file.path(out_dir, "FigureS1.pdf"),
  plot     = four_panel,
  width    = 14, height = 10, units = "in",
  device   = cairo_pdf
)
message("4-panel correlation figure saved.")

# ---------------------------------------------------------------------------
# Top-20 centrality bubble plots (lollipop charts)
# ---------------------------------------------------------------------------

plot_degr_df <- tree_df %>%
  dplyr::arrange(desc(degr_cent)) %>%
  dplyr::mutate(
    degr_rank   = dplyr::row_number(),
    IsIndSp_lbl = ifelse(IsIndSp, "Indicator", "Non-indicator"),
    Habitat     = factor(Habitat, levels = names(habitat_pal))
  ) %>%
  dplyr::filter(degr_rank <= 20) %>%
  dplyr::mutate(feature = factor(feature, levels = rev(feature)))

plot_betw_df <- tree_df %>%
  dplyr::arrange(desc(betw_cent)) %>%
  dplyr::mutate(
    betw_rank   = dplyr::row_number(),
    IsIndSp_lbl = ifelse(IsIndSp, "Indicator", "Non-indicator"),
    Habitat     = factor(Habitat, levels = names(habitat_pal))
  ) %>%
  dplyr::filter(betw_rank <= 20) %>%
  dplyr::mutate(feature = factor(feature, levels = rev(feature)))

plot_degr <- ggplot(plot_degr_df, aes(x = degr_cent, y = feature)) +
  geom_segment(aes(x = 0, xend = degr_cent, yend = feature),
               color = "grey70", linewidth = 0.5) +
  geom_point(aes(fill = Habitat, color = IsIndSp_lbl),
             shape = 21, size = 5, stroke = 2) +
  scale_fill_manual("Habitat",
                    values = habitat_pal,
                    guide  = guide_legend(order = 1, override.aes = list(size = 5, color = "grey30"))) +
  scale_color_manual("Species type",
                     values = c("Indicator" = "orange", "Non-indicator" = "black"),
                     guide  = guide_legend(order = 2, override.aes = list(size = 5, fill = "grey70"))) +
  labs(title = "Top 20 \u2014 Degree centrality",
       x = "Degree centrality", y = NULL) +
  theme_bw(base_size = 12) +
  theme(plot.title      = element_text(face = "bold"),
        legend.position = "bottom",
        legend.box      = "horizontal",
        legend.key      = element_rect(colour = NA))

plot_betw <- ggplot(plot_betw_df, aes(x = betw_cent, y = feature)) +
  geom_segment(aes(x = 0, xend = betw_cent, yend = feature),
               color = "grey70", linewidth = 0.5) +
  geom_point(aes(fill = Habitat, color = IsIndSp_lbl),
             shape = 21, size = 5, stroke = 2) +
  scale_fill_manual("Habitat",
                    values = habitat_pal,
                    guide  = guide_legend(order = 1, override.aes = list(size = 5, color = "grey30"))) +
  scale_color_manual("Species type",
                     values = c("Indicator" = "orange", "Non-indicator" = "black"),
                     guide  = guide_legend(order = 2, override.aes = list(size = 5, fill = "grey70"))) +
  labs(title = "Top 20 \u2014 Betweenness centrality",
       x = "Betweenness centrality", y = NULL) +
  theme_bw(base_size = 12) +
  theme(plot.title      = element_text(face = "bold"),
        legend.position = "bottom",
        legend.box      = "horizontal",
        legend.key      = element_rect(colour = NA))

bubble_combined <- (plot_degr | plot_betw) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.box      = "horizontal")

ggsave(
  filename = file.path(out_dir, "FigureS2.pdf"),
  plot     = bubble_combined,
  width    = 14, height = 7, units = "in",
  device   = cairo_pdf
)
message("Top-20 bubble plot saved.")