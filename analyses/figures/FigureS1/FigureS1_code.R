# =============================================================================
# Script name: FigureS1_code.R
# Authors: Estephe Kana & Edi Prifti & Eugeni Belda
# Date created: 2025-12-10
# Purpose: Build a multi-panel sensitivity figure across three ecorr thresholds:
#   Row 1 — Panels A/B/C: Zone-coloured co-presence networks at ecorr > 0.3/0.5/0.7
#   Row 2 — Panels D/E/F: Indicator vs non-indicator centrality boxplots per threshold
#   Row 3 — Panel G:      Network topology metrics summary (sensitivity analysis)
#   Row 4 — Panels H/I/J: Zone → Module → Habitat alluvial diagrams per threshold
# Inputs:  analyses/files/rdata/graph_data/graph_data_ecorr50_all_strat_*.rda
#          analyses/analysis_outputs/{bininter,terinter}_output_data/*_prev_3.Rda
# Outputs: analyses/figures/FigureS1/FigureS1.pdf
# =============================================================================

# -----------------------------------------------------------------------------
# Commands to run the script (from repo root):
#   export REPO_ROOT="/data/projects/aime/analyses/seamount/metabarcoding/"
#
#   From REPO_ROOT, run the following commands:
#
#   Rscript analyses/figures/FigureS1/FigureS1_code.R
# -----------------------------------------------------------------------------

# Check for required packages
required_pkgs <- c("igraph", "ggplot2", "dplyr", "tidyr", "ggpubr",
                   "cowplot", "ggalluvial", "gridExtra")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))

library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(cowplot)
library(ggalluvial)
library(gridExtra)

# Derive repo root from script location
script_path  <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
script_dir   <- dirname(script_path)                      # <repo>/analyses/figures/FigureS1
repo_root    <- dirname(dirname(dirname(script_dir)))     # <repo>/

# Allow override via environment variable (set on the remote server)
repo_root <- Sys.getenv("REPO_ROOT", unset = repo_root)

analyses_dir <- file.path(repo_root, "analyses")

out_dir         <- script_dir
rda_files       <- list.files(file.path(analyses_dir, "files", "rdata", "graph_data"),
                               pattern = "^graph_data_ecorr50_all_strat_.*\\.rda$", full.names = TRUE)
graph_data_path <- rda_files[which.max(file.mtime(rda_files))]

# =============================================================================
# PALETTES
# =============================================================================

zone_pal <- c(
  "Shallow" = "#5ae6ab",
  "Middle"  = "#ffe699",
  "Deep"    = "#25456B",
  "NS"      = "gray80"
)

alluvial_habitat_colors <- c(
  "Bay"              = "#5ae6ab",
  "Lagoon"           = "#88d941",
  "BackReef"   = "#2e7d00",
  "OuterSlope" = "#4e8273",
  "Seamount50"         = "#ffe699",
  "DeepSlope150"     = "#d79c3b",
  "Seamount250"        = "#4a6a94",
  "Seamount500"        = "#1a3250",
  "NS.habitat"       = "gray"
)

# Fixed cluster orderings per threshold (used for consistent y-axis alignment)
cluster_order_30 <- c("Cluster_08", "Cluster_05", "Cluster_03", "Cluster_02",
                      "Cluster_06", "Cluster_07", "Cluster_04", "Cluster_01", "Cluster_09")

cluster_order_50 <- c("Cluster_08", "Cluster_01", "Cluster_02", "Cluster_03",
                      "Cluster_04", "Cluster_05", "Cluster_09", "Cluster_06",
                      "Cluster_10", "Cluster_12", "Cluster_07", "Cluster_11",
                      "Cluster_13")

cluster_order_70 <- c("Cluster_07", "Cluster_10", "Cluster_13", "Cluster_14",
                      "Cluster_15", "Cluster_19", "Cluster_06", "Cluster_01",
                      "Cluster_02", "Cluster_05", "Cluster_09", "Cluster_11",
                      "Cluster_03", "Cluster_08", "Cluster_16", "Cluster_20",
                      "Cluster_21", "Cluster_12", "Cluster_17", "Cluster_18")

# Cluster order for ecorr > 0.6 (was 0.55; reverted) is left empty;
# the alluvial panel
# helper computes it dynamically via intersect(cluster_order,
# present_modules) U sort(extras), so passing an empty seed simply
# falls back to the sorted natural order of the modules detected at
# this threshold.
cluster_order_60 <- character(0)

# =============================================================================
# LOAD DATA
# =============================================================================

load(graph_data_path)   # -> edges.all, network, lay, nodes.annot,
#    sp.chisq_posthoc, sp.chisq_posthoc.habitat

# Pre-compute Figure 5 cluster assignments from the loaded ecorr > 0.5 network
fg50_ref    <- igraph::cluster_fast_greedy(igraph::as_undirected(network, mode = "collapse"))
fg50_ref_df <- data.frame(
  name   = igraph::V(network)$name,
  Module = paste0("Cluster_", sprintf("%02d", fg50_ref$membership)),
  stringsAsFactors = FALSE
)

# =============================================================================
# FIX CHI-SQ TABLES  (mirrors Figure 5 logic)
# =============================================================================

sp.chisq_habitat <- sp.chisq_posthoc.habitat
cols_to_rename   <- c("pval_chisq", "padj_chisq", "mean_presence_count",
                      "enriched_zones", "depleted_zones", "assigned_class")
colnames(sp.chisq_habitat)[colnames(sp.chisq_habitat) %in%
                             paste0(cols_to_rename, ".habitat")] <- cols_to_rename
sp.chisq_habitat$assigned_class[
  sp.chisq_habitat$assigned_class == "DeepSlope"] <- "DeepSlope150"
colnames(sp.chisq_habitat) <- gsub("DeepSlope", "DeepSlope150",
                                   colnames(sp.chisq_habitat))
sp.chisq_zone <- sp.chisq_posthoc

zone_df <- sp.chisq_zone %>%
  dplyr::transmute(
    name = feature,
    Zone = dplyr::case_when(
      !is.na(padj_chisq) & padj_chisq < 0.05 & assigned_class != "NS" ~ assigned_class,
      TRUE ~ "NS"
    )
  )

habitat_df <- sp.chisq_habitat %>%
  dplyr::transmute(
    name    = feature,
    Habitat = dplyr::case_when(
      !is.na(padj_chisq) & padj_chisq < 0.05 & assigned_class != "NS" ~ assigned_class,
      TRUE ~ "NS.habitat"
    )
  )

# =============================================================================
# PREDOMICS — IsIndSp flag per feature
# =============================================================================

objlist <- list()

x <- load(file.path(analyses_dir, "analysis_outputs", "bininter_output_data", "bininter_Predomics_all_analyses_overall_data_strat_group_prev_3.Rda"))
for (i in x) objlist[["bininter"]][[i]] <- get(i)
rm(list = x)

x <- load(file.path(analyses_dir, "analysis_outputs", "terinter_output_data", "terinter_Predomics_all_analyses_overall_data_strat_group_prev_3.Rda"))
for (i in x) objlist[["terinter"]][[i]] <- get(i)
rm(list = x)

feat_summary <- rbind(objlist$bininter$predout.bin.sub,
                      objlist$terinter$predout.bin.sub) %>%
  dplyr::group_by(feature) %>%
  dplyr::summarise(
    IsIndSp = any(IsIndSp == 1, na.rm = TRUE),
    .groups = "drop"
  )

# =============================================================================
# HELPER — build network + extract node metrics at a given ecorr threshold
# =============================================================================

build_threshold_network <- function(edges_all, threshold, zone_df, habitat_df,
                                    feat_summary) {
  edges_filt <- edges_all[abs(edges_all$ecorr) > threshold, ]
  if (nrow(edges_filt) == 0)
    stop(sprintf("No edges survive ecorr threshold = %.2f", threshold))
  
  node_names <- unique(c(edges_filt$from, edges_filt$to))
  
  g_dir <- igraph::graph_from_data_frame(
    d        = edges_filt,
    directed = TRUE,
    vertices = data.frame(name = node_names)
  )
  g  <- igraph::as_undirected(g_dir, mode = "collapse")
  fg <- igraph::cluster_fast_greedy(g)
  
  V(g)$module  <- paste0("Cluster_", sprintf("%02d", fg$membership))
  V(g)$Zone    <- zone_df$Zone[match(V(g)$name, zone_df$name)]
  V(g)$Habitat <- habitat_df$Habitat[match(V(g)$name, habitat_df$name)]
  V(g)$Zone[is.na(V(g)$Zone)]       <- "NS"
  V(g)$Habitat[is.na(V(g)$Habitat)] <- "NS.habitat"
  
  betAll      <- igraph::betweenness(g, directed = FALSE) /
    (((igraph::vcount(g) - 1) * (igraph::vcount(g) - 2)) / 2)
  betAll.norm <- (betAll - min(betAll)) / (max(betAll) - min(betAll))
  
  node_df <- data.frame(
    name      = V(g)$name,
    degr_cent = igraph::degree(g),
    betw_cent = betAll.norm,
    Zone      = V(g)$Zone,
    Habitat   = V(g)$Habitat,
    Module    = V(g)$module,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(
      feat_summary %>% dplyr::select(feature, IsIndSp),
      by = c("name" = "feature")
    ) %>%
    dplyr::mutate(
      IsIndSp = tidyr::replace_na(IsIndSp, FALSE)
    )
  
  list(
    graph      = g,
    fg         = fg,
    node_df    = node_df,
    n_nodes    = igraph::vcount(g),
    n_edges    = igraph::ecount(g),
    modularity = igraph::modularity(fg),
    n_modules  = length(unique(fg$membership))
  )
}

# =============================================================================
# BUILD THE THREE NETWORKS
# =============================================================================

message("Building network at ecorr > 0.3 ...")
net30 <- build_threshold_network(edges.all, 0.3, zone_df, habitat_df, feat_summary)

message("Building network at ecorr > 0.5 ...")
net50 <- build_threshold_network(edges.all, 0.5, zone_df, habitat_df, feat_summary)

# Override ecorr > 0.5 cluster assignments with Figure 5 reference
net50$node_df <- net50$node_df %>%
  dplyr::select(-Module) %>%
  dplyr::left_join(fg50_ref_df, by = "name") %>%
  dplyr::mutate(Module = tidyr::replace_na(Module, "NS"))
net50$modularity <- igraph::modularity(fg50_ref)
net50$n_modules  <- length(unique(fg50_ref_df$Module))

message("Building network at ecorr > 0.6 ...")
net60 <- build_threshold_network(edges.all, 0.6, zone_df, habitat_df, feat_summary)

message("Building network at ecorr > 0.7 ...")
net70 <- build_threshold_network(edges.all, 0.7, zone_df, habitat_df, feat_summary)

# =============================================================================
# PANELS A / B / C — Zone-coloured network plots
# =============================================================================

make_network_panel <- function(net_obj, threshold_label, panel_label,
                               zone_pal, sp.chisq_zone, layout_seed = 100) {
  g  <- net_obj$graph
  df <- net_obj$node_df
  
  node_color  <- zone_pal[df$Zone]
  node_color[is.na(node_color)] <- "gray80"
  
  padj_vals <- sp.chisq_zone$padj_chisq[match(df$name, sp.chisq_zone$feature)]
  node_size <- ifelse(!is.na(padj_vals) & padj_vals < 0.05,
                      pmin(-log10(padj_vals) * 1.5, 12), 4)
  
  frame_col   <- ifelse(df$IsIndSp, "firebrick1", "gray50")
  frame_width <- ifelse(df$IsIndSp, 2.5, 1)
  
  set.seed(layout_seed)
  lo <- igraph::layout_with_fr(g)
  
  cowplot::as_grob(function() {
    par(mar = c(0, 0.5, 3, 0.5))
    
    plot(g,
         layout             = lo,
         vertex.label       = NA,
         vertex.color       = node_color,
         vertex.size        = node_size,
         vertex.frame.color = frame_col,
         vertex.frame.width = frame_width,
         edge.color         = "gray70",
         edge.width         = 0.4,
         edge.arrow.size    = 0,
         asp                = FALSE,
         rescale            = TRUE)
    
    mtext(panel_label,
          side = 3, line = 0.5, adj = 0.02, cex = 3, font = 2)
    # mtext(
    #   sprintf("ecorr > %s  |  %d nodes  |  %d edges  |  modularity_score = %.3f | %d modules",
    #           threshold_label, net_obj$n_nodes, net_obj$n_edges, net_obj$modularity, net_obj$n_modules),
    #   side = 3, line = -1, adj = 0.5, cex = 1.25, col = "gray30"
    # )
  })
}

## Per-threshold panel rename to match the new 4-row x 4-column layout:
##   row 1 (ecorr > 0.3)  -> A network | B centrality | C alluvial | D degree
##   row 2 (ecorr > 0.5)  -> E         | F            | G          | H
##   row 3 (ecorr > 0.6) -> I         | J            | K          | L
##   row 4 (ecorr > 0.7)  -> M         | N            | O          | P
##
## 0.6 is used as the 4th threshold for the round-number sensitivity
## scan readable progression (0.3, 0.5, 0.6, 0.7). A fine-grained scan
## (0.30-0.85 step 0.05) showed that the best scale-free fit (highest
## log-log R^2 = 0.91, gamma = 2.09) is reached at ecorr ~ 0.55, just
## inside the 0.5-0.6 interval; the manuscript caption states this
## explicitly so a reader doesn't have to interpolate it from the
## panels.
panel_A <- make_network_panel(net30, "0.3",  "A", zone_pal, sp.chisq_zone)
panel_E <- make_network_panel(net50, "0.5",  "E", zone_pal, sp.chisq_zone)
panel_I <- make_network_panel(net60, "0.6", "I", zone_pal, sp.chisq_zone)
panel_M <- make_network_panel(net70, "0.7",  "M", zone_pal, sp.chisq_zone)

# Standalone zone and habitat legends
## Legends now live in a single shared bottom row. Each legend uses
## guide_legend(title.position = "top", title.hjust = 0) so the
## title 'Zone' / 'Indicator status' / 'Habitat' sits on a line of
## its own above the colour chips, rather than next to them where it
## was overlapping the chips at the chosen font sizes.
zone_legend_plot <- ggplot(
  data.frame(Zone = factor(names(zone_pal), levels = names(zone_pal))),
  aes(x = 1, y = Zone, fill = Zone)
) +
  geom_tile() +
  scale_fill_manual("Zone", values = zone_pal) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 28, face = "bold"),
    legend.text     = element_text(size = 24),
    legend.key.size = unit(1.2, "cm"),
    legend.margin   = margin(t = 10, b = 10)
  ) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0))

## Extract the legend grob directly from the plot's gtable. cowplot's
## get_legend / get_plot_component were inconsistent across versions
## here -- some returned an empty zeroGrob, others returned the chip
## strip without the title text. Walking the gtable layout for any
## name matching 'guide-box*' is the robust path.
get_legend_safe <- function(p) {
  g  <- ggplot2::ggplotGrob(p)
  ix <- grep("guide-box", g$layout$name)
  if (length(ix) == 0) return(grid::nullGrob())
  for (i in ix) {
    candidate <- g$grobs[[i]]
    if (!inherits(candidate, "zeroGrob")) return(candidate)
  }
  g$grobs[[ix[1]]]
}

zone_legend_grob <- get_legend_safe(zone_legend_plot)

habitat_legend_plot <- ggplot(
  data.frame(Habitat = factor(names(alluvial_habitat_colors), levels = names(alluvial_habitat_colors))),
  aes(x = 1, y = Habitat, fill = Habitat)
) +
  geom_tile() +
  scale_fill_manual("Habitat", values = alluvial_habitat_colors) +
  theme_void() +
  theme(
    legend.position      = "bottom",
    legend.justification = "left",  # left-align in its row-allocated slot
    legend.title         = element_text(size = 28, face = "bold"),
    legend.text          = element_text(size = 24),
    legend.key.size      = unit(1.2, "cm"),
    legend.margin        = margin(t = 10, b = 10)
  ) +
  guides(fill = guide_legend(ncol = 3, title.position = "top", title.hjust = 0))
habitat_legend_grob <- get_legend_safe(habitat_legend_plot)

# =============================================================================
# PANELS D / E / F — Centrality boxplots: indicator vs non-indicator per threshold
# =============================================================================

make_boxplot_panel <- function(net_obj, threshold_label, panel_label, show_legend = FALSE) {
  metrics_long <- net_obj$node_df %>%
    dplyr::mutate(
      IsIndSp = factor(
        ifelse(IsIndSp, "Indicator", "Non-Indicator"),
        levels = c("Non-Indicator", "Indicator")
      )
    ) %>%
    dplyr::select(name, IsIndSp, degr_cent, betw_cent) %>%
    tidyr::pivot_longer(
      cols      = c(degr_cent, betw_cent),
      names_to  = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      metric = dplyr::recode(metric,
                             degr_cent = "Degree",
                             betw_cent = "Betweenness")
    )
  
  p <- ggplot(metrics_long, aes(x = IsIndSp, y = value, fill = IsIndSp)) +
    geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA,
                 fill = "white", color = "grey30", alpha = 0.8) +
    geom_jitter(aes(color = IsIndSp), width = 0.08, size = 1.2,
                alpha = 0.5, show.legend = FALSE) +
    ggpubr::stat_compare_means(
      method      = "wilcox.test",
      label       = "p.signif",
      label.x     = 1.5,
      ## label.y.npc anchors the significance mark to 90% of the plot
      ## panel's vertical extent so it sits inside the data area (high
      ## up but clearly below the facet strip band that names
      ## 'Betweenness' / 'Degree'). The previous vjust = -0.5 pushed
      ## the mark right against the top edge, so '***' / 'ns' was
      ## visually touching or sitting under the strip text.
      label.y.npc = 0.90,
      size        = 12
    ) +
    facet_wrap(~ metric, scales = "free_y", ncol = 2) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    scale_fill_manual(values  = c("Non-Indicator" = "#7fbfff", "Indicator" = "#ff7f7f")) +
    scale_color_manual(values = c("Non-Indicator" = "#3a7fc1", "Indicator" = "#c13a3a")) +
    labs(
      x        = NULL,
      y        = "Centrality Value",
      fill     = NULL,
      title    = panel_label,
    ) +
    theme_bw(base_size = 26) +
    theme(
      strip.text       = element_text(face = "bold", size = 26),
      axis.text.x      = element_text(size = 18, angle = 30, hjust = 1),
      axis.text.y      = element_text(size = 22),
      axis.title       = element_text(size = 26),
      legend.position  = "bottom",
      legend.text      = element_text(size = 22),
      legend.title     = element_text(size = 24),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", hjust = 0, size = 32),
    )
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

## Centrality violins are column 2 of each row:
##   B (ecorr 0.3) / F (ecorr 0.5) / J (ecorr 0.6) / N (ecorr 0.7).
panel_B <- make_boxplot_panel(net30, "0.3",  "B", show_legend = FALSE)
panel_F <- make_boxplot_panel(net50, "0.5",  "F", show_legend = FALSE)
panel_J <- make_boxplot_panel(net60, "0.6", "J", show_legend = FALSE)
panel_N <- make_boxplot_panel(net70, "0.7",  "N", show_legend = FALSE)

indic_legend_plot <- ggplot(
  data.frame(IsIndSp = factor(c("Non-Indicator", "Indicator"),
                              levels = c("Non-Indicator", "Indicator"))),
  aes(x = 1, y = IsIndSp, fill = IsIndSp)
) +
  geom_tile() +
  scale_fill_manual("Indicator status",
                    values = c("Non-Indicator" = "#7fbfff", "Indicator" = "#ff7f7f")) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 28, face = "bold"),
    legend.text     = element_text(size = 24),
    legend.key.size = unit(1.2, "cm"),
    legend.margin   = margin(t = 10, b = 10)
  ) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0))
indic_legend_grob <- get_legend_safe(indic_legend_plot)

# =============================================================================
# PANEL G
# =============================================================================

make_topology_panel <- function(net_list, threshold_labels, panel_label) {
  
  topo_df <- do.call(rbind, lapply(threshold_labels, function(thr) {
    obj <- net_list[[thr]]
    data.frame(
      threshold      = thr,
      Nodes          = obj$n_nodes,
      Edges          = obj$n_edges,
      Modules        = obj$n_modules,
      Modularity = round(obj$modularity, 3),
      check.names    = FALSE,
      stringsAsFactors = FALSE
    )
  })) %>%
    dplyr::mutate(threshold = factor(threshold, levels = threshold_labels)) %>%
    tidyr::pivot_longer(
      cols      = c(Nodes, Edges, Modules, `Modularity`),
      names_to  = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      metric = factor(metric, levels = c("Nodes", "Edges", "Modules", "Modularity"))
    )
  
  ggplot(topo_df, aes(x = threshold, y = value, group = 1)) +
    geom_line(color = "grey60", linewidth = 0.9, linetype = "dashed") +
    geom_point(fill="grey70" ,shape = 21, size = 7, color = "grey30", stroke = 1) +
    geom_text(
      aes(label = ifelse(metric == "Modularity",
                         sprintf("%.3f", value),
                         as.character(as.integer(value)))),
      vjust = -1.4, size = 8, fontface = "bold", color = "grey20"
    ) +
    facet_wrap(~ metric, scales = "free_y", ncol = 4) +
    ## Top expansion bumped 0.25 -> 0.55 because the per-point value
    ## labels (vjust = -1.4) were getting clipped by the facet strip
    ## band on the highest data point (e.g. 318 nodes at ecorr > 0.3).
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.55))) +
    labs(
      x     = "ecorr_threshold",
      y     = "Topology_value",
      title = panel_label
    ) +
    theme_bw(base_size = 26) +
    theme(
      strip.text       = element_text(face = "bold", size = 26),
      strip.background = element_rect(fill = "grey92"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x      = element_text(face = "bold", size = 22),
      axis.text.y      = element_text(size = 20),
      axis.title       = element_text(size = 26),
      plot.title       = element_text(face = "bold", hjust = 0, size = 32),
      plot.margin      = margin(t = 5, r = 15, b = 10, l = 5)
    )
}

## Topology summary panel reinstated as panel Q at the bottom of the
## figure (just above the legend row). The per-threshold degree
## distributions (D / H / L / P) carry the qualitative story, but the
## summary line plot is a useful single-glance reference showing how
## the four scalar topology metrics (n_nodes / n_edges / n_modules /
## modularity) co-vary with the ecorr threshold across the four sampled
## values.
panel_Q <- make_topology_panel(
  net_list         = list("0.3"  = net30, "0.5"  = net50,
                          "0.6" = net60, "0.7"  = net70),
  threshold_labels = c("0.3", "0.5", "0.6", "0.7"),
  panel_label      = "Q"
)

# =============================================================================
# PANELS C / G / K — Zone → Module → Habitat alluvial per threshold
# =============================================================================

make_alluvial_panel <- function(net_obj, threshold_label, panel_label,
                                sp.chisq_zone, sp.chisq_habitat,
                                cluster_order, show_legend = TRUE) {
  
  zone_order    <- c("Shallow", "Middle", "Deep", "NS.zone")
  habitat_order <- c("Bay", "Lagoon", "BackReef", "OuterSlope",
                     "Seamount50", "DeepSlope150", "Seamount250", "Seamount500", "NS.habitat")
  
  present_modules <- unique(net_obj$node_df$Module)
  extra_modules   <- sort(setdiff(present_modules, cluster_order))
  module_order    <- c(intersect(cluster_order, present_modules), extra_modules)
  
  alluvial_df <- net_obj$node_df %>%
    dplyr::select(name, Module) %>%
    dplyr::left_join(
      sp.chisq_habitat %>% dplyr::select(feature, assigned_habitat = assigned_class),
      by = c("name" = "feature")
    ) %>%
    dplyr::left_join(
      sp.chisq_zone %>% dplyr::select(feature, assigned_zone = assigned_class),
      by = c("name" = "feature")
    ) %>%
    dplyr::mutate(
      Habitat = dplyr::case_when(
        !is.na(assigned_habitat) & assigned_habitat != "NS" ~ assigned_habitat,
        TRUE ~ "NS.habitat"
      ),
      Zone = dplyr::case_when(
        !is.na(assigned_zone) & assigned_zone != "NS" ~ assigned_zone,
        TRUE ~ "NS.zone"
      )
    ) %>%
    dplyr::count(Zone, Module, Habitat, name = "n_species") %>%
    dplyr::mutate(
      Zone    = factor(Zone,    levels = zone_order),
      Module  = factor(Module,  levels = module_order),
      Habitat = factor(Habitat, levels = habitat_order)
    )
  
  ggplot(
    alluvial_df,
    aes(axis1 = Zone, axis2 = Module, axis3 = Habitat, y = n_species)
  ) +
    ggalluvial::geom_alluvium(
      aes(fill = Habitat),
      width = 0.25, alpha = 0.75, knot.pos = 0.4, curve_type = "sigmoid"
    ) +
    ggalluvial::geom_stratum(
      width = 0.25, fill = "grey92", color = "grey40", linewidth = 0.3
    ) +
    ggplot2::geom_text(
      stat          = ggalluvial::StatStratum,
      aes(label     = sub("^Cluster_", "", as.character(after_stat(stratum)))),
      size          = 6,
      lineheight    = 1.3,
      fontface      = "bold",
      color         = "grey20",
      check_overlap = TRUE
    ) +
    scale_x_discrete(
      limits = c("Zone", "Module", "Habitat"),
      expand = c(0.08, 0.01)
    ) +
    scale_fill_manual(
      values   = alluvial_habitat_colors[habitat_order],
      name     = "Habitat",
      na.value = "grey70",
      drop     = FALSE
    ) +
    scale_y_continuous(name = "MOTUs", expand = c(0, 0)) +
    theme_minimal(base_size = 24) +
    theme(
      axis.text.x     = element_text(size = 24, face = "bold"),
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      axis.title.x    = element_blank(),
      axis.title.y    = element_text(size = 24, face = "bold"),
      panel.grid      = element_blank(),
      legend.title    = element_text(size = 26, face = "bold"),
      legend.text     = element_text(size = 22),
      legend.key.size = unit(1.1, "cm"),
      legend.position = if (show_legend) "right" else "none",
      plot.title      = element_text(face = "bold", hjust = 0, size = 32),
      plot.margin     = margin(t = 10, r = 10, b = 10, l = 10)
    ) +
    labs(
      title    = panel_label,
      # subtitle = sprintf("ecorr > %s  |  %d nodes  |  %d modules",
      #                    threshold_label, net_obj$n_nodes, net_obj$n_modules)
    )
}

## Alluvials are column 3 of each row:
##   C (ecorr 0.3) / G (ecorr 0.5) / K (ecorr 0.6) / O (ecorr 0.7).
panel_C <- make_alluvial_panel(net30, "0.3",  "C", sp.chisq_zone, sp.chisq_habitat,
                               cluster_order_30, show_legend = FALSE)
panel_G <- make_alluvial_panel(net50, "0.5",  "G", sp.chisq_zone, sp.chisq_habitat,
                               cluster_order_50, show_legend = FALSE)
panel_K <- make_alluvial_panel(net60, "0.6", "K", sp.chisq_zone, sp.chisq_habitat,
                               cluster_order_60, show_legend = FALSE)
panel_O <- make_alluvial_panel(net70, "0.7",  "O", sp.chisq_zone, sp.chisq_habitat,
                               cluster_order_70, show_legend = FALSE)

# =============================================================================
# PANELS D / H / L — Degree distribution (scale-free check) per threshold
# =============================================================================

make_degree_dist_panel <- function(net_obj, threshold_label, panel_label) {

  deg     <- igraph::degree(net_obj$graph)
  deg_tab <- as.data.frame(table(deg), stringsAsFactors = FALSE)
  deg_tab$k    <- as.integer(deg_tab$deg)
  deg_tab$P    <- deg_tab$Freq / sum(deg_tab$Freq)

  ## Power-law reference line: fit log(P) = -gamma * log(k) + c on the
  ## right tail (k >= 2) and overlay as a dashed grey line so a reader
  ## can eye-check how close the empirical points lie to the
  ## scale-free expectation at each threshold.
  fit_df <- deg_tab[deg_tab$k >= 2 & deg_tab$P > 0, ]
  if (nrow(fit_df) >= 3) {
    fit   <- lm(log10(P) ~ log10(k), data = fit_df)
    gamma <- -coef(fit)[["log10(k)"]]
    intercept_c <- coef(fit)[["(Intercept)"]]
    line_df <- data.frame(
      k = 10 ^ seq(log10(min(deg_tab$k[deg_tab$k > 0])),
                   log10(max(deg_tab$k)),
                   length.out = 100)
    )
    line_df$P <- 10 ^ (intercept_c + log10(line_df$k) * coef(fit)[["log10(k)"]])
    gamma_label <- sprintf("gamma ~ %.2f", gamma)
  } else {
    line_df     <- NULL
    gamma_label <- ""
  }

  p <- ggplot(deg_tab, aes(x = k, y = P)) +
    geom_point(size = 4, color = "#1F4E79", alpha = 0.8) +
    scale_x_log10() +
    scale_y_log10()

  if (!is.null(line_df)) {
    p <- p + geom_line(data = line_df, aes(x = k, y = P),
                       linetype = "dashed", color = "grey50", linewidth = 0.7)
  }

  p +
    annotate("text",
             x     = max(deg_tab$k),
             y     = max(deg_tab$P),
             label = gamma_label,
             hjust = 1, vjust = 1, size = 9, fontface = "bold",
             color = "grey25") +
    labs(
      x     = "Degree k (log scale)",
      y     = "P(k) (log scale)",
      title = panel_label
    ) +
    theme_bw(base_size = 26) +
    theme(
      panel.grid.minor   = element_blank(),
      axis.text          = element_text(size = 22),
      axis.title         = element_text(size = 26),
      plot.title         = element_text(face = "bold", hjust = 0, size = 32),
      plot.margin        = margin(t = 10, r = 15, b = 10, l = 10)
    )
}

panel_D <- make_degree_dist_panel(net30, "0.3",  "D")
panel_H <- make_degree_dist_panel(net50, "0.5",  "H")
panel_L <- make_degree_dist_panel(net60, "0.6", "L")
panel_P <- make_degree_dist_panel(net70, "0.7",  "P")

# =============================================================================
# COMPOSE & SAVE
# =============================================================================

message("Composing figure ...")

## Layout: 4 rows (one per ecorr threshold) of 4 data panels each,
## then a topology summary line plot Q spanning the full width, then
## a final bottom row with all three legends consolidated.
##
##   row 1 (ecorr > 0.3): A network | B centrality | C alluvial | D degree
##   row 2 (ecorr > 0.5): E         | F            | G          | H
##   row 3 (ecorr > 0.6): I        | J            | K          | L
##   row 4 (ecorr > 0.7): M         | N            | O          | P
##   row 5 (legends):     zone      | indicator    | habitat (2 cols wide)
##   row 6 (summary):     Q (nodes / edges / modules / modularity vs threshold)

row1 <- gridExtra::arrangeGrob(
  panel_A,
  ggplotGrob(panel_B),
  ggplotGrob(panel_C),
  ggplotGrob(panel_D),
  ncol   = 4,
  widths = c(1, 1, 1, 1)
)

row2 <- gridExtra::arrangeGrob(
  panel_E,
  ggplotGrob(panel_F),
  ggplotGrob(panel_G),
  ggplotGrob(panel_H),
  ncol   = 4,
  widths = c(1, 1, 1, 1)
)

row3 <- gridExtra::arrangeGrob(
  panel_I,
  ggplotGrob(panel_J),
  ggplotGrob(panel_K),
  ggplotGrob(panel_L),
  ncol   = 4,
  widths = c(1, 1, 1, 1)
)

row4 <- gridExtra::arrangeGrob(
  panel_M,
  ggplotGrob(panel_N),
  ggplotGrob(panel_O),
  ggplotGrob(panel_P),
  ncol   = 4,
  widths = c(1, 1, 1, 1)
)

## Bottom legend row: habitat gets two slots because it has 9 entries
## (laid out 3 cols x 3 rows in the legend) while zone (4 entries)
## and indicator (2 entries) fit comfortably in one slot each.
legend_row <- gridExtra::arrangeGrob(
  zone_legend_grob,
  indic_legend_grob,
  habitat_legend_grob,
  ncol   = 3,
  widths = c(1, 1, 2)
)

## Canvas stretched vertically to fit the 4 data rows + legend row +
## Q summary at the bottom: 28 wide x 31 tall.
pdf(file.path(out_dir, "FigureS1.pdf"),
    width = 28, height = 31)
gridExtra::grid.arrange(row1, row2, row3, row4,
                        legend_row,
                        ggplotGrob(panel_Q),
                        nrow    = 6,
                        heights = c(1, 1, 1, 1, 0.32, 0.75))
dev.off()

message("Done — FigureS1 saved to: ", out_dir)

