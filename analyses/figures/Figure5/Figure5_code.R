# =============================================================================
# Script name: Figure5_code.R
# Authors: Estephe Kana & Edi Prifti & Eugeni Belda
# Purpose: Build a 4-panel publication figure:
#   Panel A — Co-occurrence network coloured by habitat
#   Panel B — Same network coloured by fast-greedy module membership
#   Panel C — Alluvial diagram: Habitat → Module membership
#   Panel D — Module enrichment heatmap (Zone + Habitat, enrichment only)
# =============================================================================

# -----------------------------------------------------------------------------
# Commands to run the script (from repo root):
# Rscript analyses/figures/Figure5/Figure5_code.R
# -----------------------------------------------------------------------------

# Check for required packages
# Note: piano — BiocManager::install("piano")
required_pkgs <- c("igraph", "ggplot2", "dplyr", "tidyr", "cowplot",
                   "ggalluvial", "piano", "gridExtra")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))

library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggalluvial)
library(piano)
library(gridExtra)

# =============================================================================
# PATHS
# =============================================================================

# Derive repo root from script location
script_path  <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
script_dir   <- dirname(script_path)                      # <repo>/analyses/figures/Figure5
repo_root    <- dirname(dirname(dirname(script_dir)))     # <repo>/
data_dir     <- file.path(repo_root, "data")
analyses_dir <- file.path(repo_root, "analyses")
source(file.path(repo_root, "analyses", "scripts", "utils.R"))

rda_files       <- list.files(file.path(analyses_dir, "files", "rdata", "graph_data"),
                               pattern = "^graph_data_ecorr50_all_strat_", full.names = TRUE)
graph_data_path <- rda_files[which.max(file.mtime(rda_files))]
dataset_path    <- file.path(data_dir, "seamount_integrated_dataset.rda")

species_prev_rate <- 3

out_pdf <- file.path(script_dir, "Figure5.pdf")

# =============================================================================
# COLOUR PALETTES
# =============================================================================

habitat_colors <- c(
  "Bay"              = "#5ae6ab",
  "Lagoon"           = "#88d941",
  "Soft_back_reef"   = "#2e7d00",
  "Reef_outer_slope" = "#4e8273",
  "Summit50"         = "#ffe699",
  "DeepSlope150"     = "#d79c3b",
  "Summit250"        = "#4a6a94",
  "Summit500"        = "#1a3250",
  "NS"               = "gray"
)


distinct_13 <- c(
  "#E6194B", "#3CB44B", "#4363D8", "#F58231",
  "#911EB4", "#42D4F4", "#F032E6", "#BFEF45",
  "#FABED4", "#469990", "#DCBEFF", "#9A6324",
  "#800000"
)

alluvial_zone_colors <- c(
  "Shallow" = "#5ae6ab",
  "Middle"  = "#ffe699",
  "Deep"    = "#25456B",
  "NS.zone" = "gray"
)

alluvial_habitat_colors <- c(
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

# =============================================================================
# LOAD DATA
# =============================================================================

load(graph_data_path)   

# -- eDNA abundance table -----------------------------------------------------
load(dataset_path)

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

# =============================================================================
# GSEA — Zone + Habitat level
# =============================================================================

run_gsea <- function(posthoc_df, modules_list, padj_cols, resid_cols) {
  
  species_list <- unique(unlist(modules_list))
  
  gsc <- piano::loadGSC(
    data.frame(
      gene    = unlist(modules_list),
      geneset = rep(names(modules_list), lengths(modules_list))
    )
  )
  results_list <- list()
  for (grp in names(padj_cols)) {
    padj_vec <- setNames(posthoc_df[[padj_cols[[grp]]]], posthoc_df$feature)
    res_vec  <- setNames(posthoc_df[[resid_cols[[grp]]]], posthoc_df$feature)
    padj_vec <- padj_vec[names(padj_vec) %in% species_list]
    res_vec  <- res_vec[names(res_vec)  %in% species_list]
    
    gsea_res <- piano::runGSA(
      geneLevelStats = padj_vec,
      directions     = res_vec,
      gsc            = gsc,
      nPerm          = 1000,
      verbose        = FALSE
    )
    
    gsea_df <- piano::GSAsummaryTable(gsea_res) %>%
      dplyr::rename(
        cluster     = Name,
        n_species   = `Genes (tot)`,
        p_nondir    = `p (non-dir.)`,
        padj_nondir = `p adj (non-dir.)`,
        p_up        = `p (dist.dir.up)`,
        p_down      = `p (dist.dir.dn)`,
        padj_up     = `p adj (dist.dir.up)`,
        padj_down   = `p adj (dist.dir.dn)`,
        n_up        = `Genes (up)`,
        n_down      = `Genes (down)`
      ) %>%
      dplyr::select(cluster, n_species, p_nondir, padj_nondir,
                    p_up, p_down, padj_up, padj_down, n_up, n_down) %>%
      dplyr::mutate(
        cluster   = paste0("Cluster_", sprintf("%02d", as.numeric(cluster))),
        group_var = grp
      )
    results_list[[grp]] <- gsea_df
  }
  dplyr::bind_rows(results_list)
}

zone_padj_cols  <- list(Shallow = "padj_PH_Shallow",  Middle = "padj_PH_Middle",  Deep = "padj_PH_Deep")
zone_resid_cols <- list(Shallow = "residual_Shallow", Middle = "residual_Middle", Deep = "residual_Deep")

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

message("Running GSEA: Zone ...")
gsea_zone <- run_gsea(sp.chisq_zone, fast_greedy_modules,
                      zone_padj_cols, zone_resid_cols)
gsea_zone$level <- "Zone"

hab_levels     <- gsub("padj_PH_", "",
                       grep("padj_PH_", names(sp.chisq_habitat), value = TRUE))

# Define the desired order
desired_hab_order <- c("Bay", "Lagoon", "Soft_back_reef", "Reef_outer_slope",
                       "Summit50", "DeepSlope150", "Summit250", "Summit500")

hab_levels <- desired_hab_order[desired_hab_order %in% hab_levels]

hab_padj_cols  <- setNames(paste0("padj_PH_",  hab_levels), hab_levels)
hab_resid_cols <- setNames(paste0("residual_", hab_levels), hab_levels)

message("Running GSEA: Habitat ...")

gsea_habitat <- run_gsea(sp.chisq_habitat, fast_greedy_modules,
                         hab_padj_cols, hab_resid_cols)
gsea_habitat$level <- "Habitat"

# Significant clusters (enrichment direction only, padj_up < 0.05)
gsea_all            <- dplyr::bind_rows(gsea_zone, gsea_habitat)
all_significant_clusters <- sort(unique(gsea_all$cluster[!is.na(gsea_all$padj_up) &
                                                           gsea_all$padj_up < 0.05 &
                                                           gsea_all$n_up > 1]))
# get significant clusters by zone
significant_clusters.zone <- sort(unique(gsea_zone$cluster[!is.na(gsea_zone$padj_up) &
                                                             gsea_zone$padj_up < 0.05 &
                                                             gsea_zone$n_up > 1]))
# get significant clusters by habitat
significant_clusters.habitat <- sort(unique(gsea_habitat$cluster[!is.na(gsea_habitat$padj_up) &
                                                                   gsea_habitat$padj_up < 0.05 &
                                                                   gsea_habitat$n_up > 1]))

# =============================================================================
# PANEL A — Network coloured by habitat
# =============================================================================

panel_A <- cowplot::as_grob(function() {
  par(mar = c(0, 0.5, 2, 12),
      oma = c(0, 0, 0, 0),
      xpd = FALSE)
  
  plot(network,
       layout             = lay,
       vertex.label       = V(network)$label,
       vertex.color       = V(network)$color.habitat,
       vertex.shape       = V(network)$shape,
       vertex.size        = V(network)$size,
       vertex.frame.color = V(network)$frame.color,
       vertex.frame.width = V(network)$frame.width,
       edge.color         = E(network)$color,
       edge.width         = E(network)$edge_width,
       asp                = FALSE,
       rescale            = TRUE,
       edge.arrow.size    = 0.3,
       vertex.label.cex   = 0.5,
       vertex.label.dist  = V(network)$label.dist,
       vertex.label.font  = V(network)$label.font,
       vertex.label.color = V(network)$label.color,
  )
  
  mtext("A",
        side  = 3,
        line  = 0.2,
        adj   = 0.05,
        cex   = 2.5,
        font  = 2,
        outer = FALSE,
        padj  = 2.5)  
  
  par(xpd = TRUE)
  
  legend(x = 1.1, y = 1.1,
         legend = c("padj < 0.05", "padj >= 0.05"),
         col    = c("#782832", "gray80"),
         pch    = 19, pt.cex = 2.0, bty = "n", cex = 1.2,
         title  = "Node label (chi-sq)")
  
  legend(x = 1.1, y = 0.93,
         legend = c("Indicator sp.", "Other"),
         col    = c("firebrick1", "gray50"),
         pt.bg  = "gray80",
         pch    = 21, pt.lwd = 2,
         pt.cex = 2.0, bty = "n", cex = 1.2,
         title  = "Node border")
  
  legend(x = 1.1, y = 0.78,
         legend = names(habitat_colors),
         col    = "gray40",
         pt.bg  = habitat_colors,
         pch    = 21, pt.cex = 2.0, bty = "n", cex = 1.2,
         title  = "Habitat")
  
  par(xpd = FALSE)
})

# =============================================================================
# PANEL B — Same network coloured by fast-greedy module membership
# =============================================================================

panel_B <- cowplot::as_grob(function() {
  
  par(mar = c(0, 0.5, 2, 12),
      oma = c(0, 0, 0, 0),
      xpd = FALSE)
  
  plot(network,
       layout             = lay,
       vertex.label       = V(network)$label,
       vertex.color       = node_cluster_colors,
       vertex.shape       = V(network)$shape,
       vertex.size        = V(network)$size,
       vertex.frame.color = V(network)$frame.color,
       vertex.frame.width = V(network)$frame.width,
       edge.color         = E(network)$color,
       edge.width         = E(network)$edge_width,
       asp                = FALSE,
       rescale            = TRUE,
       edge.arrow.size    = 0.3,
       vertex.label.cex   = 0.5,
       vertex.label.dist  = V(network)$label.dist,
       vertex.label.font  = V(network)$label.font,
       vertex.label.color = V(network)$label.color,
  )
  mtext("B",
        side  = 3,
        line  = 0.2,
        adj   = 0.05,
        cex   = 2.5,
        font  = 2,
        outer = FALSE,
        padj  = 2.5)  
  
  par(xpd = TRUE)
  
  legend(x = 1.1, y = 1.1,
         legend = c("padj < 0.05", "padj >= 0.05"),
         col    = c("#782832", "gray80"),
         pch    = 19, pt.cex = 2.0, bty = "n", cex = 1.2,
         title  = "Node label (chi-sq)")
  
  legend(x = 1.1, y = 0.93,
         legend = c("Indicator sp.", "Other"),
         col    = c("firebrick1", "gray50"),
         pt.bg  = "gray80",
         pch    = 21, pt.lwd = 2,
         pt.cex = 2.0, bty = "n", cex = 1.2,
         title  = "Node border")
  
  cluster_node_counts <- table(species_to_cluster[V(network)$name])
  leg_colors <- cluster_palette[cluster_labels]
  
  leg_labels <- sapply(cluster_labels, function(cl) {
    n   <- as.integer(cluster_node_counts[cl])
    n   <- ifelse(is.na(n), 0L, n)
    
    rows <- gsea_zone[gsea_zone$cluster == cl, ]
    min_padj <- if (nrow(rows) == 0) NA_real_ else min(rows$padj, na.rm = TRUE)
    
    sig <- dplyr::case_when(
      !is.na(min_padj) & min_padj < 0.001 ~ "***",
      !is.na(min_padj) & min_padj < 0.01  ~ "**",
      !is.na(min_padj) & min_padj < 0.05  ~ "*",
      TRUE ~ ""
    )
    
    paste0(sig, cl, " (n=", n, ")")
  })
  
  leg_colors <- c(leg_colors, "gray80")
  leg_labels <- c(leg_labels, "NS")
  
  legend(x = 1.1, y = 0.78,
         legend = leg_labels,
         col    = "gray40",
         pt.bg  = leg_colors,
         pch    = 21, pt.cex = 2.0, bty = "n", cex = 1.2,
         title  = "Cluster membership")
  
  par(xpd = FALSE)
})

# =============================================================================
# PANEL C — Alluvial diagram: Zone → Module → Habitat
# =============================================================================

zone_order <- c("Shallow", "Middle", "Deep", "NS.zone")

module_order <- c("Cluster_08", "Cluster_01", "Cluster_02", "Cluster_03",
                  "Cluster_04", "Cluster_05", "Cluster_09", "Cluster_06",
                  "Cluster_10", "Cluster_12", "Cluster_07","Cluster_11", 
                  "Cluster_13")

habitat_order <- c("Bay", "Lagoon", "Soft_back_reef", "Reef_outer_slope",
                   "Summit50", "DeepSlope150", "Summit250", "Summit500", "NS.habitat")

alluvial_df <- modularity.df[, c("name", "fast_greedy")] %>%
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
      TRUE                                                ~ "NS.habitat"
    ),
    Zone = dplyr::case_when(
      !is.na(assigned_zone) & assigned_zone != "NS" ~ assigned_zone,
      TRUE                                          ~ "NS.zone"
    ),
    Module = paste0("Cluster_", sprintf("%02d", fast_greedy))
  ) %>%
  dplyr::count(Zone, Module, Habitat, name = "n_species") %>%
  dplyr::mutate(
    Zone    = factor(Zone,    levels = zone_order),
    Module  = factor(Module,  levels = module_order),
    Habitat = factor(Habitat, levels = habitat_order)
  )

# Habitat fill colours aligned to habitat_colors palette
alluvial_habitat_colors <- alluvial_habitat_colors[habitat_order]

panel_C <- ggplot(
  alluvial_df,
  aes(axis1 = Zone, axis2 = Module, axis3 = Habitat, y = n_species)
) +
  geom_alluvium(
    aes(fill = Habitat),
    width      = 0.25,
    alpha      = 0.75,
    knot.pos   = 0.4,
    curve_type = "sigmoid"
  ) +
  geom_stratum(
    width     = 0.25,
    fill      = "grey92",
    color     = "grey40",
    linewidth = 0.3
  ) +
  geom_text(
    stat       = "stratum",
    aes(label  = after_stat(stratum)),
    size       = 5.5,
    lineheight = 1.3,
    fontface   = "bold",
    color      = "grey20"
  ) +
  scale_x_discrete(
    limits = c("Zone", "Module", "Habitat"),
    expand = c(0.08, 0.01)
  ) +
  scale_fill_manual(
    values   = alluvial_habitat_colors,
    name     = "Habitat",
    na.value = "grey70",
    drop     = FALSE
  ) +
  scale_y_continuous(
    name   = "MOTUs",
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x         = element_text(size = 18, face = "bold"),
    axis.text.y         = element_blank(),
    axis.ticks.y        = element_blank(),
    axis.ticks.length.y = unit(0, "pt"),
    axis.title.x        = element_text(size = 16, face = "bold"),
    axis.title.y        = element_text(size = 16, face = "bold", vjust = 0.5),
    panel.grid.major    = element_blank(),
    panel.grid.minor    = element_blank(),
    legend.title        = element_text(size = 18, face = "bold", margin = margin(r = 25)),
    legend.text         = element_text(size = 16),
    legend.key.size     = unit(0.6, "cm"),
    legend.position     = "bottom",
    legend.direction    = "horizontal",
    plot.title          = element_text(face = "bold", hjust = 0, size = 32),
    plot.margin         = margin(t = 10, r = 10, b = 10, l = 30)
  ) +
  guides(
    fill = guide_legend(nrow = 3)
  ) +
  labs(title = "C")

# =============================================================================
# PANEL D — Module enrichment heatmap (Zone + Habitat, enrichment only)
# =============================================================================

# Filter GSEA results for significant enrichments (padj_up < 0.05) and n_up > 1

zone_heat <- gsea_zone %>%
  dplyr::filter(cluster %in% all_significant_clusters & n_up > 1) %>%
  dplyr::select(cluster, group_var, padj_up, n_species) %>%
  dplyr::mutate(level = "Zone")

hab_heat <- gsea_habitat %>%
  dplyr::filter(cluster %in% all_significant_clusters & n_up > 1) %>%
  dplyr::select(cluster, group_var, padj_up, n_species) %>%
  dplyr::mutate(level = "Habitat")


panel_D_data <- dplyr::bind_rows(zone_heat, hab_heat) %>%
  dplyr::mutate(
    neg_log10_padj = ifelse(is.na(padj_up) | padj_up >= 1, 0, -log10(padj_up)),
    sig_label = dplyr::case_when(
      padj_up < 0.001 ~ "***",
      padj_up < 0.01  ~ "**",
      padj_up < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    level     = factor(level, levels = c("Zone", "Habitat")),
    group_var = factor(
      group_var,
      levels = c("Shallow", "Middle", "Deep", hab_levels)
    )
  )

# Complete grid — every cluster x group combination
full_grid_D <- expand.grid(
  cluster   = all_significant_clusters,
  group_var = levels(panel_D_data$group_var),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    level = ifelse(group_var %in% c("Shallow", "Middle", "Deep"), "Zone", "Habitat"),
    level = factor(level, levels = c("Zone", "Habitat")),
    group_var = factor(group_var, levels = levels(panel_D_data$group_var))
  ) %>%
  dplyr::left_join(
    panel_D_data %>% dplyr::select(cluster, group_var, neg_log10_padj, sig_label),
    by = c("cluster", "group_var")
  ) %>%
  dplyr::mutate(
    sig_label      = tidyr::replace_na(sig_label, ""),
    neg_log10_padj = ifelse(sig_label == "", 0, neg_log10_padj) 
  )

panel_D <- ggplot(full_grid_D,
                  aes(x = group_var, y = cluster, fill = neg_log10_padj)) +
  geom_tile(color = "black", linewidth = 0.25) +
  geom_text(aes(label = sig_label),
            size = 5.5,
            vjust = 0.75, color = "white", fontface = "bold") +
  facet_grid(. ~ level, scales = "free_x", space = "free_x") +
  scale_fill_gradient(
    low      = "white",
    high     = "#782832",
    na.value = "grey90",
    name     = expression(-log[10](p[adj]))
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(
    expand = c(0, 0),
    limits = rev(intersect(module_order, all_significant_clusters))
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 16),
    axis.title.x     = element_text(size = 15, face = "bold", margin = margin(r = 8)),
    axis.text.y      = element_text(size = 16),
    axis.title.y     = element_text(size = 16, face = "bold", margin = margin(r = 20)),
    strip.text       = element_text(size = 16, face = "bold"),
    strip.background = element_rect(fill = "grey92"),
    legend.title     = element_text(size = 16, face = "bold"),
    legend.text      = element_text(size = 14),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA),
    plot.title       = element_text(face = "bold", hjust = 0, size = 32),
  ) +
  labs(
    title = "D",
    x     = NULL,
    y     = "Module (Fast greedy)"
  )

# =============================================================================
# COMPOSE & EXPORT
# =============================================================================

message("Composing final figure ...")

# Top row    : Panel A (habitat colours) | Panel B (module colours)
top_row    <- gridExtra::arrangeGrob(panel_A, panel_B, ncol = 2)

# Give Panel C more space in the layout
middle_row <- gridExtra::arrangeGrob(
  ggplotGrob(panel_C),
  ggplotGrob(panel_D),
  ncol   = 2,
  widths = c(1.6, 1)
)

pdf(out_pdf, width = 35, height = 30) 

gridExtra::grid.arrange(
  top_row,
  middle_row,
  nrow    = 2,
  heights = c(1.2, 1.2)
)

dev.off()

message("Figure saved to: ", out_pdf)