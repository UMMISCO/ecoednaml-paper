# =============================================================================
# Script name: Figure4_code.R
# Author: Estephe Kana & Edi Prifti & Eugeni Belda
# Purpose: Combine network plot (Panel A) and centrality metrics (Panel B)
#          into a single merged PDF figure
#          >> Centrality metrics restricted to Degree and Betweenness 
# =============================================================================

# -----------------------------------------------------------------------------
# Commands to run the script (from repo root):
# Rscript analyses/figures/Figure4/Figure4_code.R
# -----------------------------------------------------------------------------

# Check for required packages
required_pkgs <- c("igraph", "ggplot2", "dplyr", "tidyr", "ggpubr", "cowplot")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))

library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(cowplot)

# =============================================================================
# LOAD DATA
# =============================================================================

# Derive repo root from script location
script_path  <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
script_dir   <- dirname(script_path)                      # <repo>/analyses/figures/Figure4
repo_root    <- dirname(dirname(dirname(script_dir)))     # <repo>/
analyses_dir <- file.path(repo_root, "analyses")

rda_files       <- list.files(file.path(analyses_dir, "files", "rdata", "graph_data"),
                               pattern = "^graph_data_ecorr50_all_strat_", full.names = TRUE)
graph_data_path <- rda_files[which.max(file.mtime(rda_files))]
load(graph_data_path)

# load indicator species — terinter model
load(file.path(analyses_dir, "analysis_outputs", "terinter_output_data", "terinter_Predomics_all_analyses_overall_data_strat_group_prev_3.Rda"))
indicSp_ter.df <- predout.bin.sub
rm(adonis_pred.bin, adonis_pred.maxn, predout.bin, predout.bin.sub, predout.maxn, predout.maxn.sub)

# load indicator species — bininter model
load(file.path(analyses_dir, "analysis_outputs", "bininter_output_data", "bininter_Predomics_all_analyses_overall_data_strat_group_prev_3.Rda"))
indicSp_bin.df <- predout.bin.sub
rm(adonis_pred.bin, adonis_pred.maxn, predout.bin, predout.bin.sub, predout.maxn, predout.maxn.sub)

keySpecies_bin <- indicSp_bin.df[indicSp_bin.df$IsIndSp == 1, ]
keySpecies_ter <- indicSp_ter.df[indicSp_ter.df$IsIndSp == 1, ]

nodes.annot$IsIndSp <- ifelse(
  nodes.annot$name %in% keySpecies_bin$feature | nodes.annot$name %in% keySpecies_ter$feature,
  "Yes", "No"
)

# =============================================================================
# PANEL A — Network plot
# =============================================================================

panel_A <- cowplot::as_grob(function() {
  
  par(mar = c(0, 0, 0, 10), xpd = FALSE)  
  
  plot(network,
       layout             = lay,
       vertex.label       = V(network)$label,
       vertex.color       = V(network)$color,
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
       vertex.label.color = V(network)$label.color
  )
  
  # xpd = TRUE allows drawing outside the plot region into the margin
  par(xpd = TRUE)
  
  legend(x = 1.2, y = 1,          
         legend = c("padj < 0.05", "padj >= 0.05"),
         col    = c("#782832", "gray80"),
         pch    = 19, pt.cex = 1.0, bty = "n", cex = 0.8,
         title  = "Node label (chi-sq)")
  
  legend(x = 1.2, y = 0.7,
         legend = c("Shallow", "Middle", "Deep", "NS"),
         col    = c("#06d6a0", "#ffd166", "#25456B", "gray"),
         pch    = 19, pt.cex = 1.0, bty = "n", cex = 0.8,
         title  = "Zone (chi-sq PH)")
  
  legend(x = 1.2, y = 0.4,
         legend = c("Indicator species", "Other"),
         col    = c("firebrick1", NA),
         pch    = 21, pt.bg = "gray80", pt.lwd = 2,
         pt.cex = 1.0, bty = "n", cex = 0.8,
         title  = "Node border")
  
  par(xpd = FALSE)
})

# =============================================================================
# PANEL B — Degree & Betweenness violin/boxplot
# =============================================================================

metrics_long <- nodes.annot %>%
  dplyr::select(name, IsIndSp, degr_cent, betw_cent) %>%
  tidyr::pivot_longer(
    cols      = c(degr_cent, betw_cent),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    IsIndSp = factor(IsIndSp,
                     levels = c("No", "Yes"),
                     labels = c("Non-Indicator", "Indicator")),
    metric = recode(metric,
                    degr_cent = "Degree Centrality",
                    betw_cent = "Betweenness Centrality")
  )

panel_B <- ggplot(metrics_long, aes(x = IsIndSp, y = value, fill = IsIndSp)) +
  
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA,
               fill = "white", color = "grey30", alpha = 0.8) +
  geom_jitter(aes(color = IsIndSp), width = 0.08, size = 1.2,
              alpha = 0.5, show.legend = FALSE) +
  
  stat_compare_means(
    method  = "wilcox.test",
    label   = "p.signif",
    label.x = 1.5,
    vjust   = -0.5,
    size    = 5
  ) +
  
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  
  scale_fill_manual(values  = c("Non-Indicator" = "#7fbfff", "Indicator" = "#ff7f7f")) +
  scale_color_manual(values = c("Non-Indicator" = "#3a7fc1", "Indicator" = "#c13a3a")) +
  
  labs(x = NULL, y = "Centrality Value", fill = NULL) +
  
  theme_bw(base_size = 13) +
  theme(
    strip.text       = element_text(face = "bold", size = 12),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

# =============================================================================
# MERGE PANELS A + B and save PDF
# =============================================================================

out_pdf <- file.path(script_dir, "Figure4.pdf")

pdf(file = out_pdf, width = 20, height = 12)

cowplot::plot_grid(
  panel_A, panel_B,
  ncol           = 2,
  rel_widths     = c(1.1, 0.9),
  labels         = c("A", "B"),
  label_size     = 16,
  label_fontface = "bold"
)

dev.off()

message("Saved merged figure to: ", out_pdf)