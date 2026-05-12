# =============================================================================
# Script name: Figure3_code.R
# Authors: Estephe Kana & Edi Prifti & Eugeni Belda
# Purpose: Build a 3-panel publication figure:
#   Panel A — Feature Importance and Cliff's delta of indicator MOTUs by Comparison
#   Panel B — AUC performances in Testing of Predomics FBM models by comparison
#   Panel C — PERMANOVA of set of MOTUs from all the community or from indicator MOTUs identified by Predomics bininter and terinter
# =============================================================================

# -----------------------------------------------------------------------------
# Commands to run the script (from repo root):
# Rscript analyses/figures/Figure3/Figure3_code.R
# -----------------------------------------------------------------------------

# Check for required packages
# Note: predomics — remotes::install_github("eprifti/predomics")
required_pkgs <- c("ggplot2", "ggh4x", "patchwork", "dplyr", "ggpubr",
                   "vegan", "viridis", "predomics", "reshape2")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))

# Load required packages
library(ggplot2)
library(ggh4x)
library(patchwork)
library(dplyr)
library(ggpubr)
library(vegan)
library(viridis)
library(predomics)
library(reshape2)

# Derive repo root from script location
script_path  <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
script_dir   <- dirname(script_path)                      # <repo>/analyses/figures/Figure3
repo_root    <- dirname(dirname(dirname(script_dir)))     # <repo>/
analyses_dir <- file.path(repo_root, "analyses")

out_pdf <- file.path(script_dir, paste0("Figure3_", Sys.Date(), ".pdf"))

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

# Combine the predomics results into a single data frame
mergedf.all <- rbind(objlist$bininter$predout.bin.sub,
                     objlist$terinter$predout.bin.sub)

mergedf.all.bplot <- mergedf.all[mergedf.all$IsIndSp == 1,]
mergedf.all.bplot <- data.frame(table(mergedf.all.bplot$source, mergedf.all.bplot$data, mergedf.all.bplot$comparison))

## Visualization of predomics results

## redefine getImportanceFeaturesFBMobjects from predomics package
getImportanceFeaturesFBMobjects <- function (clf_res, X, y, verbose = TRUE, filter.cv.prev = 0, 
                                             scaled.importance = FALSE, k_penalty = 0, k_max = 0) 
{
  mode <- NULL
  if (!(predomics::isExperiment(clf_res))) {
    stop("analyzeLearningFeatures: please provide a valid experiment results!")
  }
  if (clf_res$classifier$params$objective == "cor") {
    mode <- "regression"
  }
  else {
    mode <- "classification"
  }
  if (!is.null(mode)) {
    cat(paste("... Estimating mode: ", mode, "\n"))
  }
  else {
    stop("analyzeImportanceFeatures: mode not founding stopping ...")
  }
  pop <- modelCollectionToPopulation(clf_res$classifier$models)
  if (verbose) 
    print(paste("There are", length(pop), "models in this population"))
  pop <- selectBestPopulation(pop = pop, p = 0.05, k_penalty = k_penalty, k_max = k_max)
  if (verbose) 
    print(paste("There are", length(pop), "models in this population after selection of the best"))

  pop.df <- populationToDataFrame(pop = pop)
  pop.noz <- listOfModelsToDenseCoefMatrix(clf = clf_res$classifier, 
                                           X = X, y = y, list.models = pop)
  if (verbose) 
    print(paste("Pop noz object is created with", nrow(pop.noz), 
                "features and", ncol(pop.noz), "models"))
  fa <- makeFeatureAnnot(pop = pop, X = X, y = y, clf = clf_res$classifier)
  pop.noz <- fa$pop.noz
  pop.noz <- data.frame(pop.noz)
  pop.noz$feature <- rownames(pop.noz)
  pop.noz <- reshape::melt(pop.noz)
  pop.noz$learner <- unlist(lapply(strsplit(as.character(pop.noz$variable), 
                                            split = "_"), function(x) {
                                              x[1]
                                            }))
  pop.noz$learner <- paste(pop.noz$learner, clf_res$classifier$params$language, 
                           sep = ".")
  pop.noz$model <- as.character(unlist(lapply(strsplit(as.character(pop.noz$variable), 
                                                       split = "_"), function(x) {
                                                         x[2]
                                                       })))
  pop.noz$value <- factor(pop.noz$value, levels = c(-1, 0, 
                                                    1))
  pop.noz$value <- droplevels(pop.noz$value)
  lr <- list(clf_res)
  names(lr) <- paste(clf_res$classifier$learner, clf_res$classifier$params$language, 
                     sep = ".")
  feat.import <- mergeMeltImportanceCV(list.results = lr, filter.cv.prev = filter.cv.prev, 
                                       min.kfold.nb = FALSE, learner.grep.pattern = "*", nb.top.features = nrow(X), 
                                       feature.selection = NULL, scaled.importance = scaled.importance, 
                                       make.plot = TRUE)
  if (is.null(feat.import)) {
    stop("analyzeImportanceFeatures: no feature importance data found... returning empty handed.")
  }
  if (mode == "regression") {
    featPrevPlot <- plotPrevalence(features = rownames(X), 
                                   X, y = NULL)
  }
  else {
    featPrevPlot <- plotPrevalence(features = rownames(X), 
                                   X, y)
  }
  effSizes.df <- computeEffectSizes(X = X, y = y, mode = mode)
  outlist <- list(featprevFBM = pop.noz, featImp = feat.import$summary, 
                  effectSizes = effSizes.df, featPrevGroups = featPrevPlot$data)
  class(outlist) <- "listFBMfeatures"
  return(outlist)
}

# Get the feature importance and effect sizes for each comparison and model from the Predomics results
alldatacomparison <- list()

for (i in unique(mergedf.all$comparison)){
  alldatacomparison[[i]][['terinter']] = list(presAbs=objlist$terinter$predout.bin[[i]])
  alldatacomparison[[i]][['bininter']] = list(presAbs=objlist$bininter$predout.bin[[i]])
}

model_perf <- list()

for (i in unique(mergedf.all$comparison)) {
  model_perf[[i]] <- list()
  
  for (model in unique(mergedf.all$source)) {
    pred.obj <- alldatacomparison[[i]][[model]][["presAbs"]]
    res_clf <- pred.obj$fit
    X <- pred.obj$Comp_data$X
    y <- pred.obj$Comp_data$y
    
    model_perf[[i]][[model]] <- getImportanceFeaturesFBMobjects(
      clf_res = res_clf, 
      X = t(X), 
      y = y, 
      verbose = TRUE, 
      k_penalty = 0, 
      k_max = 0, 
      filter.cv.prev = 0
    )
  }
}

# Combine the effect sizes and feature importance results into data frames for plotting

effsize.df <- do.call(rbind, lapply(names(model_perf), function(comp) {
  do.call(rbind, lapply(names(model_perf[[comp]]), function(src) {
    df         <- model_perf[[comp]][[src]][["effectSizes"]]
    df$source  <- src
    df$fdr     <- p.adjust(df$pval.wilcox, method = "BH")
    df$Comparison <- comp
    df
  }))
}))


featImp.df <- do.call(rbind, lapply(names(model_perf), function(comp) {
  do.call(rbind, lapply(names(model_perf[[comp]]), function(src) {
    df            <- model_perf[[comp]][[src]][["featImp"]]
    df$source     <- src
    df$Comparison <- comp
    df
  }))
}))

tfeats.fbm <- mergedf.all[mergedf.all$IsIndSp==1 & mergedf.all$data=="pres/abs", ]

# Create a key in both data frames
featImp.df$key  <- paste(featImp.df$feature, featImp.df$Comparison)
tfeats.fbm$key  <- paste(tfeats.fbm$feature, tfeats.fbm$comparison)

# Filter
featImp.df.filtered <- featImp.df[featImp.df$key %in% tfeats.fbm$key, ]
featImp.df.filtered$key <- NULL

featImp.df.filtered.dcast <- dcast(data = featImp.df.filtered, formula = feature~Comparison, value.var="value")

effsize.df$key  <- paste(effsize.df$feature, effsize.df$Comparison)
effsize.df.filtered <- effsize.df[effsize.df$key %in% tfeats.fbm$key, ]
effsize.df.filtered$key <- NULL
effsize.df.filtered$label <- ifelse(
  effsize.df.filtered$fdr < 0.05, "**", 
  ifelse(effsize.df.filtered$pval.wilcox < 0.05, "*",""))

# Get all unique features and set consistent order across all plots

feature_order <- featImp.df.filtered %>%
  group_by(feature) %>%  
  summarise(mean_importance = mean(abs(value), na.rm = TRUE)) %>%
  arrange(desc(mean_importance)) %>%
  pull(feature)

# Convert feature to ordered factor in all dataframes
tfeats.fbm$feature <- factor(tfeats.fbm$feature, levels = feature_order)
featImp.df.filtered$feature <- factor(featImp.df.filtered$feature, levels = feature_order)
effsize.df.filtered$feature <- factor(effsize.df.filtered$feature, levels = feature_order)

# Plot 1: Feature presence tiles
plot1 <- ggplot(tfeats.fbm, aes(x = comparison, y = feature)) + 
  geom_tile(fill = "black", colour = "white", linewidth = 0.5) + 
  xlab("Comparison") + 
  ylab("Indicator MOTUs") +
  theme_minimal() + 
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=14),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    panel.background = element_rect(fill = "grey95", colour = NA),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  scale_y_discrete(drop = FALSE)

# Plot 2: Feature importance
featImp.df.filtered$sign <- factor(featImp.df.filtered$sign, levels=c("-1","1"))

featImp.df.filtered$source <- factor(featImp.df.filtered$source, 
                                     levels=c("bininter","terinter"))

source_colours <- c("bininter" = "black",
                    "terinter" = "gray80")

plot2 <- ggplot(featImp.df.filtered, aes(x=feature, y=value)) +
  geom_hline(yintercept = min(0, featImp.df.filtered$value, na.rm = TRUE), col="gray") +
  scale_colour_manual("Source", values = source_colours) +
  ylab("Feature importance") +
  xlab("") +
  facet_grid(.~Comparison) +
  theme_bw() +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +  
  geom_errorbar(aes(ymin = value - se, ymax = value + se, colour = source), width=.1, position=position_dodge(0.7)) +
  geom_point(position = position_dodge(0.7), size=2, aes(colour = source)) +
  guides(colour = "none", fill = "none") +
  theme_bw() + 
  theme(
    strip.text.x = element_text(size=14),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 16),
    strip.background = element_rect(fill = NA)
  )

# Plot 3: Effect sizes
plot3 <- ggplot(effsize.df.filtered, aes(x=feature, y=cdelta, fill=cdelta)) + 
  geom_bar(stat = "identity") + 
  geom_text(aes(label=label)) + 
  scale_fill_gradient(low = "deepskyblue1", high = "firebrick1", limits=c(-1,1)) + 
  scale_x_discrete(drop = FALSE) +  
  ylab("Cliff's delta 1 vs. -1") + 
  ggtitle("") + 
  facet_grid(.~Comparison) + 
  coord_flip() + 
  theme_bw() + 
  theme(
    strip.text.x = element_text(size=14),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 16),
    strip.background = element_rect(fill = NA),
    legend.position = "none"
  )

# ── Plot 4: AUC boxplot 
predout_terinter.bin=objlist$terinter$predout.bin
predout_bininter.bin=objlist$bininter$predout.bin

fbmAUC_terinter.df <- do.call(
  rbind,
  lapply(names(predout_terinter.bin), function(comp) {
    df <- predout_terinter.bin[[comp]][["fbm"]]
    df$comparison <- comp
    df$model_id <- rownames(df)
    df$source     <- "terinter"
    df
  })
)

fbmAUC_bininter.df <- do.call(
  rbind,
  lapply(names(predout_bininter.bin), function(comp) {
    df <- predout_bininter.bin[[comp]][["fbm"]]
    df$comparison <- comp
    df$model_id <- rownames(df)
    df$source     <- "bininter"
    df
  })
)

fbmAUC.df <- rbind(fbmAUC_bininter.df, fbmAUC_terinter.df)
fbmAUC.df$language <- factor(fbmAUC.df$language, levels=c("bininter","terinter"))

# Automatic pairwise comparisons
comparisons_list <- levels(factor(fbmAUC.df$comparison))
comparison_pairs <- combn(comparisons_list, 2, simplify = FALSE)

# Add model counts label per comparison
model_counts <- fbmAUC.df %>%
  group_by(comparison) %>%
  summarise(n = n(), .groups = "drop")

model_counts$label <- paste0("n=", model_counts$n)

plot4 <- ggplot(data = fbmAUC.df, aes(x = comparison, y = auc_)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_point(aes(colour = language),
             position = position_jitter(width = 0.2, seed = 123),
             size = 2,
             alpha = 0.6
  ) +
  geom_text(data = model_counts,
            aes(x = comparison, y = Inf, label = label),
            vjust = 1.5, size = 6, inherit.aes = FALSE) +
  scale_colour_manual("Source", values = source_colours) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = comparison_pairs,
    label = "p.signif",
    step.increase = 0.08,
    hide.ns = TRUE,
    size = 7    
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  ylab("Test AUC") +
  xlab("FBM models") +
  guides(colour = guide_legend(title = "Model"), size = "none") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
    strip.background = element_rect(fill = NA),
    axis.title = element_text(size = 16),
    legend.position = "bottom",
    legend.text  = element_text(size = 16),
    legend.title = element_text(size = 16)
  )


# ── Plot 5 - PERMANOVA

# Initialize a list to store results of permanova
adonis_pred <- list()

# Function to perform PERMANOVA on all OTUs for a given comparison and model
permanova_all_otus <- function(comparison, source_name = "terinter") {
  data  <- objlist[[source_name]]$predout.bin[[comparison]][["Comp_data"]]
  all.X <- data$X
  all.y <- data$y
  names(all.y) <- rownames(all.X)
  
  jacc     <- vegdist(all.X, method = "jaccard", binary = TRUE)
  all.meta <- data.frame(sample = rownames(all.X), class = all.y)
  set.seed(100)
  res <- adonis2(jacc ~ class, data = all.meta)
  
  out            <- data.frame(res)[1, , drop = FALSE]
  out$comparison <- comparison
  out$source     <- "all MOTUs"
  out$features   <- ncol(all.X)
  out
}

# Function to perform PERMANOVA on indicator MOTUs for a given comparison and model
permanova_analysis_source <- function(comparison, source_name) {
  
  data  <- objlist[[source_name]]$predout.bin[[comparison]][["Comp_data"]]
  all.X <- data$X
  all.y <- data$y
  names(all.y) <- rownames(all.X)

  fbm.species <- unique(
    tfeats.fbm[tfeats.fbm$comparison == comparison &
                 tfeats.fbm$source     == source_name, "feature"]
  )
  
  if (length(fbm.species) == 0) return(NULL)
  
  subdf <- all.X[, fbm.species, drop = FALSE]
  subdf <- subdf[rowSums(subdf) > 0, , drop = FALSE]
  subdf.class <- all.y[rownames(subdf)]
  
  jacc   <- vegdist(subdf, method = "jaccard", binary = TRUE)
  meta   <- data.frame(sample = rownames(subdf), class = subdf.class)
  set.seed(100)
  adonis.res <- adonis2(jacc ~ class, data = meta)
  
  out <- data.frame(adonis.res)[1, , drop = FALSE]
  out$comparison <- comparison
  out$source     <- paste0(source_name, " MOTUs")
  out$features   <- ncol(subdf)
  out
}

# Loop through each comparison and perform PERMANOVA for all OTUs and indicator OTUs from both Predomics models
for (i in unique(tfeats.fbm$comparison)) {
  rows <- list()
  rows[["all"]] <- permanova_all_otus(i)
  
  for (src in c("bininter", "terinter")) {
    r <- permanova_analysis_source(i, src)
    if (!is.null(r)) rows[[src]] <- r
  }
  adonis_pred[[i]] <- do.call(rbind, rows)
}

# Combine results into a single data frame for plotting
adonis.alldf <- do.call(rbind, adonis_pred)
adonis.alldf$pval.cat <- dplyr::case_when(
  adonis.alldf$Pr..F. < 0.001 ~ "***",
  adonis.alldf$Pr..F. < 0.01  ~ "**",
  adonis.alldf$Pr..F. < 0.05  ~ "*",
  TRUE            ~ "")
adonis.alldf$source   <- factor(adonis.alldf$source,
                                levels = c("all MOTUs",
                                           "bininter MOTUs",
                                           "terinter MOTUs"))

# Plot 5 - PERMANOVA results
plot5 <- ggplot(adonis.alldf, aes(x = source, y = R2)) +
  geom_bar(stat = "identity", aes(fill = features)) +
  geom_text(aes(label = pval.cat), size = 7, vjust = 0.1) +  
  ylab("R2") +
  xlab("MOTUs set") +
  scale_fill_viridis() +
  facet_grid(. ~ comparison) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    strip.text.x = element_text(size = 14),
    legend.position = "bottom",
    axis.title = element_text(size = 16),
    legend.title = element_text(color = "black", size = 16, vjust = 0.9, margin = margin(r = 5))
  )


plot2_fixed <- plot2 + force_panelsizes(cols = unit(1.5, "cm") * 3)
plot3_fixed <- plot3 + force_panelsizes(cols = unit(1.5, "cm") * 3)

plot.list <- list("featFBM"=plot1,
                  "featImp"=plot2_fixed,
                  "effectSizes"=plot3_fixed)

layout <- "
      AABC
      "
# combine plots into a single figure
combined_plot <- patchwork::wrap_plots(plot.list, design = layout)+
  plot_annotation(tag_levels = NULL)

# Add panel labels directly to the plots
combined_plot_labeled <- combined_plot + 
  plot_annotation(
    title = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 16))
  )

plot4_labeled <- plot4 + 
  ggtitle("B") + 
  theme(plot.title = element_text(face = "bold", size = 16))

plot5_labeled <- plot5 + 
  ggtitle("C") + 
  theme(plot.title = element_text(face = "bold", size = 16))

# Wrap combined_plot1 as a grob so patchwork treats it as a single unit
final_figure <- wrap_elements(full = patchworkGrob(combined_plot_labeled)) /
  (plot4_labeled | plot5_labeled) +
  plot_layout(heights = c(2, 1))

# Save the final figure to a PDF
pdf(out_pdf, h = 20, w = 20)
print(final_figure)
dev.off()

message("Figure saved to: ", out_pdf)