# =============================================================================
# Script name: Figure1_code.R
# Authors: Estephe Kana & Edi Prifti & Eugeni Belda
# Date created: 2025-12-10
# Purpose: Build a 4-panel publication figure:
#   Panel A — Map of sampling sites (New Caledonia, coloured by Habitat/Zone)
#   Panel B — Heatmap of MOTU abundance across samples (ordered by clustering)
#   Panel C — Species richness barplot per sample and habitat
#   Panel D — PCoA of Jaccard beta-diversity with environmental fitting
# Inputs:  data/seamount_integrated_dataset.rda
# Outputs: analyses/figures/Figure1/Figure1.pdf
# =============================================================================

# -----------------------------------------------------------------------------
# Commands to run the script (from repo root):
# Rscript analyses/figures/Figure1/Figure1_code.R
# -----------------------------------------------------------------------------

# Check for required packages
required_pkgs <- c("ggplot2", "sf", "vegan", "reshape2", "patchwork", "ggspatial",
                   "maps", "rnaturalearth", "ggrepel", "viridis", "dplyr", "ggh4x")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))

# Load required libraries
library(ggplot2)
library(sf)
library(vegan)
library(reshape2)
library(patchwork)
library(ggspatial)
library(maps)
library(rnaturalearth)
library(ggrepel)
library(viridis)
library(dplyr)
library(ggh4x)


## Mapview and heatmap of eDNA data

# Derive repo root from script location
script_path  <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
script_dir   <- dirname(script_path)                      # <repo>/analyses/figures/Figure1
repo_root    <- dirname(dirname(dirname(script_dir)))     # <repo>/
data_dir     <- file.path(repo_root, "data")
analyses_dir <- file.path(repo_root, "analyses")
source(file.path(repo_root, "analyses", "scripts", "utils.R"))

# load eDNA dataset
load(file.path(data_dir, "seamount_integrated_dataset.rda"))

# get environmental data
env_eDNA_df <- sm$sample_info
env_eDNA_df$Habitat[env_eDNA_df$Habitat == "DeepSlope"] <- "DeepSlope150"

# env_eDNA_df$Zone <- ifelse(env_eDNA_df$Habitat %in% c("Bay","Lagoon","Reef_outer_slope", "Soft_back_reef"), "INSHORE", "OFFSHORE")
env_eDNA_df$Zone <- ifelse(env_eDNA_df$Habitat %in% c("Bay","Lagoon","Reef_outer_slope", "Soft_back_reef"), "Shallow",
                           ifelse(env_eDNA_df$Habitat %in% c("Summit50", "DeepSlope150"), "Middle",
                                  ifelse(env_eDNA_df$Habitat %in% c("Summit250", "Summit500"), "Deep",
                                         NA)))
env_eDNA_sf <- st_as_sf(env_eDNA_df, coords = c('Longitude', 'Latitude'), crs=4326)

# Use the bounding box of the sample points to set the map extent
bbox <- st_bbox(env_eDNA_sf)
buffer <- 0.045
bbox_expanded <- bbox + c(-buffer, -buffer, buffer, buffer)

# Using built-in world map data; subset to oceania and new caledonia card
world <- ne_countries(scale="medium", continent = "oceania", country = "new caledonia", returnclass = "sf")

## Add labels for regions
states <- st_as_sf(map("world", regions = "new caledonia", plot = FALSE, fill = TRUE))
states <- cbind(states, st_coordinates(st_centroid(states)))

#Crop the map object with st_crop to New caledonia island
world_cropped <- st_crop(world, bbox_expanded)

# Increase latitude range to make the map taller
lat_expansion <- 0.15
long_expansion <- 0.07

# Modify bbox_expanded manually
bbox_expanded["ymin"] <- bbox["ymin"] - lat_expansion
bbox_expanded["ymax"] <- bbox["ymax"] + lat_expansion
bbox_expanded["xmin"] <- bbox["xmin"] - long_expansion
bbox_expanded["xmax"] <- bbox["xmax"] + long_expansion

# factorize Habitat and Zone
env_eDNA_df$Habitat <- factor(env_eDNA_df$Habitat, levels= c("Bay","Lagoon", "Soft_back_reef", "Reef_outer_slope","Summit50", "DeepSlope150", "Summit250", "Summit500"))
env_eDNA_df$Zone <- factor(env_eDNA_df$Zone, levels= c("Shallow","Middle", "Deep"))

p1 <- ggplot() +
  annotation_map_tile(data = world_cropped, type = "osm", zoom = 14) +
  geom_point(
    data = env_eDNA_df,
    aes(x = Longitude, y = Latitude, color = Habitat, shape = Zone),
    size = 3.5,
    stroke = 1,
    alpha = 0.8,
    position = position_jitter(width = 0.15, height = 0.15)
  ) +
  scale_color_manual(values = c(
    "Bay" = "#5ae6ab", 
    "Lagoon" = "#88d941", 
    "Soft_back_reef" = "#2e7d00",
    "Reef_outer_slope" = "#4e8273", 
    "Summit50" = "#ffe699", 
    "DeepSlope150" = "#d79c3b", 
    "Summit250" = "#4a6a94", 
    "Summit500" = "#1a3250"
  )) +
  scale_shape_manual(values = c(16, 17, 15)) +
  coord_sf(
    xlim = c(bbox_expanded["xmin"], bbox_expanded["xmax"]),
    ylim = c(bbox_expanded["ymin"], bbox_expanded["ymax"]),
    expand = TRUE
  ) +
  theme_minimal() +
  labs(
    x = "Longitude",
    y = "Latitude",
    color = "Habitat",  
    shape = "Zone"
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),  
    shape = guide_legend(nrow = 1)                 
  ) +
  theme(
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.position = "none",
    legend.box = "vertical"  
  )

# Heatmap of co-abundance patterns across samples

#Get the raw data for the heatmap
dfaims <- sm$X
#Get the metadata for the heatmap
dfaims_meta <- sm$sample_info
dfaims_meta$Habitat[dfaims_meta$Habitat == "DeepSlope"] <- "DeepSlope150"
# dfaims_meta$Zone <- ifelse(dfaims_meta$Habitat %in% c("Bay","Lagoon","Reef_outer_slope", "Soft_back_reef"), "INSHORE", "OFFSHORE")
dfaims_meta$Zone <- ifelse(dfaims_meta$Habitat %in% c("Bay","Lagoon","Reef_outer_slope", "Soft_back_reef"), "Shallow",
                           ifelse(dfaims_meta$Habitat %in% c("Summit50", "DeepSlope150"), "Middle",
                                  ifelse(dfaims_meta$Habitat %in% c("Summit250", "Summit500"), "Deep",
                                         NA)))
rownames(dfaims_meta) <- dfaims_meta$Spygen
#Get taxo info for species
dfaims_taxo <- sm$taxonomy
dfaims_melt <- as.data.frame(dfaims)
dfaims_melt$feature <- rownames(dfaims)
dfaims_melt <- melt(dfaims_melt)
dfaims_clust.sp <- hclust(dist(dfaims, method = "euclidean"), method = "ward.D") #cluster species with ward method from euclidean distances
dfaims_clust.samples <- hclust(dist(t(dfaims), method = "euclidean"), method = "ward.D") #cluster samples with ward method from euclidean distances
#Fix the samples and species order from the clustering results in the melted object
dfaims_melt$feature <- factor(dfaims_melt$feature, levels = dfaims_clust.sp$labels[dfaims_clust.sp$order])
dfaims_melt$variable <- factor(dfaims_melt$variable, levels = dfaims_clust.samples$labels[dfaims_clust.samples$order])
dfaims_melt <- merge(dfaims_melt, dfaims_meta[,c("Habitat", "Zone", "Site")], by.x="variable", by.y=0, all.x=TRUE)
dfaims_taxo$feature <- rownames(dfaims_taxo)
dfaims_melt <- merge(dfaims_melt, dfaims_taxo[,c("order","family","genus", "tax_name")], by.x="feature", by.y=0, all.x=TRUE)
dfaims_melt$Habitat[dfaims_melt$Habitat == "DeepSlope"] <- "DeepSlope150"

# get samples with species at least 10% prevalence across samples
dfaims_3_prev_species= get_sample_by_prevalence(t(dfaims), 3)


# Add prevalence tag "prev<10" or "prev>=10" to each feature
dfaims_melt$prev_rate <- ifelse(dfaims_melt$feature %in% colnames(dfaims_3_prev_species), "prev>=3",
                                ifelse(!(dfaims_melt$feature %in% colnames(dfaims_3_prev_species)), "prev<3",
                                       NA)
)

# Factorize prev_rate
dfaims_melt$prev_rate= as.factor(dfaims_melt$prev_rate)

dfaims_melt$Habitat= factor(dfaims_melt$Habitat, levels= c("Bay","Lagoon", "Soft_back_reef", "Reef_outer_slope","Summit50", "DeepSlope150", "Summit250", "Summit500"))
dfaims_melt$Zone= factor(dfaims_melt$Zone, levels= c("Shallow","Middle", "Deep"))


# Get the order of samples based on Site variable in dfaims_melt
site_order <- dfaims_melt %>%
  select(variable, Site) %>%
  distinct() %>%
  arrange(Site) %>%
  pull(variable)

# Apply the ordered factor to dfaims_melt
dfaims_melt$variable <- factor(dfaims_melt$variable, levels = site_order)


p2 <- ggplot(dfaims_melt, aes(x=variable, y=feature, fill=value)) +
  geom_tile(colour=NA,linewidth = 0.5) +
  scale_fill_gradient(low = "white", high = viridis::viridis(100))+
  ylab("MOTUs") +
  labs(fill="MOTU abundance") +
  facet_nested(prev_rate~Zone+Habitat, scales = "free", space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.position = "right",
        legend.title = element_text(color = "black", size = 16, angle= 90, vjust = 0.9, margin = margin(b = 15)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 16, angle =0)
  )

#richness barplot

dfaims_meta$Habitat <- factor(dfaims_meta$Habitat, levels= c("Bay","Lagoon", "Soft_back_reef", "Reef_outer_slope","Summit50", "DeepSlope150", "Summit250", "Summit500"))
# Apply the same order to dfaims_meta (matching by sample name)
dfaims_meta$Spygen <- factor(dfaims_meta$Spygen, levels = site_order)

p3 <- ggplot(dfaims_meta, aes(x=Spygen, y=otu_richness, fill=Habitat)) + 
  geom_bar(stat = "identity") + 
  xlab("samples (Habitat)") +
  ylab("MOTU_richness") +
  labs(fill="Habitat") +
  facet_nested(.~Habitat, scales = "free_x", space = "free") + 
  theme_bw() +
  scale_fill_manual(values = c("Bay" = "#5ae6ab", 
                               "Lagoon" = "#88d941", 
                               "Soft_back_reef" = "#2e7d00",
                               "Reef_outer_slope" = "#4e8273", 
                               "Summit50" = "#ffe699", 
                               "DeepSlope150" = "#d79c3b", 
                               "Summit250" = "#4a6a94", 
                               "Summit500" = "#1a3250"))+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size = 16),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16, margin = margin(r = 15)),
        legend.position = "bottom")


################
# PCoA analyses
################


df.all <- sm$X

## Compute beta-diversity with jaccard  method (pairwise sample distance from presence/absence data)
X.jac <- vegdist(x = t(df.all), method = "jaccard", binary = TRUE)
#Do a PCoA from the X.jac matrix with vegan::cmdscale function
X.jac.pcoa <- cmdscale(X.jac, eig = TRUE)
# This ensures consistency between ordination and envfit
colnames(X.jac.pcoa$points) <- paste0("Dim", 1:ncol(X.jac.pcoa$points))
#Next extract the coordinate points for plotting with ggplot the two first ordination axis
X.jac.pcoa.df <- data.frame(X.jac.pcoa$points)
#Add the metadata to colour points by sample variables
X.jac.pcoa.df <- merge(X.jac.pcoa.df, sm$sample_info, by.x = 0, by.y = "Spygen", all.x=TRUE)
X.jac.pcoa.df$Habitat[X.jac.pcoa.df$Habitat == "DeepSlope"] <- "DeepSlope150"
X.jac.pcoa.df$Habitat <- factor(X.jac.pcoa.df$Habitat, levels= c("Bay","Lagoon", "Soft_back_reef", "Reef_outer_slope","Summit50", "DeepSlope150", "Summit250", "Summit500"))

X.jac.pcoa.df$Zone <- ifelse(X.jac.pcoa.df$Habitat %in% c("Bay","Lagoon","Reef_outer_slope", "Soft_back_reef"), "Shallow",
                             ifelse(X.jac.pcoa.df$Habitat %in% c("Summit50", "DeepSlope150"), "Middle",
                                    ifelse(X.jac.pcoa.df$Habitat %in% c("Summit250", "Summit500"), "Deep",
                                           NA)))

X.jac.pcoa.df$Zone <- factor(X.jac.pcoa.df$Zone, levels=c("Shallow", "Middle", "Deep"))

# Prepare environmental variables
numvars <- sapply(env_eDNA_df, is.numeric)
numvars <- numvars[numvars]

num_env_vars <- env_eDNA_df[, c("Spygen",names(numvars))]
rownames(num_env_vars) <- num_env_vars$Spygen ; num_env_vars <- num_env_vars[,-1]
num_env_vars <- num_env_vars[, colSums(is.na(num_env_vars)) == 0]

# select few variables for env fitting
num_env_vars <- num_env_vars[, c("Depth", "Salinity", "Chla", "seafloorTemp", "SSTmean", "TravelTime", "ReefMinDist.m", "otu_richness", "Latitude", "Longitude")]

# Rename otu_richness to MOTU_richness
colnames(num_env_vars)[colnames(num_env_vars) == "otu_richness"] <- "MOTU_richness"

# Get the row names from the ordination points
ord_samples_jac <- rownames(X.jac.pcoa$points)

# Match environmental variables to ordination sample order
num_env_vars_jac <- num_env_vars[match(ord_samples_jac, rownames(num_env_vars)), ]

# env fitting for Jaccard distance
ef.jac <- envfit(X.jac.pcoa, num_env_vars_jac, permutations=999)
ef.jac.df <- as.data.frame(ef.jac$vectors$arrows)
ef.jac.df$r2 <- ef.jac$vectors$r
ef.jac.df$pval <- ef.jac$vectors$pvals
ef.jac.df$Dim1 <- ef.jac.df$Dim1*sqrt(ef.jac.df$r2)
ef.jac.df$Dim2 <- ef.jac.df$Dim2*sqrt(ef.jac.df$r2)
ef.jac.df <- ef.jac.df[order(ef.jac.df$r2, decreasing=TRUE),] 
ef.jac.df$env_variable <- rownames(ef.jac.df)

# get only significant relationship between env_variable and the community patterns shown in the ordination
ef.jac.df <- ef.jac.df[ef.jac.df$pval < 0.05, ]

ef.jac.plot <- ggplot(X.jac.pcoa.df) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(mapping = aes(x = Dim1, y = Dim2, colour=Habitat, shape=Zone), size=3.5) +
  geom_segment(data = ef.jac.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "blue") +
  geom_text_repel(data = ef.jac.df, aes(x = Dim1, y = Dim2, label = rownames(ef.jac.df)), size = 5, fontface = "bold")+ 
  xlab(paste("PCo1_Jac [", signif((X.jac.pcoa$eig[1]/sum(X.jac.pcoa$eig)),3)*100,"%]", sep="")) +
  ylab(paste("PCo2_Jac [", signif((X.jac.pcoa$eig[2]/sum(X.jac.pcoa$eig)),3)*100,"%]", sep="")) +
  scale_color_manual(values = c("Bay" = "#5ae6ab", "Lagoon" = "#88d941", "Soft_back_reef" = "#2e7d00","Reef_outer_slope" = "#4e8273", "Summit50" = "#ffe699", "DeepSlope" = "#d79c3b", "Summit250" = "#4a6a94", "Summit500" = "#1a3250"))+
  stat_ellipse(aes(x=Dim1, y=Dim2, fill=Zone), geom="polygon", alpha=0.15) +
  scale_fill_manual(values = c(
    "Shallow" = "#5ae6ab",
    "Middle" = "#ffe699",
    "Deep" = "#25456B"
  )) +
  guides(
    color = "none",
    shape = guide_legend(nrow = 1, override.aes = list(alpha = 0.5))  
  ) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title= element_text(size=16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16, margin = margin(b = 15)),
        legend.position = "right",
        legend.box = "vertical",
        legend.spacing.y = unit(0.4, "cm")
  )

# set the figure label for the combined plot
p4 <- ef.jac.plot

layout ="
 AAABBB
 AAABBB
 AAABBB
 CCCCCC
 CCCCCC
 CCCCCC
 CCCCCC
 DDDDDD"

# set the path to save the figure
path <- paste0(script_dir, "/")

# Create the directory if it doesn't exist
if (!dir.exists(path)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# save the combined figure to a pdf
pdf(file=paste0(path, "Figure1.pdf"), h=20, w=22)
wrap_elements(full = p1) +  
  wrap_elements(full = p4) + p2 + p3 + 
  plot_layout(design = layout) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = "bold", size=24))
dev.off()