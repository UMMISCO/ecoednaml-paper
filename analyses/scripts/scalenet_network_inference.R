# =====================================================================================
# Script name: scalenet_network_inference.R
# Author: Estephe Kana & Edi Prifti & Eugeni Belda
# Date created: 2025-08-26
# Purpose: Reconstruct a co-presence network from eDNA data using ScaleNet, 
#     and annotate the network with species' habitat preferences and indicator status.
# Inputs: abundance_data_matrix, annotation data, predomics key species info
# Outputs: annotation_data.Rda
# =====================================================================================

# -----------------------------------------------------------------------------
# Commands to run the script (from repo root):
# Rscript analyses/scripts/scalenet_network_inference.R 50
# -----------------------------------------------------------------------------

# define arguments for the script
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: scalenet_network_inference.R <ecorr_percent>")
}

# Define threshold variables to select edges
ecorr_percent <- as.numeric(args[1])

# Check for required packages
# Note: scalenet — contact authors or install from source
#       momr     — remotes::install_github("eprifti/momr")
required_pkgs <- c("vegan", "reshape2", "plyr", "scalenet", "chisq.posthoc.test", "igraph")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))

# Load required libraries
library(vegan)
library(reshape2)
library(plyr)
library(scalenet)
library(chisq.posthoc.test)
library(igraph)

set.seed(42)

# define paths
script_path  <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
script_dir   <- dirname(script_path)         # <repo>/analyses/scripts
repo_root    <- dirname(dirname(script_dir)) # <repo>/
data_dir     <- file.path(repo_root, "data")
analyses_dir <- file.path(repo_root, "analyses")
source(file.path(script_dir, "utils.R"))

# load dataset
load(file.path(data_dir, "seamount_integrated_dataset.rda"))

# get samples filtered at 3% of prevalence of the total of samples
filtered_edna_abundance <- get_sample_by_prevalence(t(sm$X), 3)

# presence/absence table
filtered_edna_presenceAbsence <- decostand(filtered_edna_abundance, method = "pa")

# Function to compute chi-square + post-hoc
compute_chisq_post_hoc <- function(df, group_var = "Zone") {
  
  species        <- unique(df$species)
  group_levels   <- unique(na.omit(df[[group_var]]))  
  
  # Base result frame
  results <- data.frame(
    feature    = species,
    pval_chisq = NA_real_,
    padj_chisq = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (g in group_levels) results[[paste0("presence_count_", g)]] <- NA_integer_
  for (g in group_levels) results[[paste0("pval_PH_",        g)]] <- NA_real_
  for (g in group_levels) results[[paste0("padj_PH_",        g)]] <- NA_real_
  for (g in group_levels) results[[paste0("residual_",       g)]] <- NA_real_
  
  results$mean_presence_count <- NA_real_
  results$enriched_zones      <- NA_character_
  results$depleted_zones      <- NA_character_
  results$assigned_class      <- NA_character_
  
  # Loop over species
  for (i in seq_along(species)) {
    sp <- species[i]
    message("Processing: ", sp)
    
    sp_df <- df[df$species == sp, ]
    tbl   <- table(sp_df$presabs, sp_df[[group_var]])
    
    for (cls in c("0", "1")) {
      if (!cls %in% rownames(tbl)) {
        empty_row <- matrix(0L, nrow = 1, ncol = ncol(tbl),
                            dimnames = list(cls, colnames(tbl)))
        tbl <- rbind(tbl, empty_row)
      }
    }
    tbl <- tbl[c("0", "1"), , drop = FALSE]
    
    # Global chi-square
    chi <- suppressWarnings(chisq.test(tbl))
    results$pval_chisq[i] <- chi$p.value
    
    # Post-hoc test
    posthoc <- suppressWarnings(chisq.posthoc.test(tbl))
    
    p_row   <- posthoc[posthoc$Value == "p values"  & posthoc$Dimension == "1", ]
    res_row <- posthoc[posthoc$Value == "Residuals" & posthoc$Dimension == "1", ]
    
    presence_counts <- numeric(length(group_levels))
    
    for (j in seq_along(group_levels)) {
      g <- group_levels[j]
      
      results[[paste0("pval_PH_", g)]][i] <-
        as.numeric(gsub("\\*", "", as.character(p_row[[g]])))
      
      results[[paste0("residual_", g)]][i] <-
        as.numeric(as.character(res_row[[g]]))
      
      cnt <- if ("1" %in% rownames(tbl) && g %in% colnames(tbl)) tbl["1", g] else 0L
      results[[paste0("presence_count_", g)]][i] <- cnt
      presence_counts[j] <- cnt
    }
    
    results$mean_presence_count[i] <- mean(presence_counts)
  }
  
  # Adjust p-values
  results$padj_chisq <- p.adjust(results$pval_chisq, method = "BH")
  
  for (g in group_levels) {
    results[[paste0("padj_PH_", g)]] <-
      p.adjust(results[[paste0("pval_PH_", g)]], method = "BH")
  }
  
  # Classify enriched / depleted + assigned_class
  for (i in seq_len(nrow(results))) {
    enriched <- character(0)
    depleted <- character(0)
    
    for (g in group_levels) {
      padj <- results[[paste0("padj_PH_", g)]][i]
      res  <- results[[paste0("residual_", g)]][i]
      
      if (!is.na(padj) && !is.na(res) && padj < 0.05) {
        if (res > 0) enriched <- c(enriched, g)
        if (res < 0) depleted <- c(depleted, g)
      }
    }
    
    results$enriched_zones[i] <- if (length(enriched)) paste(enriched, collapse = ";") else "NS"
    results$depleted_zones[i] <- if (length(depleted)) paste(depleted, collapse = ";") else "NS"
    
    # ── Assigned class ──────────────────────────────────────────────────────
    if (!is.na(results$padj_chisq[i]) && results$padj_chisq[i] < 0.05 && length(enriched) > 0) {
      
      if (length(enriched) == 1) {
        # Only one enriched zone → assign directly
        results$assigned_class[i] <- enriched
        
      } else {
        # Multiple enriched zones → pick the one with the lowest adjusted p-value
        padj_enriched <- sapply(enriched, function(g) results[[paste0("padj_PH_", g)]][i])
        results$assigned_class[i] <- enriched[which.min(padj_enriched)]
      }
      
    } else {
      results$assigned_class[i] <- "NS"
    }
  }
  
  return(results)
}

# Add Habitat info presence/absence matrix

sample.info <- sm$sample_info
sample.info$Zone <- ifelse(sample.info$Habitat %in% c("Bay","Lagoon","Reef_outer_slope", "Soft_back_reef"), "Shallow",
                           ifelse(sample.info$Habitat %in% c("Summit50", "DeepSlope"), "Middle",
                                  ifelse(sample.info$Habitat %in% c("Summit250", "Summit500"), "Deep",
                                         NA)))

if (any(is.na(sample.info$Zone))) {
  warning(paste("Unmatched habitats:", 
                paste(unique(sample.info$Habitat[is.na(sample.info$Zone)]), collapse=", ")))
}
table(sample.info$Habitat, sample.info$Zone)

abund.df <- as.data.frame(t(filtered_edna_abundance))
abund.df$species <- rownames(abund.df)
abund.df <- reshape2::melt(abund.df, id.vars = "species")
abund.df$presabs <- ifelse(abund.df$value>0,1,0)
df.abund.sample.info <- merge(abund.df, sample.info[,c("Spygen","Zone", "Habitat")], by.x="variable", by.y="Spygen", all.x=TRUE)

# Compute chisq + PH at zone level
sp.chisq_posthoc <- compute_chisq_post_hoc(df.abund.sample.info)

# rename fishes' name by replacing " " with .
sp.chisq_posthoc$feature <- gsub(" ", ".", sp.chisq_posthoc$feature)

# Compute chisq + PH at Habitat level
sp.chisq_posthoc.habitat <- compute_chisq_post_hoc(df.abund.sample.info, group_var = "Habitat")

# rename fishes' name by replacing " " with .
sp.chisq_posthoc.habitat$feature <- gsub(" ", ".", sp.chisq_posthoc.habitat$feature)

data_path <- file.path(analyses_dir, "files", "txt", "presanceAbsence_table_prev_3.txt")
df <- read.delim(data_path)

results_path <- file.path(analyses_dir, "analysis_outputs", "scalenet_results")

#-------------------------------------------------------------------------------------------
# 2. ScaleNet NETWORK
#-------------------------------------------------------------------------------------------

# Remove recursively the content of results_path to avoid errors
if (dir.exists(results_path)) {
  # unlink(results_path, recursive = TRUE)
  message("Directory existed: ", results_path)
} else {
  message("Directory does not exist: ", results_path)
  # Create a consensus from Scalenet
  tmp <- scs(workspaceDir = results_path,
             argInData = as.data.frame(df),
             argReconsMeth = c("aracne", "bayes_hc"),
             argReconsMethInfo = list(aracne = list(ort = "n", eweight = "epresenceScore"), bayes_hc = list(ort = "y", eweight = "ecorr")),
             argEmbReconsParam = list(aracne = list(estimator="mi.mm", epsilon=0.001), bayes_hc = list(score="bde", restart=21), varPerc = 0.2),
             argPresFreqThresh = c(0.3, 0.5, 0.8),
             clean.workspace = FALSE,
             argDiscretize = TRUE,
             argVerbose = TRUE
  )
}

taxo <- sp.chisq_posthoc

fname <- paste(results_path, "consensusGraph/globalNet_presFreq0.8/edgesList.txt", sep = "/")

# load indicator species for strat comparisons for terinter model
load(file.path(analyses_dir, "analysis_outputs", "terinter_output_data", "terinter_Predomics_all_analyses_overall_data_strat_group_prev_3.Rda"))

# select results for indicator species identification on pres/abs data
indicSp_ter.df <- predout.bin.sub

# remove loaded variables to avoid override
rm(adonis_pred.bin, adonis_pred.maxn, predout.bin, predout.bin.sub, predout.maxn, predout.maxn.sub)

# load indicator species for strat comparisons for bininter model
load(file.path(analyses_dir, "analysis_outputs", "bininter_output_data", "bininter_Predomics_all_analyses_overall_data_strat_group_prev_3.Rda"))

# select results for indicator species identification on pres/abs data
indicSp_bin.df <- predout.bin.sub

# remove loaded variables
rm(adonis_pred.bin, adonis_pred.maxn, predout.bin, predout.bin.sub, predout.maxn, predout.maxn.sub)

# select key species for each model
keySpecies_bin <- indicSp_bin.df[indicSp_bin.df$IsIndSp==1, ]

keySpecies_ter <- indicSp_ter.df[indicSp_ter.df$IsIndSp==1, ]

allIndicSpecies <- rbind(keySpecies_bin, keySpecies_ter)

select_network_attributes <- function(edgesListPath, ecorr_percent, taxo) {

  # load the edge information for spectral3off2 network
  edges <- read.delim(edgesListPath, as.is = TRUE); dim(edges) # 50403 edges and 7 columns
  edges.raw <- edges
  # remove edges that do not have an orientation
  edges <- edges[!is.na(edges$eorientScore),]; dim(edges) # 1607 edges and 7 columns
  colnames(edges)[1:2] <- c("from","to")
  # edges$from <- gsub("_",":",edges$from)
  # edges$to <- gsub("_",":",edges$to)
  rownames(edges) <- paste(edges$from, edges$to, sep=" => ")
  
  # get the list of species with significance appearance across Deep, shallow and Middle
  SpList <- taxo[taxo$padj_chisq <0.05,]; dim(SpList)
  
  edges$SpeciesAppearance <- ifelse(edges$from %in% SpList$feature | edges$to %in% SpList$feature, TRUE, FALSE)
  
  edges$IsIndicSp_bin <- ifelse(edges$from %in% unique(keySpecies_bin$feature) | edges$to %in% unique(keySpecies_bin$feature), TRUE, FALSE)
  
  edges$IsIndicSp_ter <- ifelse(edges$from %in% unique(keySpecies_ter$feature) | edges$to %in% unique(keySpecies_ter$feature), TRUE, FALSE)
  
  dim(edges)
  # select links where we have a key indicator species or ecorr > ecorr_threshold
  ecorr_threshold <- ecorr_percent / 100
  edges.filt <- edges[abs(edges$ecorr) > ecorr_threshold,]
  dim(edges.filt)

  # ANNOTATION of the edges
  nodes <- unique(c(edges.filt$from, edges.filt$to))
  
  nodes.annot <- taxo[taxo$feature %in% nodes,]
  
  colnames(nodes.annot)[colnames(nodes.annot) == "feature"] <- "name"
  
  # return structured output
  return(list(
    edges.all   = edges,
    edges.filt  = edges.filt,
    nodes       = nodes,
    nodes.annot = nodes.annot
  ))
}

# get network attributes
network.attributes <- select_network_attributes(fname, ecorr_percent, taxo)
  
# Build the igraph network
# create the igraph object
network <- igraph::graph_from_data_frame(d = network.attributes$edges.filt, directed = TRUE, vertices = network.attributes$nodes.annot)

isolated_vertices <- V(network)[igraph::degree(network, mode = "all") == 0]
if (length(isolated_vertices)) {
  network <- delete_vertices(network, isolated_vertices)
}

# Calculate degree for all nodes
degAll <- igraph::degree(network, v = V(network), mode = "all")
V(network)$degree <- degAll   # add as vertex attribute

# Calculate betweenness for all nodes
betAll <- igraph::betweenness(network, v = V(network), directed = FALSE) / (((vcount(network) - 1) * (vcount(network)-2)) / 2)
betAll.norm <- (betAll - min(betAll))/(max(betAll) - min(betAll)); rm(betAll)

# Compute the closeness centraility
clos_cent <- igraph::closeness(network)

# Compute the eigenvector centrality of our network
eign_cent <- eigen_centrality(network)
eign_cent <- eign_cent$vector

nodes.annot <- network.attributes$nodes.annot

# modify first some variables in PH for habitat level
cols_to_rename <- c("pval_chisq", "padj_chisq", "mean_presence_count", "enriched_zones", "depleted_zones", "assigned_class")
colnames(sp.chisq_posthoc.habitat)[colnames(sp.chisq_posthoc.habitat) %in% cols_to_rename] <- paste0(cols_to_rename, ".habitat")

# merge nodes.annot with chisq Post Hoc results from Habitat level
nodes.annot <-  merge(nodes.annot, sp.chisq_posthoc.habitat, by.x="name", by.y="feature", all.x=TRUE)

# add centrality measures to node annotation dataframe
nodes.annot$degr_cent <- degAll
nodes.annot$betw_cent <- betAll.norm
nodes.annot$clos_cent <- clos_cent
nodes.annot$eign_cent <- eign_cent


# Align comparison with indics species
comparisons <- allIndicSpecies$comparison
names(comparisons) <- allIndicSpecies$feature

V(network)$comparison <- comparisons[V(network)$name]

# Create a vector to store the adjusted p-values for assigned classes
assigned_padj <- numeric(nrow(nodes.annot))

for (i in seq_len(nrow(nodes.annot))) {
  class <- nodes.annot$assigned_class[i]

  # If assigned_class is not "NS", get the corresponding post-hoc p-value
  if (class != "NS") {
    padj_col <- paste0("padj_PH_", class)
    assigned_padj[i] <- nodes.annot[[padj_col]][i]
  } else {
    # For NS or missing values, use the global p-value or set to NA
    assigned_padj[i] <- nodes.annot$padj_chisq[i]
  }
}

# Create a vector to store the adjusted p-values for assigned classes for Habitat
assigned_padj.habitat <- numeric(nrow(nodes.annot))

for (i in seq_len(nrow(nodes.annot))) {
  class <- nodes.annot$assigned_class.habitat[i]
  
  # If assigned_class is not "NS", get the corresponding post-hoc p-value
  if (class != "NS") {
    padj_col <- paste0("padj_PH_", class)
    assigned_padj.habitat[i] <- nodes.annot[[padj_col]][i]
  } else {
    # For NS or missing values, use the global p-value or set to NA
    assigned_padj.habitat[i] <- nodes.annot$padj_chisq.habitat[i]
  }
}

# Calculate -log10 of the adjusted p-values
nodes.annot$assigned_padj <- assigned_padj

nodes.annot$node_size <- -log10(assigned_padj)

# # Replace infinite values (from p-value = 0) with a reasonable maximum
nodes.annot$node_size[is.infinite(nodes.annot$node_size)] <- max(nodes.annot$node_size[!is.infinite(nodes.annot$node_size)], na.rm = TRUE)

# Optionally scale the sizes to a reasonable range for visualization
V(network)$size <- scales::rescale(nodes.annot$node_size, to = c(0.5, 5))

# Calculate -log10 of the adjusted p-values
nodes.annot$assigned_padj.habitat <- assigned_padj.habitat

nodes.annot$node_size.habitat <- -log10(assigned_padj.habitat)

# # Replace infinite values (from p-value = 0) with a reasonable maximum
nodes.annot$node_size.habitat[is.infinite(nodes.annot$node_size.habitat)] <- max(nodes.annot$node_size.habitat[!is.infinite(nodes.annot$node_size.habitat)], na.rm = TRUE)

# Optionally scale the sizes to a reasonable range for visualization
V(network)$size.habitat <- scales::rescale(nodes.annot$node_size.habitat, to = c(0.5, 5))

# Add new node/edge attributes based on the calculated node properties/similarities
V(network)$betweenness <- betAll.norm
V(network)$label       <- nodes.annot$name
V(network)$label.cex   <- -log10(nodes.annot$padj_chisq)
V(network)$label.color <- ifelse(nodes.annot$padj_chisq < 0.05, "#782832", "gray80")
V(network)$color       <- c("#06d6a0","#ffd166", "#25456B", "gray")[as.factor(factor(nodes.annot$assigned_class, levels=c("Shallow", "Middle", "Deep", "NS")))]

habitat_colors <- c(
  "Bay"              = "#5ae6ab",
  "Lagoon"           = "#88d941",
  "Soft_back_reef"   = "#2e7d00",
  "Reef_outer_slope" = "#4e8273",
  "Summit50"         = "#ffe699",
  "DeepSlope"        = "#d79c3b",
  "Summit250"        = "#4a6a94",
  "Summit500"        = "#1a3250",
  "NS"               = "gray"
)

V(network)$color.habitat <- habitat_colors[factor(nodes.annot$assigned_class.habitat,
                                          levels = names(habitat_colors))]

V(network)$frame.color <- ifelse(nodes.annot$name %in% keySpecies_bin$feature | nodes.annot$name %in% keySpecies_ter$feature, "firebrick1", NA)
V(network)$frame.width <- 3
V(network)$alpha       <- 0.6
V(network)$label.dist  <- 0.5
V(network)$label.font  <- ifelse(nodes.annot$padj_chisq < 0.05, 2, 1)

# Calculate edge properties and add to the network

#Calculate Dice similarities between all pairs of nodes
dsAll <- similarity(network, vids = V(network), mode = "all", method = "dice")
# The following function will transform a square matrix to an edge driven one and add values to each edge
F1 <- function(x) {data.frame(dice = dsAll[which(V(network)$name == as.character(x$from)), which(V(network)$name == as.character(x$to))])}

edges.filt.ext <- ddply(network.attributes$edges.filt, .variables=c("from", "to"), function(x) data.frame(F1(x))); dim(edges.filt.ext)

# merge edges data
edges.filt.ext <- merge(edges.filt.ext, network.attributes$edges.filt, all.x = TRUE)

# Add Edge betweness as edge attribute
edges.filt.ext$edge_betweeness <- edge_betweenness(network)

# edge width based on absolute ecorr
edge_widths <- abs(E(network)$ecorr)
edge_range <- range(edge_widths, na.rm = TRUE)
if (is.finite(diff(edge_range)) && diff(edge_range) > 0) {
  edge_widths_scaled <- 1 + 4 * (edge_widths - edge_range[1]) / diff(edge_range)
} else {
  edge_widths_scaled <- rep(1, length(edge_widths))
}

E(network)$similarity <- 0
E(network)[as.character(edges.filt.ext$from) %--% as.character(edges.filt.ext$to)]$similarity <- as.numeric(edges.filt.ext$dice)
E(network)$edge_width <- edge_widths_scaled
E(network)$color <- c("orange","darkgray")[as.factor(factor(sign(E(network)$ecorr), levels=c('-1','1')))]

# fix orientation coding
# 1 backward
# 2 forward
# 3 bidirected
E(network)$eorient[E(network)$eorient==-1] <- 1
E(network)$eorient[E(network)$eorient==1] <- 2
# E(network)$eorient[E(network)$eorient==0] <- 3

# Check the attributes
# Print number of nodes and edges
print(paste("There are",vcount(network),"nodes and",ecount(network),"edges"))

summary(network)

set.seed(100)
lay <- layout_with_fr(network)

n_nodes <- vcount(network)
n_edges <- ecount(network)
avg_degree <- mean(igraph::degree(network))
density <- edge_density(network)
n_components <- igraph::components(network)$no
diameter <- diameter(network, directed = FALSE, weights = NA)
transitivity <- transitivity(network)

# Build a multi-line title
main_title <- paste0(
  "Strat_group: ", "All_strat",
  " | Nodes: ", n_nodes,
  " | Edges: ", n_edges,
  " | Avg degree: ", round(avg_degree, 2),
  " | Density: ", round(density, 3),
  " | Components: ", n_components,
  " | Diameter: ", diameter,
  " | Transitivity: ", round(transitivity, 3)
)

pdf(file = file.path(analyses_dir, "figures", paste0("scalenet_network_ecorr", ecorr_percent, "_all_strat_group_", Sys.Date(), ".pdf")), w = 15, h = 10)

# Reduce margins to maximize plot area
par(mar = c(0, 0, 0, 10))  # Remove all margins

plot(network,
     layout = lay,
     vertex.label = V(network)$label,
     vertex.color= V(network)$color,
     vertex.shape= V(network)$shape,
     vertex.size= V(network)$size,
     vertex.frame.color= V(network)$frame.color,
     vertex.frame.width= V(network)$frame.width,
     edge.color= E(network)$color,
     edge.width= E(network)$width,
     asp= FALSE,  
     rescale=TRUE,
     edge.arrow.size = 0.3,
     vertex.label.cex = 0.5,
     vertex.label.dist= V(network)$label.dist,
     vertex.label.font = V(network)$label.font,
     vertex.label.color= V(network)$label.color
)

# Add legends (optional - comment out if you want no legends)
par(xpd=TRUE)

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

dev.off()

# Save the graph data
save_path <- file.path(analyses_dir, "files", "rdata", "graph_data")

# Create the directory if it doesn't exist
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
}

# Define the filename
filename <- paste(save_path, paste0("graph_data_ecorr", ecorr_percent,"_all_strat_",Sys.Date(),".rda"), sep = "/")

edges.all     <- network.attributes$edges.all
nodes         <- network.attributes$nodes

save(network, edges.all, nodes, nodes.annot, edges.filt.ext, lay, sp.chisq_zone=sp.chisq_posthoc, sp.chisq_habitat=sp.chisq_posthoc.habitat,
     file = filename)

message("Network saved to: ", filename)

# Record session for reproducibility
session_file <- sub("\\.rda$", "_sessionInfo.txt", filename)
writeLines(capture.output(sessionInfo()), session_file)