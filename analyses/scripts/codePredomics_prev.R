# =============================================================================
# Script name: codePredomics_prev.R
# Author: Estephe Kana & Eugeni Belda & Edi Prifti
# Date created: 2025-12-10
# Purpose: Run Predomics (bininter/terinter) pairwise habitat comparisons on
#          eDNA presence/absence data filtered by prevalence.
# Inputs:  data/smX_pres_abs_matrix.rda
# Outputs: analyses/analysis_outputs/<algorithm>_output_data/
#            <algorithm>_Predomics_all_analyses_overall_data_<group>_prev_<N>.Rda
#          analyses/analysis_outputs/sessionInfo_predomics_<algorithm>_<date>.txt
# =============================================================================

# -----------------------------------------------------------------------------
# Commands to run the script (from repo root):
#   export REPO_ROOT="/data/projects/aime/analyses/seamount/metabarcoding/"
#   export DATA_DIR="/data/projects/aime/data/seamount/metabarcoding"
#
#   From REPO_ROOT, run the following commands:
#
#   Rscript analyses/scripts/codePredomics_prev.R terinter 3 strat_group
#   Rscript analyses/scripts/codePredomics_prev.R bininter 3 strat_group
#
# To run pairwise comparisons of habitats (Shallow, Middle, Deep) instead of Inshore/Offshore:
#   Rscript analyses/scripts/codePredomics_prev.R terinter 3 hab_inoff
#   Rscript analyses/scripts/codePredomics_prev.R terinter 3 Shallow
#   Rscript analyses/scripts/codePredomics_prev.R terinter 3 Middle
#   Rscript analyses/scripts/codePredomics_prev.R terinter 3 Deep
# -----------------------------------------------------------------------------

# Check for required packages
# Note: predomics — remotes::install_github("eprifti/predomics")
required_pkgs <- c("vegan", "predomics")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))

# Load required libraries
library(vegan)
library(predomics)

set.seed(42)

# define arguments for the script
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2)
  stop("Usage: Rscript codePredomics_prev.R <language> <prevalence> [<habitat_group>]\n",
       "  language: 'terinter' or 'bininter'\n",
       "  prevalence: numeric 0-100 (e.g. 3)\n",
       "  habitat_group: 'strat_group' | 'hab_inoff' | 'Shallow' | 'Middle' | 'Deep'  (default: hab_inoff)")

algorithm_language <- args[1]
prevalence_rate    <- suppressWarnings(as.numeric(args[2]))
habitat_group      <- if (length(args) >= 3) args[3] else "hab_inoff"

if (!algorithm_language %in% c("terinter", "bininter"))
  stop("'language' must be 'terinter' or 'bininter', got: ", algorithm_language)
if (is.na(prevalence_rate) || prevalence_rate < 0 || prevalence_rate > 100)
  stop("'prevalence' must be a number between 0 and 100, got: ", args[2])

# Derive repo root from script location
script_path  <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
script_dir   <- dirname(script_path)         # <repo>/analyses/scripts
repo_root    <- dirname(dirname(script_dir)) # <repo>/
data_dir     <- file.path(repo_root, "data")

# Allow override via environment variable (set on the remote server)
repo_root <- Sys.getenv("REPO_ROOT", unset = repo_root)
data_dir  <- Sys.getenv("DATA_DIR",  unset = file.path(repo_root, "data"))
analyses_dir <- file.path(repo_root, "analyses")
source(file.path(script_dir, "utils.R"))

# load dataset: samples x MOTU presence/absence, with Habitat/Zone/hab_inoff
load(file.path(data_dir, "smX_pres_abs_matrix.rda"))
rownames(smX_pres_abs_matrix) <- smX_pres_abs_matrix$Station
base_meta_cols <- c("Station", "Habitat", "Zone", "hab_inoff")
edna_presenceAbsence <- as.matrix(smX_pres_abs_matrix[, !colnames(smX_pres_abs_matrix) %in% base_meta_cols])

# Already prevalence-filtered upstream (make_db_object_clean.R); prevalence_rate only labels the output filename.
# filtered_edna_presenceAbsence <- get_sample_by_prevalence(edna_presenceAbsence, prevalence_rate)
filtered_edna_presenceAbsence <- edna_presenceAbsence

acc <- data.frame(filtered_edna_presenceAbsence)
acc$Station <- rownames(acc)

# merge presence/absence table with sample info (Zone/hab_inoff already derived from Habitat)
acc <- merge(acc, smX_pres_abs_matrix[, base_meta_cols], by = 'Station', all.x = TRUE, sort = FALSE)
rownames(acc) <- acc$Station

# rename Habitat/Zone to match the habitat-group variable names used below
colnames(acc)[colnames(acc) == "Habitat"] <- "hab"
colnames(acc)[colnames(acc) == "Zone"]    <- "strat_group"

# define comparisons and variables for analysis based on chosen habitat group

if (habitat_group == "hab_inoff") {
  hab_type <- "hab_inoff"
  comp     <- combn(x=unique(acc$hab_inoff), m=2, simplify = FALSE)
} else if (habitat_group == "strat_group") {
  hab_type <- "strat_group"
  comp     <- combn(x=unique(acc$strat_group), m=2, simplify = FALSE)
} else if (habitat_group == "Shallow"){
  hab_type <-"hab"
  comp     <- combn(x=c('Bay', 'Lagoon', 'OuterSlope', 'BackReef'), m = 2, simplify = FALSE)
} else if (habitat_group == "Middle"){
  hab_type <- "hab"
  comp     <- combn(x=c('DeepSlope150', 'Seamount50'), m = 2, simplify = FALSE)
} else if (habitat_group == "Deep"){
  hab_type <- "hab"
  comp     <- combn(x=c('Seamount250', 'Seamount500'), m = 2, simplify = FALSE)
}else {
  stop("Invalid habitat group specified. Use 'hab_inoff' or 'group'.")
}

print(paste("Selected habitat grouping:", habitat_group))

# Function to prepare data for training

prepare_data_for_analysis <- function(sample_data, habitat_type, hab1, hab2) {

  if (hab1 == hab2)
    stop("It seems that the habitats you entered are the same.")

  meta_cols <- c('Station', 'hab', 'hab_inoff', 'strat_group')

  X <- sample_data[, !colnames(sample_data) %in% meta_cols, drop = FALSE]
  sample_X <- X[, colSums(X > 0) > 0, drop = FALSE]

  unprevalent_species <- setdiff(colnames(X), colnames(sample_X))
  sample_data <- sample_data[, !colnames(sample_data) %in% unprevalent_species, drop = FALSE]
  sample_data <- sample_data[rownames(sample_data) %in% rownames(sample_X), ]
  sample_data <- sample_data[rowSums(sample_data > 0) > 0, ]

  sample_data.class <- ifelse(sample_data[[habitat_type]] == hab1, -1,
                               ifelse(sample_data[[habitat_type]] == hab2, 1, NA))

  sample_data.X <- sample_data[, !colnames(sample_data) %in% meta_cols, drop = FALSE]

  return(list(sample.info = sample_data, sample.X = sample_data.X, sample.class = sample_data.class))
}

# Function for Predomics analysis
predomics_analysis <- function(X, y, algorithm) {

  # set the learner (nCores = 1 for full reproducibility)
  if (algorithm == "terinter") {
    clf <- terBeam(language = 'terinter', nCores = 1, objective = "auc", seed = 20, plot = TRUE)
  } else if (algorithm == "bininter") {
    clf <- terBeam(language = 'bininter', nCores = 1, objective = "auc", seed = 20, plot = TRUE)
  } else {
    stop("algorithm must be 'terinter' or 'bininter', got: ", algorithm)
  }

  #### Run the learner (training) ####
  res_clf <- predomics::fit(X = t(X), y = y, clf = clf, cross.validate = TRUE, nfolds = 5)
  
  # Build a master list to save predomics analysis results
  results_list <- list()
  
  # store data for training and testing
  results_list$Comp_data <- list(X= X, y= y)
  
  # save clf and classification results in the master list
  results_list$fit <- res_clf
  results_list$clf <- clf
  
  # save digest of results
  res_clf.dig <- digest(obj = res_clf, penalty = 0.75/100, plot = TRUE)
  results_list$digest <- res_clf.dig
  
  #### Family of Best Models (FBM) ####
  
  # get the population of models scrambled by model size
  pop <- modelCollectionToPopulation(res_clf$classifier$models)
  pop.df <- populationToDataFrame(pop)
  results_list$model_pop <- pop.df
  
  # select the best
  fbm <- selectBestPopulation(pop)
  fbm.df <- populationToDataFrame(fbm)
  results_list$fbm <- fbm.df
  
  #### Family of Best Models (FBM) - Feature Annotation ####
  
  fa <- makeFeatureAnnot(pop = fbm,
                         X = t(X),
                         y = y,
                         clf = clf)
  
  results_list$featureAnnotFBM <- fa
  
  #### Feature importance ####
  
  feat1.import <- mergeMeltImportanceCV(list.results = list(terBeam = res_clf),
                                        filter.cv.prev = 0,
                                        min.kfold.nb = FALSE,
                                        learner.grep.pattern = "*",
                                        #nb.top.features = 50,
                                        feature.selection = rownames(fa$pop.noz),
                                        scaled.importance = TRUE,
                                        make.plot = TRUE,
                                        cv.prevalence = FALSE)
  
  
  feat2.import <- mergeMeltImportanceCV(list.results = list(terBeam = res_clf),
                                        filter.cv.prev = 0,
                                        min.kfold.nb = FALSE,
                                        learner.grep.pattern = "*",
                                        nb.top.features = 148,
                                        #feature.selection = rownames(fa$pop.noz),
                                        scaled.importance = TRUE,
                                        make.plot = TRUE,
                                        cv.prevalence = FALSE)
  
  
  results_list$FI_fmbFeats <- feat1.import
  results_list$FI_allfeat <- feat2.import
  return(results_list)
}
# Function to save predomics results object
save_pred_results <- function(results_pred, X, ilevels){
  
  # get only species from FBM
  subdf <- results_pred$FI_allfeat$summary[, c('feature', 'value','sign')]
  # replace NA values by 0
  subdf$value[is.nan(subdf$value)] <- 0
  # round feature importance value
  subdf$value = round(subdf$value, digits = 3)
  # get the class vector of results
  subdf$sign = ifelse(subdf$sign == -1, ilevels[[1]], ilevels[[2]])
  # get the list of FBM species
  fbm.species = results_pred$FI_fmbFeats$summary$feature
  # add a column to indicate if a species is an indicator or not
  subdf$IsIndsp = as.integer(subdf$feature %in% fbm.species)
  colnames(subdf) <- c("feature","featureImportance", "class","IsIndSp")
  #compute prevalence of each species
  subdf$prevalence <- results_pred$FI_allfeat$fprev$value
  results_pred$pred_out_fbm <- subdf
  
  save_results <- results_pred
  
  return(save_results)
  
}
# compute permanova analysis on presence/absence data (Jaccard distances)
permanova_analysis <- function(sample.info, sample.X, sample.class, results_pred, hab_type, ilevels, algorithm_language) {

  # Get features from FBM
  fbm.species <- results_pred$FI_fmbFeats$summary$feature

  # Subset the presence/absence table to these species
  subdf.features.df <- sample.X[, fbm.species]

  # Exclude samples with no species
  subdf.features.df <- subdf.features.df[rowSums(subdf.features.df) > 0, ]

  # Get class according to rows of subdf.features.df
  subdf.features.df.class <- sample.info[rownames(subdf.features.df),]
  subdf.features.df.class <- ifelse(subdf.features.df.class[[hab_type]] == ilevels[[1]], -1,
                                    ifelse(subdf.features.df.class[[hab_type]] == ilevels[[2]], 1, NA))

  # Check for any NA classes
  na_indices <- which(is.na(subdf.features.df.class))
  if (length(na_indices) > 0) {
    warning(paste("NA classes found at indices:", paste(na_indices, collapse = ", ")))
  }

  # Compute Jaccard distances
  subdf.features.df.jaccard <- vegdist(subdf.features.df, method = "jaccard", binary = TRUE)

  # Do the PERMANOVA with subset species
  subdf.meta <- data.frame(sample = rownames(subdf.features.df), class = subdf.features.df.class)
  set.seed(100)
  subdf.features.df.jaccard.adonis <- adonis2(subdf.features.df.jaccard ~ class, data = subdf.meta)

  # Extract results
  subdf.features.df.jaccard.adonis <- data.frame(subdf.features.df.jaccard.adonis)[1, , drop = FALSE]
  subdf.features.df.jaccard.adonis$comparison <- paste(ilevels, collapse = "_")
  subdf.features.df.jaccard.adonis$data <- "pres/abs"
  subdf.features.df.jaccard.adonis$source <- paste0(algorithm_language, "Species")
  subdf.features.df.jaccard.adonis$features <- ncol(subdf.features.df)

  # Compute Jaccard distances for the entire community
  ilevels.df.jaccard <- vegdist(sample.X, method = "jaccard", binary = TRUE)
  alldf.meta <- data.frame(sample = rownames(sample.X), class = sample.class)
  set.seed(100)
  ilevels.df.jaccard.adonis <- adonis2(ilevels.df.jaccard ~ class, data = alldf.meta)

  # Extract results
  ilevels.df.jaccard.adonis <- data.frame(ilevels.df.jaccard.adonis)[1, , drop = FALSE]
  ilevels.df.jaccard.adonis$comparison <- paste(ilevels, collapse = "_")
  ilevels.df.jaccard.adonis$data <- "pres/abs"
  ilevels.df.jaccard.adonis$source <- paste("allSpecies", algorithm_language, sep = '_')
  ilevels.df.jaccard.adonis$features <- ncol(sample.X)

  # Combine results
  rbind(ilevels.df.jaccard.adonis, subdf.features.df.jaccard.adonis)
}


print("##################### starting predomics analyses on data in presence/absence ###############################")

####################
## Analyses on presence/absence data
####################

predout.bin <- list()
adonis_pred.bin <- list()

for(i in 1:length(comp))
{
  # print(i)
  ##Get the levels to compare
  ilevels <- comp[[i]]
  print(ilevels)
  ##Get the data limited to the levels to compare
  ilevels.df <- acc[acc[[hab_type]] %in% ilevels,]

  # prepare data for analysis
  data <- prepare_data_for_analysis(ilevels.df, habitat_type = hab_type, ilevels[[1]], ilevels[[2]])
  # get elements returned
  ilevels.df= data$sample.info
  ilevels.df.X= data$sample.X
  ilevels.df.class= data$sample.class

  # compute predomics analysis
  predomics_res_list= predomics_analysis(ilevels.df.X, ilevels.df.class, algorithm_language)

  # save predomics results
  predout.bin[[paste(ilevels, collapse = "_")]] <- save_pred_results(predomics_res_list, ilevels.df.X, ilevels)
  # Permanova analysis
  adonis_pred.bin[[paste(ilevels, collapse = "_")]] <- permanova_analysis(ilevels.df,ilevels.df.X, ilevels.df.class, predomics_res_list, hab_type = hab_type, ilevels, algorithm_language)
}

##Get the pred_out table
predout.bin.sub <- lapply(predout.bin, function(x){x[["pred_out_fbm"]]})
for(i in names(predout.bin.sub))
{
  predout.bin.sub[[i]][,"comparison"] <- i
}
predout.bin.sub <- do.call("rbind", predout.bin.sub)
predout.bin.sub$source <- algorithm_language
predout.bin.sub$data <- "pres/abs"


#### Save the results
save_path <- file.path(analyses_dir, "analysis_outputs", paste0(algorithm_language, "_output_data"))

# Create the directory if it doesn't exist
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
}

# Define the filename  
filename <- paste(save_path, paste(algorithm_language, "Predomics_all_analyses_overall_data", habitat_group, "prev", paste0(prevalence_rate,".Rda"), sep = "_"), sep = "/") 

save(adonis_pred.bin, #adonis results all comparisons; presence/absence data
     predout.bin, # indval results all comparisons, presence/absence data
     predout.bin.sub, file = filename)

# Print the successful message
message(paste("Analysis is successfully completed and results are saved to:", filename))

# Record session for reproducibility
session_file <- file.path(analyses_dir, "analysis_outputs", paste0("sessionInfo_predomics_", algorithm_language, "_", Sys.Date(), ".txt"))
writeLines(capture.output(sessionInfo()), session_file)
