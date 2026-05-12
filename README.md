# EcoEdnaNet: A network-based pipeline for reconstructing taxa co-occurrence patterns from environmental DNA in marine fish communities

## 🧭 Overview

***EcoEdnaNet*** is a repository that contains the code and analysis pipeline for reconstructing **co-occurrence networks of marine fish assemblages** from environmental DNA (eDNA) metabarcoding data collected around the seamounts and coastal habitats of **New Caledonia** (southwestern Pacific Ocean).

The project combines an interpretable machine learning framework (**Predomics**) for indicator MOTU identification with a network inference approach (**ScaleNet**) to characterise fish community structure across a depth and habitat gradient ranging from coastal bays to deep seamount summits (500 m). Module enrichment analysis using Gene Set Enrichment Analysis (**GSEA**) links network community structure to fine-grained habitat associations.

## 🌟 Main findings

- Predomics identified **55 indicator MOTUs** across three pairwise zone comparisons (Deep vs. Shallow, Deep vs. Middle, Middle vs. Shallow), with six taxa retained in all three comparisons, capturing the full depth gradient.
- ScaleNet reconstructed a co-presence network of **261 nodes and 579 edges** from 318 MOTUs filtered at 3% prevalence; 31 of the 55 indicator MOTUs were recovered in the network.
- Indicator taxa occupy structurally marginal positions in the network (significantly lower betweenness centrality than non-indicator taxa), acting as zone-specific diagnostic signatures rather than community connectivity hubs.
- Fast-Greedy modularity partitioned the network into **12 modules** that mirror depth zonation and fine-grained habitat structure beyond the three broad zones.
- GSEA confirmed that modules are non-randomly enriched for taxa associated with specific habitats (Bay, Lagoon, Soft\_back\_reef, Reef\_outer\_slope, Summit50, DeepSlope150, Summit250, Summit500), revealing that eDNA co-presence captures habitat-associated trophic guilds rather than stochastic co-detection.

## 🔄 Project pipeline

![Project pipeline](analyses/figures/Figure2/Figure2.svg)

The analysis proceeds in five sequential steps. All scripts derive their working directory from their own file path, so they can be run from any location. Place all three raw data files in `data/` before starting (see [data/README.md](data/README.md)).

| Step | Script | Output |
|------|--------|--------|
| 0. Build dataset | `analyses/scripts/make_db_object.R` | `data/seamount_integrated_dataset.rda` |
| 1. Prepare tables | `analyses/scripts/save_table_data.R` | Presence/absence tables for ScaleNet |
| 2. Predomics analyses | `analyses/scripts/codePredomics_prev.R` | Indicator MOTU results (`.Rda`) |
| 3. Network inference | `analyses/scripts/scalenet_network_inference.R` | `graph_data/graph_data_ecorr50_all_strat_<date>.rda` |
| 4. Figures | `analyses/figures/FigureN/FigureN_code.R` | Publication figures (PDF) |

## ⚙️ Methods

- **Predomics `bininter` / `terinter`**: Interpretable machine learning framework that identifies indicator MOTUs discriminating pairs of depth zones (Shallow, Middle, Deep). `bininter` assigns strictly positive unit coefficients (+1), making each taxon a positive presence-diagnostic indicator; `terinter` extends this to signed unit coefficients (+1 or −1), allowing absence-diagnostic indicators. Models are trained with 5-fold cross-validation optimising AUC; the Family of Best Models (FBM) selects parsimonious indicator sets. See [Prifti et al., 2020](https://academic.oup.com/gigascience/article/9/3/giaa010/5801229).

- **ScaleNet**: A network inference method combining ensemble Bayesian and ARACNE approaches (`scs()` function). Applied to a consensus presence/absence table at 3% prevalence threshold across all stations. Edge filtering at |ecorr| > 0.5 retains ecologically supported co-occurrences. See [Prifti et al., 2021](https://link.springer.com/article/10.1186/s12859-016-1308-y).

- **Chi-squared post-hoc tests**: Each MOTU is independently tested for differential occurrence across zones and habitats using a χ² test; significant taxa (Benjamini–Hochberg adjusted *p* < 0.05) are assigned to their associated zone/habitat via post-hoc pairwise χ² tests. These assignments annotate network nodes.

- **Fast-Greedy modularity**: Community detection method applied to the co-presence network using the Fast-Greedy algorithm (`igraph`). Selected over Louvain (equally high modularity, same number of modules) for its full determinism and reproducibility without random initialisation.

- **GSEA**: Gene Set Enrichment Analysis (Piano package) tests whether each Fast-Greedy module is non-randomly enriched for taxa associated with specific zones or habitats. Taxa are ranked by Benjamini–Hochberg adjusted *p*-value with the sign of the Pearson residual as direction indicator (1,000 permutations). See [Väremo et al., 2013](https://pubmed.ncbi.nlm.nih.gov/23444143/).

- **PERMANOVA**: Community-level significance of habitat separation was assessed with `adonis2` (Bray-Curtis for abundance, Jaccard for presence/absence), run both on all MOTUs and on indicator MOTUs identified by Predomics.


## 🐟 Data description

This study uses an environmental DNA dataset collected by Baletaud et al. (2023) and Mathon et al. (2025) in marine ecosystems of New Caledonia. Water samples collected at each station were filtered and amplified using fish-specific 12S rRNA primers. MOTUs (Molecular Operational Taxonomic Units) were identified by matching amplicon sequences against a curated reference database. Read counts were normalised per PCR replicate using the `normFreqTC` function from the `momr` package.

![Dataset overview — map, species richness heatmap, and PCoA ordination](analyses/figures/Figure1/Figure1.svg)

### Data overview

| **Feature** | **Description** |
|-------------|-----------------|
| **Location** | New Caledonia seamounts and coastal reefs (southwestern Pacific Ocean) |
| **Sampling method** | Environmental DNA water filtration + 12S rRNA metabarcoding |
| **Molecular marker** | 12S rRNA (fish-specific primers) |
| **Number of MOTUs** | 967 |
| **Number of habitat types** | 8 (Bay, Lagoon, Reef_outer_slope, Soft_back_reef, Summit50, DeepSlope150, Summit250, Summit500) |
| **Depth range** | 0 – 500 m |
| **Depth strata (for analyses)** | Shallow (Bay, Lagoon, Reef_outer_slope, Soft_back_reef), Middle (Summit50 / DeepSlope150), Deep (Summit250 / Summit500) |
| **Environmental covariates** | 13 (Sea Surface Temperature (SST) [SSTmax, SSTmean, SSTmin, SSTsd], SeafloorTemp, Salinity, chlorophyll-a (chla), current velocity (EastwardVelocity, NorthwardVelocity), reef Min distance (ReefMinDist.m), land distance, OTU richness, travel time to Nouméa) |
| **Abundance normalisation** | `normFreqTC` (per-PCR relative frequency normalisation) |
| **Prevalence filter (analyses)** | 3% of samples |

### Ecological context

- **Habitat gradient**: Stations span a continuous shallow–deep gradient, from sheltered bays and lagoons to exposed outer reefs and isolated deep seamount summits.
- **Community structure**: Beta-diversity analyses (Bray-Curtis PCoA) reveal clear separation of shallow coastal from deep seamount assemblages, motivating the three-stratum zonation (Shallow / Middle / Deep) used in Predomics comparisons.
- **Network structure**: The ScaleNet co-occurrence network captures habitat-specific species associations and depth-partitioned modules consistent with known biogeographic patterns in the southwestern Pacific.

## 💾 Data availability

The processed eDNA dataset used here are archived in the public repository [Zenodo](DOI: [add DOI once deposited]). The New Caledonian legislation regarding sensitive environmental data does not permit unrestricted public access. Accordingly, access to the data will require a Data Use Agreement (DUA), which will be systematically granted for reproducibility purposes.


## 📁 Repository structure

| Folder | Description |
|--------|--------------|
| `analyses/analysis_outputs/` | Contains the outputs of Predomics and ScaleNet analyses |
| `analyses/figures/` | Contains all the vizualization with their code as presented in the paper linked to this project |
| `analyses/files/` | Contains others output files of the pipeline|
| `analyses/scripts/` | Contains all the scripts for curating data, indicator taxa identification and network inference analyses |
| `data/` | Provides instructions to get the dataset used in this project |


## 🧩 Install main required packages for analyses

```r
# ── CRAN packages ──────────────────────────────────────────────────────────────
install.packages(c(
  # Data handling
  "readr", "readxl", "data.table",
  # Plotting
  "ggplot2", "patchwork", "viridis", "ggrepel", "ggpubr",
  "cowplot", "ggalluvial", "ggh4x", "ggspatial",
  "maps", "rnaturalearth", "sf", "scatterplot3d", "ggplotify",
  # Analysis
  "vegan", "igraph", "reshape2", "plyr",
  "chisq.posthoc.test", "dplyr", "tidyr"
))

# ── Bioconductor packages ──────────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("piano", "ggtree", "ggtreeExtra", "ggnewscale"))

# ── GitHub packages ────────────────────────────────────────────────────────────
install.packages("remotes")

# momr — read-count normalisation
remotes::install_github("eprifti/momr")

# Predomics — interpretable machine learning for omics data
remotes::install_github("predomics/predomicspkg", dependencies = TRUE)

# ScaleNet — network inference
remotes::install_github("UMMISCO/scalenet")
```

## 🧠 Dependencies

- **R version**: 4.4.0 or later
- **Key packages**: `predomics` (1.1.0), `scalenet` (1.2.3), `momr`, `igraph` (2.2.1), `vegan` (2.7.2), `piano` (2.22.0), `ggtree`, `ggplot2`, `patchwork`

## 💰 Acknowledgments

*Funding:* This work was supported by the **AIME (Artificial Intelligence for Marine Ecosystems)** project, co-funded by the **French National Research Agency (ANR)** and the **French Development Agency (AFD)**.

## 📚 References

- [*Prifti et al., 2020. Interpretable and accurate prediction models for metagenomics data. GigaScience.*](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giaa010/5801229)

- [*Affeldt et al., 2016. Spectral consensus strategy for accurate reconstruction of large biological networks. BMC Bioinformatics.*](https://link.springer.com/article/10.1186/s12859-016-1308-y)

- [*Väremo et al., 2013. Enriching the gene set analysis of genome-wide data by incorporating directionality of gene expression and combining statistical hypotheses and methods. Nucleic Acids Research.*](https://pubmed.ncbi.nlm.nih.gov/23444143/)

- [*Oksanen et al., 2022. vegan: Community Ecology Package. R package.*](https://CRAN.R-project.org/package=vegan)

- [*Csardi & Nepusz, 2006. The igraph software package for complex network research. InterJournal.*](https://igraph.org)
