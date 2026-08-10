# Data folder

The data files required to run the pipeline are **not publicly distributed** in this repository due to New Caledonian legislation on sensitive environmental data. Access requires a Data Use Agreement (DUA), which will be granted systematically for reproducibility purposes — contact the corresponding author. A processed version of the dataset is archived on [Zenodo](https://doi.org/10.5281/zenodo.21502975).

For reproducibility purposes, `smX_pres_abs_matrix.rda` — the presence/absence matrix used as input to the pipeline — is already available in this folder. With this file alone, all analyses are reproducible except Figure 1 and Figure 6.

To reproduce the pipeline completely, gain access via *Zenodo* and place the following files in this `data/` directory:

| File | Description |
|------|-------------|
| `eDNA_Data.csv` | Raw eDNA relative abundance data: long-format abundance table (one row per MOTU × station combination). Columns: `Station`, `Spygen`, `Site`, `Habitat`, `Depth`, `Latitude`, `Longitude`, `mean_pcr_count_reads`, `best_taxonomic_assignment`. |
| `Environmental_Variables_Data.csv` | Raw environmental variables per station. Columns: `Station`, `Latitude`, `Longitude`, `Site`, `Habitat`, `Depth`, `Salinity`, `SSTmean`, `seafloorTemp`, `Chla`, `TravelTime`, `ReefMinDist.m`. |
| `smX_pres_abs_matrix.rda` | Samples × MOTU presence/absence matrix (filtered at 3% prevalence), with `Habitat`/`Zone`/`hab_inoff` columns. Used throughout the overall pipeline (Predomics, ScaleNet, figures). |
| `Sample_Data.csv` | Per-sample metadata (station, habitat, environmental covariates, MOTU richness). Used to reproduce Figure 1. |
| `Taxonomy_Data.csv` | Per-MOTU taxonomy (kingdom → species). Used to reproduce Figures 1 and 6. |
