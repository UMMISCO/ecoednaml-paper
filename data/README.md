# Data folder

The raw data files required to run the pipeline are **not publicly distributed** in this repository due to New Caledonian legislation on sensitive environmental data. Access requires a Data Use Agreement (DUA), which will be granted systematically for reproducibility purposes — contact the corresponding author.

Once access is obtained, place the following files in this `data/` directory:

| File | Description |
|------|-------------|
| `eDNA_Data_SEAMOUNTS_REEF3.0_merged.csv` | Long-format eDNA abundance table (one row per MOTU × sample combination). Required columns: `Station`, `sequence`, `new_scientific_name_ncbi`, `SpeciesFB`, `mean_pcr_count_reads`. |
| `eDNA_SEAMOUNTS_REEF3.0_merged_Environmental_Variables_raw.csv` | Environmental metadata per station. Required columns: `Station`, `Site`, `Habitat`, `Substrate`, `Stratum`, `Date`, `Time`, `Depth`, `Sampling_depth`, `Latitude`, `Longitude`, `Project`, `Method`, plus oceanographic covariates (`EastwardVelocity`, `NorthwardVelocity`, `Salinity`, `SuspendedParticulateMatter`, `SSTmax`, `SSTmean`, `SSTmin`, `SSTsd`, `seafloorTemp`, `Chla`, `TravelTime`, `ReefMinDist.m`). |
| `DBtaxonomy.xlsx` | Taxonomy reference table. Required columns: `tax_name`, `tax_id`, `genus`, `family`, `order`, `class`. |

## Expected data format

- **Abundance CSV**: UTF-8, comma-separated. Each row corresponds to one MOTU observed at one station in one replicate. `mean_pcr_count_reads` is the normalised read count per PCR replicate.
- **Environmental CSV**: UTF-8, comma-separated. One row per station. All numeric covariates should be in SI units.
- **Taxonomy XLSX**: Standard Excel format. Species names in `tax_name` must match `new_scientific_name_ncbi` in the abundance CSV.

## Generated file

After running `Rscript analyses/scripts/make_db_object.R`, the following file is created here:

| File | Description |
|------|-------------|
| `seamount_integrated_dataset.rda` | Integrated R dataset object (`sm`) containing the wide-format abundance matrix, sample metadata, and taxonomy table. Used by all downstream scripts. |

This generated file is excluded from version control (see `.gitignore`).
