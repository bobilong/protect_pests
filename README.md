# Code and data for global analyses of forest pest-record visibility in protected forests

This repository contains the R code and processed data used to produce the analyses and figures for a manuscript "Uneven forest pest observations constrain carbon monitoring in protected forests worldwide".

## Repository contents

```text
protect_pests/
├── README.md
├── fig1.R                 # Main Figure 1 and Supplementary Figures S1-S3
├── fig2.R                 # Main Figure 2 and Supplementary Figures S4-S5
├── fig3.r                 # Main Figure 3 and Supplementary Figure S6
└── data/
    ├── *.tif              # Processed spatial raster layers
    ├── match_table.csv    # Grid-level data before matching
    ├── matchit_data.csv   # Matched grid-level analysis data
    ├── gam_fig2_model.rds # Fitted GAM used for Figure 2
    └── *.csv              # Processed and figure-ready summary data
```

## Figures and analyses

| Manuscript output | Script | Description |
| --- | --- | --- |
| Figure 1 | `fig1.R` | Global distribution of protected forests and pest records; variation in pest-record coverage and intensity across forest-carbon percentiles |
| Figure S1 | `fig1.R` | Background biodiversity-sampling coverage and intensity across carbon percentiles |
| Figure S2 | `fig1.R` | Bivariate spatial patterns of forest carbon and pest-record intensity by forest biome |
| Figure S3 | `fig1.R` | Pest-record coverage across carbon percentiles within forest biomes |
| Figure 2 | `fig2.R` | Matched comparison of pest-record probability inside and outside protected areas, including absolute and relative differences |
| Figure S4 | `fig2.R` | Covariate balance before and after matching |
| Figure S5 | `fig2.R` | Effect estimates by carbon percentile and forest biome, with matched protected-forest area |
| Figure 3 | `fig3.r` | Pest-record visibility by governance type, IUCN management category, and country |
| Figure S6 | `fig3.r` | Forest Landscape Integrity Index (FLII), pest-record probability, and background sampling effort |

## Analytical overview

### Matching protected and unprotected forest grids

The analysis in `fig2.R` defines protected forest grids as cells with at least 50% protected-area coverage and control grids as cells with no protected-area coverage. Forest grids are restricted to tree cover greater than 10%.

Protected and control grids are matched using `MatchIt` with:

- exact matching within ecoregion-by-country strata;
- 1:1 nearest-neighbour matching without replacement;
- a generalised linear model distance and a caliper of 0.2; and
- forest carbon percentile, tree cover, human footprint, climate, elevation, and forest fragmentation as matching covariates.

### Generalised additive model

A binomial generalised additive model is fitted with `mgcv::bam()`. The response indicates whether a grid contains at least one pest record. The model includes protected-area status, status-specific smooths of forest-carbon percentile, background GBIF sampling effort, climate, human footprint, elevation, tree cover, forest fragmentation, and a two-dimensional spatial smooth.

Uncertainty in the protected-versus-control contrasts is estimated from 1,000 simulations of the model coefficients using the fitted coefficient covariance matrix. Predictions are summarised by forest-carbon decile.

## Software requirements

The scripts were developed in an R 4.3 library environment. R 4.3 or later is recommended.

Required R packages can be installed from CRAN with:

```r
install.packages(c(
  "terra", "sf", "tidyterra",
  "data.table", "dplyr", "tidyr", "stringr", "purrr", "forcats",
  "ggplot2", "scales", "patchwork",
  "MatchIt", "mgcv", "MASS",
  "countrycode", "rnaturalearth", "rnaturalearthdata"
))
```

The spatial packages `terra` and `sf` may require system installations of GDAL, GEOS, and PROJ. Exact package versions are not currently locked with `renv`; adding an `renv.lock` file is recommended before creating the archived manuscript release.

## Reproducing the analyses

1. Clone or download the repository.

2. Start R from the repository root, or set the working directory explicitly:

   ```r
   setwd("path/to/protect_pests")
   ```

3. Install the packages listed above.

4. Confirm that the required files are available in `data/`. The repository includes processed model inputs, fitted model output, and figure-ready summary tables. Large grid-level files require substantial memory.

5. Update the machine-specific paths described under [Reproducibility notes](#reproducibility-notes).

6. Open the required script in R or RStudio and run its sections in order. The scripts create named `ggplot2`/`patchwork` objects in the active R session. They do not currently save all figures automatically.

   After the path and object issues below have been resolved, a complete script can be run with, for example:

   ```r
   source("fig2.R", echo = TRUE)
   ```

7. Export the required plot object at the dimensions used for the manuscript. For example:

   ```r
   ggsave(
     filename = "Figure_2.png",
     plot = fig2a + fig2b,
     width = 10,
     height = 6,
     dpi = 300
   )
   ```

## Key processed data files

| File | Description |
| --- | --- |
| `data/match_table.csv` | Grid-level forest, protection, pest-record, sampling-effort, climate, topographic, and human-pressure variables used before matching |
| `data/matchit_data.csv` | Matched protected and control grids, including matching weights and subclasses |
| `data/gam_fig2_model.rds` | Fitted binomial GAM used for the Figure 2 predictions |
| `data/fig2_df.csv` | Protected and control predictions and contrasts by forest-carbon decile |
| `data/fig2_pred_long.csv` | Figure-ready protected and control predictions with uncertainty intervals |
| `data/balance_plot_data.csv` | Standardised mean differences before and after matching |
| `data/biome_carbon_decile_pest_coverage.csv` | Pest-record coverage by forest biome and carbon decile |
| `data/fig3_attribute_pest_visibility_summary.csv` | Protected-area summaries by governance type and IUCN category |
| `data/fig3_country_pest_visibility_summary.csv` | Country-level protected-forest and high-probability-area summaries |
| `data/predicted_pest_probability_by_FLII.csv` | Modelled pest-record probability along the FLII gradient |
| `data/predicted_sampling_effort_by_FLII.csv` | Modelled background sampling effort along the FLII gradient |


## Data availability

This repository contains processed spatial layers, model inputs, fitted model output, and figure-ready summaries used by the manuscript scripts.

## Citation

Citation details for the manuscript and the archived code release will be added upon publication.

## Licence

No software licence is currently included in this repository. A licence should be added before public release; until then, no permission for reuse is granted beyond applicable statutory exceptions.

## Contact

For questions about the analysis or code, please contact the corresponding author listed in the manuscript or open an issue in this repository.
