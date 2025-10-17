# Why Ocean Protection Is Stalling—and How to Regain Momentum

## Directory layout

```
habitat_protection_protectedseas/
├── _targets.R                     # pipeline definition
├── data/                          # raw inputs + processed rasters
├── scripts/                       # modular R helpers (data, analysis, viz, tables)
├── scripts_py/                    # Python rasterisation step
├── outputs/                       # auto-generated figures and tables
└── manuscript/                    # manuscript + response letter sources
```

## Data inputs

All raw files reside under `data/raw/`:

- `protected_seas_navigator_all_sites/` – ProtectedSeas Navigator (03‑28‑25 release, available on request from [ProtectedSeas](https://www.protectedseas.net/navigator)).
- `eez_land/` – EEZ polygon boundaries (filtered to single-sovereign claims).
- `habitats_original/` – UNEP-WCMC mapped habitat layers (download from the [UNEP-WCMC Ocean Data Viewer](https://www.unep-wcmc.org/resources-and-data); see FAIR section).
- `protected_seas_raster_legacy/` – Legacy 1 km ProtectedSeas raster and yearly counts reused from the original workflow.
- `World_EEZ_v12_20231025_gpkg/` – Reference copy of the full VLIZ EEZ release.
- `Enabling_conditions_equitable_sustainable_Blue_Economy_-_Cisneros-Montemayor_et_al_2021.csv` – Original enabling-condition indicators from Cisneros-Montemayor et al. (2021).

Processed reference tables bundled with this repository live in `data/processed/reference_tables/`:

- `habitat_area.csv` – Country-level habitat area summaries exported from the legacy workflow.
- `habitat_area_PS.csv` – Legacy rasterised habitat coverage intersected with MPAs and heavily restricted MPAs.
- `percent_protected_world2.csv` – Global habitat percentage statistics.
- `blue_economy_intersected_data.RDS` – Digitised enabling-condition scores joined to EEZ territories.

## Reproducing the analysis

1. **Install required R packages**

   ```r
   install.packages(c(
     "targets", "tarchetypes", "dplyr", "sf", "readr", "stringr",
     "tidyr", "purrr", "ggplot2", "patchwork", "ggrepel", "scales",
     "terra", "htmltools", "rmarkdown", "rnaturalearth", "rnaturalearthdata"
   ))
   ```

2. **Run the pipeline**

   ```r
   library(targets)
   tar_make()
   ```

3. **Inspect outputs**

   - Figures: `outputs/figures/figure_1.png`, `figure_2.png`, `figure_2_raster.png`
   - Tables: `outputs/tables/*.csv`
   - Manuscript: `manuscript/manuscript.html`
   - Response letter: `manuscript/response_to_reviewers.docx`

4. **Query intermediate targets**

   ```r
   tar_read(yearly_stats)
   tar_read(country_protection)
   tar_read(habitat_raster_summary)
   ```

## Workflow highlights

- **Raster-first pipeline**: the legacy 1 km ProtectedSeas raster (max LFP per cell) is ingested and reused to avoid polygon double-counting while keeping results consistent with the original analysis.
- **Habitat rasterisation**: `scripts_py/02_habitat_raster_analysis.py` reproduces the 1 km UNEP-WCMC habitat rasters (mangroves, seagrasses, saltmarshes, coral reefs, cold-water corals, and seamounts/knolls) using mapped—not modelled—datasets following Kumagai et al. (2021).
- **Targets orchestration**: `scripts/02_analysis.R` ingests the raster outputs to derive cumulative trajectories, country-level statistics, and an improved protection opportunity index linked to enabling conditions.
- **Cartography + figures**: `scripts/03_visualization.R` produces a trends figure, global high-protection map, and upgraded opportunity scatterplot for the manuscript.
- **Integrated writing**: `tarchetypes::tar_render()` keeps the manuscript and reviewer response document in sync with analytical outputs.

## FAIR data & provenance

- **Findable**: Raw UNEP-WCMC habitat shapefiles can be obtained via the [UNEP-WCMC Ocean Data Viewer](https://www.unep-wcmc.org/resources-and-data) (e.g. Global Mangrove Watch, Global Seagrass, Global Saltmarsh, Coral Reef, and Cold-water coral layers). Place the downloaded `.shp` files in `data/raw/habitats_original/`. ProtectedSeas Navigator data are accessible upon request (license-free for research) via the [ProtectedSeas portal](https://www.protectedseas.net/navigator). The enabling-condition scores derive from Cisneros-Montemayor et al. (2021) and can be retrieved from the Dryad repository referenced in the manuscript.
- **Accessible**: Processed rasters are written to `data/processed/habitat_rasters/` with metadata captured in `data/processed/habitat_raster_summary.csv`. ProtectedSeas Navigator source files reside under `data/raw/protected_seas_navigator_all_sites/` (not redistributed here because of upstream licensing); the precomputed legacy raster and yearly counts required by the workflow live in `data/raw/protected_seas_raster_legacy/`; processed auxiliary tables (e.g. enabling-condition scores and legacy habitat matrices) are stored in `data/processed/reference_tables/`.
- **Interoperable**: All spatial outputs use the equal-area EPSG:6933 CRS and 1 km grid resolution; tabular outputs are CSV with UTF-8 encoding.
- **Reusable**: Scripts are version-controlled, documented, and rely only on open-source tooling (R, Python, GDAL stack). The README records source citations and download locations. All processed outputs (`data/processed/`, `outputs/`) are distributed to allow full reproduction of manuscript tables and figures without bundling proprietary raw datasets.

## Processed outputs packaged with this repository

- `data/processed/mpa_raster/lfp_max_1km.tif` and `data/processed/mpa_raster/lfp_yearly_counts.csv` – global 1 km LFP raster and yearly protection trajectories derived from the ProtectedSeas Navigator dataset.
- `data/processed/habitat_rasters/` – 1 km UNEP-WCMC habitat rasters (mangroves, seagrasses, saltmarshes, warm-water coral reefs, cold-water corals) plus `data/processed/habitat_raster_summary.csv` capturing pixel counts and habitat area (km²).
- `data/processed/reference_tables/` – Legacy processed tables (habitat coverage matrices, enabling-condition scores, global habitat percentages) reused by the current pipeline.
- `outputs/tables/` – CSV tables used in the manuscript (annual protection, annual growth, upgrade opportunities, global habitat summary).
- `outputs/figures/figure_1.png`, `figure_2.png`, `figure_2_raster.png` – publication-ready figures referenced in the manuscript and supplementary materials.
- `manuscript/manuscript.html` and `manuscript/response_to_reviewers.docx` – knitted outputs for the main text and response letter.
