# Targets pipeline --------------------------------------------------------

library(targets)
library(tarchetypes)

tar_option_set(
  packages = c(
    "dplyr",
    "sf",
    "readr",
    "stringr",
    "tidyr",
    "purrr",
    "forcats",
    "ggplot2",
    "patchwork",
    "ggrepel",
    "scales",
    "terra",
    "htmltools",
    "rmarkdown",
    "rnaturalearth",
    "rnaturalearthdata",
    "testthat"
  ),
  format = "rds"
)

invisible(lapply(
  c(
    "scripts/01_data_preparation.R",
    "scripts/02_analysis.R",
    "scripts/03_visualization.R",
    "scripts/04_summary_tables.R"
  ),
  source
))

list(
  tar_target(path_protectedseas, "data/raw/protected_seas_navigator_all_sites"),
  tar_target(eez_sf, load_eez("data/raw/eez_land/EEZ_Land_v3_202030.shp")),
  tar_target(
    raster_outputs,
    run_mpa_raster_pipeline(
      data_dir = path_protectedseas,
      output_dir = "data/processed/mpa_raster",
      end_year = 2025
    )
  ),
  tar_target(yearly_stats, tidy_yearly_stats(raster_outputs)),
  tar_target(country_protection, compute_country_protection(raster_outputs$mpas$raster, eez_sf, aggregation_factor = 1)),
  tar_target(country_protection_other, compute_country_protection(raster_outputs$other$raster, eez_sf, aggregation_factor = 1)),

  tar_target(habitat_matrix, load_habitat_matrix("data/processed/reference_tables/habitat_area_PS.csv", eez_sf)),
  tar_target(enabling_conditions, load_enabling_conditions("data/processed/reference_tables/blue_economy_intersected_data.RDS")),
  tar_target(upgrade_opportunities, derive_upgrade_opportunities(habitat_matrix, enabling_conditions)),
  tar_target(upgrade_by_country, aggregate_upgrade_by_country(upgrade_opportunities)),
  tar_target(protection_enabling_summary, combine_protection_sources(country_protection, country_protection_other, enabling_conditions)),
  tar_target(habitat_enabling_summary, summarise_habitat_by_enabling(habitat_matrix, enabling_conditions)),
  tar_target(habitat_extent_enabling, summarise_habitat_extent_by_enabling(habitat_matrix, enabling_conditions)),
  tar_target(country_upgrade_gap, summarise_country_upgrade_gap(upgrade_opportunities, enabling_conditions)),
  tar_target(high_enabling_gaps, identify_high_enabling_gaps(upgrade_opportunities)),
  tar_target(global_habitat_summary, load_global_habitat_summary("data/processed/reference_tables/percent_protected_world2.csv")),
  tar_target(habitat_raster_outputs, run_habitat_raster_pipeline()),
  tar_target(habitat_raster_summary, readr::read_csv(habitat_raster_outputs$summary, show_col_types = FALSE)),

  tar_target(fig_figure1, plot_figure1(yearly_stats, global_habitat_summary, "outputs/figures/figure_1.png"), format = "file"),
  tar_target(fig_figure2, plot_figure2(country_upgrade_gap, habitat_extent_enabling, habitat_enabling_summary, high_enabling_gaps, eez_sf, "outputs/figures/figure_2.png"), format = "file"),

  tar_target(table_annual_protection, write_annual_summary(yearly_stats, "outputs/tables/annual_protection_by_group.csv"), format = "file"),
  tar_target(table_growth, write_growth_summary(yearly_stats, "outputs/tables/annual_growth_by_group.csv"), format = "file"),
  tar_target(table_upgrade_opportunities, write_upgrade_table(upgrade_opportunities, "outputs/tables/upgrade_opportunities.csv"), format = "file"),
  tar_target(table_global_habitats, write_global_habitat_table(global_habitat_summary, "outputs/tables/global_habitat_summary.csv"), format = "file"),
  
  tarchetypes::tar_render(
    manuscript_report,
    "manuscript/ISCIENCE-D-25-20442_main_manuscript_reviewed.Rmd",
    output_file = "ISCIENCE-D-25-20442_main_manuscript_reviewed.docx",
    output_dir = "manuscript",
    quiet = TRUE
  ),
  tarchetypes::tar_render(
    response_letter,
    "manuscript/ISCIENCE-D-25-20442_response_to_reviewers.Rmd",
    output_file = "ISCIENCE-D-25-20442_response_to_reviewers.docx",
    output_dir = "manuscript",
    quiet = TRUE
  )
)
