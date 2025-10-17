# Analysis helpers --------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(terra)
  library(sf)
})

GLOBAL_OCEAN_AREA_KM2 <- 361000000

#' Run the Python rasterisation pipeline and return output paths.
run_mpa_raster_pipeline <- function(data_dir,
                                    output_dir,
                                    python = "python3",
                                    resolution = 1000,
                                    start_year = 2000,
                                    end_year = 2025) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  data_dir_abs <- normalizePath(data_dir, mustWork = TRUE)
  output_dir_abs <- normalizePath(output_dir, mustWork = TRUE)

  # MPA rasters come from the curated legacy bundle to preserve historical trends.
  legacy_dir <- file.path("data", "raw", "protected_seas_raster_legacy")
  legacy_raster <- file.path(legacy_dir, "lfp_max_1km.tif")
  legacy_counts <- file.path(legacy_dir, "lfp_yearly_counts.csv")

  if (!all(file.exists(c(legacy_raster, legacy_counts)))) {
    stop("Legacy MPA rasters not found in data/raw/protected_seas_raster_legacy/.")
  }

  mpas_dir <- file.path(output_dir_abs, "mpa")
  dir.create(mpas_dir, recursive = TRUE, showWarnings = FALSE)
  mpas_raster_dest <- file.path(mpas_dir, "mpa_lfp_max_1km.tif")
  mpas_counts_dest <- file.path(mpas_dir, "mpa_lfp_yearly_counts.csv")

  file.copy(legacy_raster, mpas_raster_dest, overwrite = TRUE)
  file.copy(legacy_counts, mpas_counts_dest, overwrite = TRUE)

  # Non-MPA measures are generated dynamically to capture alternative tools.
  script_path <- file.path("scripts_py", "01_mpa_raster_analysis.py")
  if (!file.exists(script_path)) {
    stop("Raster script not found at ", script_path)
  }
  script_abs <- normalizePath(script_path, mustWork = TRUE)

  other_dir <- file.path(output_dir_abs, "other")
  dir.create(other_dir, recursive = TRUE, showWarnings = FALSE)

  other_args <- c(
    script_abs,
    "--data-dir", data_dir_abs,
    "--output-dir", other_dir,
    "--resolution", as.character(resolution),
    "--start-year", as.character(start_year),
    "--end-year", as.character(end_year),
    "--category-mode", "other",
    "--output-prefix", "other",
    "--exclude-categories", "Fisheries Management Area"
  )

  status <- system2(python, args = other_args)
  if (!identical(status, 0L)) {
    stop("Raster pipeline failed for prefix other with exit status ", status)
  }

  list(
    mpas = list(
      raster = mpas_raster_dest,
      yearly_counts = mpas_counts_dest
    ),
    other = list(
      raster = file.path(other_dir, "other_lfp_max_1km.tif"),
      yearly_counts = file.path(other_dir, "other_lfp_yearly_counts.csv")
    )
  )
}

#' Run the habitat raster pipeline (requires UNEP-WCMC shapefiles).
run_habitat_raster_pipeline <- function(python = "python3") {
  script_path <- file.path("scripts_py", "02_habitat_raster_analysis.py")
  if (!file.exists(script_path)) {
    stop("Habitat raster script not found at ", script_path)
  }

  raw_dir <- file.path("data", "raw", "habitats_original")
  shapefiles <- list.files(raw_dir, pattern = "\\.shp$", full.names = TRUE)
  summary_path <- file.path("data", "processed", "habitat_raster_summary.csv")

  if (length(shapefiles) == 0) {
    warning(
      "No shapefiles found in data/raw/habitats_original. " ,
      "Skipping habitat rasterisation and writing an empty summary."
    )
    readr::write_csv(tibble::tibble(habitat = character(), cells = integer(), area_km2 = numeric()), summary_path)
    return(list(summary = summary_path))
  }

  status <- system2(python, args = script_path)
  if (!identical(status, 0L)) {
    stop("Habitat raster pipeline failed with exit status ", status)
  }

  list(summary = summary_path)
}

#' Convert yearly LFP counts into tidy cumulative statistics for multiple sources.
tidy_yearly_stats <- function(counts, ocean_area = GLOBAL_OCEAN_AREA_KM2) {
  if (!is.list(counts)) {
    stop("`counts` must be a list object returned by run_mpa_raster_pipeline().")
  }

  resolve_counts_path <- function(entry) {
    if (is.list(entry) && !is.null(entry$yearly_counts)) {
      return(entry$yearly_counts)
    }
    if (is.character(entry)) {
      return(entry)
    }
    stop("Unable to resolve counts path from entry.")
  }

  source_labels <- c(mpas = "MPAs", other = "Other measures")

  stats <- purrr::imap_dfr(counts, function(entry, key) {
    counts_path <- resolve_counts_path(entry)
    if (!file.exists(counts_path)) {
      stop("Counts file not found: ", counts_path)
    }

    source_label <- if (key %in% names(source_labels)) {
      source_labels[[key]]
    } else {
      stringr::str_to_title(key)
    }

    readr::read_csv(counts_path, show_col_types = FALSE) %>%
      mutate(
        year = as.integer(year),
        minimal_cells = lfp_1 + lfp_2 + lfp_3,
        high_cells = lfp_4 + lfp_5
      ) %>%
      select(year, minimal_cells, high_cells) %>%
      pivot_longer(
        cols = c(minimal_cells, high_cells),
        names_to = "portfolio",
        values_to = "cells"
      ) %>%
      mutate(
        source = source_label,
        protection_level = dplyr::case_when(
          portfolio == "minimal_cells" ~ "Minimally or Moderately Protected",
          portfolio == "high_cells" ~ "Highly/Fully Protected",
          TRUE ~ portfolio
        ),
        series = sprintf("%s – %s", source, protection_level),
        cumulative_area_km2 = cells,
        cumulative_percent = 100 * cumulative_area_km2 / ocean_area
      ) %>%
      group_by(series) %>%
      arrange(year) %>%
      mutate(
        annual_delta_km2 = cumulative_area_km2 - lag(cumulative_area_km2),
        annual_delta_percent = cumulative_percent - lag(cumulative_percent)
      ) %>%
      ungroup()
  })

  series_levels <- c(
    "MPAs – Minimally or Moderately Protected",
    "MPAs – Highly/Fully Protected",
    "Other measures – Minimally or Moderately Protected",
    "Other measures – Highly/Fully Protected"
  )

  stats %>%
    mutate(
      series = factor(series, levels = series_levels),
      protection_level = factor(
        protection_level,
        levels = c("Minimally or Moderately Protected", "Highly/Fully Protected")
      )
    )
}

#' Summarise protected area by enabling band for a single source.
summarise_protection_by_enabling <- function(country_tbl, enabling_tbl, source_label) {
  enabling_levels <- c("<50", "50-69", ">=70", "Missing")

  country_tbl %>%
    left_join(enabling_tbl, by = c("sovereign", "territory")) %>%
    mutate(
      enabling_conditions = dplyr::if_else(
        territory == "Antarctica",
        NA_real_,
        enabling_conditions
      ),
      enabling_band = case_when(
        is.na(enabling_conditions) ~ "Missing",
        enabling_conditions < 50 ~ "<50",
        enabling_conditions < 70 ~ "50-69",
        TRUE ~ ">=70"
      ),
      source = source_label
    ) %>%
    group_by(source, enabling_band) %>%
    summarise(
      `Highly/Fully Protected` = sum(high_area_km2, na.rm = TRUE),
      `Minimally or Moderately Protected` = sum(minimal_area_km2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(enabling_band != "Missing") %>%
    pivot_longer(
      cols = c(`Highly/Fully Protected`, `Minimally or Moderately Protected`),
      names_to = "protection_level",
      values_to = "area_km2"
    ) %>%
    mutate(
      enabling_band = factor(enabling_band, levels = enabling_levels),
      source = factor(source, levels = c("MPAs", "Other measures")),
      area_million_km2 = area_km2 / 1e6
    )
}

#' Combine protection summaries for multiple sources.
combine_protection_sources <- function(mpa_tbl, other_tbl, enabling_tbl) {
  bind_rows(
    summarise_protection_by_enabling(mpa_tbl, enabling_tbl, "MPAs"),
    summarise_protection_by_enabling(other_tbl, enabling_tbl, "Other measures")
  )
}

#' Summarise habitat protection by enabling band.
summarise_habitat_by_enabling <- function(habitat_tbl, enabling_tbl) {
  enabling_levels <- c("<50", "50-69", ">=70", "Missing")
  base_habitats <- c("Cold-water corals", "Warm-water corals", "Mangroves", "Saltmarshes", "Seagrasses")

  habitat_tbl %>%
    left_join(enabling_tbl, by = c("sovereign", "territory")) %>%
    filter(habitat %in% base_habitats) %>%
    mutate(
      enabling_conditions = dplyr::if_else(
        territory == "Antarctica",
        NA_real_,
        enabling_conditions
      ),
      enabling_band = case_when(
        is.na(enabling_conditions) ~ "Missing",
        enabling_conditions < 50 ~ "<50",
        enabling_conditions < 70 ~ "50-69",
        TRUE ~ ">=70"
      ),
      minimal_area_km2 = pmax(mpa - heavily_restricted, 0),
      high_area_km2 = heavily_restricted
    ) %>%
    group_by(habitat, enabling_band) %>%
    summarise(
      `Highly/Fully Protected` = sum(high_area_km2, na.rm = TRUE),
      `Minimally or Moderately Protected` = sum(minimal_area_km2, na.rm = TRUE),
      total_habitat_km2 = sum(habitat_area, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(enabling_band != "Missing") %>%
    pivot_longer(
      cols = c(`Highly/Fully Protected`, `Minimally or Moderately Protected`),
      names_to = "protection_level",
      values_to = "area_km2"
    ) %>%
    mutate(
      enabling_band = factor(enabling_band, levels = enabling_levels),
      area_thousand_km2 = area_km2 / 1e3
    )
}

#' Summarise total habitat extent by enabling band (no protection split).
summarise_habitat_extent_by_enabling <- function(habitat_tbl, enabling_tbl) {
  base_habitats <- c("Cold-water corals", "Warm-water corals", "Mangroves", "Saltmarshes", "Seagrasses")

  habitat_tbl %>%
    left_join(enabling_tbl, by = c("sovereign", "territory")) %>%
    filter(habitat %in% base_habitats) %>%
    mutate(
      enabling_conditions = dplyr::if_else(
        territory == "Antarctica",
        NA_real_,
        enabling_conditions
      ),
      enabling_band = case_when(
        is.na(enabling_conditions) ~ "Missing",
        enabling_conditions < 50 ~ "<50",
        enabling_conditions < 70 ~ "50-69",
        TRUE ~ ">=70"
      )
    ) %>%
    filter(enabling_band != "Missing") %>%
    group_by(habitat, enabling_band) %>%
    summarise(total_habitat_km2 = sum(habitat_area, na.rm = TRUE), .groups = "drop") %>%
    group_by(habitat) %>%
    mutate(
      total_global_km2 = sum(total_habitat_km2, na.rm = TRUE),
      share_percent = if_else(total_global_km2 > 0, 100 * total_habitat_km2 / total_global_km2, NA_real_)
    ) %>%
    ungroup()
}

#' Identify high-enabling jurisdictions with low fully protected coverage.
identify_high_enabling_gaps <- function(opportunities_tbl, min_high_share = 30) {
  opportunities_tbl %>%
    filter(!is.na(enabling_conditions), enabling_conditions >= 70) %>%
    group_by(sovereign, territory, habitat) %>%
    summarise(
      enabling_mean = mean(enabling_conditions, na.rm = TRUE),
      habitat_area_km2 = sum(habitat_area, na.rm = TRUE),
      minimal_area_km2 = sum(minimal_area_km2, na.rm = TRUE),
      high_area_km2 = sum(high_area_km2, na.rm = TRUE),
      global_share_percent = sum(global_share, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      high_share_percent = if_else(habitat_area_km2 > 0, 100 * high_area_km2 / habitat_area_km2, NA_real_),
      minimal_share_percent = if_else(habitat_area_km2 > 0, 100 * minimal_area_km2 / habitat_area_km2, NA_real_),
      upgrade_gap_km2 = pmax(minimal_area_km2 - high_area_km2, 0)
    ) %>%
    filter(
      !is.na(high_share_percent),
      high_share_percent < min_high_share,
      upgrade_gap_km2 > 0
    ) %>%
    arrange(desc(upgrade_gap_km2))
}

#' Summarise country-level upgrade gap by enabling capacity.
summarise_country_upgrade_gap <- function(opportunities_tbl, enabling_tbl) {
  opportunities_tbl %>%
    group_by(sovereign, territory) %>%
    summarise(
      minimal_area_km2 = sum(minimal_area_km2, na.rm = TRUE),
      high_area_km2 = sum(high_area_km2, na.rm = TRUE),
      habitat_area_km2 = sum(habitat_area, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(enabling_tbl, by = c("sovereign", "territory")) %>%
    mutate(
      enabling_conditions = dplyr::if_else(
        territory == "Antarctica",
        NA_real_,
        enabling_conditions
      ),
      enabling_band = case_when(
        is.na(enabling_conditions) ~ "Missing",
        enabling_conditions < 50 ~ "<50",
        enabling_conditions < 70 ~ "50-69",
        TRUE ~ ">=70"
      ),
      total_protected_km2 = minimal_area_km2 + high_area_km2,
      high_share_protected = if_else(total_protected_km2 > 0, 100 * high_area_km2 / total_protected_km2, NA_real_),
      high_share_habitat = if_else(habitat_area_km2 > 0, 100 * high_area_km2 / habitat_area_km2, NA_real_),
      minimal_share_habitat = if_else(habitat_area_km2 > 0, 100 * minimal_area_km2 / habitat_area_km2, NA_real_),
      gap_km2 = pmax(minimal_area_km2 - high_area_km2, 0)
    ) %>%
    filter(!is.na(enabling_conditions))
}

#' Compute country-level shares of high and minimal protection using the LFP raster.
compute_country_protection <- function(raster_path, eez_sf, aggregation_factor = 5) {
  lfp_raster <- terra::rast(raster_path)
  base_cell_area_km2 <- (terra::res(lfp_raster)[1] * terra::res(lfp_raster)[2]) / 1e6

  high_raster <- terra::ifel(lfp_raster >= 4, 1, 0)
  minimal_raster <- terra::ifel(lfp_raster > 0 & lfp_raster < 4, 1, 0)

  if (aggregation_factor > 1) {
    high_raster <- terra::aggregate(high_raster, fact = aggregation_factor, fun = sum, na.rm = TRUE)
    minimal_raster <- terra::aggregate(minimal_raster, fact = aggregation_factor, fun = sum, na.rm = TRUE)
  }

  eez_proj <- sf::st_transform(eez_sf, terra::crs(high_raster, proj = TRUE))
  eez_vect <- terra::vect(eez_proj)

  cell_area_km2 <- base_cell_area_km2

  high_sum <- terra::extract(high_raster, eez_vect, fun = sum, na.rm = TRUE) %>%
    rename(high_cells = 2)
  minimal_sum <- terra::extract(minimal_raster, eez_vect, fun = sum, na.rm = TRUE) %>%
    rename(minimal_cells = 2)

  eez_area_km2 <- as.numeric(sf::st_area(eez_proj)) / 1e6

  tibble::tibble(
    row_id = eez_proj$row_id,
    sovereign = eez_proj$SOVEREIGN1,
    territory = eez_proj$TERRITORY1,
    eez_area_km2 = eez_area_km2,
    high_cells = replace_na(high_sum$high_cells, 0),
    minimal_cells = replace_na(minimal_sum$minimal_cells, 0)
  ) %>%
    mutate(
      high_area_km2 = high_cells * cell_area_km2,
      minimal_area_km2 = minimal_cells * cell_area_km2,
      high_percent = if_else(eez_area_km2 > 0, 100 * high_area_km2 / eez_area_km2, NA_real_),
      minimal_percent = if_else(eez_area_km2 > 0, 100 * minimal_area_km2 / eez_area_km2, NA_real_)
    )
}

#' Combine habitat coverage with enabling conditions to estimate upgrade gains.
derive_upgrade_opportunities <- function(habitat_tbl, enabling_tbl) {
  habitat_totals <- habitat_tbl %>%
    mutate(habitat_area = if_else(is.na(habitat_area), 0, habitat_area)) %>%
    group_by(habitat) %>%
    summarise(global_habitat = sum(habitat_area, na.rm = TRUE), .groups = "drop")

  habitat_tbl %>%
    mutate(
      minimal_area_km2 = pmax(mpa - heavily_restricted, 0),
      high_area_km2 = heavily_restricted,
      upgrade_gain_percent = if_else(
        habitat_area > 0,
        100 * minimal_area_km2 / habitat_area,
        NA_real_
      ),
      high_protection_percent = if_else(
        habitat_area > 0,
        100 * high_area_km2 / habitat_area,
        NA_real_
      )
    ) %>%
    left_join(habitat_totals, by = "habitat") %>%
    mutate(
      global_share = if_else(
        global_habitat > 0,
        100 * habitat_area / global_habitat,
        NA_real_
      )
    ) %>%
    left_join(enabling_tbl, by = c("sovereign", "territory")) %>%
    mutate(
      enabling_band = case_when(
        is.na(enabling_conditions) ~ "Data missing",
        enabling_conditions >= 70 ~ ">=70",
        enabling_conditions >= 50 ~ "50-69",
        TRUE ~ "<50"
      ),
      enabling_band = factor(
        enabling_band,
        levels = c(">=70", "50-69", "<50", "Data missing")
      ),
      opportunity_score = minimal_area_km2 * (coalesce(enabling_conditions, 0) / 100)
    )
}

#' Aggregate upgrade opportunities at sovereign/territory level.
aggregate_upgrade_by_country <- function(opportunities_tbl) {
  opportunities_tbl %>%
    group_by(sovereign, territory) %>%
    summarise(
      upgrade_minimal_km2 = sum(minimal_area_km2, na.rm = TRUE),
      upgrade_high_km2 = sum(high_area_km2, na.rm = TRUE),
      enabling_mean = mean(enabling_conditions, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      enabling_mean = if_else(is.nan(enabling_mean), NA_real_, enabling_mean),
      upgrade_score = upgrade_minimal_km2 * (coalesce(enabling_mean, 0) / 100)
    )
}
