# Data preparation helpers -------------------------------------------------
#
# These functions load the raw datasets needed for the analysis and return
# clean, ready-to-use objects. All paths are passed in explicitly so the
# pipeline can remain portable and testable.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(sf)
  library(stringr)
  library(tidyr)
})

#' Load the ProtectedSeas Navigator spatial dataset and attributes.
#'
#' @param root Path to the folder that contains both the shapefile directory
#'   and the CSV attributes (e.g., data/raw/protected_seas_navigator_all_sites).
#' @return sf object with attributes merged and geometry column preserved.
load_protectedseas <- function(root) {
  shp_path <- file.path(root, "Navigator_AllSites_032825_shp", "Navigator_AllSites_032825.shp")
  corrected_csv <- file.path(root, "Navigator_AllSites_032825_corrected.csv")
  csv_path <- if (file.exists(corrected_csv)) {
    corrected_csv
  } else {
    file.path(root, "Navigator_AllSites_032825.csv")
  }

  shp <- sf::st_read(shp_path, quiet = TRUE)
  attrs <- readr::read_csv(csv_path, show_col_types = FALSE)

  shp %>%
    rename(site_id = SITE_ID) %>%
    left_join(attrs, by = "site_id") %>%
    filter(!is.na(geometry))
}

#' Resolve duplicated parent-child polygons by keeping the record with the
#' highest protection level.
#'
#' The ProtectedSeas database stores relationships where a parent site can
#' contain multiple sub-sites (children). To avoid double counting we retain
#' either the parent or the child with the strongest protection (LFP) score.
#'
#' @param mpa_sf sf object returned by load_protectedseas().
#' @return sf object with duplicates removed.
clean_protectedseas <- function(mpa_sf) {
  mpa_sf %>%
    mutate(
      child = coalesce(as.logical(child), FALSE),
      parent_id = na_if(parent_id, ""),
      site_id = as.character(site_id),
      lfp = as.numeric(lfp),
      year_est = suppressWarnings(as.integer(year_est))
    ) -> data

  parent_ids <- unique(na.omit(data$parent_id))
  keep_ids <- character(0)

  for (pid in parent_ids) {
    parent_row <- data %>% filter(site_id == pid)
    child_rows <- data %>% filter(parent_id == pid)

    if (nrow(child_rows) == 0) {
      next
    }

    best_child <- child_rows %>%
      slice_max(order_by = lfp, n = 1, with_ties = TRUE) %>%
      slice(1)

    parent_lfp <- if (nrow(parent_row) > 0) parent_row$lfp[1] else NA_real_

    if (!is.na(parent_lfp) && parent_lfp >= best_child$lfp) {
      keep_ids <- c(keep_ids, pid)
    } else {
      keep_ids <- c(keep_ids, best_child$site_id[1])
    }
  }

  data %>%
    mutate(
      keep = case_when(
        site_id %in% keep_ids ~ TRUE,
        !child & is.na(parent_id) ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    filter(keep) %>%
    select(-keep)
}

#' Compute area (km2) for each MPA after cleaning.
#'
#' @param mpa_clean sf object returned by clean_protectedseas().
#' @param area_crs Equal-area CRS for stable area calculations.
#' @return sf object with area_km2 column.
add_area_km2 <- function(mpa_clean, area_crs = "+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs") {
  mpa_clean %>%
    st_make_valid() %>%
    st_transform(crs = area_crs) %>%
    mutate(area_km2 = as.numeric(st_area(geometry)) / 1e6) %>%
    st_transform(crs = 4326)
}

#' Attach simple protection category labels used throughout the analysis.
#'
#' @param mpa_sf sf object with numeric `lfp`.
#' @return sf object with `protection_class` factor column.
add_protection_class <- function(mpa_sf) {
  mpa_sf %>%
    mutate(
      protection_class = case_when(
        is.na(lfp) ~ "Unclassified",
        lfp >= 4 ~ "Highly or Fully Protected",
        lfp == 3 ~ "Moderately Protected",
        lfp > 0 ~ "Minimally Protected",
        TRUE ~ "Unclassified"
      ),
      protection_class = factor(
        protection_class,
        levels = c("Highly or Fully Protected", "Moderately Protected", "Minimally Protected", "Unclassified")
      )
    )
}

#' Load EEZ boundaries and keep primary sovereign/territory names.
#'
#' @param path Path to EEZ shapefile (e.g., data/raw/eez_land/EEZ_Land_v3_202030.shp).
#' @return sf object filtered to unique sovereign territories.
load_eez <- function(path) {
  sf::st_read(path, quiet = TRUE) %>%
    filter(!str_detect(UNION, "Joint regime"), !str_detect(UNION, "Overlapping claim")) %>%
    mutate(row_id = dplyr::row_number())
}

#' Load habitat coverage matrices exported from the original workflow and
#' attach sovereign / territory names from the EEZ polygons.
#'
#' @param csv_path Path to the processed `habitat_area_PS.csv` table (e.g. `data/processed/reference_tables/habitat_area_PS.csv`).
#' @param eez_sf sf object from load_eez().
#' @return Tibble with area metrics by sovereign, territory, and habitat.
load_habitat_matrix <- function(csv_path, eez_sf) {
  raw <- readr::read_csv(csv_path, show_col_types = FALSE) %>%
    rename(row_id = 1)

  if (nrow(raw) != nrow(eez_sf)) {
    stop("EEZ geometry count and habitat matrix rows do not match; join is unsafe.")
  }

  raw %>%
    mutate(row_id = dplyr::row_number()) %>%
    pivot_longer(
      cols = -row_id,
      names_to = "metric",
      values_to = "value"
    ) %>%
    separate(metric, into = c("habitat", "descriptor"), sep = "_(?=[^_]+$)", remove = FALSE) %>%
    mutate(
      descriptor = case_when(
        str_detect(metric, "PS_MPA$") ~ "heavily_restricted",
        str_detect(metric, "All_MPA$") ~ "mpa",
        str_detect(metric, "All_MPA_with_PS") ~ "mpa_with_ps",
        TRUE ~ "habitat_area"
      ),
      habitat = str_remove_all(habitat, "_All|_PS")
    ) %>%
    filter(!str_detect(habitat, "Knoll")) %>%
    select(row_id, habitat, descriptor, value) %>%
    left_join(
      eez_sf %>%
        st_drop_geometry() %>%
        select(row_id, sovereign = SOVEREIGN1, territory = TERRITORY1),
      by = "row_id"
    ) %>%
    group_by(sovereign, territory, habitat, descriptor) %>%
    summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = descriptor, values_from = value) %>%
    select(-dplyr::any_of("mpa_with_ps")) %>%
    mutate(
      habitat = recode(
        habitat,
        "mangroves" = "Mangroves",
        "Seagrasses" = "Seagrasses",
        "Saltmarshes" = "Saltmarshes",
        "CoralReefs" = "Warm-water corals",
        "ColdCorals" = "Cold-water corals",
        .default = habitat
      ),
      habitat_area = replace_na(habitat_area, 0),
      mpa = replace_na(mpa, 0),
      heavily_restricted = replace_na(heavily_restricted, 0)
    )
}

#' Load global habitat summary table.
#'
#' @param csv_path Path to the processed global summary (e.g. `data/processed/reference_tables/percent_protected_world2.csv`).
load_global_habitat_summary <- function(csv_path) {
  readr::read_csv(csv_path, show_col_types = FALSE) %>%
    mutate(
      percent_all = readr::parse_number(percent_all),
      percent_PA = readr::parse_number(percent_PA),
      habitat = recode(
        habitat,
        "ColdCorals" = "Cold-water corals",
        "CoralReefs" = "Warm-water corals",
        "KnollSeamounts" = "Seamounts",
        .default = habitat
      )
    )
}

#' Prepare MPA geometries and attributes for raster workflows.
prepare_mpa_for_raster <- function(mpa_sf, target_crs = "EPSG:6933") {
  mpa_sf %>%
    mutate(
      lfp = suppressWarnings(as.numeric(lfp)),
      lfp = replace_na(lfp, 0),
      year_est = suppressWarnings(as.integer(year_est)),
      year_est = replace_na(year_est, 1999L)
    ) %>%
    st_make_valid() %>%
    st_transform(target_crs)
}

#' Create template metadata for rasterisation at a specified resolution (m).
create_mpa_template_definition <- function(mpa_sf_transformed, resolution = 1000) {
  ext <- st_bbox(mpa_sf_transformed)
  list(
    extent = c(
      xmin = ext[["xmin"]],
      xmax = ext[["xmax"]],
      ymin = ext[["ymin"]],
      ymax = ext[["ymax"]]
    ),
    resolution = resolution,
    crs = st_crs(mpa_sf_transformed)$wkt
  )
}

#' Read enabling condition scores (already spatial) and drop geometry.
#'
#' @param rds_path Path to the processed enabling-condition table (e.g. `data/processed/reference_tables/blue_economy_intersected_data.RDS`).
load_enabling_conditions <- function(rds_path) {
  readRDS(rds_path) %>%
    sf::st_drop_geometry() %>%
    rename_with(~ tolower(.x)) %>%
    mutate(across(where(is.character), stringr::str_trim)) %>%
    group_by(sovereign, territory) %>%
    summarise(
      enabling_conditions = suppressWarnings(mean(enabling_conditions, na.rm = TRUE)),
      blue_economy = suppressWarnings(mean(blue_economy, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      enabling_conditions = ifelse(is.nan(enabling_conditions), NA_real_, enabling_conditions),
      blue_economy = ifelse(is.nan(blue_economy), NA_real_, blue_economy)
    )
}
