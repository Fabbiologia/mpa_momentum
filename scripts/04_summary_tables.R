# Table writers -----------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

write_annual_summary <- function(yearly_tbl, path) {
  yearly_tbl %>%
    arrange(series, year) %>%
    mutate(
      cumulative_area_km2 = round(cumulative_area_km2, digits = 1),
      cumulative_percent = round(cumulative_percent, digits = 4)
    ) %>%
    readr::write_csv(path)

  path
}

write_growth_summary <- function(yearly_tbl, path) {
  yearly_tbl %>%
    arrange(series, year) %>%
    mutate(
      annual_delta_km2 = round(annual_delta_km2, 1),
      annual_delta_percent = round(annual_delta_percent, 4)
    ) %>%
    readr::write_csv(path)

  path
}

write_upgrade_table <- function(opportunities_tbl, path) {
  opportunities_tbl %>%
    select(
      sovereign,
      territory,
      habitat,
      habitat_area,
      mpa,
      heavily_restricted,
      minimal_area_km2,
      high_area_km2,
      upgrade_gain_percent,
      opportunity_score,
      high_protection_percent,
      enabling_conditions,
      enabling_band,
      global_share
    ) %>%
    arrange(desc(upgrade_gain_percent)) %>%
    readr::write_csv(path)

  path
}

write_global_habitat_table <- function(global_tbl, path) {
  global_tbl %>%
    mutate(
      percent_all = round(percent_all, 2),
      percent_PA = round(percent_PA, 2)
    ) %>%
    readr::write_csv(path)

  path
}
