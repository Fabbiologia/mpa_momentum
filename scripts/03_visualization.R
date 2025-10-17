# Visualisation helpers ---------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
  library(tidyr)
  library(sf)
  library(rnaturalearth)
  library(forcats)
  library(stringr)
})

plot_figure1 <- function(yearly_stats, global_habitat_summary, output_path) {
  palette <- c(
    "Highly/Fully Protected" = "#1f78b4",
    "Minimally or Moderately Protected" = "#33a02c"
  )

  mpas_stats <- yearly_stats %>%
    filter(source == "MPAs") %>%
    mutate(
      year = as.integer(year),
      protection_group = forcats::fct_relevel(
        protection_level,
        "Minimally or Moderately Protected",
        "Highly/Fully Protected"
      )
    )

  y_max <- max(10, max(mpas_stats$cumulative_percent, na.rm = TRUE) * 1.05)

  line_plot <- mpas_stats %>%
    ggplot(aes(x = year, y = cumulative_percent, color = protection_group)) +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = palette, name = "") +
    scale_x_continuous(breaks = seq(min(mpas_stats$year), max(mpas_stats$year), by = 2)) +
    scale_y_continuous(labels = label_percent(scale = 1), limits = c(0, y_max)) +
    labs(
      x = "Year",
      y = "Cumulative ocean coverage (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")

  growth_data <- mpas_stats %>%
    mutate(
      annual_delta_percent = replace_na(annual_delta_percent, 0),
      protection_group = forcats::fct_rev(protection_group)
    )

 growth_plot <- growth_data %>%
    filter(year >= min(year) + 1) %>%
    ggplot(aes(x = year, y = annual_delta_percent, fill = protection_group)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = palette, name = "") +
    scale_x_continuous(breaks = seq(min(mpas_stats$year), max(mpas_stats$year), by = 2)) +
    scale_y_continuous(labels = label_number(accuracy = 0.01)) +
    labs(
      x = "Year",
      y = "Annual change (percentage points)"
    ) +
    geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")

  habitat_long <- global_habitat_summary %>%
    mutate(
      habitat = factor(habitat, levels = rev(global_habitat_summary$habitat)),
      `Highly protected` = percent_PA,
      `Minimal protection` = pmax(percent_all - percent_PA, 0),
      `Unprotected` = pmax(100 - percent_all, 0)
    ) %>%
    select(habitat, `Highly protected`, `Minimal protection`, `Unprotected`) %>%
    pivot_longer(-habitat, names_to = "protection", values_to = "percent")

  bar_plot <- habitat_long %>%
    ggplot(aes(x = percent, y = habitat, fill = protection)) +
    geom_col(colour = "white", width = 0.75) +
    scale_fill_manual(
      values = c(
        "Highly protected" = "#1f78b4",
        "Minimal protection" = "#33a02c",
        "Unprotected" = "#bdbdbd"
      ),
      name = ""
    ) +
    scale_x_continuous(labels = label_percent(scale = 1), expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Share of habitat protected (%)", y = NULL, fill = "") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")

  time_data <- tibble(time = seq(0, 20, length.out = 200))

  social_costs <- time_data %>%
    mutate(
      fp_cost = 100 / (time + 1),
      mu_cost = 20 + 2 * time + 0.1 * time^2,
      upgrade_cost = if_else(
        time < 10,
        20 + 2 * time + 0.1 * time^2,
        80 * exp(-0.098 * (time - 10))
      )
    )

  concept_long <- social_costs %>%
    select(time, fp_cost, mu_cost, upgrade_cost) %>%
    pivot_longer(-time, names_to = "scenario", values_to = "cost") %>%
    mutate(
      scenario = factor(
        scenario,
        levels = c("mu_cost", "upgrade_cost", "fp_cost"),
        labels = c("Multi-use baseline", "Upgrade (from year 10)", "Fully protected")
      )
    )

  concept_plot <- ggplot(concept_long, aes(time, cost, color = scenario, linetype = scenario)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c("#33a02c", "#1f78b4", "#1f78b4"), name = "") +
    scale_linetype_manual(values = c("solid", "solid", "dashed"), name = "") +
    annotate("segment", x = 10, xend = 10, y = 80, yend = 40, arrow = arrow(length = unit(0.08, "inches"))) +
    annotate("text", x = 10.2, y = 82, label = "Upgrade decision", hjust = 0, size = 3) +
    labs(x = "Time (years)", y = "Social cost (arbitrary units)") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  output_path_abs <- file.path(getwd(), output_path)

  combined <- (line_plot | growth_plot) / (bar_plot | concept_plot) +
    plot_annotation(
      tag_levels = "A"
        )

  ggplot2::ggsave(filename = output_path_abs, plot = combined, width = 14, height = 10, dpi = 300)
  output_path_abs
}

plot_figure2 <- function(country_upgrade_gap,
                         habitat_extent_summary,
                         habitat_enabling_summary,
                         high_enabling_gaps,
                         eez_sf,
                         output_path) {
  robinson_proj <- "+proj=robin +datum=WGS84"
  global_habitat_total <- sum(country_upgrade_gap$habitat_area_km2, na.rm = TRUE)
  total_is_valid <- is.finite(global_habitat_total) && global_habitat_total > 0

  map_data <- eez_sf %>%
    mutate(
      sovereign = SOVEREIGN1,
      territory = TERRITORY1
    ) %>%
    select(sovereign, territory, geometry) %>%
    left_join(country_upgrade_gap, by = c("sovereign", "territory")) %>%
    filter(!is.na(enabling_conditions), enabling_conditions >= 70) %>%
    mutate(
      gap_share_percent = if (total_is_valid) {
        100 * gap_km2 / global_habitat_total
      } else {
        rep(NA_real_, length(gap_km2))
      },
      high_share_habitat = dplyr::coalesce(high_share_habitat, 0)
    ) %>%
    st_transform(robinson_proj)

  world_outline <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
    st_transform(robinson_proj)
  land_union <- sf::st_make_valid(st_union(world_outline))

  map_data <- map_data %>%
    mutate(
      geometry = sf::st_make_valid(geometry),
      geometry = sf::st_difference(geometry, land_union),
      geometry = sf::st_collection_extract(geometry, "POLYGON")
    ) %>%
    filter(!sf::st_is_empty(geometry))

  compute_percent_breaks <- function(values) {
    positives <- values[is.finite(values) & values > 0]
    if (length(positives) == 0) {
      return(NULL)
    }
    candidate <- c(0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10, 20, 30, 40, 50)
    breaks <- candidate[candidate <= max(positives)]
    if (length(breaks) < 3) {
      extra <- scales::extended_breaks()(positives)
      breaks <- sort(unique(c(breaks, extra)))
    }
    breaks
  }

  legend_breaks <- compute_percent_breaks(map_data$gap_share_percent)

  map_plot <- ggplot() +
    geom_sf(
      data = map_data,
      aes(fill = gap_share_percent),
      colour = "white",
      size = 0.1
    ) +
    geom_sf(data = world_outline, fill = "grey93", colour = "white", size = 0.1) +
    scale_fill_viridis_c(
      option = "plasma",
      trans = scales::pseudo_log_trans(sigma = 0.01),
      breaks = legend_breaks,
      labels = label_number(accuracy = 0.01, suffix = "%"),
      name = "Upgrade gap (share of global habitat)",
      na.value = "grey90",
      guide = guide_colorbar( 
        direction = "horizontal",
        title.position = "top",
        barheight = grid::unit(0.4, "cm"),
        barwidth = grid::unit(6, "cm"),
        label.position = "bottom",
        label.theme = element_text(angle = 90, hjust = 1, vjust = 0.5),
        ticks = FALSE
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )

  extent_data <- habitat_extent_summary %>%
    mutate(
      habitat = stringr::str_to_title(habitat),
      enabling_band = factor(enabling_band, levels = c("<50", "50-69", ">=70"))
    )

  extent_plot <- extent_data %>%
    ggplot(aes(x = habitat, y = share_percent, fill = enabling_band)) +
    geom_col(width = 0.75, colour = "white") +
    coord_flip() +
    scale_fill_manual(
      values = c("<50" = "#fb8500", "50-69" = "#219ebc", ">=70" = "#023047"),
      name = "Enabling band"
    ) +
    labs(
      x = NULL,
      y = "Share of global habitat (%)") +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      panel.grid.major.y = element_blank()
    )

  gap_rank_data <- high_enabling_gaps %>%
    mutate(
      gap_share_percent = if (total_is_valid) {
        100 * upgrade_gap_km2 / global_habitat_total
      } else {
        rep(NA_real_, length(upgrade_gap_km2))
      },
      label = paste0(territory, " — ", habitat)
    ) %>%
    slice_max(order_by = upgrade_gap_km2, n = 12) %>%
    arrange(gap_share_percent)

  gap_plot <- gap_rank_data %>%
    ggplot(aes(x = gap_share_percent, y = factor(label, levels = label))) +
    geom_col(fill = "#1f78b4", width = 0.7) +
    labs(
      x = "Upgrade gap (% of global habitat)",
      y = NULL
      ) +
    scale_x_continuous(labels = label_number(accuracy = 0.01, suffix = "%")) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank()
    )

  combined <- map_plot / (extent_plot | gap_plot) + plot_annotation(
    tag_levels = "A"
  )

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  output_path_abs <- file.path(getwd(), output_path)
  ggplot2::ggsave(output_path_abs, combined, width = 14, height = 14, dpi = 300)
  output_path_abs
}
