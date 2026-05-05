# data clenaing 

library(tidyverse)
library(tsibble)
library(fpp3)


# Helpers -----------------------------------------------------------------

canon_metric <- function(x) {
  x |>
    stringr::str_remove("^\\d+\\.\\s*") |>
    stringr::str_squish()
}

agg_value <- function(x, method) {
  x_non_missing <- x[!is.na(x)]

  if (length(x_non_missing) == 0) {
    return(NA_real_)
  }

  case_when(
    method == "mean" ~ mean(x_non_missing),
    method == "sum"  ~ sum(x_non_missing),
    method == "max"  ~ max(x_non_missing),
    method == "none" ~ dplyr::last(x_non_missing),
    TRUE             ~ NA_real_
  )
}

# Read data ---------------------------------------------------------------

zip_path <- "data/turingAI_forecasting_challenge_dataset.zip"

csv_name <- unzip(zip_path, list = TRUE) |>
  filter(grepl("\\.csv$", Name)) |>
  pull(Name) |>
  first()

nhs_ead <- read_csv(unz(zip_path, csv_name), show_col_types = FALSE)
aggregation_guide <- read_csv("data/aggregation_guide.csv", show_col_types = FALSE)
metric_metadata <- read_csv("data/metric_metadata.csv", show_col_types = FALSE)

# Standardise guide and metadata keys -------------------------------------

aggregation_guide_clean <- aggregation_guide |>
  mutate(metric_key = canon_metric(metric_name)) |>
  arrange(metric_key) |>
  distinct(metric_key, .keep_all = TRUE)

metric_metadata_clean <- metric_metadata |>
  mutate(metric_key = canon_metric(metric_name)) |>
  group_by(metric_key) |>
  summarise(
    frequency = first(frequency),
    freq_label = first(freq_label),
    n_coverages_expected = n_distinct(coverage_label),
    coverages_expected = paste(sort(unique(coverage_label)), collapse = ", "),
    .groups = "drop"
  )

# Outcome -----------------------------------------------------------------

outcome <- nhs_ead |>
  filter(variable_type == "outcome", value >= 0) |>
  mutate(dt = as_date(dt)) |>
  as_tsibble(index = dt)

has_gaps(outcome)

outcome |>
  autoplot(value) +
  labs(
    x = "Date",
    y = "Outcome"
  ) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "lightgrey", fill = NA))

# Features: join guide, align to midday forecast day ----------------------

features_raw <- nhs_ead |>
  filter(variable_type != "outcome", value != -9999) |>
  mutate(
    metric_key = canon_metric(metric_name),
    date = as_date(dt),
    time = format(dt, "%H:%M:%S"),
    forecast_date = if_else(time <= "12:00:00", date, date + days(1))
  ) |>
  left_join(
    aggregation_guide_clean |>
      select(metric_key, temporal_agg, coverage_agg, guide_metric_name = metric_name),
    by = "metric_key"
  ) |>
  left_join(
    metric_metadata_clean,
    by = "metric_key"
  )

# Quality checks ----------------------------------------------------------

# missing_guide_metrics <- features_raw |>
#   filter(is.na(temporal_agg) | is.na(coverage_agg)) |>
#   distinct(metric_name, metric_key)

# if (nrow(missing_guide_metrics) > 0) {
#   warning("Some metrics did not match the aggregation guide after standardisation.")
#   print(missing_guide_metrics)
# }

# Stage 1: temporal aggregation within coverage ---------------------------

daily_by_coverage <- features_raw |>
  group_by(
    forecast_date,
    metric_key,
    metric_name,
    guide_metric_name,
    coverage,
    coverage_label,
    variable_type,
    frequency,
    freq_label,
    n_coverages_expected,
    coverages_expected,
    temporal_agg,
    coverage_agg
  ) |>
  summarise(
    value_daily_coverage = agg_value(value, first(temporal_agg)),
    n_obs_in_window = sum(!is.na(value)),
    window_start = min(dt, na.rm = TRUE),
    window_end = max(dt, na.rm = TRUE),
    .groups = "drop"
  )

# Stage 2: coverage aggregation to one system-level row -------------------

daily_system_long <- daily_by_coverage |>
  group_by(
    forecast_date,
    metric_key,
    metric_name,
    guide_metric_name,
    variable_type,
    frequency,
    freq_label,
    n_coverages_expected,
    coverages_expected,
    temporal_agg,
    coverage_agg
  ) |>
  summarise(
    value = agg_value(value_daily_coverage, first(coverage_agg)),
    n_coverages_observed = n_distinct(coverage_label[!is.na(value_daily_coverage)]),
    n_obs_in_window = sum(n_obs_in_window, na.rm = TRUE),
    window_start = min(window_start, na.rm = TRUE),
    window_end = max(window_end, na.rm = TRUE),
    coverage_missing_flag = n_coverages_observed < first(n_coverages_expected),
    .groups = "drop"
  ) |>
  arrange(forecast_date, metric_key)

# Save long table ---------------------------------------------------------

write_rds(daily_system_long, "data/daily_features_long.rds")

# Build wide modeling table -----------------------------------------------

daily_features_wide <- daily_system_long |>
  select(forecast_date, metric_key, value) |>
  pivot_wider(
    id_cols = forecast_date,
    names_from = metric_key,
    values_from = value
  ) |>
  arrange(forecast_date)

outcome_daily <- nhs_ead |>
  filter(variable_type == "outcome", value >= 0) |>
  transmute(
    forecast_date = as_date(dt),
    outcome = value
  )

nhs_ead_tidy <- outcome_daily |>
  full_join(daily_features_wide, by = "forecast_date") |>
  arrange(forecast_date) |> 
  rename(estimated_avoidable_deaths = outcome)

write_rds(nhs_ead_tidy, "data/nhs_ead_tidy.rds")

