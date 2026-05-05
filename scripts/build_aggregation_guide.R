# build_aggregation_guide.R
# Creates an aggregation guide for NHS EAD metrics, specifying how to:
#   1. Aggregate sub-daily observations to daily (temporal aggregation)
#   2. Aggregate across coverage areas/sites (coverage aggregation)
#
# Inputs:
#   - data/metric_metadata.csv   (metric names, coverages, frequencies)
#   - data/metric_details.xlsx   (metric descriptions)
#   - `feat` data frame in session (unique metric names)
#
# Output:
#   - data/aggregation_guide.csv

library(dplyr)
library(stringr)
library(readxl)

# -- Load metadata -----------------------------------------------------------

meta_csv  <- read.csv("data/metric_metadata.csv")
meta_xlsx <- read_excel("data/metric_details.xlsx")

# Unique metric list — either from session or derived from metadata
if (exists("feat")) {
  feat_df <- feat
} else {
  feat_df <- meta_csv |>
    distinct(metric_name) |>
    as_tibble()
}

# -- Summarise metadata per metric --------------------------------------------

meta_summary <- meta_csv |>
  group_by(metric_name) |>
  summarise(
    coverages   = paste(unique(coverage_label), collapse = ", "),
    n_coverages = n_distinct(coverage_label),
    frequency   = first(frequency),
    freq_label  = first(freq_label),
    .groups = "drop"
  )

joined <- feat_df |>
  left_join(meta_summary, by = "metric_name") |>
  left_join(meta_xlsx, by = c("metric_name" = "Metric Name"))

# -- Helper: regex detector (case-insensitive) --------------------------------

rx <- function(pattern) function(x) str_detect(x, regex(pattern, ignore_case = TRUE))

# -- Classify temporal aggregation (sub-daily -> daily) -----------------------
#
# Priority order matters — earlier rules take precedence.
# "Since Midnight" cumulative metrics must be caught before other patterns.

agg_guide <- joined |>
  mutate(
    temporal_agg = case_when(
      # Already daily or unknown frequency
      frequency >= 1440 | is.na(frequency) ~ "none",

      # --- Cumulative "Since Midnight" counters: take end-of-day value -------
      rx("since midnight")(metric_name) ~ "max",

      # --- Longest wait: daily peak ------------------------------------------
      rx("longest wait")(metric_name) ~ "max",

      # --- OPEL / escalation scores: worst case ------------------------------
      rx("OPEL|escalation level")(metric_name) ~ "max",

      # --- 90th percentile response times: daily worst case ------------------
      # These are NOT counts — summing percentiles is meaningless
      rx("90th response")(metric_name) ~ "max",

      # --- Discharge counts per interval: sum for daily total ----------------
      rx("P0 discharges|complex discharges")(metric_name) ~ "sum",

      # --- Flow counts in rolling windows: sum for daily total ---------------
      rx("last hour|last 15 min|in last")(metric_name) ~ "sum",

      # --- Rates, percentages, averages: daily mean --------------------------
      rx("^average|^mean|performance|%|proportion")(metric_name) ~ "mean",

      # --- Counts reported per interval (5-hourly): sum ----------------------
      rx("number of home visits|number of treatment centre")(metric_name) ~ "sum",

      # --- Snapshot levels (queues, patient counts, capacity, DTAs) ----------
      rx(str_c(
        "queue|patient count|patients in|en route|active calls|",
        "waiting calls|cases on|beds available|capacity|resuscitation|",
        "DTA|cohorting|slots"
      ))(metric_name) ~ "mean",

      # Fallback
      TRUE ~ "mean"
    ),

    temporal_reason = case_when(
      temporal_agg == "none" ~
        "Already at daily granularity; no temporal aggregation needed.",

      temporal_agg == "max" ~
        "Cumulative counter resetting at midnight; the last (max) value of the day equals the daily total.",

      temporal_agg == "max" & rx("longest wait")(metric_name) ~
        "Peak wait time; daily maximum captures the worst case.",
      temporal_agg == "max" & rx("OPEL|escalation")(metric_name) ~
        "Escalation/pressure score; daily maximum reflects the highest pressure point.",
      temporal_agg == "max" & rx("90th response")(metric_name) ~
        "Percentile-based response time; daily maximum captures worst hourly performance.",

      temporal_agg == "sum" & rx("P0 discharges|complex discharges")(metric_name) ~
        "Discharge count per interval; summing gives the daily total.",
      temporal_agg == "sum" & rx("home visits|treatment centre")(metric_name) ~
        "Activity count per reporting interval; summing approximates daily total.",
      temporal_agg == "sum" ~
        "Flow count over a rolling window; summing intervals gives the daily total.",

      temporal_agg == "mean" & rx("average|mean|time to")(metric_name) ~
        "Time or average metric; daily mean represents typical performance across the day.",
      temporal_agg == "mean" & rx("%|performance|proportion")(metric_name) ~
        "Rate or percentage; daily mean reflects average performance over the day.",
      temporal_agg == "mean" & rx("queue|patient|DTA|en route|calls|cases|beds|capacity|resuscitation|cohorting|slots")(metric_name) ~
        "Point-in-time snapshot; daily mean captures average load throughout the day.",

      temporal_agg == "mean" ~
        "Point-in-time or rate metric; daily mean summarises typical intra-day values.",
      TRUE ~ ""
    ),

    # -- Classify coverage aggregation (across sites -> system-wide) ----------

    coverage_agg = case_when(
      n_coverages <= 1 | is.na(n_coverages) ~ "none",

      # OPEL / escalation: worst site
      rx("OPEL|escalation")(metric_name) ~ "max",

      # Rates, percentages, averages: mean across sites
      rx(str_c(
        "performance|%|proportion|average|mean|",
        "time to|occupancy|speed|handling|median"
      ))(metric_name) ~ "mean",

      # Counts / volumes: sum across sites
      rx(str_c(
        "patient|breach|DTA|handover|conveyed|queue|arrival|discharge|",
        "calls|referral|NCtR|schedule|beds|home visit|appointment|cases|",
        "en route|admission|attend|time lost|capacity|resuscitation|cohorting"
      ))(metric_name) ~ "sum",

      # Fallback for multi-site
      TRUE ~ "sum"
    ),

    coverage_reason = case_when(
      coverage_agg == "none" ~
        "Single coverage area; no cross-site aggregation needed.",
      coverage_agg == "max" ~
        "Escalation score; system-wide pressure is defined by the most stressed site.",
      coverage_agg == "mean" ~
        "Rate, average, or percentage; mean across sites approximates the system-wide value (weighted mean preferred if volumes available).",
      coverage_agg == "sum" ~
        "Count or volume metric; summing across sites gives the system-wide total.",
      TRUE ~ ""
    )
  )

# -- Build final table --------------------------------------------------------

aggregation_guide <- agg_guide |>
  select(
    metric_name, Description,
    frequency, freq_label, n_coverages, coverages,
    temporal_agg, temporal_reason,
    coverage_agg, coverage_reason
  )

# -- Save ---------------------------------------------------------------------

write.csv(aggregation_guide, "data/aggregation_guide.csv", row.names = FALSE)

cat("Aggregation guide saved to data/aggregation_guide.csv\n")
cat(sprintf("  Total metrics: %d\n", nrow(aggregation_guide)))
cat("\n--- Temporal aggregation summary ---\n")
print(count(aggregation_guide, temporal_agg))
cat("\n--- Coverage aggregation summary ---\n")
print(count(aggregation_guide, coverage_agg))
