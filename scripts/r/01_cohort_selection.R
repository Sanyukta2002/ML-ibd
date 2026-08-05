# ============================================================
# 01_cohort_selection.R
# Pull sampleMetadata, filter to our 4 cohorts + stool + IBD/control,
# apply NielsenHB_2014's ESP restriction, baseline-filter to one row
# per subject. Produces the checkpoint everything downstream depends on.
# ============================================================

library(here)
library(yaml)
library(curatedMetagenomicData)
library(dplyr)
library(tidyr)

cfg <- yaml::read_yaml(here("config", "config.yaml"))

dir.create(here(cfg$paths$checkpoints_dir), recursive = TRUE, showWarnings = FALSE)
dir.create(here(cfg$paths$tables_dir),      recursive = TRUE, showWarnings = FALSE)

writeLines(
  capture.output(sessionInfo()),
  here(cfg$paths$checkpoints_dir, "sessionInfo_01_cohort_selection.txt")
)

# ---- Cohort selection ----
selected_studies <- cfg$cohorts$selected_studies
esp_restriction_study <- "NielsenHB_2014"
esp_country <- cfg$cohorts$nielsen_country_restriction

full_metadata <- sampleMetadata

cohort_metadata <- full_metadata %>%
  filter(study_name %in% selected_studies) %>%
  filter(body_site == "stool")

cohort_metadata <- cohort_metadata %>%
  filter(!(study_name == esp_restriction_study & country != esp_country))

n_by_cohort_raw <- cohort_metadata %>% count(study_name, name = "n_samples_raw")
print(n_by_cohort_raw)
write.csv(n_by_cohort_raw, here(cfg$paths$tables_dir, "step2_cohort_raw_counts.csv"), row.names = FALSE)

# ---- Baseline filtering ----
timepoint_availability <- cohort_metadata %>%
  group_by(study_name) %>%
  summarise(
    has_visit_number = any(!is.na(visit_number)),
    has_days_col      = any(!is.na(days_from_first_collection)),
    .groups = "drop"
  )
print(timepoint_availability)
write.csv(timepoint_availability, here(cfg$paths$tables_dir, "step3_timepoint_availability.csv"), row.names = FALSE)

baseline_metadata <- cohort_metadata %>%
  group_by(study_name, subject_id) %>%
  arrange(study_name, subject_id, dplyr::coalesce(visit_number, days_from_first_collection, 0)) %>%
  slice(1) %>%
  ungroup()

baseline_metadata <- baseline_metadata %>%
  mutate(group = case_when(
    disease == "healthy" ~ "control",
    grepl("IBD", disease, ignore.case = TRUE) ~ "IBD",
    TRUE ~ disease
  ))

n_by_cohort_group <- baseline_metadata %>%
  count(study_name, group) %>%
  pivot_wider(names_from = group, values_from = n, values_fill = 0)
print(n_by_cohort_group)
write.csv(n_by_cohort_group, here(cfg$paths$tables_dir, "step3_baseline_counts_by_group.csv"), row.names = FALSE)

# ---- Sanity checks ----
one_row_per_subject_check <- baseline_metadata %>%
  group_by(study_name) %>%
  summarise(n_rows = n(), n_unique_subjects = n_distinct(subject_id), .groups = "drop")
stopifnot(all(one_row_per_subject_check$n_rows == one_row_per_subject_check$n_unique_subjects))

cat("Total baseline samples:", nrow(baseline_metadata), "\n")
stopifnot(nrow(baseline_metadata) == 353)   # known-good count from the original run

saveRDS(baseline_metadata, here(cfg$paths$checkpoints_dir, "baseline_metadata.rds"))
cat("Checkpoint written: baseline_metadata.rds\n")
