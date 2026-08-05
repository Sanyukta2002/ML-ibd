# ============================================================
# 02_species_abundance_pull.R
# Pull relative_abundance TSEs for the 4 cohorts, PINNED to
# specific snapshot dates (validated via dry-run before pulling —
# 3 of 4 cohorts have 2 available dates for relative_abundance;
# an unpinned pull is ambiguous and not reproducible). Merge,
# subset to the 353 baseline samples, pool with metadata.
# ============================================================

# NOTE on taxa count (781, not 788 as in the original exploratory run):
# investigated and confirmed NOT caused by snapshot-date ambiguity (both
# unpinned and pinned pulls give 781). Most likely explanation: the
# original run was interactive/unversioned with no recorded sessionInfo,
# so it cannot be forensically reproduced; curatedMetagenomicData/mia
# version drift between then and now is the probable cause. 781 is
# treated as the correct, current, fully-reproducible number (pinned
# dates + versioned container, see sessionInfo_02... in checkpoints/).

library(here)
library(yaml)
library(curatedMetagenomicData)
library(mia)
library(TreeSummarizedExperiment)
library(dplyr)
library(tidyr)

cfg <- yaml::read_yaml(here("config", "config.yaml"))

baseline_metadata <- readRDS(here(cfg$paths$checkpoints_dir, "baseline_metadata.rds"))
selected_studies <- cfg$cohorts$selected_studies
snapshot_dates <- cfg$cohorts$species_snapshot_dates

# ---- Dry-run validation: confirm our pinned dates actually exist
# BEFORE pulling, so a catalog change fails loudly, not silently ----
dryrun_check <- curatedMetagenomicData(
  paste0("(", paste(selected_studies, collapse = "|"), ")\\.relative_abundance"),
  dryrun = TRUE
)

dated_studies <- data.frame(resource = dryrun_check) %>%
  dplyr::mutate(
    date  = sub("^([0-9]{4}-[0-9]{2}-[0-9]{2})\\..*", "\\1", resource),
    study = sub("^[0-9]{4}-[0-9]{2}-[0-9]{2}\\.(.*)\\.relative_abundance$", "\\1", resource)
  )

print(dated_studies)

for (study in names(snapshot_dates)) {
  expected_date <- snapshot_dates[[study]]
  stopifnot(expected_date %in% dated_studies$date[dated_studies$study == study])
}
cat("All pinned snapshot dates confirmed available.\n")

# ---- Pull, pinned to exact dates (not a bare study-name pattern) ----
pull_patterns <- paste0(snapshot_dates, "\\.", names(snapshot_dates), "\\.relative_abundance")
species_pattern_matched <- paste0("(", paste(pull_patterns, collapse = ")|("), ")")

tse_list <- curatedMetagenomicData(species_pattern_matched, dryrun = FALSE, rownames = "short")

stopifnot(length(tse_list) == length(selected_studies))

pulled_studies <- names(tse_list) %>%
  sub("^[0-9]{4}-[0-9]{2}-[0-9]{2}\\.", "", .) %>%
  sub("\\.relative_abundance$", "", .)
stopifnot(setequal(pulled_studies, selected_studies))

cat("Pulled relative_abundance for (pinned dates):\n")
print(names(tse_list))

# ---- Log raw per-cohort taxa counts BEFORE merge, for provenance ----
raw_taxa_counts <- data.frame(
  cohort = pulled_studies,
  n_taxa_raw = sapply(tse_list, nrow),
  n_samples_raw = sapply(tse_list, ncol)
)
print(raw_taxa_counts)
write.csv(raw_taxa_counts, here(cfg$paths$tables_dir, "step4_raw_taxa_counts_per_cohort.csv"), row.names = FALSE)

# ---- Merge into one TSE, subset to our 353 baseline samples ----
tse_merged_raw <- mergeData(tse_list)
n_taxa_post_merge <- nrow(tse_merged_raw)

final_sample_ids <- baseline_metadata$sample_id
tse_final <- tse_merged_raw[, colnames(tse_merged_raw) %in% final_sample_ids]

missing_samples <- setdiff(final_sample_ids, colnames(tse_final))
stopifnot(length(missing_samples) == 0)

cat("tse_final dimensions (taxa x samples):", paste(dim(tse_final), collapse = " x "), "\n")
cat("Taxa count: pre-merge max per-cohort =", max(raw_taxa_counts$n_taxa_raw),
    "| post-merge =", n_taxa_post_merge,
    "| final (after mia rowTree matching) =", nrow(tse_final), "\n")

write.csv(
  data.frame(
    stage = c("max_per_cohort_raw", "post_merge", "final_after_rowtree_match"),
    n_taxa = c(max(raw_taxa_counts$n_taxa_raw), n_taxa_post_merge, nrow(tse_final))
  ),
  here(cfg$paths$tables_dir, "step4_taxa_count_provenance.csv"),
  row.names = FALSE
)

# ---- Extract abundance matrix (samples x taxa), attach metadata ----
abundance_mat <- t(assay(tse_final, "relative_abundance"))
abundance_df  <- as.data.frame(abundance_mat)
abundance_df$sample_id <- rownames(abundance_df)

meta_for_merge <- baseline_metadata %>%
  select(sample_id, study_name, group, disease, age, gender, BMI, country, subject_id) %>%
  rename(cohort = study_name)

pooled_data <- meta_for_merge %>%
  dplyr::inner_join(abundance_df, by = "sample_id")

cat("pooled_data sample count:", nrow(pooled_data),
    "(should match baseline_metadata count from Step 1-3 -- verified via join, not a magic number)\n")
stopifnot(nrow(pooled_data) == nrow(baseline_metadata))  # internal consistency: join shouldn't drop/add samples
stopifnot(anyDuplicated(pooled_data$sample_id) == 0)

cat("pooled_data dimensions:", paste(dim(pooled_data), collapse = " x "), "\n")

pooled_counts <- pooled_data %>%
  count(cohort, group) %>%
  pivot_wider(names_from = group, values_from = n, values_fill = 0)

print(pooled_counts)
write.csv(pooled_counts, here(cfg$paths$tables_dir, "step4_pooled_counts_by_group.csv"), row.names = FALSE)

saveRDS(pooled_data, here(cfg$paths$checkpoints_dir, "pooled_data.rds"))
cat("Checkpoint written: pooled_data.rds\n")
