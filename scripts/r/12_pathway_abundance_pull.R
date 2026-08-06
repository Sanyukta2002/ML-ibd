# 12_pathway_abundance_pull.R
# Pull pathway_abundance (pinned dates), filter to unstratified
# pathways, merge, subset to baseline samples, drop zero-total sample.

library(here); library(yaml); library(curatedMetagenomicData)
library(mia); library(TreeSummarizedExperiment); library(dplyr); library(tidyr)

cfg <- yaml::read_yaml(here("config", "config.yaml"))
baseline_metadata <- readRDS(here(cfg$paths$checkpoints_dir, "baseline_metadata.rds"))
selected_studies <- cfg$cohorts$selected_studies
snapshot_dates <- cfg$cohorts$pathway_snapshot_dates

dryrun_check <- curatedMetagenomicData(
  paste0("(", paste(selected_studies, collapse = "|"), ")\\.pathway_abundance"),
  dryrun = TRUE
)
dated_studies <- data.frame(resource = dryrun_check) %>%
  mutate(date = sub("^([0-9]{4}-[0-9]{2}-[0-9]{2})\\..*", "\\1", resource),
         study = sub("^[0-9]{4}-[0-9]{2}-[0-9]{2}\\.(.*)\\.pathway_abundance$", "\\1", resource))
for (study in names(snapshot_dates)) {
  stopifnot(snapshot_dates[[study]] %in% dated_studies$date[dated_studies$study == study])
}

pull_patterns <- paste0(snapshot_dates, "\\.", names(snapshot_dates), "\\.pathway_abundance")
pattern_matched <- paste0("(", paste(pull_patterns, collapse = ")|("), ")")
tse_list <- curatedMetagenomicData(pattern_matched, dryrun = FALSE, rownames = "short")
stopifnot(length(tse_list) == length(selected_studies))

filter_to_pathways <- function(tse) {
  rn <- rownames(tse)
  tse[!grepl("\\|", rn) & rn != "UNMAPPED" & rn != "UNINTEGRATED", ]
}
tse_list_filtered <- lapply(tse_list, filter_to_pathways)

tse_merged_raw <- mergeData(tse_list_filtered)
final_sample_ids <- baseline_metadata$sample_id
tse_final <- tse_merged_raw[, colnames(tse_merged_raw) %in% final_sample_ids]
stopifnot(length(setdiff(final_sample_ids, colnames(tse_final))) == 0)

pathway_abund_mat <- t(assay(tse_final, "pathway_abundance"))
pathway_abund_df <- as.data.frame(pathway_abund_mat)
pathway_abund_df$sample_id <- rownames(pathway_abund_df)

meta_for_merge <- baseline_metadata %>%
  select(sample_id, study_name, group, disease, age, gender, BMI, country, subject_id) %>%
  rename(cohort = study_name)
pooled_data_pathway <- meta_for_merge %>% inner_join(pathway_abund_df, by = "sample_id")
stopifnot(nrow(pooled_data_pathway) == nrow(baseline_metadata))
stopifnot(anyDuplicated(pooled_data_pathway$sample_id) == 0)

# drop zero-total sample(s), if any
meta_cols <- c("sample_id","cohort","group","disease","age","gender","BMI","country","subject_id")
pathway_cols <- setdiff(colnames(pooled_data_pathway), meta_cols)
pathway_matrix <- as.matrix(pooled_data_pathway[, pathway_cols])
rownames(pathway_matrix) <- pooled_data_pathway$sample_id
zero_total_samples <- rownames(pathway_matrix)[rowSums(pathway_matrix) == 0]
cat("Zero-total pathway samples dropped:", length(zero_total_samples), "\n")
pooled_data_pathway <- pooled_data_pathway %>% filter(!sample_id %in% zero_total_samples)

pooled_counts <- pooled_data_pathway %>% count(cohort, group) %>%
  pivot_wider(names_from = group, values_from = n, values_fill = 0)
print(pooled_counts)
write.csv(pooled_counts, here(cfg$paths$tables_dir, "step_pathway_pooled_counts.csv"), row.names = FALSE)

saveRDS(pooled_data_pathway, here(cfg$paths$checkpoints_dir, "pooled_data_pathway.rds"))
cat("n samples:", nrow(pooled_data_pathway), "| n pathways:", length(pathway_cols), "\n")
