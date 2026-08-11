# 16_pathway_maaslin2_tier2_age.R
library(here); library(yaml); library(dplyr); library(Maaslin2)

cfg <- yaml::read_yaml(here("config", "config.yaml"))
dir.create(here(cfg$paths$maaslin2_dir), recursive = TRUE, showWarnings = FALSE)
pooled_final <- readRDS(here(cfg$paths$checkpoints_dir, "pooled_data_pathway_final.rds"))
final_pathways <- readRDS(here(cfg$paths$checkpoints_dir, "final_pathways.rds"))
covariate_availability <- readRDS(here(cfg$paths$checkpoints_dir, "covariate_availability.rds"))

age_cohorts <- covariate_availability$age
data_tier2 <- pooled_final %>% filter(cohort %in% age_cohorts, !is.na(age))
cat("Tier2 pathway subset n =", nrow(data_tier2), "\n")
stopifnot(nrow(data_tier2) > 0)

input_data_t2 <- as.data.frame(data_tier2[, final_pathways])
rownames(input_data_t2) <- data_tier2$sample_id
input_metadata_t2 <- as.data.frame(data_tier2[, c("sample_id","cohort","group","age")])
rownames(input_metadata_t2) <- input_metadata_t2$sample_id
input_metadata_t2$sample_id <- NULL
stopifnot(identical(rownames(input_data_t2), rownames(input_metadata_t2)))

input_metadata_t2$group  <- relevel(factor(input_metadata_t2$group), ref = "control")
input_metadata_t2$cohort <- factor(input_metadata_t2$cohort)

fit_tier2 <- Maaslin2(
  input_data = input_data_t2, input_metadata = input_metadata_t2,
  output = here(cfg$paths$maaslin2_dir, "pathway_tier2_age"),
  fixed_effects = c("group","age"), random_effects = c("cohort"),
  normalization = cfg$maaslin2$normalization, transform = cfg$maaslin2$transform,
  analysis_method = cfg$maaslin2$analysis_method, correction = cfg$maaslin2$correction,
  min_prevalence = 0, min_abundance = 0, standardize = TRUE, plot_scatter = FALSE
)

tier2_group_results <- fit_tier2$results %>% filter(metadata == "group")
hit_summary <- data.frame(
  qval_threshold = c(cfg$maaslin2$qval_threshold_primary, cfg$maaslin2$qval_threshold_secondary),
  n_significant = c(sum(tier2_group_results$qval < cfg$maaslin2$qval_threshold_primary, na.rm=TRUE),
                     sum(tier2_group_results$qval < cfg$maaslin2$qval_threshold_secondary, na.rm=TRUE)),
  n_total_features_tested = nrow(tier2_group_results), n_samples = nrow(data_tier2)
)
print(hit_summary)
write.csv(hit_summary, here(cfg$paths$tables_dir, "step_pathway_maaslin2_tier2_hit_summary.csv"), row.names = FALSE)
saveRDS(fit_tier2, here(cfg$paths$checkpoints_dir, "fit_tier2_pathway.rds"))
saveRDS(tier2_group_results, here(cfg$paths$checkpoints_dir, "tier2_group_results_pathway.rds"))
