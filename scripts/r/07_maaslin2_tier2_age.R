# ============================================================
# 07_maaslin2_tier2_age.R
# STEP 9: MaAsLin2 Tier 2 -- age-adjusted robustness check.
# Cohort subset derived from actual data (covariate_availability.rds),
# NOT the original script's hardcoded 3-cohort/n=314 assumption --
# see 05_covariate_availability_check.R for why that drifted.
# ============================================================

library(here)
library(yaml)
library(dplyr)
library(Maaslin2)

cfg <- yaml::read_yaml(here("config", "config.yaml"))
dir.create(here(cfg$paths$maaslin2_dir), recursive = TRUE, showWarnings = FALSE)

pooled_data_filtered <- readRDS(here(cfg$paths$checkpoints_dir, "pooled_data_filtered.rds"))
final_species <- readRDS(here(cfg$paths$checkpoints_dir, "final_species.rds"))
covariate_availability <- readRDS(here(cfg$paths$checkpoints_dir, "covariate_availability.rds"))

age_cohorts <- covariate_availability$age
cat("Tier 2 age-adjusted cohorts (data-derived):", paste(age_cohorts, collapse = ", "), "\n")

data_tier2 <- pooled_data_filtered %>%
  filter(cohort %in% age_cohorts, !is.na(age))

cat("Tier 2 subset n =", nrow(data_tier2), "\n")
stopifnot(nrow(data_tier2) > 0)

input_data_t2 <- as.data.frame(data_tier2[, final_species])
rownames(input_data_t2) <- data_tier2$sample_id

input_metadata_t2 <- as.data.frame(data_tier2[, c("sample_id", "cohort", "group", "age")])
rownames(input_metadata_t2) <- input_metadata_t2$sample_id
input_metadata_t2$sample_id <- NULL

stopifnot(identical(rownames(input_data_t2), rownames(input_metadata_t2)))

input_metadata_t2$group  <- relevel(factor(input_metadata_t2$group), ref = "control")
input_metadata_t2$cohort <- factor(input_metadata_t2$cohort)
# age stays numeric -- MaAsLin2's standardize=TRUE z-scores it internally

fit_tier2 <- Maaslin2(
  input_data      = input_data_t2,
  input_metadata  = input_metadata_t2,
  output          = here(cfg$paths$maaslin2_dir, "species_tier2_age"),
  fixed_effects   = c("group", "age"),
  random_effects  = c("cohort"),
  normalization   = cfg$maaslin2$normalization,
  transform       = cfg$maaslin2$transform,
  analysis_method = cfg$maaslin2$analysis_method,
  correction      = cfg$maaslin2$correction,
  min_prevalence  = 0,
  min_abundance   = 0,
  standardize     = TRUE,
  plot_scatter    = FALSE   # MaAsLin2's continuous-covariate scatter plot code
                             # uses ggplot2 syntax incompatible with the
                             # installed ggplot2 (>=3.4.0, unnamed labs()
                             # args) -- crashes on the age scatter plot
                             # specifically. Model fitting is unaffected;
                             # we only consume $results downstream, not
                             # the diagnostic plots. Confirmed via log:
                             # fit completed (boxplots for group already
                             # generated) before the age-scatter step failed.
)

# ---- Extract results for the 'group' term specifically ----
tier2_group_results <- fit_tier2$results %>% filter(metadata == "group")

tier2_hit_summary <- data.frame(
  qval_threshold = c(cfg$maaslin2$qval_threshold_primary, cfg$maaslin2$qval_threshold_secondary),
  n_significant  = c(
    sum(tier2_group_results$qval < cfg$maaslin2$qval_threshold_primary, na.rm = TRUE),
    sum(tier2_group_results$qval < cfg$maaslin2$qval_threshold_secondary, na.rm = TRUE)
  ),
  n_total_features_tested = nrow(tier2_group_results),
  n_samples = nrow(data_tier2)
)

print(tier2_hit_summary)
write.csv(tier2_hit_summary, here(cfg$paths$tables_dir, "step9_maaslin2_tier2_hit_summary.csv"), row.names = FALSE)

saveRDS(fit_tier2, here(cfg$paths$checkpoints_dir, "fit_tier2.rds"))
saveRDS(tier2_group_results, here(cfg$paths$checkpoints_dir, "tier2_group_results.rds"))
cat("Checkpoints written: fit_tier2.rds, tier2_group_results.rds\n")
