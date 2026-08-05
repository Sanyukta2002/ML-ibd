# ============================================================
# 06_maaslin2_tier1.R
# STEP 8: MaAsLin2 differential abundance, Tier 1 -- all 353
# samples, group as fixed effect, cohort as random effect.
# ============================================================

library(here)
library(yaml)
library(Maaslin2)

cfg <- yaml::read_yaml(here("config", "config.yaml"))

pooled_data_filtered <- readRDS(here(cfg$paths$checkpoints_dir, "pooled_data_filtered.rds"))
final_species <- readRDS(here(cfg$paths$checkpoints_dir, "final_species.rds"))

# ---- Prepare input_data: species abundance, samples as rownames ----
input_data <- as.data.frame(pooled_data_filtered[, final_species])
rownames(input_data) <- pooled_data_filtered$sample_id

# ---- Prepare input_metadata ----
input_metadata <- as.data.frame(pooled_data_filtered[, c("sample_id", "cohort", "group")])
rownames(input_metadata) <- input_metadata$sample_id
input_metadata$sample_id <- NULL

stopifnot(identical(rownames(input_data), rownames(input_metadata)))

# ---- "control" as reference level: positive coefficient = enriched in IBD ----
input_metadata$group  <- relevel(factor(input_metadata$group), ref = "control")
input_metadata$cohort <- factor(input_metadata$cohort)
stopifnot(levels(input_metadata$group)[1] == "control")

fit_tier1 <- Maaslin2(
  input_data      = input_data,
  input_metadata  = input_metadata,
  output          = here(cfg$paths$maaslin2_dir, "species_tier1"),
  fixed_effects   = c("group"),
  random_effects  = c("cohort"),
  normalization   = cfg$maaslin2$normalization,
  transform       = cfg$maaslin2$transform,
  analysis_method = cfg$maaslin2$analysis_method,
  correction      = cfg$maaslin2$correction,
  min_prevalence  = 0,      # already filtered upstream
  min_abundance   = 0,      # already filtered upstream
  standardize     = TRUE
)

tier1_hit_summary <- data.frame(
  qval_threshold = c(cfg$maaslin2$qval_threshold_primary, cfg$maaslin2$qval_threshold_secondary),
  n_significant  = c(
    sum(fit_tier1$results$qval < cfg$maaslin2$qval_threshold_primary, na.rm = TRUE),
    sum(fit_tier1$results$qval < cfg$maaslin2$qval_threshold_secondary, na.rm = TRUE)
  ),
  n_total_features_tested = nrow(fit_tier1$results)
)

print(tier1_hit_summary)
write.csv(tier1_hit_summary, here(cfg$paths$tables_dir, "step8_maaslin2_tier1_hit_summary.csv"), row.names = FALSE)

saveRDS(fit_tier1, here(cfg$paths$checkpoints_dir, "fit_tier1.rds"))
cat("Checkpoint written: fit_tier1.rds\n")
cat("MaAsLin2 output written to:", here(cfg$paths$maaslin2_dir, "species_tier1"), "\n")
