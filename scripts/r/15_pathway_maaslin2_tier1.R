# 15_pathway_maaslin2_tier1.R
library(here); library(yaml); library(Maaslin2)

cfg <- yaml::read_yaml(here("config", "config.yaml"))
dir.create(here(cfg$paths$maaslin2_dir), recursive = TRUE, showWarnings = FALSE)
pooled_final <- readRDS(here(cfg$paths$checkpoints_dir, "pooled_data_pathway_final.rds"))
final_pathways <- readRDS(here(cfg$paths$checkpoints_dir, "final_pathways.rds"))

input_data <- as.data.frame(pooled_final[, final_pathways])
rownames(input_data) <- pooled_final$sample_id
input_metadata <- as.data.frame(pooled_final[, c("sample_id","cohort","group")])
rownames(input_metadata) <- input_metadata$sample_id
input_metadata$sample_id <- NULL
stopifnot(identical(rownames(input_data), rownames(input_metadata)))

input_metadata$group  <- relevel(factor(input_metadata$group), ref = "control")
input_metadata$cohort <- factor(input_metadata$cohort)

fit_tier1 <- Maaslin2(
  input_data = input_data, input_metadata = input_metadata,
  output = here(cfg$paths$maaslin2_dir, "pathway_tier1"),
  fixed_effects = c("group"), random_effects = c("cohort"),
  normalization = cfg$maaslin2$normalization, transform = cfg$maaslin2$transform,
  analysis_method = cfg$maaslin2$analysis_method, correction = cfg$maaslin2$correction,
  min_prevalence = 0, min_abundance = 0, standardize = TRUE, plot_scatter = FALSE
)

hit_summary <- data.frame(
  qval_threshold = c(cfg$maaslin2$qval_threshold_primary, cfg$maaslin2$qval_threshold_secondary),
  n_significant = c(sum(fit_tier1$results$qval < cfg$maaslin2$qval_threshold_primary, na.rm=TRUE),
                     sum(fit_tier1$results$qval < cfg$maaslin2$qval_threshold_secondary, na.rm=TRUE)),
  n_total_features_tested = nrow(fit_tier1$results)
)
print(hit_summary)
write.csv(hit_summary, here(cfg$paths$tables_dir, "step_pathway_maaslin2_tier1_hit_summary.csv"), row.names = FALSE)
saveRDS(fit_tier1, here(cfg$paths$checkpoints_dir, "fit_tier1_pathway.rds"))
