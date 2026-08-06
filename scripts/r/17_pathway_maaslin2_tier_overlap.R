# 17_pathway_maaslin2_tier_overlap.R
# Tier1 vs Tier2 overlap, q<0.25 threshold (matches original pathway
# script's reasoning: Tier2 has too few hits at q<0.05 to be
# meaningful there). Reporting/robustness only -- Tier2 does NOT
# gate candidate selection for pathway (unlike species), per the
# original script's design: selection = Tier1(q<0.25) intersect
# Wilcoxon(q<0.05) only. See corrections notes.

library(here); library(yaml); library(dplyr)

cfg <- yaml::read_yaml(here("config", "config.yaml"))
qval_secondary <- cfg$maaslin2$qval_threshold_secondary

fit_tier1 <- readRDS(here(cfg$paths$checkpoints_dir, "fit_tier1_pathway.rds"))
tier2_group_results <- readRDS(here(cfg$paths$checkpoints_dir, "tier2_group_results_pathway.rds"))

tier1_sig <- fit_tier1$results %>% filter(qval < qval_secondary) %>% pull(feature)
tier2_sig <- tier2_group_results %>% filter(qval < qval_secondary) %>% pull(feature)

overlap <- intersect(tier1_sig, tier2_sig)
comparison_df <- data.frame(feature = union(tier1_sig, tier2_sig)) %>%
  mutate(sig_tier1 = feature %in% tier1_sig, sig_tier2 = feature %in% tier2_sig, sig_both = sig_tier1 & sig_tier2)
write.csv(comparison_df, here(cfg$paths$tables_dir, "step_pathway_tier1_tier2_comparison.csv"), row.names = FALSE)

robust_core <- comparison_df %>% filter(sig_both) %>%
  left_join(fit_tier1$results %>% select(feature, coef_tier1 = coef, qval_tier1 = qval), by = "feature") %>%
  left_join(tier2_group_results %>% select(feature, coef_tier2 = coef, qval_tier2 = qval), by = "feature") %>%
  arrange(qval_tier1)
write.csv(robust_core, here(cfg$paths$tables_dir, "step_pathway_robust_core_tier1_tier2.csv"), row.names = FALSE)

sign_flips <- robust_core %>% filter(sign(coef_tier1) != sign(coef_tier2))
if (nrow(sign_flips) > 0) { cat("Sign flips found:\n"); print(sign_flips) }
stopifnot(nrow(sign_flips) == 0)

overlap_summary <- data.frame(
  metric = c("tier1_sig_q25","tier2_sig_q25","overlap","sign_flips"),
  n = c(length(tier1_sig), length(tier2_sig), length(overlap), nrow(sign_flips))
)
print(overlap_summary)
write.csv(overlap_summary, here(cfg$paths$tables_dir, "step_pathway_tier_overlap_summary.csv"), row.names = FALSE)
