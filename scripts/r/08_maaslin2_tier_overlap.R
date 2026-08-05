# ============================================================
# 08_maaslin2_tier_overlap.R
# STEP 10: Tier 1 vs Tier 2 overlap -- does the disease signal
# survive age-adjustment and cohort restriction? Sign-flip check
# on the overlapping set (hard stop if any found).
# ============================================================

library(here)
library(yaml)
library(dplyr)

cfg <- yaml::read_yaml(here("config", "config.yaml"))
qval_primary <- cfg$maaslin2$qval_threshold_primary

fit_tier1 <- readRDS(here(cfg$paths$checkpoints_dir, "fit_tier1.rds"))
tier2_group_results <- readRDS(here(cfg$paths$checkpoints_dir, "tier2_group_results.rds"))

tier1_sig <- fit_tier1$results %>% dplyr::filter(qval < qval_primary) %>% dplyr::pull(feature)
tier2_sig <- tier2_group_results %>% dplyr::filter(qval < qval_primary) %>% dplyr::pull(feature)

overlap    <- intersect(tier1_sig, tier2_sig)
only_tier1 <- setdiff(tier1_sig, tier2_sig)
only_tier2 <- setdiff(tier2_sig, tier1_sig)

comparison_df <- data.frame(feature = union(tier1_sig, tier2_sig)) %>%
  dplyr::mutate(
    sig_tier1 = feature %in% tier1_sig,
    sig_tier2 = feature %in% tier2_sig,
    sig_both  = sig_tier1 & sig_tier2
  )

write.csv(comparison_df, here(cfg$paths$tables_dir, "step10_tier1_tier2_comparison.csv"), row.names = FALSE)

robust_core <- comparison_df %>%
  dplyr::filter(sig_both) %>%
  dplyr::left_join(fit_tier1$results %>% dplyr::select(feature, coef_tier1 = coef, qval_tier1 = qval), by = "feature") %>%
  dplyr::left_join(tier2_group_results %>% dplyr::select(feature, coef_tier2 = coef, qval_tier2 = qval), by = "feature") %>%
  dplyr::arrange(qval_tier1)

write.csv(robust_core, here(cfg$paths$tables_dir, "step10_robust_core_tier1_tier2.csv"), row.names = FALSE)

# ---- Sign-flip check: hard stop if MaAsLin2 disagrees with itself
# across tiers on direction (would indicate a genuine model instability,
# not just a robustness edge case) ----
sign_flips <- robust_core %>% dplyr::filter(sign(coef_tier1) != sign(coef_tier2))

if (nrow(sign_flips) > 0) {
  cat("Sign flips found between Tier1 and Tier2:\n")
  print(sign_flips)
}
stopifnot(nrow(sign_flips) == 0)

tier_overlap_summary <- data.frame(
  metric = c("tier1_significant", "tier2_significant", "overlap_both_tiers",
             "only_tier1", "only_tier2", "sign_flips_in_overlap"),
  n = c(length(tier1_sig), length(tier2_sig), length(overlap),
        length(only_tier1), length(only_tier2), nrow(sign_flips))
)

print(tier_overlap_summary)
write.csv(tier_overlap_summary, here(cfg$paths$tables_dir, "step10_tier_overlap_summary.csv"), row.names = FALSE)

saveRDS(robust_core, here(cfg$paths$checkpoints_dir, "robust_core_tier1_tier2.rds"))
cat("Checkpoint written: robust_core_tier1_tier2.rds\n")
