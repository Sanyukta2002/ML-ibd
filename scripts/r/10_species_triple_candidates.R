# ============================================================
# 10_species_triple_candidates.R
# STEP 13 (corrected): Triple-candidate species selection.
# Tier1 (q<0.05) intersect Tier2 (q<0.05) intersect blocked
# Wilcoxon (q<0.05), THEN a genuine 3-way direction check:
# Wilcoxon-corrected vs. MaAsLin2 BEFORE computing gFC, then gFC
# vs. both -- mirroring the pathway pipeline's rigor, which
# species never had (see corrections notes: species previously
# only checked gFC-vs-MaAsLin2, never verified Wilcoxon's
# direction against anything).
# ============================================================

library(here)
library(yaml)
library(dplyr)

cfg <- yaml::read_yaml(here("config", "config.yaml"))
qval_primary <- cfg$maaslin2$qval_threshold_primary

pooled_data_filtered <- readRDS(here(cfg$paths$checkpoints_dir, "pooled_data_filtered.rds"))
fit_tier1 <- readRDS(here(cfg$paths$checkpoints_dir, "fit_tier1.rds"))
tier2_group_results <- readRDS(here(cfg$paths$checkpoints_dir, "tier2_group_results.rds"))
blocked_wilcox_df <- readRDS(here(cfg$paths$checkpoints_dir, "blocked_wilcox_df.rds"))
name_map2 <- readRDS(here(cfg$paths$checkpoints_dir, "name_map2.rds"))

strip_all <- function(x) tolower(gsub("[^a-zA-Z0-9]", "", x))
strip_maaslin_feature <- function(x) {
  x <- sub("^X(?=[^a-zA-Z0-9])", "", x, perl = TRUE)
  strip_all(x)
}

tier1_results_stripped <- fit_tier1$results
tier1_results_stripped$feature_stripped <- strip_maaslin_feature(tier1_results_stripped$feature)

tier1_sig_stripped  <- strip_maaslin_feature(fit_tier1$results$feature[fit_tier1$results$qval < qval_primary])
tier2_sig_stripped  <- strip_maaslin_feature(tier2_group_results$feature[tier2_group_results$qval < qval_primary])
wilcox_sig_stripped <- strip_all(blocked_wilcox_df$feature[blocked_wilcox_df$qval < qval_primary])

candidate_stripped_raw <- Reduce(intersect, list(tier1_sig_stripped, tier2_sig_stripped, wilcox_sig_stripped))
cat("Raw intersection (Tier1 q<", qval_primary, " ∩ Tier2 q<", qval_primary, " ∩ Wilcoxon q<", qval_primary, "), before direction check:",
    length(candidate_stripped_raw), "\n")

candidate_names <- name_map2$original[match(candidate_stripped_raw, name_map2$stripped)]
stopifnot(sum(is.na(candidate_names)) == 0)

# ============================================================
# 3-way direction check: Wilcoxon-corrected vs MaAsLin2 FIRST
# ============================================================

candidate_pool <- data.frame(
  feature_stripped = candidate_stripped_raw,
  species = candidate_names,
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    tier1_results_stripped %>% dplyr::select(feature_stripped, coef_maaslin_tier1 = coef, qval_maaslin_tier1 = qval),
    by = "feature_stripped"
  ) %>%
  dplyr::left_join(
    blocked_wilcox_df %>% dplyr::mutate(feature_stripped = strip_all(feature)) %>%
      dplyr::select(feature_stripped, statistic_wilcox_corrected = statistic_corrected,
                    pval_wilcox = pval, qval_wilcox = qval),
    by = "feature_stripped"
  ) %>%
  dplyr::mutate(sign_agree_maaslin_wilcox = sign(coef_maaslin_tier1) == sign(statistic_wilcox_corrected)) %>%
  dplyr::arrange(qval_maaslin_tier1)

cat("Direction agreement (MaAsLin2 vs corrected Wilcoxon) within raw intersection:\n")
print(table(candidate_pool$sign_agree_maaslin_wilcox))

if (any(!candidate_pool$sign_agree_maaslin_wilcox)) {
  cat("\nSpecies excluded for MaAsLin2/Wilcoxon direction disagreement:\n")
  print(candidate_pool %>% dplyr::filter(!sign_agree_maaslin_wilcox) %>%
          dplyr::select(species, coef_maaslin_tier1, statistic_wilcox_corrected))
}

triple_candidates <- candidate_pool %>% dplyr::filter(sign_agree_maaslin_wilcox)
cat("\nAfter Wilcoxon/MaAsLin2 direction-agreement filter:", nrow(triple_candidates), "candidates\n")

# ============================================================
# gFC as third independent effect-size check
# ============================================================

q_range_cfg <- cfg$gfc$quantile_range
q_range <- seq(q_range_cfg[1], q_range_cfg[2], q_range_cfg[3])

compute_gfc <- function(x, group, q_range) {
  x_ibd     <- x[group == "IBD"]
  x_control <- x[group == "control"]
  pc <- min(x[x > 0], na.rm = TRUE) / 2
  q_ibd     <- quantile(log10(x_ibd + pc), probs = q_range, na.rm = TRUE)
  q_control <- quantile(log10(x_control + pc), probs = q_range, na.rm = TRUE)
  mean(q_ibd - q_control)
}

triple_candidates$gFC <- sapply(triple_candidates$species, function(sp) {
  compute_gfc(pooled_data_filtered[[sp]], pooled_data_filtered$group, q_range)
})

triple_candidates$gfc_maaslin_agree <- sign(triple_candidates$gFC) == sign(triple_candidates$coef_maaslin_tier1)
triple_candidates$gfc_wilcox_agree  <- sign(triple_candidates$gFC) == sign(triple_candidates$statistic_wilcox_corrected)

disagreements <- triple_candidates %>% dplyr::filter(!gfc_maaslin_agree | !gfc_wilcox_agree)

if (nrow(disagreements) > 0) {
  cat("\n", nrow(disagreements), "species show gFC disagreement -- investigate before finalizing:\n")
  print(disagreements %>% dplyr::select(species, gFC, coef_maaslin_tier1, statistic_wilcox_corrected))
  write.csv(disagreements, here(cfg$paths$tables_dir, "step13_gfc_disagreements_for_investigation.csv"), row.names = FALSE)
}

triple_candidates_clean <- triple_candidates %>% dplyr::filter(gfc_maaslin_agree, gfc_wilcox_agree)

cat("\nFinal triple-candidate species (post 3-way agreement):", nrow(triple_candidates_clean), "\n")
print(triple_candidates_clean$species)

selection_funnel <- data.frame(
  stage = c("raw_intersection_tier1_tier2_wilcoxon", "after_maaslin_wilcoxon_direction_filter",
            "after_gfc_direction_filter"),
  n_species = c(length(candidate_stripped_raw), nrow(triple_candidates), nrow(triple_candidates_clean))
)
print(selection_funnel)
write.csv(selection_funnel, here(cfg$paths$tables_dir, "step13_selection_funnel.csv"), row.names = FALSE)

# ---- Post-selection correlation re-check (confirmatory -- Step 5b
# already deduped upstream, this confirms nothing slipped through
# to this late stage; ported from pathway's equivalent check) ----
if (nrow(triple_candidates_clean) >= 2) {
  cor_matrix_final <- cor(as.matrix(pooled_data_filtered[, triple_candidates_clean$species]), method = cfg$dedup$method)
  high_cor_idx <- which(cor_matrix_final > cfg$dedup$correlation_flag_threshold & upper.tri(cor_matrix_final), arr.ind = TRUE)
  high_cor_df <- data.frame(f1 = rownames(cor_matrix_final)[high_cor_idx[,1]], f2 = colnames(cor_matrix_final)[high_cor_idx[,2]])
  if (nrow(high_cor_df) > 0) {
    high_cor_df$identical <- mapply(function(a,b) identical(pooled_data_filtered[[a]], pooled_data_filtered[[b]]), high_cor_df$f1, high_cor_df$f2)
  } else high_cor_df$identical <- logical(0)
  cat("Post-selection correlation re-check:", nrow(high_cor_df), "high-correlation pairs found (",
      sum(high_cor_df$identical), "truly identical -- would indicate a dedup miss upstream )\n")
  write.csv(high_cor_df, here(cfg$paths$tables_dir, "step13_final_candidates_correlation_recheck.csv"), row.names = FALSE)
  stopifnot(all(!high_cor_df$identical))
}

write.csv(triple_candidates_clean, here(cfg$paths$tables_dir, "step13_triple_candidates.csv"), row.names = FALSE)
saveRDS(triple_candidates_clean, here(cfg$paths$checkpoints_dir, "triple_candidates.rds"))
cat("Checkpoint written: triple_candidates.rds\n")
