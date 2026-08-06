# 19_pathway_triple_candidates.R
# Pathway candidate selection: Tier1(q<0.25) INTERSECT Wilcoxon(q<0.05),
# THEN direction-agreement (Wilcoxon-corrected vs MaAsLin2), THEN gFC
# as a third independent check. Tier2 does NOT gate (reporting-only,
# per original pathway script's design -- underpowered in this rebuild
# too, see tier_overlap rule). Post-selection correlation re-check
# included (item 8 -- confirmatory, Step 6/13 already deduped upstream).

library(here); library(yaml); library(dplyr)

cfg <- yaml::read_yaml(here("config", "config.yaml"))
qval_secondary <- cfg$maaslin2$qval_threshold_secondary
qval_primary <- cfg$maaslin2$qval_threshold_primary

pooled_final <- readRDS(here(cfg$paths$checkpoints_dir, "pooled_data_pathway_final.rds"))
fit_tier1 <- readRDS(here(cfg$paths$checkpoints_dir, "fit_tier1_pathway.rds"))
blocked_wilcox_df <- readRDS(here(cfg$paths$checkpoints_dir, "blocked_wilcox_df_pathway.rds"))
name_map_pw <- readRDS(here(cfg$paths$checkpoints_dir, "name_map_pathway.rds"))

strip_all <- function(x) tolower(gsub("[^a-zA-Z0-9]", "", x))
strip_maaslin_feature <- function(x) {
  x <- sub("^X(?=[^a-zA-Z0-9])", "", x, perl = TRUE)
  strip_all(x)
}

tier1_results_stripped <- fit_tier1$results
tier1_results_stripped$feature_stripped <- strip_maaslin_feature(tier1_results_stripped$feature)

tier1_sig_stripped    <- strip_maaslin_feature(fit_tier1$results$feature[fit_tier1$results$qval < qval_secondary])
wilcox_sig_stripped   <- strip_all(blocked_wilcox_df$feature[blocked_wilcox_df$qval < qval_primary])
candidate_stripped_raw <- intersect(tier1_sig_stripped, wilcox_sig_stripped)
cat("Raw intersection (Tier1 q<", qval_secondary, " AND Wilcoxon q<", qval_primary, "):", length(candidate_stripped_raw), "\n")

candidate_names <- name_map_pw$original[match(candidate_stripped_raw, name_map_pw$stripped)]
stopifnot(sum(is.na(candidate_names)) == 0)

candidate_pool <- data.frame(feature_stripped = candidate_stripped_raw, pathway = candidate_names, stringsAsFactors = FALSE) %>%
  left_join(tier1_results_stripped %>% select(feature_stripped, coef_maaslin_tier1 = coef, qval_maaslin_tier1 = qval), by = "feature_stripped") %>%
  left_join(blocked_wilcox_df %>% mutate(feature_stripped = strip_all(feature)) %>%
              select(feature_stripped, statistic_wilcox_corrected = statistic_corrected, pval_wilcox = pval, qval_wilcox = qval), by = "feature_stripped") %>%
  mutate(sign_agree_maaslin_wilcox = sign(coef_maaslin_tier1) == sign(statistic_wilcox_corrected)) %>%
  arrange(qval_maaslin_tier1)

cat("Direction agreement (MaAsLin2 vs corrected Wilcoxon):\n"); print(table(candidate_pool$sign_agree_maaslin_wilcox))
if (any(!candidate_pool$sign_agree_maaslin_wilcox)) {
  cat("Excluded for direction disagreement:\n")
  print(candidate_pool %>% filter(!sign_agree_maaslin_wilcox) %>% select(pathway, coef_maaslin_tier1, statistic_wilcox_corrected))
}

triple_candidates <- candidate_pool %>% filter(sign_agree_maaslin_wilcox)
cat("After direction filter:", nrow(triple_candidates), "\n")

q_range_cfg <- cfg$gfc$quantile_range
q_range <- seq(q_range_cfg[1], q_range_cfg[2], q_range_cfg[3])
compute_gfc <- function(x, group, q_range) {
  x_ibd <- x[group == "IBD"]; x_control <- x[group == "control"]
  pc <- min(x[x > 0], na.rm = TRUE) / 2
  mean(quantile(log10(x_ibd + pc), probs = q_range, na.rm = TRUE) - quantile(log10(x_control + pc), probs = q_range, na.rm = TRUE))
}
triple_candidates$gFC <- sapply(triple_candidates$pathway, function(pw) compute_gfc(pooled_final[[pw]], pooled_final$group, q_range))
triple_candidates$gfc_maaslin_agree <- sign(triple_candidates$gFC) == sign(triple_candidates$coef_maaslin_tier1)
triple_candidates$gfc_wilcox_agree  <- sign(triple_candidates$gFC) == sign(triple_candidates$statistic_wilcox_corrected)

disagreements <- triple_candidates %>% filter(!gfc_maaslin_agree | !gfc_wilcox_agree)
if (nrow(disagreements) > 0) {
  cat(nrow(disagreements), "pathways show gFC disagreement -- investigate:\n")
  print(disagreements %>% select(pathway, gFC, coef_maaslin_tier1, statistic_wilcox_corrected))
  write.csv(disagreements, here(cfg$paths$tables_dir, "step_pathway_gfc_disagreements.csv"), row.names = FALSE)
}

triple_candidates_clean <- triple_candidates %>% filter(gfc_maaslin_agree, gfc_wilcox_agree)
cat("Final triple-candidate pathways:", nrow(triple_candidates_clean), "\n")
print(triple_candidates_clean$pathway)

# ---- Item 8: post-selection correlation re-check (confirmatory) ----
if (nrow(triple_candidates_clean) >= 2) {
  cor_matrix_final <- cor(as.matrix(pooled_final[, triple_candidates_clean$pathway]), method = cfg$dedup$method)
  high_cor_idx <- which(cor_matrix_final > cfg$dedup$correlation_flag_threshold & upper.tri(cor_matrix_final), arr.ind = TRUE)
  high_cor_df <- data.frame(f1 = rownames(cor_matrix_final)[high_cor_idx[,1]], f2 = colnames(cor_matrix_final)[high_cor_idx[,2]])
  if (nrow(high_cor_df) > 0) {
    high_cor_df$identical <- mapply(function(a,b) identical(pooled_final[[a]], pooled_final[[b]]), high_cor_df$f1, high_cor_df$f2)
  } else high_cor_df$identical <- logical(0)
  cat("Post-selection correlation re-check:", nrow(high_cor_df), "high-correlation pairs found (",
      sum(high_cor_df$identical), "truly identical -- would indicate a dedup miss upstream )\n")
  write.csv(high_cor_df, here(cfg$paths$tables_dir, "step_pathway_final_candidates_correlation_recheck.csv"), row.names = FALSE)
  stopifnot(all(!high_cor_df$identical))
}

write.csv(triple_candidates_clean, here(cfg$paths$tables_dir, "step_pathway_triple_candidates.csv"), row.names = FALSE)
saveRDS(triple_candidates_clean, here(cfg$paths$checkpoints_dir, "triple_candidates_pathway.rds"))
