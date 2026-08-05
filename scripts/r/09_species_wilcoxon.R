# ============================================================
# 09_species_wilcoxon.R
# STEP 11: Species name mapping utility (strip_all lookup, to
#          safely cross-reference MaAsLin2's sanitized feature
#          names against the original species names).
# STEP 12 (corrected): Blocked Wilcoxon test (coin package,
#          cohort as block), WITH empirically verified + corrected
#          sign convention -- coin::wilcox_test's raw statistic
#          does not match MaAsLin2's "positive = enriched in IBD"
#          convention (same issue the pathway pipeline found and
#          fixed; species never had this correction -- see
#          corrections notes). Verified here directly against real
#          group medians, not assumed to transfer from pathway.
# ============================================================

library(here)
library(yaml)
library(dplyr)
library(coin)

cfg <- yaml::read_yaml(here("config", "config.yaml"))

pooled_data_filtered <- readRDS(here(cfg$paths$checkpoints_dir, "pooled_data_filtered.rds"))
final_species <- readRDS(here(cfg$paths$checkpoints_dir, "final_species.rds"))
fit_tier1 <- readRDS(here(cfg$paths$checkpoints_dir, "fit_tier1.rds"))

# ============================================================
# STEP 11: Name mapping (strip everything non-alphanumeric ->
# robust lookup key between original names and MaAsLin2's
# sanitized feature names)
# ============================================================

# MaAsLin2 prepends "X" to sanitized feature names that would otherwise
# start with a non-alphanumeric character after its internal cleanup --
# this includes names starting with a digit (pathway pipeline's known
# case) AND bracket-notation species names like "[Ruminococcus] torques"
# (MaAsLin2 converts "[" -> ".", then prepends X since the result starts
# with a dot). Confirmed via direct investigation 2026-08-05: 8/65 Tier1
# significant species use bracket notation and were silently failing to
# map without this fix.
strip_maaslin_feature <- function(x) {
  x <- sub("^X(?=[^a-zA-Z0-9])", "", x, perl = TRUE)
  strip_all(x)
}

strip_all <- function(x) tolower(gsub("[^a-zA-Z0-9]", "", x))

name_map2 <- data.frame(
  original = final_species,
  stripped = strip_all(final_species),
  stringsAsFactors = FALSE
)

# ---- Confirm the stripped mapping is injective (no collisions) ----
stopifnot(anyDuplicated(name_map2$stripped) == 0)

# ---- Sanity check: does this successfully map Tier 1's significant species? ----
sig_species_raw <- fit_tier1$results$feature[fit_tier1$results$qval < cfg$maaslin2$qval_threshold_primary]
sig_stripped <- strip_maaslin_feature(sig_species_raw)
sig_species_original <- name_map2$original[match(sig_stripped, name_map2$stripped)]

stopifnot(sum(is.na(sig_species_original)) == 0)
cat("name_map2 built and verified: all", length(sig_species_raw), "Tier 1 significant species mapped successfully\n")

saveRDS(name_map2, here(cfg$paths$checkpoints_dir, "name_map2.rds"))

# ============================================================
# STEP 12: Blocked Wilcoxon test, with sign-convention correction
# ============================================================

wilcox_data <- pooled_data_filtered %>%
  dplyr::select(sample_id, cohort, group, dplyr::all_of(final_species))

wilcox_data$cohort <- factor(wilcox_data$cohort)
wilcox_data$group  <- factor(wilcox_data$group, levels = c("control", "IBD"))

blocked_wilcox_results <- lapply(final_species, function(sp) {
  df <- wilcox_data[, c("group", "cohort", sp)]
  colnames(df) <- c("group", "cohort", "abundance")
  tryCatch({
    wt <- coin::wilcox_test(abundance ~ group | cohort, data = df, distribution = "asymptotic")
    data.frame(
      feature             = sp,
      statistic           = as.numeric(coin::statistic(wt)),
      statistic_corrected = -as.numeric(coin::statistic(wt)),
      pval                = as.numeric(coin::pvalue(wt))
    )
  }, error = function(e) {
    data.frame(feature = sp, statistic = NA, statistic_corrected = NA, pval = NA)
  })
})

blocked_wilcox_df <- do.call(rbind, blocked_wilcox_results)
blocked_wilcox_df$qval <- p.adjust(blocked_wilcox_df$pval, method = cfg$maaslin2$correction)
blocked_wilcox_df <- blocked_wilcox_df[order(blocked_wilcox_df$pval), ]

n_failed <- sum(is.na(blocked_wilcox_df$pval))
if (n_failed > 0) {
  cat("WARNING:", n_failed, "species failed the Wilcoxon test:\n")
  print(blocked_wilcox_df[is.na(blocked_wilcox_df$pval), ])
}

# ---- Empirical verification of the sign correction: pick the
# species with the strongest signal, check its ACTUAL group medians,
# confirm the corrected statistic's sign matches reality ----
check_sp <- blocked_wilcox_df$feature[which.min(blocked_wilcox_df$pval)]
medians_check <- wilcox_data %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(median_abund = median(.data[[check_sp]]), .groups = "drop")

cat("Sign-convention empirical check on:", check_sp, "\n")
print(medians_check)

raw_stat <- blocked_wilcox_df$statistic[blocked_wilcox_df$feature == check_sp]
corrected_stat <- blocked_wilcox_df$statistic_corrected[blocked_wilcox_df$feature == check_sp]
cat("Raw statistic:", raw_stat, "| Corrected statistic:", corrected_stat, "\n")

ibd_higher <- medians_check$median_abund[medians_check$group == "IBD"] >
              medians_check$median_abund[medians_check$group == "control"]
corrected_positive <- corrected_stat > 0

stopifnot(ibd_higher == corrected_positive)  # hard stop: correction direction must match reality
cat("Verified: corrected statistic sign matches actual group medians (IBD higher =", ibd_higher, ").\n")

wilcox_hit_summary <- data.frame(
  qval_threshold = c(cfg$maaslin2$qval_threshold_primary, cfg$maaslin2$qval_threshold_secondary),
  n_significant  = c(
    sum(blocked_wilcox_df$qval < cfg$maaslin2$qval_threshold_primary, na.rm = TRUE),
    sum(blocked_wilcox_df$qval < cfg$maaslin2$qval_threshold_secondary, na.rm = TRUE)
  ),
  n_total_features_tested = nrow(blocked_wilcox_df),
  n_failed_tests = n_failed
)

print(wilcox_hit_summary)
write.csv(wilcox_hit_summary, here(cfg$paths$tables_dir, "step12_wilcoxon_hit_summary.csv"), row.names = FALSE)
write.csv(blocked_wilcox_df, here(cfg$paths$tables_dir, "step12_wilcoxon_full_results.csv"), row.names = FALSE)

saveRDS(blocked_wilcox_df, here(cfg$paths$checkpoints_dir, "blocked_wilcox_df.rds"))
cat("Checkpoint written: blocked_wilcox_df.rds\n")
