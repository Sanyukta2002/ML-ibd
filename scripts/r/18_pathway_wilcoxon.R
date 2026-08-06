# 18_pathway_wilcoxon.R
# Blocked Wilcoxon + sign-convention correction, re-verified empirically
# against THIS run's data (original pathway script already had this
# correction -- re-confirming it still holds, not assuming it transfers).

library(here); library(yaml); library(dplyr); library(coin)

cfg <- yaml::read_yaml(here("config", "config.yaml"))
pooled_final <- readRDS(here(cfg$paths$checkpoints_dir, "pooled_data_pathway_final.rds"))
final_pathways <- readRDS(here(cfg$paths$checkpoints_dir, "final_pathways.rds"))
fit_tier1 <- readRDS(here(cfg$paths$checkpoints_dir, "fit_tier1_pathway.rds"))

strip_all <- function(x) tolower(gsub("[^a-zA-Z0-9]", "", x))
strip_maaslin_feature <- function(x) {
  x <- sub("^X(?=[^a-zA-Z0-9])", "", x, perl = TRUE)
  strip_all(x)
}

name_map_pw <- data.frame(original = final_pathways, stripped = strip_all(final_pathways), stringsAsFactors = FALSE)
stopifnot(anyDuplicated(name_map_pw$stripped) == 0)

wilcox_data <- pooled_final %>% select(sample_id, cohort, group, all_of(final_pathways))
wilcox_data$cohort <- factor(wilcox_data$cohort)
wilcox_data$group  <- factor(wilcox_data$group, levels = c("control", "IBD"))

blocked_results <- lapply(final_pathways, function(pw) {
  df <- wilcox_data[, c("group", "cohort", pw)]; colnames(df) <- c("group", "cohort", "abundance")
  tryCatch({
    wt <- coin::wilcox_test(abundance ~ group | cohort, data = df, distribution = "asymptotic")
    data.frame(feature = pw, statistic = as.numeric(coin::statistic(wt)),
               statistic_corrected = -as.numeric(coin::statistic(wt)), pval = as.numeric(coin::pvalue(wt)))
  }, error = function(e) data.frame(feature = pw, statistic = NA, statistic_corrected = NA, pval = NA))
})
blocked_wilcox_df <- do.call(rbind, blocked_results)
blocked_wilcox_df$qval <- p.adjust(blocked_wilcox_df$pval, method = cfg$maaslin2$correction)
blocked_wilcox_df <- blocked_wilcox_df[order(blocked_wilcox_df$pval), ]

n_failed <- sum(is.na(blocked_wilcox_df$pval))
if (n_failed > 0) cat("WARNING:", n_failed, "pathways failed the test\n")

# empirical sign-convention re-verification
check_pw <- blocked_wilcox_df$feature[which.min(blocked_wilcox_df$pval)]
medians_check <- wilcox_data %>% group_by(group) %>% summarise(median_abund = median(.data[[check_pw]]), .groups = "drop")
cat("Sign-convention empirical check on:", check_pw, "\n"); print(medians_check)
raw_stat <- blocked_wilcox_df$statistic[blocked_wilcox_df$feature == check_pw]
corrected_stat <- blocked_wilcox_df$statistic_corrected[blocked_wilcox_df$feature == check_pw]
cat("Raw:", raw_stat, "| Corrected:", corrected_stat, "\n")
ibd_higher <- medians_check$median_abund[medians_check$group == "IBD"] > medians_check$median_abund[medians_check$group == "control"]
stopifnot(ibd_higher == (corrected_stat > 0))
cat("Verified: corrected statistic sign matches actual group medians.\n")

hit_summary <- data.frame(
  qval_threshold = c(cfg$maaslin2$qval_threshold_primary, cfg$maaslin2$qval_threshold_secondary),
  n_significant = c(sum(blocked_wilcox_df$qval < cfg$maaslin2$qval_threshold_primary, na.rm=TRUE),
                     sum(blocked_wilcox_df$qval < cfg$maaslin2$qval_threshold_secondary, na.rm=TRUE)),
  n_total_features_tested = nrow(blocked_wilcox_df), n_failed_tests = n_failed
)
print(hit_summary)
write.csv(hit_summary, here(cfg$paths$tables_dir, "step_pathway_wilcoxon_hit_summary.csv"), row.names = FALSE)
write.csv(blocked_wilcox_df, here(cfg$paths$tables_dir, "step_pathway_wilcoxon_full_results.csv"), row.names = FALSE)
saveRDS(blocked_wilcox_df, here(cfg$paths$checkpoints_dir, "blocked_wilcox_df_pathway.rds"))
saveRDS(name_map_pw, here(cfg$paths$checkpoints_dir, "name_map_pathway.rds"))
