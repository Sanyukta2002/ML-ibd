# ============================================================
# 03_species_filter_dedup.R
# STEP 5: Abundance/prevalence filter intersected with majority-
#         cohort rule.
# STEP 5b: Systematic duplicate/near-duplicate scan (ported from
#          pathway pipeline's Step 6) — runs BEFORE any differential
#          abundance testing, so MaAsLin2/Wilcoxon never see
#          duplicate features as independent.
# ============================================================

# NOTE: the original exploratory run found an "Adlercreutzia equolifaciens_1"
# mergeData() name-collision duplicate requiring special-case collapsing.
# Investigated in this rebuild (2026-08-05): NOT present in this pull —
# only two distinct species (A. equolifaciens, A. caecimuris) exist, no
# collision. Consistent with the known taxa-count drift (781 vs 788,
# see 02_species_abundance_pull.R) from package/database version
# differences vs. the unversioned original run. No special-case handling
# needed; the systematic dedup scan below remains as a general safety
# net for whatever duplicate patterns (if any) exist in THIS pull.

library(here)
library(yaml)
library(dplyr)
library(igraph)

cfg <- yaml::read_yaml(here("config", "config.yaml"))

pooled_data <- readRDS(here(cfg$paths$checkpoints_dir, "pooled_data.rds"))

meta_cols <- c("sample_id", "cohort", "group", "disease", "age", "gender", "BMI", "country", "subject_id")
taxa_cols <- setdiff(colnames(pooled_data), meta_cols)

cat("Total taxa columns:", length(taxa_cols), "\n")

abund_matrix <- as.matrix(pooled_data[, taxa_cols])
rownames(abund_matrix) <- pooled_data$sample_id

# ---- Scale sanity check: cMD relative_abundance is 0-100 (percentage) ----
row_sum_summary <- summary(rowSums(abund_matrix, na.rm = TRUE))
print(row_sum_summary)
stopifnot(row_sum_summary["Median"] > 50)

# ============================================================
# STEP 5: Filter 1 (abundance + prevalence) intersected with
# Filter 2 (majority-cohort rule)
# ============================================================

abund_threshold <- cfg$species_filter$abund_threshold
prev_threshold  <- cfg$species_filter$prev_threshold
majority_prev   <- cfg$species_filter$majority_cohort_prevalence
majority_min_n  <- cfg$species_filter$majority_cohort_min_n

mean_rel_abund <- colMeans(abund_matrix, na.rm = TRUE)
prevalence     <- colMeans(abund_matrix > 0, na.rm = TRUE)

keep_species <- names(mean_rel_abund)[mean_rel_abund >= abund_threshold & prevalence >= prev_threshold]

full_species_by_cohort <- sapply(taxa_cols, function(sp) {
  tapply(abund_matrix[, sp] > 0, pooled_data$cohort, mean)
})
full_species_by_cohort_df <- as.data.frame(t(full_species_by_cohort))

n_cohorts_passing <- rowSums(full_species_by_cohort_df >= majority_prev)
majority_species <- rownames(full_species_by_cohort_df)[n_cohorts_passing >= majority_min_n]

final_species <- intersect(keep_species, majority_species)

filtering_summary <- data.frame(
  stage = c("total_taxa", "pass_abundance_prevalence_filter", "pass_majority_cohort_rule", "final_intersection"),
  n_species = c(length(taxa_cols), length(keep_species), length(majority_species), length(final_species))
)
print(filtering_summary)
write.csv(filtering_summary, here(cfg$paths$tables_dir, "step5_species_filtering_funnel.csv"), row.names = FALSE)

species_filter_detail <- data.frame(
  species = taxa_cols,
  mean_rel_abund = mean_rel_abund,
  overall_prevalence = prevalence,
  n_cohorts_passing_threshold = n_cohorts_passing,
  passed_abundance_prevalence_filter = taxa_cols %in% keep_species,
  passed_majority_cohort_rule = taxa_cols %in% majority_species,
  kept_final = taxa_cols %in% final_species
)
write.csv(species_filter_detail, here(cfg$paths$tables_dir, "step5_species_filter_detail.csv"), row.names = FALSE)

pooled_data_filtered <- pooled_data[, c(meta_cols, final_species)]
stopifnot(nrow(pooled_data_filtered) == 353)
cat("Final filtered species set:", length(final_species), "species\n")

# ============================================================
# STEP 5b: Systematic duplicate/near-duplicate scan
# (runs BEFORE MaAsLin2/Wilcoxon — see rationale in corrections notes)
# ============================================================

corr_threshold <- cfg$dedup$correlation_flag_threshold
corr_method <- cfg$dedup$method

species_matrix_for_dedup <- as.matrix(pooled_data_filtered[, final_species])

cor_matrix_species <- cor(species_matrix_for_dedup, method = corr_method)

near_perfect_pairs_sp <- which(cor_matrix_species > corr_threshold & upper.tri(cor_matrix_species), arr.ind = TRUE)

near_perfect_pairs_sp_df <- data.frame(
  species_1 = rownames(cor_matrix_species)[near_perfect_pairs_sp[, 1]],
  species_2 = colnames(cor_matrix_species)[near_perfect_pairs_sp[, 2]]
)

if (nrow(near_perfect_pairs_sp_df) > 0) {
  near_perfect_pairs_sp_df$truly_identical <- mapply(function(s1, s2) {
    identical(pooled_data_filtered[[s1]], pooled_data_filtered[[s2]])
  }, near_perfect_pairs_sp_df$species_1, near_perfect_pairs_sp_df$species_2)

  near_perfect_pairs_sp_df$max_abs_diff <- mapply(function(s1, s2) {
    max(abs(pooled_data_filtered[[s1]] - pooled_data_filtered[[s2]]))
  }, near_perfect_pairs_sp_df$species_1, near_perfect_pairs_sp_df$species_2)
} else {
  near_perfect_pairs_sp_df$truly_identical <- logical(0)
  near_perfect_pairs_sp_df$max_abs_diff <- numeric(0)
}

print(near_perfect_pairs_sp_df)
write.csv(near_perfect_pairs_sp_df, here(cfg$paths$tables_dir, "step5b_species_duplicate_investigation.csv"), row.names = FALSE)

# ---- Near-perfect-but-distinct pairs: COLLAPSED into a combined
# feature (preserves the signal as one feature rather than splitting
# importance credit across two near-identical columns downstream) ----
near_dup_not_identical <- near_perfect_pairs_sp_df %>% dplyr::filter(!truly_identical)

collapsed_features <- list()
species_replaced <- character(0)

if (nrow(near_dup_not_identical) > 0) {
  for (i in seq_len(nrow(near_dup_not_identical))) {
    s1 <- near_dup_not_identical$species_1[i]
    s2 <- near_dup_not_identical$species_2[i]
    combined_name <- paste0(s1, "_combined")
    pooled_data_filtered[[combined_name]] <- pooled_data_filtered[[s1]] + pooled_data_filtered[[s2]]
    collapsed_features[[combined_name]] <- c(s1, s2)
    species_replaced <- c(species_replaced, s1, s2)
  }
}

# ---- Truly identical pairs: cluster-based keep-first-drop-rest
# (handles A=B=C chains, not just simple pairs) ----
identical_pairs_sp_df <- near_perfect_pairs_sp_df %>% dplyr::filter(truly_identical)
species_to_drop <- character(0)

if (nrow(identical_pairs_sp_df) > 0) {
  g_sp <- igraph::graph_from_data_frame(identical_pairs_sp_df[, c("species_1", "species_2")], directed = FALSE)
  clusters_sp <- igraph::components(g_sp)

  cluster_membership_sp <- data.frame(
    species = names(clusters_sp$membership),
    cluster_id = clusters_sp$membership
  )
  write.csv(cluster_membership_sp, here(cfg$paths$tables_dir, "step5b_species_duplicate_clusters.csv"), row.names = FALSE)

  species_to_drop <- cluster_membership_sp %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::arrange(species, .by_group = TRUE) %>%
    dplyr::filter(dplyr::row_number() > 1) %>%
    dplyr::pull(species)
}

cat("Near-perfect-correlation pairs found:", nrow(near_perfect_pairs_sp_df),
    "| truly identical (deduped):", nrow(identical_pairs_sp_df),
    "| distinct-but-collapsed:", length(collapsed_features), "\n")

final_species_deduped <- setdiff(final_species, c(species_to_drop, species_replaced))
final_species_deduped <- c(final_species_deduped, names(collapsed_features))

dedup_summary_sp <- data.frame(
  metric = c("pre_dedup_species_count", "n_dropped_as_exact_duplicates",
             "n_collapsed_as_near_duplicates", "post_dedup_species_count"),
  n = c(length(final_species), length(species_to_drop), length(collapsed_features), length(final_species_deduped))
)
print(dedup_summary_sp)
write.csv(dedup_summary_sp, here(cfg$paths$tables_dir, "step5b_species_dedup_summary.csv"), row.names = FALSE)

pooled_data_filtered <- pooled_data_filtered[, c(meta_cols, final_species_deduped)]
final_species <- final_species_deduped

stopifnot(nrow(pooled_data_filtered) == 353)
cat("Final deduped species set:", length(final_species), "species,", nrow(pooled_data_filtered), "samples\n")

saveRDS(pooled_data_filtered, here(cfg$paths$checkpoints_dir, "pooled_data_filtered.rds"))
saveRDS(final_species, here(cfg$paths$checkpoints_dir, "final_species.rds"))
cat("Checkpoints written: pooled_data_filtered.rds, final_species.rds\n")
