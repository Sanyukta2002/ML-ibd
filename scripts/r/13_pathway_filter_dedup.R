# 13_pathway_filter_dedup.R
library(here); library(yaml); library(dplyr); library(igraph)

cfg <- yaml::read_yaml(here("config", "config.yaml"))
pooled_data_pathway <- readRDS(here(cfg$paths$checkpoints_dir, "pooled_data_pathway.rds"))

meta_cols <- c("sample_id","cohort","group","disease","age","gender","BMI","country","subject_id")
pathway_cols <- setdiff(colnames(pooled_data_pathway), meta_cols)
pw_matrix <- as.matrix(pooled_data_pathway[, pathway_cols])
rownames(pw_matrix) <- pooled_data_pathway$sample_id

row_sum_summary <- summary(rowSums(pw_matrix, na.rm = TRUE))
print(row_sum_summary)
# pathway_abundance is 0-1 scale (unlike species' 0-100) -- confirm median is
# well below 50 as the scale sanity check here, mirroring species' check
stopifnot(row_sum_summary["Median"] < 50)

abund_threshold <- cfg$pathway_filter$abund_threshold
prev_threshold  <- cfg$pathway_filter$prev_threshold
majority_prev   <- cfg$pathway_filter$majority_cohort_prevalence
majority_min_n  <- cfg$pathway_filter$majority_cohort_min_n

mean_abund <- colMeans(pw_matrix, na.rm = TRUE)
prevalence <- colMeans(pw_matrix > 0, na.rm = TRUE)
keep_pathways <- names(mean_abund)[mean_abund >= abund_threshold & prevalence >= prev_threshold]

by_cohort <- sapply(pathway_cols, function(pw) tapply(pw_matrix[, pw] > 0, pooled_data_pathway$cohort, mean))
by_cohort_df <- as.data.frame(t(by_cohort))
n_cohorts_passing <- rowSums(by_cohort_df >= majority_prev)
majority_pathways <- rownames(by_cohort_df)[n_cohorts_passing >= majority_min_n]

final_pathways <- intersect(keep_pathways, majority_pathways)
funnel <- data.frame(
  stage = c("total","pass_abund_prev","pass_majority_cohort","final_intersection"),
  n = c(length(pathway_cols), length(keep_pathways), length(majority_pathways), length(final_pathways))
)
print(funnel)
write.csv(funnel, here(cfg$paths$tables_dir, "step_pathway_filtering_funnel.csv"), row.names = FALSE)

pooled_filtered <- pooled_data_pathway[, c(meta_cols, final_pathways)]
stopifnot(nrow(pooled_filtered) == nrow(pooled_data_pathway))

# systematic dedup (same logic as species Step 5b)
corr_threshold <- cfg$dedup$correlation_flag_threshold
cor_matrix <- cor(as.matrix(pooled_filtered[, final_pathways]), method = cfg$dedup$method)
pairs_idx <- which(cor_matrix > corr_threshold & upper.tri(cor_matrix), arr.ind = TRUE)
pairs_df <- data.frame(f1 = rownames(cor_matrix)[pairs_idx[,1]], f2 = colnames(cor_matrix)[pairs_idx[,2]])
if (nrow(pairs_df) > 0) {
  pairs_df$identical <- mapply(function(a,b) identical(pooled_filtered[[a]], pooled_filtered[[b]]), pairs_df$f1, pairs_df$f2)
} else pairs_df$identical <- logical(0)
write.csv(pairs_df, here(cfg$paths$tables_dir, "step_pathway_duplicate_investigation.csv"), row.names = FALSE)

collapsed <- list(); replaced <- character(0)
near_dup <- pairs_df %>% filter(!identical)
if (nrow(near_dup) > 0) for (i in seq_len(nrow(near_dup))) {
  a <- near_dup$f1[i]; b <- near_dup$f2[i]; nm <- paste0(a, "_combined")
  pooled_filtered[[nm]] <- pooled_filtered[[a]] + pooled_filtered[[b]]
  collapsed[[nm]] <- c(a,b); replaced <- c(replaced, a, b)
}
identical_pairs <- pairs_df %>% filter(identical)
to_drop <- character(0)
if (nrow(identical_pairs) > 0) {
  g <- graph_from_data_frame(identical_pairs[, c("f1","f2")], directed = FALSE)
  cl <- data.frame(pw = names(components(g)$membership), cluster = components(g)$membership)
  to_drop <- cl %>% group_by(cluster) %>% arrange(pw, .by_group = TRUE) %>% filter(row_number() > 1) %>% pull(pw)
}

final_deduped <- c(setdiff(final_pathways, c(to_drop, replaced)), names(collapsed))
dedup_summary <- data.frame(metric = c("pre","dropped_identical","collapsed_near_dup","post"),
                             n = c(length(final_pathways), length(to_drop), length(collapsed), length(final_deduped)))
print(dedup_summary)
write.csv(dedup_summary, here(cfg$paths$tables_dir, "step_pathway_dedup_summary.csv"), row.names = FALSE)

pooled_final <- pooled_filtered[, c(meta_cols, final_deduped)]
stopifnot(nrow(pooled_final) == nrow(pooled_data_pathway))

saveRDS(pooled_final, here(cfg$paths$checkpoints_dir, "pooled_data_pathway_final.rds"))
saveRDS(final_deduped, here(cfg$paths$checkpoints_dir, "final_pathways.rds"))
cat("Final deduped pathway set:", length(final_deduped), "\n")
