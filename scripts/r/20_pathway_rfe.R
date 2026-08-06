# 20_pathway_rfe.R
# RFE on the 18 pathway triple-candidates -> 3 panels derived from
# THIS run's curve (same design as species: optimum / parsimony /
# full pool). No Adlercreutzia-equivalent issue here -- dedup already
# resolved upstream (Step 6/13).

library(here); library(yaml); library(caret); library(randomForest); library(ggplot2)
source(here("scripts", "r", "utils", "plot_theme.R"))

cfg <- yaml::read_yaml(here("config", "config.yaml"))
triple_candidates <- readRDS(here(cfg$paths$checkpoints_dir, "triple_candidates_pathway.rds"))
pooled_final <- readRDS(here(cfg$paths$checkpoints_dir, "pooled_data_pathway_final.rds"))

candidate_pathways <- triple_candidates$pathway
n_candidates <- length(candidate_pathways)
cat("RFE candidate pool:", n_candidates, "pathways\n")

feat_matrix <- as.data.frame(pooled_final[, candidate_pathways])
colnames(feat_matrix) <- make.names(colnames(feat_matrix))
y <- factor(pooled_final$group, levels = c("control", "IBD"))

name_map_rfe <- data.frame(original = candidate_pathways, sanitized = make.names(candidate_pathways), stringsAsFactors = FALSE)

rfFuncs_auc <- rfFuncs
rfFuncs_auc$summary <- twoClassSummary
rfFuncs_auc$fit <- function(x, y, first, last, ...) randomForest(x, y, importance = TRUE, ...)
rfFuncs_auc$pred <- function(object, x) {
  tmp <- predict(object, x, type = "prob"); out <- predict(object, x)
  data.frame(pred = out, obs = NA, tmp, check.names = FALSE)
}
ctrl_auc <- rfeControl(functions = rfFuncs_auc, method = "cv", number = cfg$rfe$cv_folds, saveDetails = TRUE)

set.seed(cfg$rfe$seed)
rfe_result <- rfe(x = feat_matrix, y = y, sizes = 2:n_candidates, rfeControl = ctrl_auc, metric = "ROC")
print(rfe_result)
saveRDS(rfe_result, here(cfg$paths$checkpoints_dir, "rfe_result_pathway.rds"))

results_df <- rfe_result$results
elbow_sizes <- data.frame(
  tolerance_pct = c(1, 2, 3),
  selected_size = c(
    caret::pickSizeTolerance(results_df, metric = "ROC", tol = 1, maximize = TRUE),
    caret::pickSizeTolerance(results_df, metric = "ROC", tol = 2, maximize = TRUE),
    caret::pickSizeTolerance(results_df, metric = "ROC", tol = 3, maximize = TRUE)
  )
)
print(elbow_sizes)
write.csv(results_df, here(cfg$paths$tables_dir, "step_pathway_rfe_results_by_size.csv"), row.names = FALSE)
write.csv(elbow_sizes, here(cfg$paths$tables_dir, "step_pathway_rfe_elbow_sizes.csv"), row.names = FALSE)

optimum_size <- results_df$Variables[which.max(results_df$ROC)]
parsimony_size <- elbow_sizes$selected_size[elbow_sizes$tolerance_pct == 1]
cat("Optimum:", optimum_size, "| Parsimony:", parsimony_size, "| Full pool:", n_candidates, "\n")

rfe_plot <- ggplot(results_df, aes(x = Variables, y = ROC)) +
  geom_ribbon(aes(ymin = ROC - ROCSD / sqrt(cfg$rfe$cv_folds), ymax = ROC + ROCSD / sqrt(cfg$rfe$cv_folds)), alpha = 0.15) +
  geom_line(color = "#0072B2") + geom_point(size = 1.6, color = "#0072B2") +
  geom_vline(xintercept = optimum_size, linetype = "solid", color = "#D55E00") +
  geom_vline(xintercept = parsimony_size, linetype = "dashed", color = "#009E73") +
  geom_vline(xintercept = n_candidates, linetype = "dotted", color = "grey40") +
  labs(title = "RFE (pathways): cross-validated AUC vs. panel size",
       subtitle = paste0(n_candidates, "-pathway candidate set | ", cfg$rfe$cv_folds, "-fold CV | ",
                          "solid: optimum | dashed: parsimonious (1% tol.) | dotted: full pool"),
       x = "Number of pathways", y = "Cross-validated ROC (AUC)") + theme_pipeline()
save_figure(rfe_plot, "step_pathway_rfe_auc_vs_panelsize", width = 8, height = 5)

varImp_rfe <- varImp(rfe_result)
build_panel <- function(n) {
  top_sanitized <- rownames(varImp_rfe)[order(-varImp_rfe$Overall)][1:n]
  top_readable <- name_map_rfe$original[match(top_sanitized, name_map_rfe$sanitized)]
  stopifnot(sum(is.na(top_readable)) == 0)
  top_readable
}

panel_definitions <- list(rfe_optimum = optimum_size, rfe_parsimony = parsimony_size, rfe_full_pool = n_candidates)
panels <- list(); size_to_name <- list()
for (nm in names(panel_definitions)) {
  n <- panel_definitions[[nm]]; key <- as.character(n)
  if (!is.null(size_to_name[[key]])) {
    cat("NOTE:", nm, "(n=", n, ") identical to", size_to_name[[key]], "-- not rebuilding\n")
  } else {
    panels[[nm]] <- build_panel(n); size_to_name[[key]] <- nm
  }
}

panel_definition_summary <- data.frame(
  panel_name = names(panel_definitions), panel_size = unlist(panel_definitions),
  definition = c("RFE global optimum (best CV ROC)", "Tolerance-based parsimonious size (1% of optimum)", "Full triple-candidate pool (unranked)")
)
print(panel_definition_summary)
write.csv(panel_definition_summary, here(cfg$paths$tables_dir, "step_pathway_panel_definitions.csv"), row.names = FALSE)

panel_summary <- data.frame(panel_name = names(panels), panel_size = sapply(panels, length), pathways = sapply(panels, paste, collapse = "; "))
print(panel_summary[, c("panel_name", "panel_size")])
write.csv(panel_summary, here(cfg$paths$tables_dir, "step_pathway_final_panels.csv"), row.names = FALSE)

for (nm in names(panels)) saveRDS(panels[[nm]], here(cfg$paths$checkpoints_dir, paste0("panel_pathway_", nm, ".rds")))
saveRDS(panels, here(cfg$paths$checkpoints_dir, "all_panels_pathway.rds"))
cat("Checkpoints written for panels:", paste(names(panels), collapse = ", "), "\n")
