# ============================================================
# 11_species_rfe.R
# STEP 14: RFE (Random Forest, AUC-optimized) on the triple-
# candidate species -> 3 panels derived from THIS run's RFE curve.
#
# NOTE: Adlercreutzia-specific investigation from the original
# script REMOVED -- fully subsumed by the systematic dedup in
# Step 5b, which runs before differential abundance testing.
#
# NOTE: candidate pool is 19 species (not the original run's 35 --
# see Step 13 commit message for why).
#
# Panel design (matches Zheng et al. 2024 Nat Med, Ning et al. 2023
# Nat Commun, Zhou et al. 2024 Front Mol Biosci -- iterative feature
# elimination to a smaller, stable panel rather than assuming the
# full set is optimal):
#   1. rfe_optimum   -- RFE's own best cross-validated ROC size
#   2. rfe_parsimony -- tolerance-based elbow (1%): smallest size
#                        within noise-tolerance of the optimum
#   3. rfe_full_pool -- all triple-candidates, unranked -- upper
#                        bound reference
# ============================================================

library(here)
library(yaml)
library(caret)
library(randomForest)
library(ggplot2)
source(here("scripts", "r", "utils", "plot_theme.R"))

cfg <- yaml::read_yaml(here("config", "config.yaml"))

triple_candidates <- readRDS(here(cfg$paths$checkpoints_dir, "triple_candidates.rds"))
pooled_data_filtered <- readRDS(here(cfg$paths$checkpoints_dir, "pooled_data_filtered.rds"))

candidate_species <- triple_candidates$species
n_candidates <- length(candidate_species)
cat("RFE candidate pool:", n_candidates, "species\n")

feat_matrix <- as.data.frame(pooled_data_filtered[, candidate_species])
colnames(feat_matrix) <- make.names(colnames(feat_matrix))
y <- factor(pooled_data_filtered$group, levels = c("control", "IBD"))

name_map_rfe <- data.frame(
  original  = candidate_species,
  sanitized = make.names(candidate_species),
  stringsAsFactors = FALSE
)

rfFuncs_auc <- rfFuncs
rfFuncs_auc$summary <- twoClassSummary
rfFuncs_auc$fit <- function(x, y, first, last, ...) {
  randomForest(x, y, importance = TRUE, ...)
}
rfFuncs_auc$pred <- function(object, x) {
  tmp <- predict(object, x, type = "prob")
  out <- predict(object, x)
  data.frame(pred = out, obs = NA, tmp, check.names = FALSE)
}

ctrl_auc <- rfeControl(
  functions   = rfFuncs_auc,
  method      = "cv",
  number      = cfg$rfe$cv_folds,
  saveDetails = TRUE
)

set.seed(cfg$rfe$seed)
rfe_result <- rfe(
  x          = feat_matrix,
  y          = y,
  sizes      = 2:n_candidates,
  rfeControl = ctrl_auc,
  metric     = "ROC"
)

print(rfe_result)
saveRDS(rfe_result, here(cfg$paths$checkpoints_dir, "rfe_result.rds"))

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
write.csv(results_df,  here(cfg$paths$tables_dir, "step14_rfe_results_by_size.csv"), row.names = FALSE)
write.csv(elbow_sizes, here(cfg$paths$tables_dir, "step14_rfe_elbow_sizes.csv"), row.names = FALSE)

# ---- Derive the 3 panel sizes from this run's actual curve ----
optimum_size <- results_df$Variables[which.max(results_df$ROC)]
parsimony_size <- elbow_sizes$selected_size[elbow_sizes$tolerance_pct == 1]

cat("RFE optimum size:", optimum_size, "(ROC =", max(results_df$ROC), ")\n")
cat("Parsimonious (1% tolerance) size:", parsimony_size, "\n")
cat("Full candidate pool size:", n_candidates, "\n")

# ---- Figure ----
rfe_plot <- ggplot(results_df, aes(x = Variables, y = ROC)) +
  geom_ribbon(aes(ymin = ROC - ROCSD / sqrt(cfg$rfe$cv_folds), ymax = ROC + ROCSD / sqrt(cfg$rfe$cv_folds)), alpha = 0.15) +
  geom_line(color = "#0072B2") +
  geom_point(size = 1.6, color = "#0072B2") +
  geom_vline(xintercept = optimum_size, linetype = "solid", color = "#D55E00") +
  geom_vline(xintercept = parsimony_size, linetype = "dashed", color = "#009E73") +
  geom_vline(xintercept = n_candidates, linetype = "dotted", color = "grey40") +
  labs(
    title = "RFE: cross-validated AUC vs. panel size",
    subtitle = paste0(n_candidates, "-species candidate set | ", cfg$rfe$cv_folds, "-fold CV | ",
                       "solid: RFE optimum | dashed: parsimonious (1% tol.) | dotted: full pool"),
    x = "Number of species", y = "Cross-validated ROC (AUC)"
  ) +
  theme_pipeline()

save_figure(rfe_plot, "step14_rfe_auc_vs_panelsize", width = 8, height = 5)

# ---- Build the 3 panels ----
varImp_rfe <- varImp(rfe_result)

build_panel <- function(n) {
  top_sanitized <- rownames(varImp_rfe)[order(-varImp_rfe$Overall)][1:n]
  top_readable  <- name_map_rfe$original[match(top_sanitized, name_map_rfe$sanitized)]
  stopifnot(sum(is.na(top_readable)) == 0)
  top_readable
}

panel_definitions <- list(
  rfe_optimum   = optimum_size,
  rfe_parsimony = parsimony_size,
  rfe_full_pool = n_candidates
)

panels <- list()
size_to_panel_name <- list()

for (panel_name in names(panel_definitions)) {
  n <- panel_definitions[[panel_name]]
  size_key <- as.character(n)
  if (!is.null(size_to_panel_name[[size_key]])) {
    cat("NOTE:", panel_name, "(n=", n, ") is identical to", size_to_panel_name[[size_key]],
        "-- not rebuilding, same species set\n")
  } else {
    panels[[panel_name]] <- build_panel(n)
    size_to_panel_name[[size_key]] <- panel_name
  }
}

panel_definition_summary <- data.frame(
  panel_name = names(panel_definitions),
  panel_size = unlist(panel_definitions),
  definition = c("RFE global optimum (best CV ROC)",
                 "Tolerance-based parsimonious size (1% of optimum)",
                 "Full triple-candidate pool (unranked)")
)
print(panel_definition_summary)
write.csv(panel_definition_summary, here(cfg$paths$tables_dir, "step14_panel_definitions.csv"), row.names = FALSE)

panel_summary <- data.frame(
  panel_name = names(panels),
  panel_size = sapply(panels, length),
  species = sapply(panels, paste, collapse = "; ")
)

print(panel_summary[, c("panel_name", "panel_size")])
write.csv(panel_summary, here(cfg$paths$tables_dir, "step14_final_panels.csv"), row.names = FALSE)

for (name in names(panels)) {
  saveRDS(panels[[name]], here(cfg$paths$checkpoints_dir, paste0("panel_", name, ".rds")))
}
saveRDS(panels, here(cfg$paths$checkpoints_dir, "all_panels.rds"))

cat("Checkpoints written for panels:", paste(names(panels), collapse = ", "), "\n")
