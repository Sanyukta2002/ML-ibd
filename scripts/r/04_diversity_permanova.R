# ============================================================
# 04_diversity_permanova.R
# STEP 6: PERMANOVA (Bray-Curtis) — cohort vs. group, plus
#         covariate submodels (age/gender/BMI).
# STEP 7: Alpha diversity (Shannon/Simpson/richness), computed on
#         the FULL raw abundance (pre-filter), intentionally.
# ============================================================

library(here)
library(yaml)
library(dplyr)
library(vegan)
library(nlme)
library(ggpubr)
library(patchwork)
source(here("scripts", "r", "utils", "plot_theme.R"))

cfg <- yaml::read_yaml(here("config", "config.yaml"))
dir.create(here(cfg$paths$figures_dir), recursive = TRUE, showWarnings = FALSE)

pooled_data <- readRDS(here(cfg$paths$checkpoints_dir, "pooled_data.rds"))
meta_cols <- c("sample_id", "cohort", "group", "disease", "age", "gender", "BMI", "country", "subject_id")
taxa_cols <- setdiff(colnames(pooled_data), meta_cols)

abund_matrix <- as.matrix(pooled_data[, taxa_cols])
rownames(abund_matrix) <- pooled_data$sample_id

# ============================================================
# STEP 6: PERMANOVA
# ============================================================

bray_dist <- vegdist(abund_matrix, method = "bray")
meta_for_permanova <- pooled_data[, meta_cols]

stopifnot(sum(is.na(meta_for_permanova$cohort)) == 0)
stopifnot(sum(is.na(meta_for_permanova$group)) == 0)

set.seed(cfg$rfe$seed)
permanova_basic <- adonis2(bray_dist ~ cohort + group, data = meta_for_permanova, permutations = 999, by = "margin")
print(permanova_basic)

run_covariate_permanova <- function(covariate, cohorts_with_covariate) {
  d <- pooled_data %>% filter(cohort %in% cohorts_with_covariate, !is.na(.data[[covariate]]))
  abund_sub <- as.matrix(d[, taxa_cols])
  rownames(abund_sub) <- d$sample_id
  bray_sub <- vegdist(abund_sub, method = "bray")
  formula_str <- paste("bray_sub ~ cohort + group +", covariate)
  set.seed(cfg$rfe$seed)
  result <- adonis2(as.formula(formula_str), data = d[, meta_cols], permutations = 999, by = "margin")
  list(n = nrow(d), result = result)
}

covariate_availability <- readRDS(here(cfg$paths$checkpoints_dir, "covariate_availability.rds"))
age_check    <- run_covariate_permanova("age", covariate_availability$age)
gender_check <- run_covariate_permanova("gender", covariate_availability$gender)
bmi_check    <- run_covariate_permanova("BMI", covariate_availability$bmi)

extract_permanova_row <- function(res, model_name, n) {
  df <- as.data.frame(res)
  df$term <- rownames(df)
  df$model <- model_name
  df$n_samples <- n
  df
}

permanova_summary <- bind_rows(
  extract_permanova_row(permanova_basic, "base_cohort_group", nrow(meta_for_permanova)),
  extract_permanova_row(age_check$result, "plus_age", age_check$n),
  extract_permanova_row(gender_check$result, "plus_gender", gender_check$n),
  extract_permanova_row(bmi_check$result, "plus_BMI", bmi_check$n)
) %>%
  filter(term != "Residual", term != "Total") %>%
  select(model, n_samples, term, R2, `Pr(>F)`)

print(permanova_summary)
write.csv(permanova_summary, here(cfg$paths$tables_dir, "step6_permanova_summary.csv"), row.names = FALSE)

# ============================================================
# STEP 7: Alpha diversity (on FULL raw abundance, pre-filter)
# ============================================================

shannon_idx <- diversity(abund_matrix, index = "shannon")
simpson_idx <- diversity(abund_matrix, index = "simpson")
observed_richness <- rowSums(abund_matrix > 0)

diversity_df <- pooled_data %>%
  select(sample_id, cohort, group, disease, age, gender, BMI) %>%
  mutate(
    shannon  = shannon_idx[match(sample_id, rownames(abund_matrix))],
    simpson  = simpson_idx[match(sample_id, rownames(abund_matrix))],
    richness = observed_richness[match(sample_id, rownames(abund_matrix))]
  )

diversity_group_summary <- diversity_df %>%
  group_by(group) %>%
  summarise(n = n(), mean_shannon = mean(shannon, na.rm=TRUE), mean_simpson = mean(simpson, na.rm=TRUE),
            mean_richness = mean(richness, na.rm=TRUE), .groups = "drop")
print(diversity_group_summary)
write.csv(diversity_group_summary, here(cfg$paths$tables_dir, "step7_diversity_group_summary.csv"), row.names = FALSE)

wilcox_shannon  <- wilcox.test(shannon ~ group, data = diversity_df)
wilcox_simpson  <- wilcox.test(simpson ~ group, data = diversity_df)
wilcox_richness <- wilcox.test(richness ~ group, data = diversity_df)

lme_shannon  <- lme(shannon  ~ group, random = ~1|cohort, data = diversity_df)
lme_simpson  <- lme(simpson  ~ group, random = ~1|cohort, data = diversity_df)
lme_richness <- lme(richness ~ group, random = ~1|cohort, data = diversity_df)

extract_lme_group_pval <- function(lme_fit) anova(lme_fit)["group", "p-value"]

diversity_test_summary <- data.frame(
  metric = c("shannon", "simpson", "richness"),
  wilcoxon_p_unadjusted = c(wilcox_shannon$p.value, wilcox_simpson$p.value, wilcox_richness$p.value),
  lme_p_cohort_adjusted = c(extract_lme_group_pval(lme_shannon), extract_lme_group_pval(lme_simpson), extract_lme_group_pval(lme_richness))
)
print(diversity_test_summary)
write.csv(diversity_test_summary, here(cfg$paths$tables_dir, "step7_diversity_test_summary.csv"), row.names = FALSE)

# ---- Figure: written to disk, headless-safe, no display needed ----
n_control <- sum(diversity_df$group == "control")
n_ibd     <- sum(diversity_df$group == "IBD")

make_diversity_panel <- function(metric_col, y_label, title) {
  y_vals <- diversity_df[[metric_col]]
  label_y <- max(y_vals, na.rm = TRUE) * 1.08
  ggplot(diversity_df, aes(x = group, y = .data[[metric_col]], fill = group)) +
    geom_boxplot(outlier.size = 0.8, alpha = 0.85, width = 0.6) +
    scale_fill_manual(values = group_palette) +
    stat_compare_means(method = "wilcox.test", label = "p.format", label.y = label_y) +
    expand_limits(y = label_y * 1.05) +
    labs(title = title, x = "Group", y = y_label) +
    theme_pipeline()
}

p_richness <- make_diversity_panel("richness", "Observed richness (# species)", "Richness")
p_shannon  <- make_diversity_panel("shannon",  "Shannon diversity index",       "Shannon")
p_simpson  <- make_diversity_panel("simpson",  "Simpson diversity index (1\u2212D)", "Simpson")

combined_diversity_plot <- p_richness + p_shannon + p_simpson +
  plot_annotation(title = "Alpha diversity: IBD vs. control",
                  subtitle = paste0("Wilcoxon rank-sum test | n = ", n_control, " control, ", n_ibd, " IBD"))

save_figure(combined_diversity_plot, "step7_diversity_combined", width = 12, height = 5)

saveRDS(diversity_df, here(cfg$paths$checkpoints_dir, "diversity_df.rds"))
cat("Checkpoint written: diversity_df.rds\n")
cat("Figures written: results/figures/step7_diversity_combined.{pdf,png}\n")