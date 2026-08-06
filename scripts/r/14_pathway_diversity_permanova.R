# 14_pathway_diversity_permanova.R
library(here); library(yaml); library(dplyr); library(vegan); library(nlme)
library(ggpubr); library(patchwork)
source(here("scripts", "r", "utils", "plot_theme.R"))

cfg <- yaml::read_yaml(here("config", "config.yaml"))
pooled_data_pathway <- readRDS(here(cfg$paths$checkpoints_dir, "pooled_data_pathway.rds"))
covariate_availability <- readRDS(here(cfg$paths$checkpoints_dir, "covariate_availability.rds"))

meta_cols <- c("sample_id","cohort","group","disease","age","gender","BMI","country","subject_id")
pathway_cols <- setdiff(colnames(pooled_data_pathway), meta_cols)
abund_matrix <- as.matrix(pooled_data_pathway[, pathway_cols])
rownames(abund_matrix) <- pooled_data_pathway$sample_id

# ---- PERMANOVA ----
bray_dist <- vegdist(abund_matrix, method = "bray")
meta_for_permanova <- pooled_data_pathway[, meta_cols]
set.seed(cfg$rfe$seed)
permanova_basic <- adonis2(bray_dist ~ cohort + group, data = meta_for_permanova, permutations = 999, by = "margin")
print(permanova_basic)

run_covariate_permanova <- function(covariate, cohorts) {
  d <- pooled_data_pathway %>% filter(cohort %in% cohorts, !is.na(.data[[covariate]]))
  abund_sub <- as.matrix(d[, pathway_cols]); rownames(abund_sub) <- d$sample_id
  bray_sub <- vegdist(abund_sub, method = "bray")
  set.seed(cfg$rfe$seed)
  list(n = nrow(d), result = adonis2(as.formula(paste("bray_sub ~ cohort + group +", covariate)),
                                       data = d[, meta_cols], permutations = 999, by = "margin"))
}

age_check    <- run_covariate_permanova("age", covariate_availability$age)
gender_check <- run_covariate_permanova("gender", covariate_availability$gender)
bmi_check    <- run_covariate_permanova("BMI", covariate_availability$bmi)

extract_row <- function(res, model_name, n) {
  df <- as.data.frame(res); df$term <- rownames(df); df$model <- model_name; df$n_samples <- n; df
}
permanova_summary <- bind_rows(
  extract_row(permanova_basic, "base_cohort_group", nrow(meta_for_permanova)),
  extract_row(age_check$result, "plus_age", age_check$n),
  extract_row(gender_check$result, "plus_gender", gender_check$n),
  extract_row(bmi_check$result, "plus_BMI", bmi_check$n)
) %>% filter(term != "Residual", term != "Total") %>% select(model, n_samples, term, R2, `Pr(>F)`)
print(permanova_summary)
write.csv(permanova_summary, here(cfg$paths$tables_dir, "step_pathway_permanova_summary.csv"), row.names = FALSE)

# ---- Alpha diversity ----
shannon_idx <- diversity(abund_matrix, index = "shannon")
simpson_idx <- diversity(abund_matrix, index = "simpson")
richness <- rowSums(abund_matrix > 0)

diversity_df <- pooled_data_pathway %>%
  select(sample_id, cohort, group, disease, age, gender, BMI) %>%
  mutate(shannon = shannon_idx[match(sample_id, rownames(abund_matrix))],
         simpson = simpson_idx[match(sample_id, rownames(abund_matrix))],
         richness = richness[match(sample_id, rownames(abund_matrix))])

wilcox_shannon <- wilcox.test(shannon ~ group, data = diversity_df)
wilcox_simpson <- wilcox.test(simpson ~ group, data = diversity_df)
wilcox_richness <- wilcox.test(richness ~ group, data = diversity_df)
lme_shannon <- lme(shannon ~ group, random = ~1|cohort, data = diversity_df)
lme_simpson <- lme(simpson ~ group, random = ~1|cohort, data = diversity_df)
lme_richness <- lme(richness ~ group, random = ~1|cohort, data = diversity_df)
extract_lme_p <- function(f) anova(f)["group", "p-value"]

diversity_test_summary <- data.frame(
  metric = c("shannon","simpson","richness"),
  wilcoxon_p_unadjusted = c(wilcox_shannon$p.value, wilcox_simpson$p.value, wilcox_richness$p.value),
  lme_p_cohort_adjusted = c(extract_lme_p(lme_shannon), extract_lme_p(lme_simpson), extract_lme_p(lme_richness))
)
print(diversity_test_summary)
write.csv(diversity_test_summary, here(cfg$paths$tables_dir, "step_pathway_diversity_test_summary.csv"), row.names = FALSE)

n_control <- sum(diversity_df$group == "control"); n_ibd <- sum(diversity_df$group == "IBD")
make_panel <- function(col, ylab, title) {
  y_vals <- diversity_df[[col]]; label_y <- max(y_vals, na.rm = TRUE) * 1.08
  ggplot(diversity_df, aes(x = group, y = .data[[col]], fill = group)) +
    geom_boxplot(outlier.size = 0.8, alpha = 0.85, width = 0.6) +
    scale_fill_manual(values = group_palette) +
    stat_compare_means(method = "wilcox.test", label = "p.format", label.y = label_y) +
    expand_limits(y = label_y * 1.05) + labs(title = title, x = "Group", y = ylab) + theme_pipeline()
}
combined_plot <- make_panel("richness", "Observed richness (# pathways)", "Richness") +
  make_panel("shannon", "Shannon diversity index", "Shannon") +
  make_panel("simpson", "Simpson diversity index (1-D)", "Simpson") +
  plot_annotation(title = "Pathway alpha diversity: IBD vs. control",
                   subtitle = paste0("Wilcoxon rank-sum test | n = ", n_control, " control, ", n_ibd, " IBD"))
save_figure(combined_plot, "step_pathway_diversity_combined", width = 12, height = 5)

saveRDS(diversity_df, here(cfg$paths$checkpoints_dir, "diversity_df_pathway.rds"))
cat("Checkpoint written: diversity_df_pathway.rds\n")
