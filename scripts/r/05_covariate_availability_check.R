# ============================================================
# 05_covariate_availability_check.R
# Derives, from actual current data, which cohorts have usable
# (non-fully-missing) age/gender/BMI -- rather than trusting
# hardcoded cohort lists carried over from the original run.
# ============================================================

library(here)
library(yaml)
library(dplyr)

cfg <- yaml::read_yaml(here("config", "config.yaml"))
pooled_data <- readRDS(here(cfg$paths$checkpoints_dir, "pooled_data.rds"))

covariate_availability <- pooled_data %>%
  group_by(cohort) %>%
  summarise(
    n_total = n(),
    age_n_present = sum(!is.na(age)),    age_pct_present = round(100 * mean(!is.na(age)), 1),
    gender_n_present = sum(!is.na(gender)), gender_pct_present = round(100 * mean(!is.na(gender)), 1),
    BMI_n_present = sum(!is.na(BMI)),    BMI_pct_present = round(100 * mean(!is.na(BMI)), 1),
    .groups = "drop"
  )

print(covariate_availability)
write.csv(covariate_availability, here(cfg$paths$tables_dir, "step_covariate_availability_by_cohort.csv"), row.names = FALSE)

# ---- Derive usable-cohort lists (>0% present -- a cohort with true
# 100% missingness contributes nothing and should not be listed as
# "available", even if it was previously) ----
min_usable_pct <- cfg$covariate_availability$min_usable_pct

derive_usable_cohorts <- function(pct_col) {
  # threshold reverse-engineered from the original script's undocumented
  # hardcoded cohort choices -- see config comment for justification
  covariate_availability$cohort[covariate_availability[[pct_col]] >= min_usable_pct]
}

age_cohorts_actual    <- derive_usable_cohorts("age_pct_present")
gender_cohorts_actual <- derive_usable_cohorts("gender_pct_present")
bmi_cohorts_actual    <- derive_usable_cohorts("BMI_pct_present")

cat("\nAge-available cohorts (derived from data):", paste(age_cohorts_actual, collapse=", "), "\n")
cat("Gender-available cohorts (derived from data):", paste(gender_cohorts_actual, collapse=", "), "\n")
cat("BMI-available cohorts (derived from data):", paste(bmi_cohorts_actual, collapse=", "), "\n")

# ---- Compare against config's hardcoded assumptions, flag drift ----
config_age_cohorts <- cfg$maaslin2$tier2_age_cohorts
age_drift <- !setequal(age_cohorts_actual, config_age_cohorts)

if (age_drift) {
  cat("\nWARNING: config's tier2_age_cohorts", paste(config_age_cohorts, collapse=", "),
      "does NOT match data-derived list", paste(age_cohorts_actual, collapse=", "),
      "-- config should be updated, downstream rules should use the derived list.\n")
} else {
  cat("\nConfig's tier2_age_cohorts matches data-derived list. No drift.\n")
}

saveRDS(
  list(age = age_cohorts_actual, gender = gender_cohorts_actual, bmi = bmi_cohorts_actual),
  here(cfg$paths$checkpoints_dir, "covariate_availability.rds")
)
cat("Checkpoint written: covariate_availability.rds\n")
