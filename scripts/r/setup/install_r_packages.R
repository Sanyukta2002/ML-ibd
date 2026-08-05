# ============================================================
# install_r_packages.R
# Installs the full R/Bioconductor package set for this pipeline
# into an external library (R_LIBS_USER, expected to be bind-mounted
# from outside the container — see envs/r_bioconductor.def).
# Idempotent: safe to re-run, BiocManager/install.packages skip
# already-installed packages by default.
#
# Known issue (documented, not a bug in this script): the CRAN
# package 'Deriv' fails to compile against this image's R 4.4.1
# (uses R_ClosureFormals, an API function added after 4.4.1 — see
# R's NEWS). This is a transitive, non-required dependency pulled
# in by the Bioconductor install below; all 4 target Bioconductor
# packages install and load successfully despite it. Left
# unresolved deliberately rather than pinning an old Deriv version
# we don't otherwise need.
# ============================================================

lib_path <- Sys.getenv("R_LIBS_USER")
if (lib_path == "") stop("R_LIBS_USER is not set — this script expects an external, bind-mounted library path.")

cat("Installing into:", lib_path, "\n")

# ---- Bioconductor packages ----
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = lib_path, repos = "https://cloud.r-project.org")
}

BiocManager::install(
  c("curatedMetagenomicData", "mia", "TreeSummarizedExperiment", "Maaslin2"),
  lib = lib_path,
  update = FALSE,
  ask = FALSE
)

# ---- CRAN packages ----
install.packages(
  c("vegan", "coin", "caret", "randomForest", "ggplot2",
    "ggpubr", "patchwork", "nlme", "igraph", "here", "yaml"),
  lib = lib_path,
  repos = "https://cloud.r-project.org"
)

# ---- Verification ----
required_pkgs <- c(
  "curatedMetagenomicData", "mia", "TreeSummarizedExperiment", "Maaslin2",
  "vegan", "coin", "caret", "randomForest", "ggplot2",
  "ggpubr", "patchwork", "nlme", "igraph", "here", "yaml"
)

results <- sapply(required_pkgs, function(p) {
  tryCatch({
    suppressPackageStartupMessages(library(p, character.only = TRUE, lib.loc = lib_path))
    TRUE
  }, error = function(e) {
    message(p, ": FAILED — ", conditionMessage(e))
    FALSE
  })
})

cat("\n---- Install verification ----\n")
print(results)

if (!all(results)) {
  stop("One or more required packages failed to install/load. See messages above.")
} else {
  cat("\nAll required packages installed and loadable.\n")
}
