# 21_export_for_python.R
# Exports FULL filtered tables (not pre-cut panels) -- CLR must be
# computed on the full filtered composition, then subset to panels
# in Python, per sub-compositional coherence (CLR is not sub-
# compositionally coherent; computing it on a pre-cut panel would
# use the wrong geometric-mean reference). Panel membership exported
# separately so Python knows which columns belong to which panel.

library(here); library(yaml); library(dplyr)

cfg <- yaml::read_yaml(here("config", "config.yaml"))
export_dir <- here(cfg$paths$python_export_dir)
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Species: full filtered+deduped table (151 species) ----
pooled_species <- readRDS(here(cfg$paths$checkpoints_dir, "pooled_data_filtered.rds"))
final_species <- readRDS(here(cfg$paths$checkpoints_dir, "final_species.rds"))
stopifnot(all(final_species %in% colnames(pooled_species)))

species_export <- pooled_species[, c("sample_id", "cohort", "group", final_species)]
write.csv(species_export, file.path(export_dir, "species_full_filtered.csv"), row.names = FALSE)
cat("Species export:", nrow(species_export), "samples x", length(final_species), "features\n")

# ---- Pathway: full filtered+deduped table (118 pathways) ----
pooled_pathway <- readRDS(here(cfg$paths$checkpoints_dir, "pooled_data_pathway_final.rds"))
final_pathways <- readRDS(here(cfg$paths$checkpoints_dir, "final_pathways.rds"))
stopifnot(all(final_pathways %in% colnames(pooled_pathway)))

pathway_export <- pooled_pathway[, c("sample_id", "cohort", "group", final_pathways)]
write.csv(pathway_export, file.path(export_dir, "pathway_full_filtered.csv"), row.names = FALSE)
cat("Pathway export:", nrow(pathway_export), "samples x", length(final_pathways), "features\n")

# ---- Panel membership (long format: panel_name, feature) ----
all_panels_species <- readRDS(here(cfg$paths$checkpoints_dir, "all_panels.rds"))
species_membership <- bind_rows(lapply(names(all_panels_species), function(nm) {
  data.frame(panel_name = nm, feature = all_panels_species[[nm]], stringsAsFactors = FALSE)
}))
write.csv(species_membership, file.path(export_dir, "species_panel_membership.csv"), row.names = FALSE)

all_panels_pathway <- readRDS(here(cfg$paths$checkpoints_dir, "all_panels_pathway.rds"))
pathway_membership <- bind_rows(lapply(names(all_panels_pathway), function(nm) {
  data.frame(panel_name = nm, feature = all_panels_pathway[[nm]], stringsAsFactors = FALSE)
}))
write.csv(pathway_membership, file.path(export_dir, "pathway_panel_membership.csv"), row.names = FALSE)

# ---- Sanity checks ----
for (nm in names(all_panels_species)) {
  stopifnot(all(all_panels_species[[nm]] %in% final_species))
}
for (nm in names(all_panels_pathway)) {
  stopifnot(all(all_panels_pathway[[nm]] %in% final_pathways))
}
cat("All panel features confirmed present in their respective full filtered tables.\n")
cat("Export complete:", export_dir, "\n")
