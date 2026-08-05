# ============================================================
# IBD cross-cohort microbiome classification pipeline
# Orchestrates R (cohort curation -> stats -> panel selection)
# and Python (nested-CV model training) stages.
# ============================================================

configfile: "config/config.yaml"

# Rule modules are added incrementally, one at a time, as each
# stage is built and verified (see rules/ — one .smk file per
# logical group). Nothing is wired into `rule all` until its
# upstream rule has been run and sanity-checked.

rule all:
    input:
        []
