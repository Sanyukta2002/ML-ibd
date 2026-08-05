rule cohort_selection:
    input:
        "config/config.yaml"
    output:
        checkpoint = f"{config['paths']['checkpoints_dir']}/baseline_metadata.rds",
        raw_counts = f"{config['paths']['tables_dir']}/step2_cohort_raw_counts.csv",
        timepoint_avail = f"{config['paths']['tables_dir']}/step3_timepoint_availability.csv",
        baseline_counts = f"{config['paths']['tables_dir']}/step3_baseline_counts_by_group.csv"
    log:
        "logs/cohort_selection.log"
    container:
        "envs/bioconductor.sif"
    shell:
        "Rscript scripts/r/01_cohort_selection.R > {log} 2>&1"

rule species_abundance_pull:
    input:
        checkpoint = f"{config['paths']['checkpoints_dir']}/baseline_metadata.rds",
        config = "config/config.yaml"
    output:
        checkpoint = f"{config['paths']['checkpoints_dir']}/pooled_data.rds",
        pooled_counts = f"{config['paths']['tables_dir']}/step4_pooled_counts_by_group.csv",
        raw_taxa_counts = f"{config['paths']['tables_dir']}/step4_raw_taxa_counts_per_cohort.csv",
        taxa_provenance = f"{config['paths']['tables_dir']}/step4_taxa_count_provenance.csv"
    log:
        "logs/species_abundance_pull.log"
    container:
        "envs/bioconductor.sif"
    resources:
        mem_mb = 4000,
        runtime = 15
    shell:
        "Rscript scripts/r/02_species_abundance_pull.R > {log} 2>&1"

rule species_filter_dedup:
    input:
        checkpoint = f"{config['paths']['checkpoints_dir']}/pooled_data.rds",
        config = "config/config.yaml"
    output:
        checkpoint = f"{config['paths']['checkpoints_dir']}/pooled_data_filtered.rds",
        final_species = f"{config['paths']['checkpoints_dir']}/final_species.rds",
        funnel = f"{config['paths']['tables_dir']}/step5_species_filtering_funnel.csv",
        filter_detail = f"{config['paths']['tables_dir']}/step5_species_filter_detail.csv",
        dedup_investigation = f"{config['paths']['tables_dir']}/step5b_species_duplicate_investigation.csv",
        dedup_summary = f"{config['paths']['tables_dir']}/step5b_species_dedup_summary.csv"
    log:
        "logs/species_filter_dedup.log"
    container:
        "envs/bioconductor.sif"
    resources:
        mem_mb = 4000,
        runtime = 15
    shell:
        "Rscript scripts/r/03_species_filter_dedup.R > {log} 2>&1"
