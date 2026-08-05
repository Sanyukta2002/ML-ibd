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
