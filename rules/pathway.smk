rule pathway_abundance_pull:
    input:
        checkpoint = f"{config['paths']['checkpoints_dir']}/baseline_metadata.rds",
        config = "config/config.yaml"
    output:
        checkpoint = f"{config['paths']['checkpoints_dir']}/pooled_data_pathway.rds",
        counts = f"{config['paths']['tables_dir']}/step_pathway_pooled_counts.csv"
    log: "logs/pathway_abundance_pull.log"
    container: "envs/bioconductor.sif"
    resources: mem_mb = 4000, runtime = 15
    shell: "Rscript scripts/r/12_pathway_abundance_pull.R > {log} 2>&1"

rule pathway_filter_dedup:
    input:
        checkpoint = f"{config['paths']['checkpoints_dir']}/pooled_data_pathway.rds",
        config = "config/config.yaml"
    output:
        checkpoint = f"{config['paths']['checkpoints_dir']}/pooled_data_pathway_final.rds",
        final_pathways = f"{config['paths']['checkpoints_dir']}/final_pathways.rds",
        funnel = f"{config['paths']['tables_dir']}/step_pathway_filtering_funnel.csv",
        dedup_summary = f"{config['paths']['tables_dir']}/step_pathway_dedup_summary.csv"
    log: "logs/pathway_filter_dedup.log"
    container: "envs/bioconductor.sif"
    resources: mem_mb = 4000, runtime = 15
    shell: "Rscript scripts/r/13_pathway_filter_dedup.R > {log} 2>&1"
