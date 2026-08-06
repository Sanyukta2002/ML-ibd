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

rule pathway_diversity_permanova:
    input:
        pooled = f"{config['paths']['checkpoints_dir']}/pooled_data_pathway.rds",
        covariate_availability = f"{config['paths']['checkpoints_dir']}/covariate_availability.rds",
        config = "config/config.yaml"
    output:
        checkpoint = f"{config['paths']['checkpoints_dir']}/diversity_df_pathway.rds",
        permanova = f"{config['paths']['tables_dir']}/step_pathway_permanova_summary.csv",
        diversity_test = f"{config['paths']['tables_dir']}/step_pathway_diversity_test_summary.csv",
        figure_pdf = f"{config['paths']['figures_dir']}/step_pathway_diversity_combined.pdf",
        figure_png = f"{config['paths']['figures_dir']}/step_pathway_diversity_combined.png"
    log: "logs/pathway_diversity_permanova.log"
    container: "envs/bioconductor.sif"
    resources: mem_mb = 4000, runtime = 20
    shell: "Rscript scripts/r/14_pathway_diversity_permanova.R > {log} 2>&1"

rule pathway_maaslin2_tier1:
    input:
        pooled = f"{config['paths']['checkpoints_dir']}/pooled_data_pathway_final.rds",
        final_pathways = f"{config['paths']['checkpoints_dir']}/final_pathways.rds",
        config = "config/config.yaml"
    output:
        checkpoint = f"{config['paths']['checkpoints_dir']}/fit_tier1_pathway.rds",
        hit_summary = f"{config['paths']['tables_dir']}/step_pathway_maaslin2_tier1_hit_summary.csv"
    log: "logs/pathway_maaslin2_tier1.log"
    container: "envs/bioconductor.sif"
    resources: mem_mb = 4000, runtime = 20
    shell: "Rscript scripts/r/15_pathway_maaslin2_tier1.R > {log} 2>&1"

rule pathway_maaslin2_tier2_age:
    input:
        pooled = f"{config['paths']['checkpoints_dir']}/pooled_data_pathway_final.rds",
        final_pathways = f"{config['paths']['checkpoints_dir']}/final_pathways.rds",
        covariate_availability = f"{config['paths']['checkpoints_dir']}/covariate_availability.rds",
        config = "config/config.yaml"
    output:
        checkpoint = f"{config['paths']['checkpoints_dir']}/fit_tier2_pathway.rds",
        group_results = f"{config['paths']['checkpoints_dir']}/tier2_group_results_pathway.rds",
        hit_summary = f"{config['paths']['tables_dir']}/step_pathway_maaslin2_tier2_hit_summary.csv"
    log: "logs/pathway_maaslin2_tier2_age.log"
    container: "envs/bioconductor.sif"
    resources: mem_mb = 4000, runtime = 20
    shell: "Rscript scripts/r/16_pathway_maaslin2_tier2_age.R > {log} 2>&1"

rule pathway_maaslin2_tier_overlap:
    input:
        fit_tier1 = f"{config['paths']['checkpoints_dir']}/fit_tier1_pathway.rds",
        tier2_group_results = f"{config['paths']['checkpoints_dir']}/tier2_group_results_pathway.rds",
        config = "config/config.yaml"
    output:
        summary = f"{config['paths']['tables_dir']}/step_pathway_tier_overlap_summary.csv",
        comparison = f"{config['paths']['tables_dir']}/step_pathway_tier1_tier2_comparison.csv",
        robust_core = f"{config['paths']['tables_dir']}/step_pathway_robust_core_tier1_tier2.csv"
    log: "logs/pathway_maaslin2_tier_overlap.log"
    container: "envs/bioconductor.sif"
    resources: mem_mb = 2000, runtime = 10
    shell: "Rscript scripts/r/17_pathway_maaslin2_tier_overlap.R > {log} 2>&1"

rule pathway_wilcoxon:
    input:
        pooled = f"{config['paths']['checkpoints_dir']}/pooled_data_pathway_final.rds",
        final_pathways = f"{config['paths']['checkpoints_dir']}/final_pathways.rds",
        fit_tier1 = f"{config['paths']['checkpoints_dir']}/fit_tier1_pathway.rds",
        config = "config/config.yaml"
    output:
        checkpoint = f"{config['paths']['checkpoints_dir']}/blocked_wilcox_df_pathway.rds",
        name_map = f"{config['paths']['checkpoints_dir']}/name_map_pathway.rds",
        hit_summary = f"{config['paths']['tables_dir']}/step_pathway_wilcoxon_hit_summary.csv",
        full_results = f"{config['paths']['tables_dir']}/step_pathway_wilcoxon_full_results.csv"
    log: "logs/pathway_wilcoxon.log"
    container: "envs/bioconductor.sif"
    resources: mem_mb = 4000, runtime = 20
    shell: "Rscript scripts/r/18_pathway_wilcoxon.R > {log} 2>&1"

rule pathway_triple_candidates:
    input:
        pooled = f"{config['paths']['checkpoints_dir']}/pooled_data_pathway_final.rds",
        fit_tier1 = f"{config['paths']['checkpoints_dir']}/fit_tier1_pathway.rds",
        blocked_wilcox_df = f"{config['paths']['checkpoints_dir']}/blocked_wilcox_df_pathway.rds",
        name_map = f"{config['paths']['checkpoints_dir']}/name_map_pathway.rds",
        config = "config/config.yaml"
    output:
        checkpoint = f"{config['paths']['checkpoints_dir']}/triple_candidates_pathway.rds",
        candidates = f"{config['paths']['tables_dir']}/step_pathway_triple_candidates.csv"
    log: "logs/pathway_triple_candidates.log"
    container: "envs/bioconductor.sif"
    resources: mem_mb = 2000, runtime = 10
    shell: "Rscript scripts/r/19_pathway_triple_candidates.R > {log} 2>&1"
