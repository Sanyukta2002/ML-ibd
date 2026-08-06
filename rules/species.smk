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

rule species_diversity_permanova:
    input:
        checkpoint = f"{config['paths']['checkpoints_dir']}/pooled_data.rds",
        covariate_availability = f"{config['paths']['checkpoints_dir']}/covariate_availability.rds",
        config = "config/config.yaml"
    output:
        checkpoint = f"{config['paths']['checkpoints_dir']}/diversity_df.rds",
        permanova_summary = f"{config['paths']['tables_dir']}/step6_permanova_summary.csv",
        diversity_group_summary = f"{config['paths']['tables_dir']}/step7_diversity_group_summary.csv",
        diversity_test_summary = f"{config['paths']['tables_dir']}/step7_diversity_test_summary.csv",
        figure_pdf = f"{config['paths']['figures_dir']}/step7_diversity_combined.pdf",
        figure_png = f"{config['paths']['figures_dir']}/step7_diversity_combined.png"
    log:
        "logs/species_diversity_permanova.log"
    container:
        "envs/bioconductor.sif"
    resources:
        mem_mb = 4000,
        runtime = 20
    shell:
        "Rscript scripts/r/04_diversity_permanova.R > {log} 2>&1"

rule covariate_availability_check:
    input:
        checkpoint = f"{config['paths']['checkpoints_dir']}/pooled_data.rds",
        config = "config/config.yaml"
    output:
        checkpoint = f"{config['paths']['checkpoints_dir']}/covariate_availability.rds",
        table = f"{config['paths']['tables_dir']}/step_covariate_availability_by_cohort.csv"
    log:
        "logs/covariate_availability_check.log"
    container:
        "envs/bioconductor.sif"
    resources:
        mem_mb = 2000,
        runtime = 10
    shell:
        "Rscript scripts/r/05_covariate_availability_check.R > {log} 2>&1"

rule species_maaslin2_tier1:
    input:
        checkpoint = f"{config['paths']['checkpoints_dir']}/pooled_data_filtered.rds",
        final_species = f"{config['paths']['checkpoints_dir']}/final_species.rds",
        config = "config/config.yaml"
    output:
        checkpoint = f"{config['paths']['checkpoints_dir']}/fit_tier1.rds",
        hit_summary = f"{config['paths']['tables_dir']}/step8_maaslin2_tier1_hit_summary.csv"
    log:
        "logs/species_maaslin2_tier1.log"
    container:
        "envs/bioconductor.sif"
    resources:
        mem_mb = 4000,
        runtime = 20
    shell:
        "Rscript scripts/r/06_maaslin2_tier1.R > {log} 2>&1"

rule species_maaslin2_tier2_age:
    input:
        checkpoint = f"{config['paths']['checkpoints_dir']}/pooled_data_filtered.rds",
        final_species = f"{config['paths']['checkpoints_dir']}/final_species.rds",
        covariate_availability = f"{config['paths']['checkpoints_dir']}/covariate_availability.rds",
        config = "config/config.yaml"
    output:
        checkpoint = f"{config['paths']['checkpoints_dir']}/fit_tier2.rds",
        group_results = f"{config['paths']['checkpoints_dir']}/tier2_group_results.rds",
        hit_summary = f"{config['paths']['tables_dir']}/step9_maaslin2_tier2_hit_summary.csv"
    log:
        "logs/species_maaslin2_tier2_age.log"
    container:
        "envs/bioconductor.sif"
    resources:
        mem_mb = 4000,
        runtime = 20
    shell:
        "Rscript scripts/r/07_maaslin2_tier2_age.R > {log} 2>&1"

rule species_maaslin2_tier_overlap:
    input:
        fit_tier1 = f"{config['paths']['checkpoints_dir']}/fit_tier1.rds",
        tier2_group_results = f"{config['paths']['checkpoints_dir']}/tier2_group_results.rds",
        config = "config/config.yaml"
    output:
        checkpoint = f"{config['paths']['checkpoints_dir']}/robust_core_tier1_tier2.rds",
        comparison = f"{config['paths']['tables_dir']}/step10_tier1_tier2_comparison.csv",
        robust_core = f"{config['paths']['tables_dir']}/step10_robust_core_tier1_tier2.csv",
        overlap_summary = f"{config['paths']['tables_dir']}/step10_tier_overlap_summary.csv"
    log:
        "logs/species_maaslin2_tier_overlap.log"
    container:
        "envs/bioconductor.sif"
    resources:
        mem_mb = 2000,
        runtime = 10
    shell:
        "Rscript scripts/r/08_maaslin2_tier_overlap.R > {log} 2>&1"

rule species_wilcoxon:
    input:
        checkpoint = f"{config['paths']['checkpoints_dir']}/pooled_data_filtered.rds",
        final_species = f"{config['paths']['checkpoints_dir']}/final_species.rds",
        fit_tier1 = f"{config['paths']['checkpoints_dir']}/fit_tier1.rds",
        config = "config/config.yaml"
    output:
        checkpoint = f"{config['paths']['checkpoints_dir']}/blocked_wilcox_df.rds",
        name_map = f"{config['paths']['checkpoints_dir']}/name_map2.rds",
        hit_summary = f"{config['paths']['tables_dir']}/step12_wilcoxon_hit_summary.csv",
        full_results = f"{config['paths']['tables_dir']}/step12_wilcoxon_full_results.csv"
    log:
        "logs/species_wilcoxon.log"
    container:
        "envs/bioconductor.sif"
    resources:
        mem_mb = 4000,
        runtime = 20
    shell:
        "Rscript scripts/r/09_species_wilcoxon.R > {log} 2>&1"

rule species_triple_candidates:
    input:
        checkpoint = f"{config['paths']['checkpoints_dir']}/pooled_data_filtered.rds",
        fit_tier1 = f"{config['paths']['checkpoints_dir']}/fit_tier1.rds",
        tier2_group_results = f"{config['paths']['checkpoints_dir']}/tier2_group_results.rds",
        blocked_wilcox_df = f"{config['paths']['checkpoints_dir']}/blocked_wilcox_df.rds",
        name_map = f"{config['paths']['checkpoints_dir']}/name_map2.rds",
        config = "config/config.yaml"
    output:
        checkpoint = f"{config['paths']['checkpoints_dir']}/triple_candidates.rds",
        funnel = f"{config['paths']['tables_dir']}/step13_selection_funnel.csv",
        candidates = f"{config['paths']['tables_dir']}/step13_triple_candidates.csv"
    log:
        "logs/species_triple_candidates.log"
    container:
        "envs/bioconductor.sif"
    resources:
        mem_mb = 2000,
        runtime = 10
    shell:
        "Rscript scripts/r/10_species_triple_candidates.R > {log} 2>&1"

rule species_rfe:
    input:
        checkpoint = f"{config['paths']['checkpoints_dir']}/triple_candidates.rds",
        pooled = f"{config['paths']['checkpoints_dir']}/pooled_data_filtered.rds",
        config = "config/config.yaml"
    output:
        checkpoint = f"{config['paths']['checkpoints_dir']}/all_panels.rds",
        rfe_result = f"{config['paths']['checkpoints_dir']}/rfe_result.rds",
        results_by_size = f"{config['paths']['tables_dir']}/step14_rfe_results_by_size.csv",
        elbow_sizes = f"{config['paths']['tables_dir']}/step14_rfe_elbow_sizes.csv",
        panel_summary = f"{config['paths']['tables_dir']}/step14_final_panels.csv",
        figure_pdf = f"{config['paths']['figures_dir']}/step14_rfe_auc_vs_panelsize.pdf",
        figure_png = f"{config['paths']['figures_dir']}/step14_rfe_auc_vs_panelsize.png"
    log:
        "logs/species_rfe.log"
    container:
        "envs/bioconductor.sif"
    resources:
        mem_mb = 8000,
        runtime = 60
    shell:
        "Rscript scripts/r/11_species_rfe.R > {log} 2>&1"
