configfile: "config/config.yaml"

include: "rules/species.smk"

rule all:
    input:
        f"{config['paths']['checkpoints_dir']}/pooled_data_filtered.rds"
