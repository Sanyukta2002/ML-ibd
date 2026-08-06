configfile: "config/config.yaml"

include: "rules/species.smk"
include: "rules/pathway.smk"

rule all:
    input:
        f"{config['paths']['checkpoints_dir']}/diversity_df.rds"
