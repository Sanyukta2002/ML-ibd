configfile: "config/config.yaml"

include: "rules/species.smk"
include: "rules/pathway.smk"
include: "rules/mlflow_test.smk"
include: "rules/ml.smk"

rule all:
    input:
        f"{config['paths']['checkpoints_dir']}/diversity_df.rds"
