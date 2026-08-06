import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))

import pandas as pd
from cv_utils import clr_transform

REPO_ROOT = Path(__file__).resolve().parents[3]
export_dir = REPO_ROOT / "results" / "python_export"

species_df = pd.read_csv(export_dir / "species_full_filtered.csv")
print("Loaded species table:", species_df.shape)

meta_cols = ["sample_id", "cohort", "group"]
feature_cols = [c for c in species_df.columns if c not in meta_cols]
print("Feature columns:", len(feature_cols))

features_only = species_df[feature_cols]
print("Any NaN in raw features?", features_only.isna().any().any())
print("Min value:", features_only.values.min(), "| Max value:", features_only.values.max())

clr_df, pc_used = clr_transform(features_only)
print("\nCLR transform succeeded.")
print("Pseudocount used:", pc_used)
print("CLR output shape:", clr_df.shape)
print("Any NaN in CLR output?", clr_df.isna().any().any())
print("Row sums (should be ~0):")
print(clr_df.sum(axis=1).describe())

# quick sanity on a real species: Adlercreutzia equolifaciens (we know
# this one well from the R-side investigation)
if "Adlercreutzia equolifaciens" in clr_df.columns:
    print("\nAdlercreutzia equolifaciens -- raw vs CLR (first 5 samples):")
    print(pd.DataFrame({
        "raw": features_only["Adlercreutzia equolifaciens"].head(),
        "clr": clr_df["Adlercreutzia equolifaciens"].head()
    }))

print("\nSMOKE TEST PASSED")
