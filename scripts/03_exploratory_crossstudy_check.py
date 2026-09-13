"""
Exploratory cross-study comparison against ChickpeaOmicsR (Alsamman et al., 2026).

Reproduces the descriptive comparison in manuscript Section 6. This is NOT a
validation of the genotype-dependent (Delta_log2FC) signal -- see the caveat
below and the manuscript text for why.

Usage:
    python 03_exploratory_crossstudy_check.py
"""
import os
import pandas as pd

RESULTS_DIR = os.path.join(os.path.dirname(__file__), "..", "results")
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")

ours = pd.read_csv(os.path.join(RESULTS_DIR, "derived_702_genotype_dependent_genes.csv"))
theirs = pd.read_csv(os.path.join(DATA_DIR, "chickpeaomicsr_query_result.csv"))

merged = pd.merge(theirs, ours, left_on="geneid", right_on="gene", how="inner")
with_drought = merged[merged["Drought"].notna()].copy()

print(f"Queried genes present in ChickpeaOmicsR: {len(merged)}")
print(f"Of these, with a non-NA aggregated Drought value: {len(with_drought)}")

# NOTE: ChickpeaOmicsR's "Drought" column is a single value aggregated across
# multiple source drought experiments. It is NOT genotype-resolved, and
# therefore CANNOT be used to test the genotype-dependent Delta_log2FC
# signal from the present study. The comparison below is purely descriptive.

with_drought["closer_to_283"] = (
    (with_drought["Drought"] - with_drought["log2FC_283"]).abs()
    < (with_drought["Drought"] - with_drought["log2FC_8261"]).abs()
)

cols = ["geneid", "log2FC_283", "log2FC_8261", "Drought", "closer_to_283", "direction_class"]
print(with_drought[cols].to_string(index=False))
print()
print(f"Closer to ICC283 (sensitive): {with_drought['closer_to_283'].sum()} / {len(with_drought)}")