"""
Derived genotype-dependent drought-response statistic.

Reproduces the analysis in Section 3.3-4.4 of the manuscript from the four
deposited processed files (GSE104609, BioProject PRJNA413294):
  - GSE104609_ICC283_DESeq2.csv   (drought-sensitive genotype)
  - GSE104609_ICC8261_DESeq2.csv  (drought-tolerant genotype)

Place the two DESeq2 CSV files in ../data/ (see data/README.md for source)
before running this script.

Usage:
    python 01_derived_genotype_response.py
"""
import os
import pandas as pd
import numpy as np
from scipy import stats

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "..", "results")

icc283 = pd.read_csv(os.path.join(DATA_DIR, "GSE104609_ICC283_DESeq2.csv"))
icc283 = icc283.rename(columns={icc283.columns[0]: "idx"}) if icc283.columns[0] == "" else icc283
if "Row.names" not in icc283.columns:
    icc283 = icc283.rename(columns={icc283.columns[1]: "Row.names"})

icc8261 = pd.read_csv(os.path.join(DATA_DIR, "GSE104609_ICC8261_DESeq2.csv"))

print("ICC283 columns:", list(icc283.columns))
print("ICC8261 columns:", list(icc8261.columns))
print("ICC283 rows:", len(icc283))
print("ICC8261 rows:", len(icc8261))

icc283_sub = icc283[["Row.names", "log2FoldChange", "lfcSE", "pvalue", "padj"]].copy()
icc283_sub.columns = ["gene", "log2FC_283", "lfcSE_283", "pvalue_283", "padj_283"]

icc8261_sub = icc8261[["Row.names", "log2FoldChange", "lfcSE", "pvalue", "padj"]].copy()
icc8261_sub.columns = ["gene", "log2FC_8261", "lfcSE_8261", "pvalue_8261", "padj_8261"]

print("\nICC283 unique genes:", icc283_sub["gene"].nunique())
print("ICC8261 unique genes:", icc8261_sub["gene"].nunique())

def within_genotype_summary(df, fc_col, p_col, label):
    sig_p = df[df[p_col] < 0.05]
    sig_both = df[(df[p_col] < 0.05) & (df[fc_col].abs() >= 1)]
    up = sig_both[sig_both[fc_col] > 0]
    down = sig_both[sig_both[fc_col] < 0]
    print(f"\n--- {label} ---")
    print(f"Adjusted P<0.05: {len(sig_p)}")
    print(f"Adjusted P<0.05 & |log2FC|>=1: {len(sig_both)}")
    print(f"Up-regulated: {len(up)}")
    print(f"Down-regulated: {len(down)}")
    return sig_both

sig283 = within_genotype_summary(icc283_sub, "log2FC_283", "padj_283", "ICC283")
sig8261 = within_genotype_summary(icc8261_sub, "log2FC_8261", "padj_8261", "ICC8261")

shared = pd.merge(sig283[["gene", "log2FC_283"]], sig8261[["gene", "log2FC_8261"]], on="gene", how="inner")
concordant = shared[np.sign(shared["log2FC_283"]) == np.sign(shared["log2FC_8261"])]
print(f"\n--- Shared (4.3) ---")
print(f"Shared genes (both sig, |FC|>=1): {len(shared)}")
print(f"Concordant direction: {len(concordant)} ({100*len(concordant)/len(shared):.1f}%)")

merged_all = pd.merge(icc283_sub, icc8261_sub, on="gene", how="inner")
print(f"\nMerged genes (present in both files): {len(merged_all)}")
print(f"NaN log2FC_283 (DESeq2 independent filtering / outliers): {merged_all['log2FC_283'].isna().sum()}")
print(f"NaN log2FC_8261: {merged_all['log2FC_8261'].isna().sum()}")

merged = merged_all.dropna(subset=["log2FC_283", "lfcSE_283", "log2FC_8261", "lfcSE_8261"]).copy()
print(f"Genes with valid estimates in both genotypes (testable): {len(merged)}")

merged["delta_log2FC"] = merged["log2FC_8261"] - merged["log2FC_283"]
merged["se_delta"] = np.sqrt(merged["lfcSE_283"]**2 + merged["lfcSE_8261"]**2)
merged["z"] = merged["delta_log2FC"] / merged["se_delta"]
merged["p_delta"] = 2 * (1 - stats.norm.cdf(np.abs(merged["z"])))

def bh_fdr(pvals):
    pvals = np.asarray(pvals)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    fdr = ranked * n / (np.arange(n) + 1)
    fdr = np.minimum.accumulate(fdr[::-1])[::-1]
    fdr = np.clip(fdr, 0, 1)
    out = np.empty(n)
    out[order] = fdr
    return out

merged["fdr_delta"] = bh_fdr(merged["p_delta"].values)

sig_fdr = merged[merged["fdr_delta"] < 0.05]
print(f"\nFDR<0.05 (before effect-size filter): {len(sig_fdr)}")

sig_final = merged[(merged["fdr_delta"] < 0.05) & (merged["delta_log2FC"].abs() >= 1)].copy()
print(f"FDR<0.05 & |delta_log2FC|>=1 (FINAL 702-gene set): {len(sig_final)}")

def classify(row):
    a, b = row["log2FC_283"], row["log2FC_8261"]
    if a < 0 and b > 0:
        return "sensitive-down_tolerant-up"
    elif a > 0 and b < 0:
        return "sensitive-up_tolerant-down"
    elif a > 0 and b > 0:
        return "both-up_diff-magnitude"
    elif a < 0 and b < 0:
        return "both-down_diff-magnitude"
    else:
        return "other"

sig_final["direction_class"] = sig_final.apply(classify, axis=1)
print("\nDirectional class counts:")
print(sig_final["direction_class"].value_counts())

out_cols = ["gene", "log2FC_283", "lfcSE_283", "padj_283",
            "log2FC_8261", "lfcSE_8261", "padj_8261",
            "delta_log2FC", "se_delta", "z", "p_delta", "fdr_delta", "direction_class"]
sig_final_out = sig_final[out_cols].sort_values("fdr_delta")
out_path = os.path.join(RESULTS_DIR, "derived_702_genotype_dependent_genes.csv")
sig_final_out.to_csv(out_path, index=False)
print(f"\nSaved: {out_path}")

print("\nTop 20 genes by FDR:")
print(sig_final_out.head(20)[["gene", "log2FC_283", "log2FC_8261", "delta_log2FC", "fdr_delta", "direction_class"]].to_string(index=False))