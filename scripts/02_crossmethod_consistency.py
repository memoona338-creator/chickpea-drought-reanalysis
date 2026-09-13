"""
Cross-method (DESeq2 vs edgeR) consistency check.

Reproduces manuscript Section 4.6. Place all four deposited processed
files from GSE104609 in ../data/ before running (see data/README.md).

Usage:
    python 02_crossmethod_consistency.py
"""
import os
import pandas as pd
import numpy as np

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")

edger283 = pd.read_csv(os.path.join(DATA_DIR, "GSE104609_ICC283_edgeR.csv"))
edger283.columns = ["gene", "logFC", "logCPM", "PValue", "FDR"]

edger8261 = pd.read_csv(os.path.join(DATA_DIR, "GSE104609_ICC8261_edgeR.csv"))
edger8261.columns = ["gene", "logFC", "logCPM", "PValue", "FDR"]

print("edgeR ICC283 rows:", len(edger283))
print("edgeR ICC8261 rows:", len(edger8261))

edger283_sig = edger283[(edger283["FDR"] < 0.05) & (edger283["logFC"].abs() >= 1)]
edger8261_sig = edger8261[(edger8261["FDR"] < 0.05) & (edger8261["logFC"].abs() >= 1)]

print(f"\nedgeR ICC283 significant (FDR<0.05, |logFC|>=1): {len(edger283_sig)}")
print(f"edgeR ICC8261 significant (FDR<0.05, |logFC|>=1): {len(edger8261_sig)}")

icc283 = pd.read_csv(os.path.join(DATA_DIR, "GSE104609_ICC283_DESeq2.csv"))
icc8261 = pd.read_csv(os.path.join(DATA_DIR, "GSE104609_ICC8261_DESeq2.csv"))

deseq283_sig = icc283[(icc283["padj"] < 0.05) & (icc283["log2FoldChange"].abs() >= 1)]
deseq8261_sig = icc8261[(icc8261["padj"] < 0.05) & (icc8261["log2FoldChange"].abs() >= 1)]

print(f"\nDESeq2 ICC283 significant: {len(deseq283_sig)}")
print(f"DESeq2 ICC8261 significant: {len(deseq8261_sig)}")

overlap283 = set(edger283_sig["gene"]) & set(deseq283_sig["Row.names"])
overlap8261 = set(edger8261_sig["gene"]) & set(deseq8261_sig["Row.names"])

print(f"\nOverlap ICC283 (edgeR ∩ DESeq2): {len(overlap283)}")
print(f"Overlap ICC8261 (edgeR ∩ DESeq2): {len(overlap8261)}")