# Chickpea Drought Transcriptomic Re-analysis

Analysis code and derived results supporting:

> Asmatullah M. "Cross-Study Identification and Validation of Transcriptional
> Signatures Associated with Drought Tolerance in Chickpea (Cicer arietinum L.)."
> Manuscript draft based on independent re-analysis of the public GSE104609 dataset.

This repository accompanies a reproducibility-oriented secondary re-analysis
of the publicly available leaf RNA-seq dataset GSE104609 (BioProject PRJNA413294),
comparing drought responses of drought-tolerant ICC8261 and drought-sensitive
ICC283 chickpea genotypes.

## Repository structure

- data/ : Input data (see data/README.md)
- scripts/ : Analysis scripts (run in numerical order)
- results/ : Output tables produced by the scripts

## How to reproduce

1. Download the four deposited processed files from NCBI GEO accession
   GSE104609 and place them in data/.
2. Run:
   python scripts/01_derived_genotype_response.py
   python scripts/02_crossmethod_consistency.py
3. Outputs are written to results/.

## License

MIT License (see LICENSE file).

## Citation

NCBI Gene Expression Omnibus. GSE104609: RNA sequencing of leaf tissues
from two contrasting chickpea genotypes reveals mechanisms for drought
tolerance.