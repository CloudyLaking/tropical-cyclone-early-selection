# Tropical cyclone ensemble post-processing

Code and derived data for **Revisiting Early-Error Selection and Weighting in Tropical Cyclone Ensembles**, Yunze Li and Jie Feng, Fudan University. Manuscript in preparation; no journal DOI or arXiv identifier is assigned here.

## Current manuscript: V30.5

[Reproduction instructions and scientific scope](v30_5/README.md)

```bash
git clone https://github.com/CloudyLaking/tropical-cyclone-early-selection.git
cd tropical-cyclone-early-selection/v30_5
python reproduce.py --verify
python -m pip install -r requirements.txt
python reproduce.py --checks --out results
```

The current workflow redraws four main figures and reconstructs principal weight/order, member-range/coverage and circle-margin estimates using frozen derived coordinates, optimized-weight certificates and paired storm draws. Every optimized weighting family uses later truth and is a diagnostic limit, not an operational forecast method. The fixed 120-h and seven-day samples are different.

This is derived-data reproduction, not raw-archive reconstruction or external replication. The [detailed README](v30_5/README.md) lists what is and is not rerun. No manuscript, private review, login material or licensed font is included.

## Earlier materials

The original V25 workflow remains at the repository root with its original five-main-figure and fifteen-supporting-figure numbering. See [README_V25.md](README_V25.md). Its figure numbers and title must not be confused with the current manuscript.

## Citation and terms

See [CITATION.cff](CITATION.cff). The current code/data release is version 2.0.0; this version number does not indicate publication of the manuscript. No software license has yet been assigned. Third-party data retain their source conditions, described in the [current](v30_5/README.md) and [earlier](data/README.md) data documentation. Public visibility does not replace those conditions.
