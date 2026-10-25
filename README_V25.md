# Revisiting early-error selection in tropical cyclone ensembles

Figure-reconstruction code and derived data for:

**Revisiting Early-Error Selection in Tropical Cyclone Ensembles: Limited Compensation for Mean-Track Lag**  
Yunze Li and Jie Feng — Department of Atmospheric and Oceanic Sciences, Fudan University.

The study examines how early member selection changes the later ensemble-mean track error and empirical member-range coverage after positional alignment. These materials correspond to manuscript V25. Publication and preprint identifiers will be added when available.

## Reproduce the figures

Tested with Python 3.12 and the versions in `requirements.txt`. Use a virtual environment for installation.

```bash
git clone https://github.com/CloudyLaking/tropical-cyclone-early-selection.git
cd tropical-cyclone-early-selection
python -m pip install -r requirements.txt
python reproduce.py
```

The default command produces the five main figures. To redraw all twenty figures:

```bash
python reproduce.py --all
```

To check the input files without installing plotting dependencies:

```bash
python reproduce.py --verify
```

Outputs are written to `results/` as vector PDF, SVG and 600-dpi PNG. Main figures are named `Figure1`–`Figure5`; supporting figures are `FigureS01`–`FigureS15`. Existing outputs are replaced when the same command is run again. Arial is used when installed, with DejaVu Sans as the fallback; text placement can differ between systems.

## Repository contents

| Location | Contents |
|---|---|
| `src/figure1.py`–`figure5.py` | Main-figure renderers |
| `src/supplementary.py`, `src/spatial.py` | Supporting-figure renderers |
| `data/` | Frozen summaries, storm scores, member ranks, trajectory examples and checksums |
| `assets/natural_earth/` | Coastline for an offline supporting map |

[Data documentation](data/README.md) describes the files, units, cohorts and original sources.

## Scope

This repository reconstructs figures from saved derived data. It does not download or process the original TIGGE archive, refit forecast models, rerun the hindsight member search, or generate new bootstrap samples. It is not an end-to-end rerun of the original study.

The 120-hour motion cohort and the smaller 168-hour cohort are distinct. Hindsight member lists and retrospective center adjustments use future observations and are diagnostic comparisons, not operational forecast methods. The saved best-found lists are not certified global optima.

## Citation and terms

Author and version information is provided in [CITATION.cff](CITATION.cff). The manuscript title above identifies the associated study; no journal DOI or arXiv identifier is currently assigned here.

No software license has yet been specified. Third-party data retain their source conditions, documented in [data/README.md](data/README.md); the Natural Earth coastline is public domain. Repository visibility does not replace those conditions.
