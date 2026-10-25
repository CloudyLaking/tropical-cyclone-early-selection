# V30.6 derived-data reproduction

**Revisiting Early-Error Selection and Weighting in Tropical Cyclone Ensembles**  
Yunze Li and Jie Feng, Fudan University. Manuscript in preparation, 25 October 2026.

This package reconstructs the four current main figures and the principal weight/order, range/coverage and circle-margin estimates from frozen derived data. It is separate from the older V25 workflow at the repository root. No manuscript PDF, private review, credential, raw TIGGE archive or licensed font is included.

## Run

Use Python 3.12 with NumPy, pandas, SciPy and Matplotlib (the tested versions are recorded in `requirements.txt`). The hash check uses only Python's standard library:

```bash
python reproduce.py --verify
python -m pip install -r requirements.txt
python reproduce.py --checks --out results
```

Omit `--checks` to redraw the four main figures and supplementary displays S4–S8/S19. Existing output folders are never replaced. The output contains vector PDF, SVG and 600-dpi PNG files. DejaVu Sans is the portable default. To reproduce the manuscript's Arial typography, supply `--font-dir /path/to/your/licensed/fonts`, containing `arial.ttf` and `arialbd.ttf`. Font files are not redistributed. Font substitution can change line lengths and must be visually checked.

## What the checks establish

* `verify_weights.py` checks feasibility, endpoints and convex optimality of all 838,240 saved monotone-weight certificates over 806 initializations and four leads, plus free nonnegative weights. It independently reconstructs fixed-split and best-split two-level limits and the paired storm means/intervals shown in Figure 3. Saved random permutations and 5,000 paired storm draws are reused.
* `verify_protocol.py` rescoring starts from 6,885 frozen member-coordinate records. It reconstructs center distance, ATE/CTE, spread, circle radius, circle/ellipse coverage, MW/CT equal-weight mean identity, and the corresponding fixed/available-cohort summaries. Ellipse scoring uses centered-data SVD.
* `verify_margin.py` reconstructs all four 120-h circle-margin comparisons from paired case records and saved storm draws.
* Plotting reads the frozen tables without fitting a new model or changing sample membership. Figure 1 is an idealized vector schematic, not an observation or a physical mechanism.

These are internal reproducibility checks, not an external replication. They do not redownload TIGGE/IBTrACS, reconstruct source-file matching or storm eligibility, rerun the earlier best-found subset searches, or regenerate the other historical SI displays. The seven-day forecast curves are redrawn from saved performance summaries; their original resampling is not rerun here. The older V25 figure workflow remains available separately, with its original numbering.

## Scientific boundaries

Members are ranked by raw 24-h error before memberwise translation; retained identities stay fixed. The primary 120-h cohort contains 300 ECMWF initializations/64 storms and 506 NCEP initializations/100 storms. The seven-day cohort is different: 100/21 and 195/41. ATE is positive forward and CTE positive to the right of observed motion.

Optimized weights use verifying observations, including the early-order optimization, and therefore describe conditional geometric limits rather than deployable forecast skill. `D` restricts weights to decrease along the order; `I` restricts them to increase; `M` takes the better of D and I. Fixed split uses k=floor(n/2); Best split also optimizes k. Both optimize total group mass alpha, with per-member weights alpha/k and (1-alpha)/(n-k), permitting reversal between groups. The two-group segment splits at alpha=k/n into D and I branches; M is invariant to complete order reversal. MW and CT have the same center for a fixed equal-weight member list but different member clouds. Empirical 70th-percentile circles and ellipses are not calibrated 70% forecast-probability regions. Intervals are pointwise, exploratory paired storm-percentile intervals.

## Files and provenance

| Folder | Contents |
|---|---|
| `data/fig2` | Seven-day performance summaries, 120-h storm scores and paired gain inputs |
| `data/fig3` | Paired weight contrasts and matched Hagibis geometry |
| `data/fig4` | MW/CT clouds, boundaries, center/radius, margin and coverage summaries |
| `data/ranking/cases` | 806 NumPy archives: member errors/identities, orders, weights, endpoints and distances |
| `data/protocol/frozen` | Frozen coordinates, fixed random subsets and saved bootstrap indices |
| `data/protocol/results` | Reference case scores and cohort summaries |
| `data/margin` | Paired score records and 120-h storm resampling indices |
| `code` | Figure sources and deterministic checks |
| `MANIFEST.json` | SHA-256 inventory; checked before every run |

Coordinates and errors are in km unless a file field explicitly indicates otherwise. Coverage fields are fractions; Figure 4 converts contrasts to percentage points. `A` denotes memberwise translation (MW), and `B` common translation (CT). Selected means use floor(n/2) members. Random lists have the same size and are averaged over 128 lists per initialization.

Original source archives are [NSF NCAR TIGGE, d330003](https://gdex.ucar.edu/datasets/d330003/) and [IBTrACS v04r01](https://doi.org/10.25921/82ty-9e16). Their source conditions continue to apply to derived materials. Pickle inputs are legacy frozen research records; load only files matching this trusted release's manifest. No software license has yet been assigned; public visibility does not grant an additional license. Publication identifiers will be added when available.

## Historical correction intervals

`code/rebootstrap_history.py` directly resamples the 45 ECMWF / 67 NCEP prior-year storms, with 5,000 paired draws (PCG64 seed 20261025). Predictions and point means are unchanged. The frozen input storm scores, old and revised summaries, ordered storm IDs and draw arrays are in `data/history/`. These conditional intervals do not refit correction models. `--checks` reconstructs them alongside the existing checks.
