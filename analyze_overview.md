# `analyze.m` — codebase overview

Study: double-blind, placebo-controlled ayahuasca trial. Two groups — Healthy (H) and treatment-resistant Depressed (D) — each split into Ayahuasca/Placebo arms, with 25 blood/saliva biomarkers measured before and after dosing. `analyze.m` is the single script that cleans the data, then runs four analyses: (1) which markers separate H from D, (2/3) which markers move toward the healthy distribution after treatment, in D and in H, (4) which marker restorations correlate with self-reported depression improvement (ΔMADRS).

External dependency: `classify_normals`, `plot_boundary`, `colorbarpzn`, `best_linear_classifier` are **not in this repo** — they come from the IntClassNorm MATLAB toolbox (a sibling Geisler-lab repo) and must be on the MATLAB path.

## Core statistical idea

Every analysis reduces to the same operation: fit Gaussian class-conditional models to two groups on one or more (z-scored) markers, take the log-likelihood ratio
$$\mathrm{LLR}(\mathbf m) = \log p(\mathbf m \mid H) - \log p(\mathbf m \mid D),$$
and turn it into a posterior $P(H\mid \mathbf m) = \sigma(\mathrm{LLR})$. A single marker gives a quadratic (two-root) decision boundary when H and D differ in variance as well as mean (e.g. `corti_sal`); combinations of markers give the same idea in higher dimensions. Every accuracy or effect number reported is a **5-fold, repeated, cross-validated** estimate, and every claim of significance is checked against a **1000-permutation label-shuffle null**.

## Pipeline (in execution order)

1. **Load & clean** (`1–32`): read `biomarkers.xlsx`; drop columns that are redundant/linear combinations of others (documented in the Notion notes' pre-processing log); NaN-out implausible `corti_sal` (>1500) and `alt` (>100); log-transform `corti_sal`.
2. **Split & z-score** (`34–64`): split by `timepoint` (before/after) × `group` (H/D) × `treatment` (Ayahuasca/Placebo) into `H0, D0, Ha, Hp, Da, Dp`. Every group is z-scored against the **healthy-baseline (H0)** mean/SD, so H0 is always $N(0,1)$ per marker.
3. **Greedy marker accumulation for H-vs-D separation** (`66–193`): 100 runs of greedy forward selection (each step adds the marker that most improves single-shot 5-fold CV accuracy), averaged into a stable mean-rank ordering (`stable_markers`), then a cumulative-accuracy curve with permutation-null bands. Produces the "mean rank" + "combined separation" figure.
4. **Top-2 2D visualization** (`195–233`): contour plot of $P(H\mid\mathbf m)$ over the best 2-marker panel from step 3, with QDA and best-linear boundaries overlaid.
5. **Per-marker H-vs-D classification** (`235–527`): for each marker, fits the QDA boundary, decomposes it into its two roots and how much each root contributes (`bds_alpha`), computes CV accuracy (QDA and linear) and a permutation null. Markers are then **re-sorted by CV accuracy** into `ind_markers` — this becomes the canonical marker order for the rest of the script. Produces the big "marker values + P(H) gradient" figure with the joint top-2 column prepended.
6. **Restoration analyses** (`529–600`): for D and then for H, greedily finds the marker combination whose **Δ P(H) after Ayahuasca exceeds Δ P(H) after Placebo** (`excess_restore_cv`), then plots per-marker restoration (`restoration_figure`). The Depressed reference distribution is the fixed Healthy baseline (and vice versa for the Healthy analysis).
7. **Treatment-trajectory animation** (`602–818`): animates each subject's (crp, creatinine) point moving from baseline to post-treatment, for all 4 group×arm combinations, exported as `treatment_trajectories.gif`.
8. **ΔMADRS correlations** (`820–1136`): per marker, cross-validated (LOOCV) $R^2$ of $\Delta P(H) \to \Delta \mathrm{MADRS}$ (Ayahuasca arm), with a permutation null; repeated using the raw z-score delta instead of $\Delta P(H)$ for comparison.
9. **Nested elastic-net** (`1137–1246`): predicts $\Delta \mathrm{MADRS}$ from all markers' $\Delta$LLR jointly (pooling Ayahuasca+Placebo), outer LOOCV + inner 5-fold CV for $\lambda$, reports out-of-fold $R^2$ and per-marker selection frequency.
10. **Exploratory plots** (`1248–1335`): pooled per-marker distributions (normality check), full Pearson correlation heatmap, and a rotating 3D boundary GIF for the fixed (`corti_sal`,`crp`,`creatinine`) triplet.
11. **Local functions** (`1337–1840`): `sigmoid`, `null_ci`/`draw_ci_bands` (permutation-null plotting), `perm_null`, `greedy_restoration`, `restoration_figure`, `cv_classify_error` (the core CV engine, QDA or linear), `excess_restore_cv` (the restoration-specific CV engine).

## Bugs found

### 1. Non-positive `corti_sal` values are silently floored, not treated as missing (medium confidence)

At `19–28`, any `corti_sal` reading $\le 0$ is replaced with the smallest *positive* value in that column, then logged — rather than set to NaN like the `>1500` exclusion just above it. This fabricates a data point instead of marking it missing, with no count/warning of how many rows were affected, and no documented justification (unlike the `>1500` rule, which the Notion notes explain).

**Fix**: unless there's a specific reason to floor rather than exclude, change line `26` to set these to NaN instead of `min_vals(c)`, and print how many rows were affected per marker so it's auditable.

## Possible improvements

1. **Duplicated subsetting logic** — FIXED: the Depressed restoration block and the Healthy restoration block each recomputed the same `intersect`-based baseline/post-treatment alignment that was already computed a few lines earlier for the greedy search. Both blocks now reuse the greedy search's `dep_preA/dep_postA/dep_preP/dep_postP` and `heal_preA/heal_postA/heal_preP/heal_postP` directly, and the Healthy greedy-search alignment itself was simplified to intersect on `H0_full`/`Ha_full`/`Hp_full` directly (matching the pattern already used for Depressed) instead of pre-filtering by treatment arm first.
2. **RNG leakage**: `rng(fold)` inside the nested elastic-net loop (`~1175`) reseeds the *global* RNG stream and is never restored. Every `rand()` call after that point (e.g. the jitter in the pooled-distribution plot, `1248+`) becomes silently deterministic across runs. Harmless here since it's cosmetic, but worth wrapping in `s = rng; ... rng(s);` if reproducibility of later figures is ever load-bearing.
3. **Magic numbers**: `corti_sal` threshold (1500), `alt` threshold (100), ridge terms (`1e-5`, `1e-6`), fold/rep counts (5, 50, 100, 1000) are scattered inline. Centralizing as named constants near the top would make the QC assumptions easier to audit — especially since the Notion notes flag several of these (units mismatch on `corti_sal`, unexplained subject-count discrepancies) as still-open questions with Fernanda/Inácio.
