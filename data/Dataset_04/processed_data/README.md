# Dataset_04 — processed_data

Source: [GSE282026](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282026) — Human kidney CosMx (~6k gene panel)

## Files

| File | Description |
| --- | --- |
| `../raw_data/sdata.zarr` | SpatialData store: cell table, boundaries, transcript points. |
| `celltyping_config.yaml` | **Single config file** — QC, spatial crop, scanpy settings, rule panels, curation sweep. |
| `celltyping_helpers.py` | Config loader, panel curation (dropout curves), rule-based label assignment. |
| `01_celltypeing_compare.ipynb` | QC, scanpy clustering (approach i), rule-based typing (approach ii), comparison plots. **Kernel:** `spatialdata`. |

## Notebook overview

`01_celltypeing_compare.ipynb` implements two independent cell-typing strategies:

| | Approach i — scanpy | Approach ii — rule-based |
| --- | --- | --- |
| **Method** | Leiden clustering + manual cluster labels | Marker panels + count thresholds |
| **Config** | `cluster_annotations` | `rule_panels`, `rule`, `curation` |
| **Output column** | `cell_type_scanpy` | `cell_type_rule` |
| **Best for** | Exploratory discovery, unknown populations | Known lineages with curated marker sets |

Both approaches write to the same `adata` object so you can compare them in §4 of the notebook.

---

## Approach ii — rule-based typing (detailed)

Rule-based typing assigns each cell a label when it expresses enough markers from a **panel** (e.g. Proximal tubule, Endothelial). Everything is driven by `celltyping_config.yaml` — you should not need to edit Python code to swap gene sets or thresholds.

### Concepts

**Panel** — a named cell type with a list of candidate marker genes.

**Detection** — a gene is detected in a cell when its raw count in `.X` is `> count_threshold` (default `0`, i.e. count ≥ 1).

**Hits** — number of panel genes detected in a cell.

**Assignment** — for each cell, count hits per panel. The panel with the most hits wins **if** hits ≥ that panel's `min_hits`. Ties → `Ambiguous`. No panel qualifies → `Unassigned`.

```text
cell counts:  ALDOB=2  LRP2=0  PECAM1=3  ENG=1
Proximal tubule panel [ALDOB, SLC13A3, APOA1, LRP2]  →  1 hit
Endothelial panel     [ENG, PECAM1, EPAS1, ESAM]      →  2 hits  ← wins (if min_hits ≤ 2)
```

### Two-step workflow

#### Step 1 — Curate panels (§3a in notebook)

You provide **candidate** gene lists per panel. Curation finds a practical subset + `min_hits` by sweeping how many genes must co-occur.

1. Genes are ranked by **detection rate** in the sample (fraction of cells with count > threshold).
2. For `k = 1, 2, …, max_k`, take the top-`k` genes and count cells where **all k** are detected (panel-only, no competition).
3. **Dropout** at step `k` = cells lost vs `k−1`.
4. **Recommendation:** pick `k* = k_cliff − 1`, where `k_cliff` is the step with the **largest dropout** (the cliff right before assignment collapses).

Example (Endothelial):

```text
k=1  [EPAS1]                         →  2,096 cells
k=2  [EPAS1, ENG]                   →    400 cells   dropout 1,696  ← cliff
k=3  [EPAS1, ENG, PECAM1]           →     12 cells
recommended: k*=1, markers=[EPAS1], min_hits=1
```

Run §3a cells → review summary table and dropout plots → copy the printed YAML snippet.

#### Step 2 — Assign labels (§3b in notebook)

Paste curated panels into `rule_panels` as fixed `markers` + per-panel `min_hits`:

```yaml
rule_panels:
  Endothelial:
    markers: [EPAS1]
    min_hits: 1
  Proximal tubule:
    markers: [ALDOB]
    min_hits: 1
```

Re-run the config cell and assignment cells. Labels land in `adata.obs["cell_type_rule"]`.

If `markers` is omitted, the notebook falls back to auto-selecting the top `rule.top_n_markers` genes by detection rate and uses the global `rule.min_hits`.

### Config reference (approach ii)

| Section | Key | Meaning |
| --- | --- | --- |
| `curation` | `max_k` | Max genes to sweep in dropout curve |
| `curation` | `max_candidates` | Top-N candidates considered per panel |
| `curation` | `count_threshold` | Count must be **>** this value to count as detected |
| `curation` | `min_assigned_cells` | If recommendation assigns fewer cells, fall back to smaller k |
| `rule_panels.*` | `candidates` | Input gene pool for curation |
| `rule_panels.*` | `markers` | Fixed curated set (skips auto-selection) |
| `rule_panels.*` | `min_hits` | Required hits for that panel (overrides global) |
| `rule` | `top_n_markers` | Auto-select fallback when `markers` not set |
| `rule` | `min_hits` | Default min hits when panel has no override |
| `rule` | `count_threshold` | Used at assignment time (should match curation) |
| `rule` | `ambiguous_label` / `unassigned_label` | Label names for edge cases |

### Helper functions (`celltyping_helpers.py`)

| Function | Role |
| --- | --- |
| `curate_rule_panels()` | Run dropout sweep for all panels |
| `recommend_panel_markers()` | Pick k* before max dropout |
| `resolve_panel_markers()` | Resolve markers + per-panel min_hits from config |
| `assign_rule_labels_from_panels()` | Apply multi-panel rules to all cells |
| `format_curated_yaml_snippet()` | Print YAML to paste into config |

### Practical tips

- **Sparse panels:** requiring all 4 top genes (`min_hits: 4`) often assigns almost nothing on CosMx kidney data. Use curation curves to find the elbow.
- **Missing genes:** curation reports `missing_genes` when candidates are absent from `adata.var_names` (common for tissue-specific markers not in the run panel).
- **Module reload:** cell 1 runs `importlib.reload(celltyping_helpers)` so edits to the helper file are picked up without a kernel restart.
- **Competition:** curation counts each panel in isolation; final assignment can differ slightly when panels overlap.

---

## Next steps: gene preference (not yet implemented)

Today, curation ranks candidates **only by detection rate**. The top gene is whichever is most often detected in the sample — not necessarily the gene you trust most biologically.

### The problem

Detection rate and biological specificity pull in different directions.

| Gene | Detection rate | Biological note |
| --- | --- | --- |
| `CASR` | Higher in this sample | Weaker TAL marker |
| `UMOD` | Lower or absent | Canonical thick ascending limb marker |

Auto-curation will prefer `CASR` because it assigns more cells, even when you know `UMOD` is the marker you want for TAL identity. There is currently **no way to express "prefer UMOD"** in the config — only `candidates` (pool), `markers` (fixed override after manual curation), and detection-driven ranking.

### What you can do today (workaround)

**Manually pin markers after reviewing the curve.** Curation tells you how many genes you can require before the cliff; you choose *which* genes:

```yaml
  Thick ascending limb:
    markers: [UMOD, SLC12A1]   # your choice, not auto top-k
    min_hits: 2                 # from dropout curve: stop before the cliff
```

Workflow:

1. Run curation → read `k_cliff` and recommended `min_hits` from the dropout plot.
2. Ignore the auto-recommended gene list if you disagree.
3. Pick your preferred genes from `candidates` (check `panel_detection_summary` for which are present).
4. Set `min_hits` to the largest value that still gives acceptable coverage **for your chosen genes** (re-count manually or temporarily set `markers` and re-run assignment).

This keeps biology in your hands while using the tooling for **stringency** (how many genes to require), not **gene choice**.

### Proposed extension: `preferred` / `required` genes

A natural next improvement is to add optional fields under each panel in `celltyping_config.yaml`:

```yaml
  Thick ascending limb:
    candidates: [UMOD, SLC12A1, KCNJ1, ...]
    required: [UMOD]          # must always be in the curated set (if detected in panel)
    preferred: [UMOD, SLC12A1]  # rank ahead of detection rate when building top-k
```

**How `required` would work**

- Always include these genes in `markers` when present in `adata.var_names`.
- Fill remaining slots up to `k*` using detection rate among other candidates.
- Dropout sweep runs on the **combined** set (required + data-driven additions).
- If a required gene is missing from the panel entirely, warn in the curation summary.

**How `preferred` would work (softer than required)**

- Re-rank candidates before the sweep: preferred genes first (in listed order), then the rest by detection rate.
- Dropout curve then answers: "given that I always try UMOD first, how many *additional* genes can I require?"
- Does not force a gene that is undetected in the sample, but stops CASR from jumping ahead of UMOD when both are present.

**Tie-breaking at assignment (optional future)**

When two panels tie on hits, prefer the panel whose **required** genes are detected, or whose preferred genes have higher total counts. Today ties become `Ambiguous` with no biology-aware tie-break.

### Implementation sketch (for developers)

Changes would live in [`celltyping_helpers.py`](celltyping_helpers.py):

1. `_rank_panel_candidates(candidates, preferred, required, detection_rates)` — custom sort order.
2. `panel_dropout_curve()` — seed curve with `required`, append preferred/data-driven genes for k > len(required).
3. `recommend_panel_markers()` — unchanged logic on the new curve.
4. `curate_rule_panels()` — surface warnings for missing required genes.
5. Config schema docs in this README + comments in `celltyping_config.yaml`.

No notebook changes beyond displaying the new rank order and warnings in the curation summary table.

### Why this matters

The dropout curve solves **"how strict should the rule be?"** (how many genes must co-occur). Gene preference solves **"which genes define the rule?"** (biology vs abundance). Both are needed for reproducible, interpretable rule-based typing on sparse spatial data.

---

## Kernel

Use **`spatialdata`** (Python ≥3.11, spatialdata ≥0.7.2) for `01_celltypeing_compare.ipynb`.

```bash
/mnt/scratch1/miniconda3/envs/spatialdata/bin/python -c "import spatialdata; print(spatialdata.__version__)"
```
