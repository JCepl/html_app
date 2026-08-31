# Restart notes — 2026-08-19, CHELSA fixed25/monotone litmus layer

Written because the Mac hit a fork/process-table exhaustion (every shell
command, even `echo`, failed with "Resource temporarily unavailable") and a
restart was needed to clear it. This is a handoff note so a fresh session
(or you) can pick up exactly where this left off, safely.

## Likely cause of the freeze
No single runaway-CPU process — looked like too many accumulated
processes/threads against macOS's per-user limit (`kern.maxprocperuid`).
Contributing: **multiple Claude Code sessions running at once** (this one,
plus at least one under Positron and one under VS Code) alongside VS Code,
ChatGPT/Codex, and Positron all running their own helper processes. If this
recurs, check with (in a real Terminal, not through Claude):
```bash
ulimit -u
sysctl kern.maxprocperuid
ps -eo pid,ppid,%cpu,%mem,nlwp,comm | sort -rn -k5 | head -20
```
and consider not running 2-3 Claude Code sessions simultaneously.

## State of `~/html_app` right now (safe, nothing destructive happened)

- Branch checked out: **`chelsa-fixed25-litmus`** (new, branched off `main`
  at commit `155b099`, identical content to main at branch time).
- Archive tag on `main`'s pre-change state: **`archive/pre-chelsa-fixed25-litmus-2026-08-19`**
  — if anything ever needs rolling back, `main` itself was never touched.
- **Stashed** (not lost, not discarded): your in-progress, uncommitted
  "Catalogue events" button WIP from before this session
  (`git stash list` → `WIP: catalogue events button (uncommitted, unrelated
  to CHELSA litmus work)`). That work is unrelated to this task — recover it
  later with `git stash pop` (only after committing/finishing the CHELSA
  branch work, to avoid mixing the two).
- **New, untracked, uncommitted file**: `webmap_data/export_chelsa_fixed25_monotone.R`
  — written and ready, but **never successfully executed** (every attempt
  died mid-R-startup on the same fork exhaustion). This is the next actual
  step.

## What this task is (context)

Goal: add ONE new litmus-test layer to the app — the newly validated
CHELSA-native bark-beetle model (fixed25/CATALOGUE25 label,
monotone-constrained XGBoost; full provenance in
`~/REENFOCE_LOCAL_MODEL_WIEN/CHELSA_NATIVE_TRAINING/PROJECT_REFERENCE_bark_beetle_model.md`
§5-9) — as a **new**, separate model entry, without touching or replacing
anything currently live. This is deliberately NOT the same as the existing
`chelsa_mpi_ssp585` layer already in the app, which is the *old*,
pre-monotone, no-held-out-CV model (see that doc's §9 for why the new one
is better-validated).

Source GeoTIFFs (already built, real, on disk, nothing to regenerate):
```
~/REENFOCE_LOCAL_MODEL_WIEN/CHELSA_NATIVE_TRAINING/anna_deliverable/
  barkbeetle_CZU_fixed25_monotone_reference_continuous_1983_2010.tif
  barkbeetle_CZU_fixed25_monotone_future_continuous_2071_2100.tif
```

## Exact next steps after restart

1. `cd ~/html_app && git status` — confirm branch is `chelsa-fixed25-litmus`,
   stash is present, working tree otherwise clean.
2. Run the export script (R + terra; source raster is CHELSA-native
   ~26M px, so the script downsamples to max_dim=768 before reprojecting —
   should take a minute or two, not longer):
   ```bash
   Rscript webmap_data/export_chelsa_fixed25_monotone.R
   ```
   Expect it to write into
   `webmap_data/forecast_bundle/chelsa_fixed25_monotone/{1983_2010,2071_2100}/`:
   `probability_ensemble_mean.tif/.png/_meta.json`,
   `risk_visible_t50.tif/.png/_meta.json`, plus
   `chelsa_fixed25_monotone/metrics_confusion.json` (with **real** CV numbers
   this time — AUC 0.663 forward, F1@0.5=0.408, precision/recall 0.368/0.580
   at the model's real calibrated threshold 0.25 — unlike the old
   `chelsa_mpi_ssp585` layer's all-`NA` metrics).
3. Add a `chelsa_fixed25_monotone` entry to
   `webmap_data/forecast_bundle/manifest.json`'s top-level `"models"` dict
   (copy the shape of the existing `"chelsa_mpi_ssp585"` entry — label e.g.
   `"CHELSA-native fixed25 (monotone XGB)"`, periods
   `["1983_2010","2071_2100"]`, scenario `ssp585`), plus two `"entries"`
   keys `"chelsa_fixed25_monotone::1983_2010"` and
   `"chelsa_fixed25_monotone::2071_2100"` (copy the shape of the existing
   `chelsa_mpi_ssp585::*` entries — `model`/`period`/`scenario`/`gcms`/`note`/
   `ensemble_probability_tif`/`threshold_visible_tifs`).
4. Bump `ASSET_VERSION` in `index.html` (search for that constant near the
   top of the `<script>` block) so browsers don't serve cached PNGs.
5. **No other `index.html` JS changes needed** — the model dropdown is
   populated dynamically from `Object.keys(manifest.models)` (confirmed by
   reading the code this session), so a correct manifest entry alone makes
   the new layer selectable.
6. Serve locally and actually look at it before calling it done:
   ```bash
   python3 -m http.server 8765
   ```
   open `http://localhost:8765/index.html`, select the new model in the
   Predictions layer card, confirm it renders and the legend/metrics panel
   isn't broken.
7. Only after it visibly works: decide with the user whether to commit on
   this branch and merge to `main` (which auto-deploys to GitHub Pages via
   `.github/workflows/` on push) — do not push without checking in, since
   that's a visible/shared action.

## Also from this session (already done, unrelated to the above, no action needed)

`~/REENFOCE_LOCAL_MODEL_WIEN/CHELSA_NATIVE_TRAINING/PROJECT_REFERENCE_bark_beetle_model.md`
got a new §9 added: closed the validation gap for rf/gam (they now have the
same spatial-block-tile + forward/backward CV as xgb/xgb_monotone already
had), on both label sources (fixed25 and efda). New files:
`eval_rf_gam_cv.py`, `merge_4model_cv_summary.py`,
`outputs/rf_gam_cv_results.csv`, `outputs/all_4learner_cv_results.csv`. This
is complete and saved to disk. Note: if `CHELSA_NATIVE_TRAINING`'s parent
repo is git-tracked, these new files are likely still uncommitted — worth a
`git status` check there too after restart.
