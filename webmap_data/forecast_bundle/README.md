# Forecast Bundle README

## 1) What this bundle is

This folder contains an app-ready, model-selectable climate risk forecast bundle for two future periods:

- 2041-2070
- 2071-2100

Supported models:

- xgb
- rf
- gam
- svm

Each model/period package includes:

- ensemble mean probability raster (continuous risk)
- thresholded rasters (below-threshold values masked to NA; suitable for transparent rendering)
- per-period metadata (`entry.json`)
- model-level metrics/confusion payload (`metrics_confusion.json`)

The top-level `manifest.json` is the entry point for app loading.

---

## 2) Where this came from (data genesis)

### 2.1 Source workspace roots

- Modeling root:
  - `/Users/jaroslavcepl/REENFOCE_LOCAL_MODEL_WIEN/FINAL_MODEL`
- App root:
  - `/Users/jaroslavcepl/Documents/D_MAC_2026/REENFORCE_postSopron/html_app`
- Bundle output root:
  - `/Users/jaroslavcepl/Documents/D_MAC_2026/REENFORCE_postSopron/html_app/webmap_data/forecast_bundle`

### 2.2 Assembly script

Bundle was assembled by:

- `/Users/jaroslavcepl/REENFOCE_LOCAL_MODEL_WIEN/FINAL_MODEL/assemble_app_forecast_data.R`

This script is the canonical provenance document for the generated files.

### 2.3 Forecast raster source experiments

Model-to-forecast source mapping:

- xgb:
  - `/Users/jaroslavcepl/REENFOCE_LOCAL_MODEL_WIEN/FINAL_MODEL/EXP_future_forecasts_from_pretrained_4gcm_v3`
- rf:
  - `/Users/jaroslavcepl/REENFOCE_LOCAL_MODEL_WIEN/FINAL_MODEL/EXP_future_forecasts_rf_from_pretrained_eval_20260505_094005`
- gam:
  - `/Users/jaroslavcepl/REENFOCE_LOCAL_MODEL_WIEN/FINAL_MODEL/EXP_future_forecasts_gam_svm_spruceclip060_20260504_223950`
- svm:
  - `/Users/jaroslavcepl/REENFOCE_LOCAL_MODEL_WIEN/FINAL_MODEL/EXP_future_forecasts_gam_svm_spruceclip060_20260504_223950`

From each valid model/GCM/period folder, the script reads:

- `map_whole_test_pred_probability.tif`
- `pixel_ever_attacked_confusion.txt`

### 2.4 Metrics source experiments

Metrics are intentionally sourced from historical evaluation/validation runs (not future truth) to provide meaningful confusion/skill tables.

- xgb metrics source:
  - `/Users/jaroslavcepl/REENFOCE_LOCAL_MODEL_WIEN/FINAL_MODEL/EXP_eval_1985_2005_vs_2006_2020_nolonlat_20260504_105347/cutoff_metrics_all_runs.tsv`
- rf metrics source:
  - `/Users/jaroslavcepl/REENFOCE_LOCAL_MODEL_WIEN/FINAL_MODEL/EXP_eval_1985_2005_vs_2006_2020_rf_nolonlat_20260505_022721/cutoff_metrics_all_runs.tsv`
- gam metrics source:
  - `/Users/jaroslavcepl/REENFOCE_LOCAL_MODEL_WIEN/FINAL_MODEL/EXP_validate_gam_svm_nolonlat_windExcluded_20260504_155444/overall_thresholds_grid_0p05.tsv` (filtered to `model=gam`)
- svm metrics source:
  - `/Users/jaroslavcepl/REENFOCE_LOCAL_MODEL_WIEN/FINAL_MODEL/EXP_validate_gam_svm_nolonlat_windExcluded_20260504_155444/overall_thresholds_grid_0p05.tsv` (filtered to `model=svm`)

---

## 3) Processing logic used to create this bundle

For each model and period (2041_2070, 2071_2100):

1. Discover valid period folders in the configured forecast source.
2. Keep only folders that contain `map_whole_test_pred_probability.tif`.
3. Stack all valid probability rasters across available GCMs.
4. Compute ensemble mean probability raster:
   - `ens_mean = mean(stack, na.rm=TRUE)`
5. Write `probability_ensemble_mean.tif` with compression:
   - `COMPRESS=DEFLATE`, `PREDICTOR=2`, `ZLEVEL=6`
6. For thresholds in `{0.50, 0.55, ..., 0.90}`:
   - write `risk_visible_tXX.tif` where values below threshold are set to `NA`
7. Parse per-GCM `pixel_ever_attacked_confusion.txt` into `forecast_confusion_by_gcm` payload.
8. Write per-period `entry.json`.

For each model:

1. Load configured metrics table.
2. Normalize schema (supports both XGB/RF and GAM/SVM table formats).
3. Build:
   - `by_gcm` metrics by cutoff
   - `aggregate_by_cutoff` metrics (summing TP/FP/FN/TN across GCMs, then recomputing precision/recall/specificity/F1/F05/red_fraction)
4. Write `metrics_confusion.json`.

Finally:

- write top-level `manifest.json` with model catalog, period catalog, threshold catalog, and all `entry` objects.

---

## 4) Folder and file structure

Top-level:

- `manifest.json`
- `xgb/`
- `rf/`
- `gam/`
- `svm/`

Per model:

- `metrics_confusion.json`
- `2041_2070/`
- `2071_2100/`

Per model-period:

- `entry.json`
- `probability_ensemble_mean.tif`
- `risk_visible_t50.tif`
- `risk_visible_t55.tif`
- `risk_visible_t60.tif`
- `risk_visible_t65.tif`
- `risk_visible_t70.tif`
- `risk_visible_t75.tif`
- `risk_visible_t80.tif`
- `risk_visible_t85.tif`
- `risk_visible_t90.tif`

---

## 5) Data contract for app usage

### 5.1 `manifest.json`

Purpose:

- discover available models
- discover available periods
- discover available threshold presets
- retrieve model metrics file path
- retrieve period entry metadata

Key fields:

- `created_utc`
- `output_root`
- `periods`
- `thresholds`
- `models.{model}.metrics_file`
- `models.{model}.periods`
- `entries.{model::period}`

### 5.2 `entry.json`

Purpose:

- load raster files for one model-period
- expose source provenance and per-GCM forecast confusion metadata

Key fields:

- `model`
- `period`
- `gcms`
- `source_probability_files`
- `ensemble_probability_tif`
- `threshold_visible_tifs.{threshold}`
- `forecast_confusion_by_gcm.{gcm}`

### 5.3 `metrics_confusion.json`

Purpose:

- render confusion/skill tables and charts by model
- allow threshold-dependent metric display in UI

Key fields:

- `model`
- `source_eval_metrics_file`
- `notes`
- `by_gcm.{gcm}[]` rows with:
  - `cutoff`, `TP`, `FP`, `FN`, `TN`, `precision`, `recall`, `specificity`, `red_fraction`, `F1`, `F05`
- `aggregate_by_cutoff[]` rows with same metrics aggregated over GCMs

---

## 6) Threshold behavior and transparency semantics

Threshold set:

- 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90

Semantics of `risk_visible_tXX.tif`:

- if `probability >= threshold`: retain probability value
- if `probability < threshold`: set to `NA`

This directly supports map rendering where sub-threshold pixels are transparent.

---

## 7) Coverage and known limitations

### 7.1 Model-period availability

Bundle includes all four models for both periods:

- xgb: 2041_2070, 2071_2100
- rf: 2041_2070, 2071_2100
- gam: 2041_2070, 2071_2100
- svm: 2041_2070, 2071_2100

### 7.2 GCM membership in ensembles

- xgb: 3 GCM contributors in current source (ICHEC, MOHC, REMO)
- rf: 3 GCM contributors in current source (ICHEC, MOHC, REMO)
- gam: 2 GCM contributors in current source (ICHEC, MOHC)
- svm: 2 GCM contributors in current source (ICHEC, MOHC)

Reason:

- some GAM/SVM source runs for REMO were failed/empty in the selected forecast source folder.

### 7.3 Forecast confusion values

In many future forecast folders, `pixel_ever_attacked_confusion.txt` has NA counts (TP/FP/FN/TN), which is expected without future observed truth labels. Those files still carry useful metadata such as selected threshold and red fraction.

### 7.4 Metrics provenance mismatch by design

- Forecast maps are future projections.
- Metrics/confusions are from historical evaluation/validation experiments.

This is intentional so UI can show model behavior statistics from known truth conditions.

---

## 8) Reproducibility

To regenerate this bundle:

1. Ensure source experiment folders and metrics TSV files listed above exist.
2. Run:

```bash
Rscript /Users/jaroslavcepl/REENFOCE_LOCAL_MODEL_WIEN/FINAL_MODEL/assemble_app_forecast_data.R
```

3. Verify:

- `manifest.json` exists and lists 4 models
- each model has both periods
- each period has `probability_ensemble_mean.tif` and 9 threshold rasters

---

## 9) Suggested app loading sequence

1. Load `manifest.json`
2. Populate model selector from `manifest.models`
3. Populate period selector from `manifest.periods`
4. Populate threshold selector from `manifest.thresholds`
5. For selected model-period:
   - use `entries.{model::period}.ensemble_probability_tif` for continuous mode
   - use `entries.{model::period}.threshold_visible_tifs[{threshold}]` for threshold-transparent mode
6. Load confusion/metrics from `models.{model}.metrics_file`

---

## 10) Contact points for future edits

Primary file to edit when data sourcing rules change:

- `/Users/jaroslavcepl/REENFOCE_LOCAL_MODEL_WIEN/FINAL_MODEL/assemble_app_forecast_data.R`

Primary output contract to keep stable for frontend:

- `manifest.json`
- `entry.json`
- `metrics_confusion.json`
