# chelsa_v41_* models (v4.1, added 2026-10-09)

12 manifest models = 4 model families x 3 scenarios (SSP1-2.6, SSP3-7.0, SSP5-8.5), each with the 1981-2010 reference and the
2041-2070 and 2071-2100 normals. CHELSA V2.1 climate, MPI-ESM1-2-HR for the futures. SSP1-2.6 and SSP3-7.0 stand in for RCP2.6 and RCP4.5.

| Key prefix | Model |
|---|---|
| `chelsa_v41_b41lat_*` | B41 (catalogue x EFDA v3 label) + latitude, monotone non-increasing and additive (damps the north by design) |
| `chelsa_v41_b41_*` | B41, climate and topography only |
| `chelsa_v41_a41lat_*` | A41 (fixed25 label) + latitude, free |
| `chelsa_v41_a41_*` | A41, climate and topography only |

Layer value = calibrated relative hazard index on spruce pixels (monotone XGBoost, isotonic calibration), **not an outbreak probability**.

**Metrics panel = testing-period metrics**, the same convention as the other CHELSA-native models: trained on 1983(85)-1997 and 80% of the 200 km
tiles, tested on 2000-2020 in the other 20% of tiles (forward direction), mean of 5 spatial splits; the calibration is fitted on the training data
only, and the cutoffs 0.25 / 0.5 apply to the calibrated layer values. The backward check (train 2000-2020, test 1983(85)-1997) is in each model's note
and in `testing_period_summary`. These are case-control samples, so precision and recall do not describe real prevalence.
Against EFDA v3 disturbance hits the models' mean yearly AUC is 0.57-0.62, no better than a fixed map: read the layers as long-term
susceptibility, not yearly forecasts. Futures are single climatological-normal snapshots from one GCM.

Files per model and period follow the existing bundle: `probability_ensemble_mean`, `risk_visible_t25`, `risk_visible_t50` (.tif, .png, _meta.json).
Built by the RE-ENFORCE BB model node (`CHELSA_NATIVE_TRAINING/v4_spei_rebuild/pipeline/stage_dss_layers.py`).
