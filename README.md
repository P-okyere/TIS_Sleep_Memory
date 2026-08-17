# MYSTI Analysis Code

Analysis code for the MYSTI study: Temporal interference
stimulation (TIS) combined with targeted memory reactivation (TMR) during
NREM sleep, EEG-based sleep microstructure analysis, and behavioural
memory outcomes.

## Pipeline order

Run in this order — each stage's output feeds the next.

### 1. `preprocessing/`
Raw EEG to artifact-corrected, analysis-ready data.
1. `A_MYSTI_FILTER.m` — filter, resample, export scoring EDF
2. `B_MYSTI_CON2EDF_SCORING.m` — bipolar montage export for manual sleep scoring
3. `C_SLEEPMARKERS2EEG_DATA.m` — import manual sleep-stage markers
4. `D_MYST_ART_CORR.m` — bad-channel interpolation, re-referencing, artifact segment removal (interactive)

### 2. `sleep_detection/`
`YASA_SO_Spindle_Detection.py` — YASA-based slow oscillation and spindle
detection (slow/fast/all bands) across all channels, all subjects.

### 3. `coupling_analysis/`
`MYSTI_SO_Spindle_Coupling.m` — trial-based local SO-spindle coupling and
phase, per channel (uses its own SO detection, independent of the YASA
detection above).

### 4. `statistics_figures/`
Cluster-permutation statistics and figures for each measure:
- Fast/slow spindle amplitude, SO amplitude (PTP)
- Fast/slow spindle density, SO density
- Coupling rate vs forgetting (brain-behaviour correlation)

Each measure has a `*_cluster_permutation_stats.m` (or
`*_TOPOplot_gamma_lme.m`-style) script that runs the stats and saves
per-cluster significant-channel CSVs, and where applicable a matching
`*_barplots.m` script that loads those CSVs (never hardcodes channel
lists) and produces the bar chart + pairwise stats.

### 5. `efield/`
E-field simulation results (Sim4Life output).
- `MYSTI_EField_Subregions.Rmd` — field strength/selectivity/collateral across hippocampal head/body/tail
- `MYSTI_EField_Behaviour.Rmd` — field strength vs recollection forgetting

### 6. `behavioural/`
- `MYSTI_Behavioural_Analysis.Rmd` — main behavioural GLMMs (familiarity, recollection, fidelity) with DHARMa diagnostics
- `MYSTI_Behavioural_PrePost.Rmd` — supplementary pre vs post panels


## Requirements

- MATLAB: EEGLAB, FieldTrip, Statistics and Machine Learning Toolbox
- Python: `mne`, `yasa`, `pandas`, `matplotlib`
- R: see `required_packages` at the top of each `.Rmd` (auto-installs on first run)

## Data availability


