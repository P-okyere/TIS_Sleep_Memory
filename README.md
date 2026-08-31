# MYSTI Analysis Code

Analysis code for the MYSTI study: Temporal interference stimulation (TIS) combined with targeted memory reactivation (TMR) during NREM sleep, EEG-based sleep microstructure analysis, and behavioural memory outcomes.

## Pipeline order

Run in this order — each stage's output feeds the next.

### 1. `preprocessing/`
Raw EEG to artifact-corrected, analysis-ready data.

1. `A_MYSTI_FILTER.m` — filter, resample, export scoring EDF
2. `B_MYSTI_CON2EDF_SCORING.m` — bipolar montage export for manual sleep scoring
3. `C_SLEEPMARKERS2EEG_DATA.m` — import manual sleep-stage markers
4. `D_MYST_ART_CORR.m` — bad-channel interpolation, re-referencing, artifact segment removal (interactive)

### 2. `sleep_oscillation_detection/`
`YASA_SO_Spindle_Detection.py` — YASA-based slow oscillation and spindle detection (slow/fast/all bands) across all channels, all subjects.

### 3. `coupling_analysis/`
`MYSTI_SO_Spindle_Coupling.m` — trial-based local SO-spindle coupling and phase, per channel (uses its own SO detection, independent of the YASA detection above).

### 4. `statistics_figures/`
Cluster-permutation statistics and figures for each measure:

- Fast/slow spindle **amplitude**, SO amplitude (PTP) — topoplot + bar-plot scripts for each
- Fast/slow spindle **density**, SO density — topoplot (Gamma GLME) scripts for each; a bar-plot script currently exists for SO density only (`SO_density_barplots.m`) — fast/slow spindle density bar-plot scripts not yet built
- Coupling rate vs forgetting (brain-behaviour correlation) — topoplot + scatter-plot scripts

Each measure has a `*_topoplot.m` script that runs the cluster-permutation stats and saves per-cluster significant-channel CSVs, and where a bar-plot script exists, it loads those CSVs and produces the bar chart + pairwise stats. **Run order matters**: bar-plot/scatter scripts depend on their topoplot script's CSV output and will fail with a clear error if run first.

### 5. `efield/`
E-field simulation results (Sim4Life output). Filtered to `Left-HP_head`/`Left-HP_body`/`Left-HP_tail` only.

- `MYSTI_EField_Subregions.Rmd` — field strength/selectivity/collateral across hippocampal head/body/tail
- `MYSTI_EField_Behaviour.Rmd` — field strength (HP-Head) vs recollection forgetting, no-covariate and covariate (Cue-only-forgetting-controlled) models, simple slopes, interaction plots

### 6. `behavioural/`
- `MYSTI_Behavioural_Analysis.Rmd` — main behavioural GLMMs (familiarity, recollection, fidelity forgetting index, TMR efficacy) with DHARMa diagnostics
- `MYSTI_Behavioural_PrePost.Rmd` — supplementary pre vs post panels (word recognition, associative recall, memory fidelity)

## Requirements

- **MATLAB:** EEGLAB, FieldTrip, Statistics and Machine Learning Toolbox
- **Python:** `mne`, `yasa`, `pandas`, `matplotlib`
- **R:** see `required_packages` at the top of each `.Rmd`/`.R` script (auto-installs on first run)

## Data availability

