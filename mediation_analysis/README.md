# Trial-level EEG feature mediation analysis — analysis code

Code accompanying the trial-level mediation analysis of TIS-TMR effects on
post-nap associative memory (recollection) accuracy.

Author: Junheng Li

---

## Repository contents

```
EEGBeh CleanedFinalSubmission Aug26/
│
├── README.md
│
├── s01_load_nap_data.m                 MATLAB   epoching + behaviour linking
├── s02_extract_eeg_features.m          MATLAB   51 trial-level EEG features
├── s03_select_features.m               MATLAB   collinearity reduction 51 -> 16
├── s04_path_a_analysis.R               R        path-a models
├── s05_mediation_joint_surrogate.R     R        joint mediation + permutation test
│
├── Plots/                              figure scripts (run after s04 / s05)
│   ├── s06_plot_path_a_heatmap.m           path-a t heatmap        (Supp. figure)
│   ├── s07_plot_normalisedACME_violin.m    normalised ACME violins (Fig. b)
│   ├── s08_plot_normalisedACME_by_region.m normalised ACME by region (Fig. c)
│   ├── s09_plot_normalisedACME_brainmap.m  scalp bubble map        (Fig. d)
│   └── Output/                             s06 figure (regenerated on run)
│       ├── path_a_heatmap.svg
│       └── path_a_heatmap.png
│
├── Data for reproducing final plots/    published results tables (s04 / s05 output)
│   ├── path_a_results.csv                  48 rows: 16 features x 3 contrasts
│   └── multivariate_mediation_results.csv  48 rows: 16 features x 3 contrasts
│
└── Functions/                           helpers that must be on the MATLAB path
    ├── chanlocs.mat                        channel locations (used by s02)
    ├── Features/
    │   ├── extract_aperiodic_component.m   FOOOF wrapper (called by s02)
    │   └── find_Spectral_Slope.m           multitaper PSD -> aperiodic exponent
    │                                       (companion helper, not called by s02)
    └── Plotting/
        └── Violinplot-Matlab-master/       bundled Violinplot-Matlab (used by s07)
```

`Functions/` is bundled, so nothing has to be supplied separately apart from
Python/FOOOF (see *Requirements*). Add it to the path before running anything:

```matlab
addpath(genpath('Functions'))
```

`Plots/Output/` ships with the current s06 figure and is overwritten when s06
is re-run. Folders written at run time and **not** shipped: `Results/Plots/`
(s07–s09 figures), and `PathA_RelBase/` + `Mediation_JointModel_RelBase/`
from the R scripts.

---

## Pipeline overview

Scripts run in numbered order. MATLAB produces the trial-level feature table,
R fits the models, MATLAB produces the figures.

```
  raw EEG (.set) + events (.tsv) + behaviour (.xlsx)
        │
  [s01] s01_load_nap_data.m ................. epoch conditions, link behaviour
        └─> NapData_AllSubjects.mat
        │
  [s02] s02_extract_eeg_features.m .......... 51 trial-level EEG features
        └─> EEG_Features_Behaviour_RelBase.csv
        │
  [s03] s03_select_features.m ............... collinearity reduction (51 → 16)
        ├─> Results_SelectedFeatures_RelBase.csv
        ├─> Feature_Clustering_Results.mat
        └─> Feature_Clustering_RelBase.svg          (Supplementary figure)
        │
        ├── [s04] s04_path_a_analysis.R ...... path-a: condition → feature
        │         └─> PathA_RelBase/path_a_results.csv, path_a_summary.txt
        │              └── [s06] Plots/s06_plot_path_a_heatmap.m   (Supp. figure)
        │
        └── [s05] s05_mediation_joint_surrogate.R
                  joint mediation model + permutation surrogate test
                  └─> Mediation_JointModel_RelBase/multivariate_mediation_results.csv
                       ├── [s07] Plots/s07_plot_normalisedACME_violin.m      (Fig. b)
                       ├── [s08] Plots/s08_plot_normalisedACME_by_region.m   (Fig. c)
                       └── [s09] Plots/s09_plot_normalisedACME_brainmap.m    (Fig. d)
```

The two results tables produced by s04 and s05 are the ones shipped in
`Data for reproducing final plots/`, so s06–s09 can be run on their own.

---

## Quick start: reproducing the figures without re-running the models

The figure scripts read the shipped tables in
`Data for reproducing final plots/`. **Start MATLAB with this folder as the
working directory** and call the scripts by name — do not use `run(...)`, which
changes the working directory to `Plots/` and breaks the relative input paths
in s07–s09.

```matlab
addpath(genpath('Functions'))
addpath('Plots')

s06_plot_path_a_heatmap             % Supplementary path-a heatmap
s07_plot_normalisedACME_violin      % Figure b
s08_plot_normalisedACME_by_region   % Figure c
s09_plot_normalisedACME_brainmap    % Figure d
```

s06 resolves its own paths from the script location, so it is safe to run from
any working directory. s07–s09 resolve theirs relative to the working
directory, which must therefore be the repository root.

This takes seconds. Re-running the full pipeline from s01 requires the raw data
and several hours of compute for the permutation test in s05 (see *Runtime*).

---

## Scripts

| # | Script | Language | What it does |
|---|--------|----------|--------------|
| s01 | `s01_load_nap_data.m` | MATLAB | Reads per-subject EEGLAB `.set` files, artefact intervals and task event tables. Epochs the three stimulation conditions (SoundOnly/Cue, TICue/TIS+Cue, TIB4Cue/TIS+Cue_Delayed), assigns each event its scored sleep stage, and matches each trial to the post-nap behavioural response by delivered word. Artefact-overlapping trials are flagged and set to NaN in full. Output: `NapData_AllSubjects.mat`. |
| s02 | `s02_extract_eeg_features.m` | MATLAB | Computes 8 EEG features × 6 regions + 3 frontal→central cross-channel coupling features = **51 features per trial**, each as post-stimulus (0.2–3.0 s) minus baseline (−1.0 to −0.2 s). Applies the behavioural coding used in the paper: trials without a valid post-nap word response are removed; trials with an incorrect Old/New judgement are coded as memory failures (0). Output: `EEG_Features_Behaviour_RelBase.csv`. |
| s03 | `s03_select_features.m` | MATLAB | Spearman correlation across the 51 features → distance matrix (1 − \|r\|) → average-linkage hierarchical clustering, cut at 0.5. Retains the centroid-closest member of each cluster, yielding **16 representative features**. Outputs the reduced trial table and the supplementary clustering figure. |
| s04 | `s04_path_a_analysis.R` | R | Standalone path-a: `feature_z ~ treat + (1\|subject)` per feature per contrast (`lmerTest`, Satterthwaite df), FDR-corrected (BH) within contrast. Outcome variable is not involved. Output: `PathA_RelBase/`. |
| s05 | `s05_mediation_joint_surrogate.R` | R | Joint mediation model. Path-b is a single binomial GLMM with **all 16 mediators plus treatment simultaneously**, so each per-feature ACME is the unique mediated effect. ACME by quasi-Bayesian Monte Carlo (1000 sims, `mediation`). Significance additionally assessed by within-subject treatment-label permutation (1000 permutations, full pipeline refitted each time, 200 sims per refit), giving a one-sided right-tailed `Perm_p`. Output: `Mediation_JointModel_RelBase/`. |
| s06 | `Plots/s06_plot_path_a_heatmap.m` | MATLAB | Path-a *t*-value heatmap, features × contrasts, grouped by region. Significance markers use the **FDR-adjusted** *p* (`path_a_p_adj`), matching the manuscript: boxed cell + asterisks for adjusted *p* < 0.05. Three cells qualify. |
| s07 | `Plots/s07_plot_normalisedACME_violin.m` | MATLAB | Distribution of the per-feature normalised ACME across the 16 features, one violin per contrast, with a two-sided one-sample *t*-test against zero. |
| s08 | `Plots/s08_plot_normalisedACME_by_region.m` | MATLAB | Per-feature normalised ACME as a horizontal bar chart, grouped by region, coloured by sign. |
| s09 | `Plots/s09_plot_normalisedACME_brainmap.m` | MATLAB | Region-averaged normalised ACME rendered as fixed-radius bubbles on a schematic scalp, diverging colour scale. |

### Terminology

**Normalised ACME** in the manuscript = the `ACME_CoV` column in the code
(ACME posterior mean ÷ SD across quasi-Bayesian simulation draws). This is the
quantity plotted in s07, s08 and s09.

---

## Conditions and contrasts

| Code name | Manuscript name |
|---|---|
| `SoundOnly` | Cue |
| `TICue` | TIS+Cue |
| `TIB4Cue` | TIS+Cue_Delayed |

Three pairwise contrasts are fitted throughout: TICue vs SoundOnly (reported in
the main text), TIB4Cue vs SoundOnly, and TICue vs TIB4Cue.

---

## Data dictionary

### `Data for reproducing final plots/path_a_results.csv`

One row per feature × contrast (16 × 3 = 48). Produced by s04, read by s06.

| Column | Description |
|---|---|
| `Contrast` | `TICue_vs_SoundOnly` \| `TIB4Cue_vs_SoundOnly` \| `TICue_vs_TIB4Cue` |
| `Feature` | `<region>_<feature>` name of the mediator |
| `N`, `N_subjects`, `N_excluded` | Trials fitted, subjects contributing, trials dropped for missing values |
| `singular` | Whether the LMM fit was singular |
| `path_a_coef`, `path_a_se` | Fixed-effect estimate of treatment on the z-scored feature, and its SE |
| `path_a_t`, `path_a_df` | *t* statistic and Satterthwaite df — `path_a_t` is the quantity plotted by s06 |
| `path_a_p` | Uncorrected *p* |
| `path_a_p_adj` | **BH/FDR-adjusted *p* within contrast. This is the *p* value reported in the manuscript and used for the significance markers in s06.** |
| `path_a_CI_lo`, `path_a_CI_hi` | 95% CI on `path_a_coef` |
| `path_a_d` | Standardised effect size |

Cells passing `path_a_p_adj < 0.05` (the three boxed cells in the s06 figure):

| Contrast | Feature | *t* | *p* | *p*<sub>adj</sub> |
|---|---|---|---|---|
| TICue vs SoundOnly | `central_aperiodic` | +3.473 | 5.23e-04 | 8.36e-03 |
| TICue vs SoundOnly | `posterior_spindle_peak_freq` | −2.859 | 4.27e-03 | 3.42e-02 |
| TIB4Cue vs SoundOnly | `central_aperiodic` | +5.856 | 5.26e-09 | 8.41e-08 |

### `Data for reproducing final plots/multivariate_mediation_results.csv`

One row per outcome × feature × contrast (48). Produced by s05, read by s07–s09.

| Column | Description |
|---|---|
| `Outcome` | `Image Accuracy` (the only outcome analysed) |
| `Contrast`, `Feature` | As above |
| `N_trials`, `N_subjects`, `Joint_singular` | Fit diagnostics for the joint model |
| `PathB_estimate`, `PathB_p`, `PathB_p_FDR` | Mediator → outcome effect from the joint binomial GLMM |
| `ACME_mean`, `ACME_sd` | Quasi-Bayesian ACME posterior mean and SD |
| `ACME_CoV` | **Normalised ACME** = `ACME_mean / ACME_sd`. Plotted in s07–s09 |
| `QB_p`, `QB_p_FDR` | Quasi-Bayesian ACME *p*, uncorrected and FDR-adjusted |
| `Perm_p`, `Perm_p_FDR` | Permutation (surrogate) ACME *p*, uncorrected and FDR-adjusted |
| `Direct_effect`, `Prop_mediated` | ADE and proportion mediated |

### Trial-level tables (produced by s02 / s03, not shipped)

`EEG_Features_Behaviour_RelBase.csv` (51 features) and
`Results_SelectedFeatures_RelBase.csv` (16 representative features) share the
same metadata columns; one row per analysed trial.

| Column | Description |
|---|---|
| `global_trial_idx` | Running trial index across the whole dataset |
| `trial_idx` | Trial index within subject × condition |
| `subject` | Subject identifier |
| `condition` | `SoundOnly` \| `TICue` \| `TIB4Cue` |
| `word_acc_post` | Old/New recognition judgement correct (1) or incorrect (0) |
| `imag_acc_post` | **Outcome variable.** Image-category (face vs object) response correct (1) or incorrect (0). Set to 0 where `word_acc_post == 0`, since recollection cannot be expressed for an unrecognised item |
| `old_new` | Trial type; all analysed trials are `Old` (foils discarded) |
| `word` | Word delivered on this trial, used to link EEG to behaviour |
| `block_post` | Post-nap test block |
| `NewManscore1` | Binarised manual recollection score (threshold ≥ 2), zeroed where `word_acc_post == 0`. Not used in the reported models |
| `<region>_<feature>` | Trial-level EEG feature, post-stimulus minus baseline |

**Regions:** `frontal`, `central`, `left_temporal`, `right_temporal`,
`posterior` (parietal), `occipital`, plus `frontal_central` for the
cross-channel coupling features.

**Features:** `spindle_power`, `spindle_peak_freq`, `theta_power`,
`sigma_burst_amp`, `so_spindle`, `delta_spindle`, `so_delta_diff`, `aperiodic`.
The `frontal_central` set contains only the three coupling features.

Channel membership per region is defined at the top of `s02`. Trials whose
features could not be estimated (artefact or failed fit) carry NaN and are
dropped by the model-fitting functions.

The 16 features retained by s03 are: `frontal_delta_spindle`,
`central_so_spindle`, `central_delta_spindle`, `central_aperiodic`,
`left_temporal_spindle_peak_freq`, `left_temporal_so_spindle`,
`left_temporal_delta_spindle`, `right_temporal_spindle_peak_freq`,
`right_temporal_so_spindle`, `right_temporal_delta_spindle`,
`posterior_spindle_peak_freq`, `posterior_theta_power`,
`posterior_sigma_burst_amp`, `posterior_delta_spindle`,
`occipital_delta_spindle`, `frontal_central_so_spindle`.

---

## Requirements

**MATLAB R2024b** (tested up to R2025b)
Signal Processing Toolbox (`pmtm`, `butter`, `filtfilt`, `hilbert`) and
Statistics and Machine Learning Toolbox (`linkage`, `cluster`,
`optimalleaforder`, `corr`, `ttest`).

**R 4.6.1** (x86_64-pc-linux-gnu, Ubuntu 22.04.5 LTS), with:

| Package | Version |
|---------|---------|
| `lme4` | 2.0.1 |
| `lmerTest` | 3.2.1 |
| `mediation` | 4.5.1 |
| `dplyr` | 1.2.1 |

**Python 3.10** with [FOOOF](https://fooof-tools.github.io/fooof/), called from
MATLAB for the aperiodic exponent. Note that FOOOF was renamed `specparam` at
version 2.0 with API changes; the version above is required to reproduce the
aperiodic feature exactly. **Please follow their instructions to install this
toolbox properly to call their functions in MATLAB.** s02 obtains the
interpreter from `pyenv` and passes it to `extract_aperiodic_component`; set
`iflinux` at the top of s02 to match your platform.

**Bundled in `Functions/`** (no separate download needed):

- `Functions/Features/extract_aperiodic_component.m` — FOOOF wrapper called by s02
- `Functions/Features/find_Spectral_Slope.m` — multitaper PSD + aperiodic fit
  over a channel set. Companion helper; the shipped s02 calls
  `extract_aperiodic_component` directly and does not use it
- `Functions/chanlocs.mat` — channel locations, must match the montage used in s01
- `Functions/Plotting/Violinplot-Matlab-master/` —
  [Violinplot-Matlab](https://github.com/bastibe/Violinplot-Matlab), used by s07

**Raw data**, sourced from the BIDS dataset for s01:
- `<BIDS_ROOT>/sub-XX/eeg/sub-XX_task-nap_eeg.set`
- `<BIDS_ROOT>/derivatives/artifact_rejection/sub-XX/eeg/sub-XX_task-nap_desc-badsamp_eeg.mat`
- `<BIDS_ROOT>/sub-XX/eeg/sub-XX_task-nap_events.tsv`
- `<BIDS_ROOT>/derivatives/behaviour/desc-postNapBehaviour_trialdata.xlsx`

Set `BIDS_ROOT` at the top of `s01_load_nap_data.m` to your local BIDS dataset
path. **Note that all the EEG data (the `.set` files) has been pre-processed
as specified in the manuscript, at this stage.**

---

## Reproducibility

Random seeds are fixed in both R scripts (`set.seed(2026)`; per-permutation
seeds `PERM_SEED + p` in s05), so the quasi-Bayesian ACME estimates and the
permutation null distributions are exactly reproducible on the versions listed
above.

---

## Runtime note

s05 is the expensive step: 1000 permutations × 3 contrasts × 16 features, with
the full path-a + joint path-b + mediation pipeline refitted per permutation.
Expect several hours on a single core. The permutation loop is embarrassingly
parallel and can be replaced with `parallel::mclapply` for a large speedup.
Setting `RUN_SURROGATE <- FALSE` skips it and returns quasi-Bayesian p-values
only, in minutes.
