"""
YASA_SO_Spindle_Detection.py

Detect slow oscillations (SOs) and sleep spindles (slow, fast, and 10-16 Hz
bands) via YASA, per channel, across the full continuous recording, for
all subjects.

REQUIRED INPUT (per subject):
    <bids_root>/<participant>/eeg/<participant>_task-nap_eeg.set

OUTPUT:
    <results_root>/<participant>/<channel>_SOs.csv
    <results_root>/<participant>/<channel>_spindles_<band>.csv
    <results_root>/<participant>/<participant>_all_channels_SOs.csv
    <results_root>/<participant>/<participant>_all_channels_spindles.csv
    <results_root>/all_participants_SOs.csv
    <results_root>/all_participants_spindles.csv
    plus one average-event plot per channel/band

MYSTI pipeline - Prince Okyere, NSN Lab / KCL
Last updated: 2026-08-13
"""

import os
import mne
import yasa
import pandas as pd
import matplotlib.pyplot as plt

# ==== EDIT FOR YOUR ENVIRONMENT ====
bids_root = "/path/to/MYSTI/NAP_BIDS"
results_root = os.path.join(bids_root, "Results_SO_Spindle")

# n=28
participants = [
    'sub-02', 'sub-03', 'sub-04', 'sub-05', 'sub-06', 'sub-08', 'sub-12',
    'sub-13', 'sub-14', 'sub-18', 'sub-19', 'sub-20', 'sub-21', 'sub-22',
    'sub-23', 'sub-24', 'sub-25', 'sub-26', 'sub-27', 'sub-29', 'sub-30',
    'sub-33', 'sub-36', 'sub-37', 'sub-38', 'sub-39', 'sub-40', 'sub-41'
]
# ====================================

# SO detection parameters
SO_PARAMS = dict(coupling=True, amp_neg=[15, 200], amp_pos=[15, 150], amp_ptp=[30, 350])

# Spindle detection parameters
SPINDLE_DURATION = [0.3, 3]  # seconds
SPINDLE_BANDS = {
    "slow": [10, 12],
    "fast": [12, 16],
    "all": [10, 16],
}

os.makedirs(results_root, exist_ok=True)

master_so_df = pd.DataFrame()
master_spindle_df = pd.DataFrame()

for participant in participants:
    print(f"\n=== {participant} ===")
    eeg_dir = os.path.join(bids_root, participant, "eeg")
    set_file = os.path.join(eeg_dir, f"{participant}_task-nap_eeg.set")

    if not os.path.exists(set_file):
        print(f"  Missing: {set_file} -- skipping")
        continue

    subj_results_dir = os.path.join(results_root, participant)
    os.makedirs(subj_results_dir, exist_ok=True)

    raw = mne.io.read_raw_eeglab(set_file, preload=True)
    sampling_rate = raw.info["sfreq"]
    channels = raw.ch_names

    subject_sos = pd.DataFrame()
    subject_spindles = pd.DataFrame()

    for ch in channels:
        print(f"  Channel: {ch}")
        data = raw.copy().pick([ch]).get_data()[0] * 1e6  # convert to uV

        # --- SO detection ---
        sw = yasa.sw_detect(data, sf=sampling_rate, **SO_PARAMS)
        if sw is not None and not sw.summary().empty:
            swdf = sw.summary()
            swdf["Channel"] = ch
            subject_sos = pd.concat([subject_sos, swdf], ignore_index=True)

            swdf.to_csv(os.path.join(subj_results_dir, f"{ch}_SOs.csv"), index=False)

            sw.plot_average(palette="Reds_r", ci=95)
            plt.title(f"Average SO - {participant} - {ch}")
            plt.savefig(os.path.join(subj_results_dir, f"Average_SO_{ch}.png"))
            plt.close()
        else:
            print(f"    No SOs detected on {ch}")

        # --- Spindle detection (slow + fast) ---
        for band_name, freq_range in SPINDLE_BANDS.items():
            sp = yasa.spindles_detect(
                data, sf=sampling_rate, ch_names=[ch],
                duration=SPINDLE_DURATION, freq_sp=freq_range
            )
            if sp is None or sp.summary().empty:
                print(f"    No {band_name} spindles on {ch}")
                continue

            spdf = sp.summary()
            spdf["Channel"] = ch
            spdf["SpindleType"] = band_name
            subject_spindles = pd.concat([subject_spindles, spdf], ignore_index=True)

            spdf.to_csv(os.path.join(subj_results_dir, f"{ch}_spindles_{band_name}.csv"), index=False)

            sp.plot_average(ci=None, palette=['tab:grey'])
            plt.title(f"{participant} - Avg {band_name.capitalize()} Spindle - {ch}")
            plt.savefig(os.path.join(subj_results_dir, f"{participant}_Average_{band_name}_Spindle_{ch}.png"))
            plt.close()

    # --- Save per-subject combined files ---
    if not subject_sos.empty:
        subject_sos.to_csv(os.path.join(subj_results_dir, f"{participant}_all_channels_SOs.csv"), index=False)
        master_so_df = pd.concat([master_so_df, subject_sos], ignore_index=True)
        print(f"  Saved combined SOs for {participant}")
    else:
        print(f"  No SOs detected for {participant}")

    if not subject_spindles.empty:
        subject_spindles["Participant"] = participant
        subject_spindles.to_csv(os.path.join(subj_results_dir, f"{participant}_all_channels_spindles.csv"), index=False)
        master_spindle_df = pd.concat([master_spindle_df, subject_spindles], ignore_index=True)
        print(f"  Saved combined spindles for {participant}")
    else:
        print(f"  No spindles detected for {participant}")

# --- Save master files across all participants ---
if not master_so_df.empty:
    master_so_df.to_csv(os.path.join(results_root, "all_participants_SOs.csv"), index=False)
    print(f"\nSaved master SO file: {os.path.join(results_root, 'all_participants_SOs.csv')}")
else:
    print("\nNo SOs detected for any participant.")

if not master_spindle_df.empty:
    master_spindle_df.to_csv(os.path.join(results_root, "all_participants_spindles.csv"), index=False)
    print(f"Saved master spindle file: {os.path.join(results_root, 'all_participants_spindles.csv')}")
else:
    print("No spindles detected for any participant.")

print("\nFinished processing all participants.")
