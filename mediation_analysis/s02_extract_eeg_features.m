%% s02 — Trial-level EEG feature extraction (baseline-referenced)
%
% EEG segments: 5 s total = 1 s pre-stimulus + 4 s post-stimulus
%
% Timing (relative to sound onset at t = 0):
%   Baseline  : -1.0 s  to  -0.2 s
%   Post-stim :  0.2 s  to   3.0 s
%
% Eight features per region, each expressed as post-stimulus minus baseline:
%   spindle_power, spindle_peak_freq, theta_power,
%   sigma_burst_amp (mean spindle-band analytic amplitude),
%   so_spindle (MVL), delta_spindle (MVL), so_delta_diff,
%   aperiodic (1/f exponent, [1-30] Hz)
% plus three cross-channel (frontal phase -> central amplitude) coupling
% features, giving 8 x 6 regions + 3 = 51 features per trial.
%
% Filtering strategy: bandpass + Hilbert applied to the FULL 5 s segment,
% phase/amplitude then windowed, to avoid edge artefacts.
%
% Behavioural coding applied here (see Methods):
%   - trials without a valid post-nap word response (foils / New) removed
%   - imag_acc_post set to 0 where the Old/New judgement was incorrect
%
% Input : NapData_AllSubjects.mat (from s01), chanlocs.mat
% Output: EEG_Features_Behaviour_RelBase.csv
%
% Requires: extract_aperiodic_component.m (FOOOF wrapper), Python 3.10 + fooof
%
% Author: Junheng Li

clear all
clc
load NapData_AllSubjects.mat
load chanlocs.mat

%% -----------------------------------------------------------------------
% User-defined parameters
% ------------------------------------------------------------------------
spindle_band      = [12 16];    % Hz
so_band           = [0.5 2];    % Hz
delta_band        = [2 4];      % Hz
theta_band        = [4 8];      % Hz
aperiodic_f_range = [1 30];     % Hz  (for 1/f exponent)

% Timing relative to sound onset (seconds)
pre_stim_dur      = 1.0;        % pre-stimulus period in trial (sound onset at 1 s)
baseline_start_t  = -1.0;       % baseline window start  (re: sound onset)
baseline_end_t    = -0.2;       % baseline window end    (re: sound onset)
poststim_start_t  =  0.2;       % post-stim window start (re: sound onset)
poststim_end_t    =  3.0;       % post-stim window end   (re: sound onset)

% Python / FOOOF settings for aperiodic component extraction
python_directory  = pyenv;
iflinux           = 1;          % 1 if running on Linux, 0 otherwise

%% -----------------------------------------------------------------------
% Define channel regions
% ------------------------------------------------------------------------
channel_regions = struct();

frontal_chs        = {'F1','F2','Fz','Fp1','Fp2','F3','F4','AF7','AF8', ...
                       'F7','F8','F5','F6','AF3','AF4','AFz'};
central_chs        = {'Cz','C1','C2','C3','C4','FC3','FC4','FCz','FC1', ...
                       'FC2','FC5','FC6','C5','C6','CP1','CP2','CP3','CP4', ...
                       'CPz','CP5','CP6'};
left_temporal_chs  = {'FT7','FT9','T7','TP7','TP9'};
right_temporal_chs = {'TP10','FT8','FT10','T8','TP8'};
parietal_chs       = {'P1','P2','P3','P4','P5','P6','Pz','P7','P8'};
occipital_chs      = {'PO3','PO4','POz','O1','O2','Oz','PO7','PO8'};

channel_regions.frontal        = find(ismember({chanlocs.labels}, frontal_chs));
channel_regions.central        = find(ismember({chanlocs.labels}, central_chs));
channel_regions.left_temporal  = find(ismember({chanlocs.labels}, left_temporal_chs));
channel_regions.right_temporal = find(ismember({chanlocs.labels}, right_temporal_chs));
channel_regions.posterior      = find(ismember({chanlocs.labels}, parietal_chs));
channel_regions.occipital      = find(ismember({chanlocs.labels}, occipital_chs));

% Cross-channel coupling: frontal phase --> central spindle amplitude
crossfrontal_chs = frontal_chs;
crosscentral_chs = central_chs;
cross_channel_frontal = find(ismember({chanlocs.labels}, crossfrontal_chs));
cross_channel_central = find(ismember({chanlocs.labels}, crosscentral_chs));

%% -----------------------------------------------------------------------
% Main analysis loop
% ------------------------------------------------------------------------
results_rel_base = table();   % Baseline-subtracted features
trial_counter    = 0;

fprintf('Starting EEG regional analysis for %d subjects...\n', height(tbl_storage));

conditions   = {'SoundOnly', 'TICue', 'TIB4Cue'};
eeg_fields   = {'EEG_soundonly',    'EEG_TICue',    'EEG_TIB4Cue'};
art_fields   = {'Artflag_soundonly','Artflag_TICue', 'Artflag_TIB4Cue'};
tinfo_fields = {'Tinfo_SO',         'Tinfo_TICue',   'Tinfo_TIB4Cue'};

for subj_idx = 1:height(tbl_storage)

    fprintf('Processing subject %d/%d (ID: %d)\n', subj_idx, height(tbl_storage), ...
        tbl_storage.Sbjidx(subj_idx));
    current_subj = tbl_storage.Sbjidx(subj_idx);

    for cond_idx = 1:length(conditions)

        cond_name = conditions{cond_idx};
        eeg_data  = tbl_storage.(eeg_fields{cond_idx}){subj_idx};
        artflag   = tbl_storage.(art_fields{cond_idx}){subj_idx};
        tinfo     = tbl_storage.(tinfo_fields{cond_idx}){subj_idx};

        if isempty(eeg_data)
            continue;
        end

        % TIB4Cue: extract the last 5 s (1s pre-cue + 4s cue period) from 8s trials
        if strcmp(cond_name, 'TIB4Cue')
            total_samps = size(eeg_data, 3);
            samps_5s    = 5 * Fs;
            eeg_data    = eeg_data(:, :, total_samps - samps_5s + 1 : end);
        end

        base_res = analyze_eeg_trials_regional( ...
            eeg_data, artflag, ...
            channel_regions, cross_channel_frontal, cross_channel_central, ...
            Fs, spindle_band, so_band, delta_band, theta_band, ...
            pre_stim_dur, baseline_start_t, baseline_end_t, ...
            poststim_start_t, poststim_end_t, ...
            aperiodic_f_range, python_directory, iflinux, trial_counter);

        % Attach metadata
        n_trials = height(base_res);
        base_res.subject   = repmat(current_subj, n_trials, 1);
        base_res.condition = repmat({cond_name},  n_trials, 1);

        % Behavioural variables
        if ~isempty(tinfo) && height(tinfo) == n_trials
            base_res.word_acc_post = tinfo.WordAcc_Post;
            base_res.imag_acc_post = tinfo.ImagAcc_Post;
            base_res.old_new       = tinfo.old_new;
            base_res.word          = tinfo.word;
            base_res.block_post    = tinfo.Block_Post;
            base_res.Manscore      = tinfo.Manual_Score;
        else
            base_res.word_acc_post = NaN(n_trials, 1);
            base_res.imag_acc_post = NaN(n_trials, 1);
            base_res.old_new       = repmat({''}, n_trials, 1);
            base_res.word          = repmat({''}, n_trials, 1);
            base_res.block_post    = NaN(n_trials, 1);
            base_res.Manscore      = NaN(n_trials, 1);
        end

        results_rel_base = [results_rel_base; base_res];
        trial_counter    = trial_counter + n_trials;
    end
end

%% -----------------------------------------------------------------------
% Recode behavioural scores and select analysed trials
% ------------------------------------------------------------------------

% Binarised manual recollection score, zeroed on failed word retrieval
manscorenow = results_rel_base.Manscore;
manscorenow(manscorenow < 2)  = 0;
manscorenow(manscorenow >= 2) = 1;
manscorenow(results_rel_base.word_acc_post == 0) = 0;
results_rel_base = addvars(results_rel_base, manscorenow, 'NewVariableNames', 'NewManscore1');

% Keep only trials with a valid post-nap word response (drops foil / New trials)
n_before = height(results_rel_base);
results_rel_base(isnan(results_rel_base.word_acc_post), :) = [];
fprintf('Removed %d trials without a valid post-nap response -> %d trials.\n', ...
    n_before - height(results_rel_base), height(results_rel_base));

% Incorrect Old/New judgement -> associative memory failure (0), not missing
n_zeroed = sum(results_rel_base.word_acc_post == 0);
results_rel_base.imag_acc_post(results_rel_base.word_acc_post == 0) = 0;
fprintf('Coded %d trials with an incorrect Old/New judgement as memory failures.\n', n_zeroed);

results_rel_base = removevars(results_rel_base, {'Manscore'});

%% Saving data

writetable(results_rel_base, 'EEG_Features_Behaviour_RelBase.csv');

fprintf('Done. Saved EEG_Features_Behaviour_RelBase.csv (%d rows).\n', ...
    height(results_rel_base));


%% ========================================================================
%                        FUNCTION DEFINITIONS
%% ========================================================================

% -------------------------------------------------------------------------
function base_results = analyze_eeg_trials_regional( ...
    eeg_data, artifact_flags, ...
    channel_regions, cross_frontal_channels, cross_central_channels, ...
    Fs, spindle_band, so_band, delta_band, theta_band, ...
    pre_stim_dur, baseline_start_t, baseline_end_t, ...
    poststim_start_t, poststim_end_t, ...
    aperiodic_f_range, python_directory, iflinux, trial_offset)
%
% Returns base_results: relative features using baseline subtraction.
%
% Timing (all relative to sound onset):
%   Baseline  : baseline_start_t  to  baseline_end_t   (e.g. -1.0 to -0.2 s)
%   Post-stim : poststim_start_t  to  poststim_end_t   (e.g.  0.2 to  3.0 s)
%
% Filtering strategy:
%   - Bandpass filter and Hilbert transform applied to FULL 5s segment
%   - Phase and amplitude then windowed into baseline / post-stim
%   - This avoids edge artifacts from filtering short segments

    num_trials = size(eeg_data, 1);
    regions    = fieldnames(channel_regions);

    % ---- Derive sample indices from timing parameters -------------------
    % Sound onset sits at sample:  snd = round(pre_stim_dur * Fs) + 1
    snd = round(pre_stim_dur * Fs) + 1;

    base_start = snd + round(baseline_start_t  * Fs);     % inclusive
    base_end   = snd + round(baseline_end_t    * Fs) - 1; % inclusive
    post_start = snd + round(poststim_start_t  * Fs);     % inclusive
    post_end   = snd + round(poststim_end_t    * Fs) - 1; % inclusive

    base_idx = base_start : base_end;
    post_idx = post_start : post_end;

    base_results = init_results_table(num_trials, regions, trial_offset);

    for trial = 1:num_trials

        if artifact_flags(trial) == 1
            continue;
        end

        % Extract [channels x time] for this trial
        trial_data = squeeze(eeg_data(trial, :, :));
        if size(eeg_data, 2) == 1
            trial_data = trial_data(:)';
        end

        total_samples = size(trial_data, 2);

        % Guard: require the full post-stim window to be present
        if total_samples < post_end
            continue;
        end

        % ------------------------------------------------------------------
        % Pre-filter FULL 5s segment for all needed bands
        % ------------------------------------------------------------------
        nyquist = Fs / 2;
        [b_so,  a_so]  = butter(4, so_band      / nyquist, 'bandpass');
        [b_del, a_del] = butter(4, delta_band   / nyquist, 'bandpass');
        [b_sp,  a_sp]  = butter(4, spindle_band / nyquist, 'bandpass');

        % Pre-compute filtered + analytic signals for ALL channels on full 5s
        n_channels = size(trial_data, 1);
        so_phase_full    = NaN(n_channels, total_samples);
        delta_phase_full = NaN(n_channels, total_samples);
        sp_amp_full      = NaN(n_channels, total_samples);

        for ch = 1:n_channels
            sig = double(trial_data(ch, :));
            if any(isnan(sig)); continue; end
            try
                so_filt    = filtfilt(b_so,  a_so,  sig);
                delta_filt = filtfilt(b_del, a_del, sig);
                sp_filt    = filtfilt(b_sp,  a_sp,  sig);

                so_phase_full(ch, :)    = angle(hilbert(so_filt));
                delta_phase_full(ch, :) = angle(hilbert(delta_filt));
                sp_amp_full(ch, :)      = abs(hilbert(sp_filt));
            catch
                % leave as NaN
            end
        end

        % ------------------------------------------------------------------
        % Regional analysis
        % ------------------------------------------------------------------
        for r = 1:length(regions)

            region     = regions{r};
            region_chs = channel_regions.(region);
            if isempty(region_chs); continue; end

            valid_chs = region_chs(region_chs <= n_channels);
            if isempty(valid_chs); continue; end

            num_chs = length(valid_chs);

            % --- Spectral features (computed independently per window) ---
            post_sp_power  = NaN(num_chs, 1);
            post_peak_freq = NaN(num_chs, 1);
            base_sp_power  = NaN(num_chs, 1);
            base_peak_freq = NaN(num_chs, 1);
            post_theta_pow = NaN(num_chs, 1);
            base_theta_pow = NaN(num_chs, 1);
            post_aperiodic = NaN(num_chs, 1);
            base_aperiodic = NaN(num_chs, 1);

            % Sigma burst amplitude: mean of spindle-band analytic envelope
            % Uses pre-filtered full-5s sp_amp_full, windowed
            post_sigma_amp = NaN(num_chs, 1);
            base_sigma_amp = NaN(num_chs, 1);

            for ch = 1:num_chs
                ch_idx = valid_chs(ch);
                p_sig  = double(trial_data(ch_idx, post_idx));
                b_sig  = double(trial_data(ch_idx, base_idx));
                if any(isnan(p_sig)) || any(isnan(b_sig)); continue; end

                try
                    [post_sp_power(ch), post_peak_freq(ch)] = ...
                        analyze_spindle_power(p_sig, Fs, spindle_band);
                catch; end
                try
                    [base_sp_power(ch), base_peak_freq(ch)] = ...
                        analyze_spindle_power(b_sig, Fs, spindle_band);
                catch; end
                try
                    post_theta_pow(ch) = analyze_band_power(p_sig, Fs, theta_band);
                catch; end
                try
                    base_theta_pow(ch) = analyze_band_power(b_sig, Fs, theta_band);
                catch; end
                try
                    post_aperiodic(ch) = compute_aperiodic_exponent( ...
                        p_sig, Fs, aperiodic_f_range, python_directory, iflinux);
                catch; end
                try
                    base_aperiodic(ch) = compute_aperiodic_exponent( ...
                        b_sig, Fs, aperiodic_f_range, python_directory, iflinux);
                catch; end

                % Sigma burst amplitude from pre-filtered envelope
                if ~any(isnan(sp_amp_full(ch_idx, :)))
                    post_sigma_amp(ch) = mean(sp_amp_full(ch_idx, post_idx));
                    base_sigma_amp(ch) = mean(sp_amp_full(ch_idx, base_idx));
                end
            end

            % --- PAC features (from pre-filtered full-5s signals) ---
            post_so_mvl    = NaN(num_chs, 1);
            post_delta_mvl = NaN(num_chs, 1);
            base_so_mvl    = NaN(num_chs, 1);
            base_delta_mvl = NaN(num_chs, 1);

            for ch = 1:num_chs
                ch_idx = valid_chs(ch);
                if any(isnan(so_phase_full(ch_idx, :))); continue; end

                % Extract windowed phase & amplitude from full-5s signals
                so_ph_post    = so_phase_full(ch_idx, post_idx);
                delta_ph_post = delta_phase_full(ch_idx, post_idx);
                sp_amp_post   = sp_amp_full(ch_idx, post_idx);

                so_ph_base    = so_phase_full(ch_idx, base_idx);
                delta_ph_base = delta_phase_full(ch_idx, base_idx);
                sp_amp_base   = sp_amp_full(ch_idx, base_idx);

                % Observed MVL on post-stim and baseline windows
                try
                    post_so_mvl(ch)    = canolty_MVL(so_ph_post,    sp_amp_post);
                    post_delta_mvl(ch) = canolty_MVL(delta_ph_post, sp_amp_post);
                    base_so_mvl(ch)    = canolty_MVL(so_ph_base,    sp_amp_base);
                    base_delta_mvl(ch) = canolty_MVL(delta_ph_base, sp_amp_base);
                catch; end
            end

            % --- Pairwise-complete check for baseline subtraction ---
            sp_power_valid  = ~isnan(post_sp_power)  & ~isnan(base_sp_power);
            peak_freq_valid = ~isnan(post_peak_freq) & ~isnan(base_peak_freq);
            theta_valid     = ~isnan(post_theta_pow) & ~isnan(base_theta_pow);
            sigma_valid     = ~isnan(post_sigma_amp) & ~isnan(base_sigma_amp);
            so_mvl_valid    = ~isnan(post_so_mvl)    & ~isnan(base_so_mvl);
            delta_mvl_valid = ~isnan(post_delta_mvl) & ~isnan(base_delta_mvl);

            if any(sp_power_valid)
                base_results.([region '_spindle_power'])(trial) = ...
                    mean(post_sp_power(sp_power_valid) - base_sp_power(sp_power_valid));
            end
            if any(peak_freq_valid)
                base_results.([region '_spindle_peak_freq'])(trial) = ...
                    mean(post_peak_freq(peak_freq_valid) - base_peak_freq(peak_freq_valid));
            end
            if any(theta_valid)
                base_results.([region '_theta_power'])(trial) = ...
                    mean(post_theta_pow(theta_valid) - base_theta_pow(theta_valid));
            end
            if any(sigma_valid)
                base_results.([region '_sigma_burst_amp'])(trial) = ...
                    mean(post_sigma_amp(sigma_valid) - base_sigma_amp(sigma_valid));
            end
            if any(so_mvl_valid)
                base_results.([region '_so_spindle'])(trial) = ...
                    mean(post_so_mvl(so_mvl_valid) - base_so_mvl(so_mvl_valid));
            end
            if any(delta_mvl_valid)
                base_results.([region '_delta_spindle'])(trial) = ...
                    mean(post_delta_mvl(delta_mvl_valid) - base_delta_mvl(delta_mvl_valid));
            end
            if any(so_mvl_valid & delta_mvl_valid)
                valid_both = so_mvl_valid & delta_mvl_valid;
                post_diff  = post_so_mvl(valid_both) - post_delta_mvl(valid_both);
                base_diff  = base_so_mvl(valid_both) - base_delta_mvl(valid_both);
                base_results.([region '_so_delta_diff'])(trial) = mean(post_diff - base_diff);
            end
            base_results.([region '_aperiodic'])(trial) = ...
                mean(post_aperiodic - base_aperiodic, 'omitnan');
        end

        % ------------------------------------------------------------------
        % Cross-channel coupling: frontal phase --> central spindle amplitude
        % ------------------------------------------------------------------
        try
            max_ch        = n_channels;
            valid_frontal = cross_frontal_channels(cross_frontal_channels <= max_ch);
            valid_central = cross_central_channels(cross_central_channels <= max_ch);

            if ~isempty(valid_frontal) && ~isempty(valid_central)

                [post_so_xc, post_delta_xc] = compute_cross_channel_mvl( ...
                    so_phase_full, delta_phase_full, sp_amp_full, ...
                    valid_frontal, valid_central, post_idx);
                [base_so_xc, base_delta_xc] = compute_cross_channel_mvl( ...
                    so_phase_full, delta_phase_full, sp_amp_full, ...
                    valid_frontal, valid_central, base_idx);

                base_results.frontal_central_so_spindle(trial)    = post_so_xc    - base_so_xc;
                base_results.frontal_central_delta_spindle(trial) = post_delta_xc - base_delta_xc;
                base_results.frontal_central_so_delta_diff(trial) = ...
                    (post_so_xc - post_delta_xc) - (base_so_xc - base_delta_xc);
            end
        catch ME
            fprintf('Warning: Cross-channel coupling failed for trial %d: %s\n', ...
                trial, ME.message);
        end
    end
end


% -------------------------------------------------------------------------
function tbl = init_results_table(num_trials, regions, trial_offset)
    tbl = table();
    tbl.global_trial_idx = trial_offset + (1:num_trials)';
    tbl.trial_idx        = (1:num_trials)';

    for r = 1:length(regions)
        rg = regions{r};
        tbl.([rg '_spindle_power'])    = NaN(num_trials, 1);
        tbl.([rg '_spindle_peak_freq'])= NaN(num_trials, 1);
        tbl.([rg '_theta_power'])      = NaN(num_trials, 1);
        tbl.([rg '_sigma_burst_amp'])  = NaN(num_trials, 1);
        tbl.([rg '_so_spindle'])       = NaN(num_trials, 1);
        tbl.([rg '_delta_spindle'])    = NaN(num_trials, 1);
        tbl.([rg '_so_delta_diff'])    = NaN(num_trials, 1);
        tbl.([rg '_aperiodic'])        = NaN(num_trials, 1);
    end

    tbl.frontal_central_so_spindle    = NaN(num_trials, 1);
    tbl.frontal_central_delta_spindle = NaN(num_trials, 1);
    tbl.frontal_central_so_delta_diff = NaN(num_trials, 1);
end


function [spindle_power, peak_freq] = analyze_spindle_power(signal, Fs, spindle_band)
    signal   = [diff(signal), 0];
    nfft     = 2^nextpow2(length(signal));
    [psd, f] = pmtm(signal, 4, nfft, Fs);
    psd_norm = psd / sum(psd);

    idx = f >= spindle_band(1) & f <= spindle_band(2);
    if ~any(idx)
        spindle_power = NaN; peak_freq = NaN; return;
    end

    spindle_power = 10*log10(sum(psd_norm(idx)));
    [~, mi]       = max(psd_norm(idx));
    freqs_in_band = f(idx);
    peak_freq     = freqs_in_band(mi);
end


function band_power = analyze_band_power(signal, Fs, freq_band)
    signal   = [diff(signal), 0];
    nfft     = 2^nextpow2(length(signal));
    [psd, f] = pmtm(signal, 4, nfft, Fs);
    psd_norm = psd / sum(psd);

    idx = f >= freq_band(1) & f <= freq_band(2);
    if ~any(idx)
        band_power = NaN; return;
    end

    band_power = 10*log10(sum(psd_norm(idx)));
end


function exponent = compute_aperiodic_exponent(signal, Fs, f_range, ...
    python_directory, iflinux)
    warning off
    signal   = detrend(signal);
    nfft     = 2^nextpow2(length(signal));
    [psd, f] = pmtm(signal, 4, nfft, Fs);
    exponent = extract_aperiodic_component(psd, f, f_range, 0, ...
        python_directory, iflinux);
end


function MVL = canolty_MVL(phase_low, amp_high)
    z   = amp_high .* exp(1i * phase_low);
    MVL = abs(mean(z)) / mean(amp_high);
end


function [so_mvl, delta_mvl] = compute_cross_channel_mvl( ...
    so_phase_full, delta_phase_full, sp_amp_full, ...
    frontal_chs, central_chs, time_idx)

    n_pairs     = length(frontal_chs) * length(central_chs);
    so_pairs    = NaN(n_pairs, 1);
    delta_pairs = NaN(n_pairs, 1);
    pair_idx    = 0;

    for fi = 1:length(frontal_chs)
        f_so_ph    = so_phase_full(frontal_chs(fi), time_idx);
        f_delta_ph = delta_phase_full(frontal_chs(fi), time_idx);
        if any(isnan(f_so_ph)); continue; end

        for ci = 1:length(central_chs)
            c_amp = sp_amp_full(central_chs(ci), time_idx);
            if any(isnan(c_amp)); continue; end

            pair_idx = pair_idx + 1;
            so_pairs(pair_idx)    = canolty_MVL(f_so_ph,    c_amp);
            delta_pairs(pair_idx) = canolty_MVL(f_delta_ph, c_amp);
        end
    end

    so_mvl    = mean(so_pairs,    'omitnan');
    delta_mvl = mean(delta_pairs, 'omitnan');
end
