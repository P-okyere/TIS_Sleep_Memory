function result = MYSTI_SO_Spindle_Coupling(base_path, participant, marker_type)
% MYSTI_SO_SPINDLE_COUPLING  Trial-based local SO-spindle coupling (per channel).
%
%   result = MYSTI_SO_Spindle_Coupling(base_path, participant, marker_type)
%
% For each channel: detects slow oscillations (SOs), takes the spindles
% already detected on that same channel, couples each spindle peak to the
% nearest SO trough within COUPLING.win seconds, and extracts SO phase
% at the spindle peak. Results are saved one .mat file per channel.
%
% INPUTS
%   base_path    - BIDS root folder (e.g. 'E:\BIDS')
%   participant  - subject ID, e.g. 'sub-03'
%   marker_type  - trial/condition label to process, e.g. 'Cue'
%
% REQUIRED INPUT FILES
%   <base_path>/<participant>/eeg/<participant>_task-nap_eeg.set
%   <base_path>/derivatives/spindle_so_detection/<participant>/eeg/
%       <participant>_task-nap_desc-fastSpindleMarkers_events.mat
%       (spindles already detected — this script does not detect them)
%
% OUTPUT
%   <base_path>/derivatives/SO_spindle_trialphase/<marker_type>/ALLCHANNELS_LOCAL/
%       <participant>_<marker_type>_<channel>_LOCAL_SOspindle_TRIALS_PHASE.mat
%
% MYSTI pipeline - Prince Okyere, NSN Lab / KCL
% Last updated: 2026-08-31

fprintf('\nSO-spindle coupling | %s | %s\n', participant, marker_type);

%% Paths
eeg_dir      = fullfile(base_path, participant, 'eeg');
eeg_file     = fullfile(eeg_dir, [participant '_task-nap_eeg.set']);

spindle_dir  = fullfile(base_path, 'derivatives', 'spindle_so_detection', participant, 'eeg');
spindle_file = fullfile(spindle_dir, [participant '_task-nap_desc-fastSpindleMarkers_events.mat']);

if ~exist(eeg_file, 'file'),     error('Missing EEG file: %s', eeg_file); end
if ~exist(spindle_file, 'file'), error('Missing spindle file: %s', spindle_file); end

out_dir = fullfile(base_path, 'derivatives', 'SO_spindle_trialphase', marker_type, 'ALLCHANNELS_LOCAL');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% Parameters
SO.band    = [0.5 2];
SO.order   = 4;
SO.min_dur = 0.8;
SO.max_dur = 2.0;
SO.amp_sd  = 1.25;          % alternative threshold option, see detect_SO_singleChannel
SO.window  = [0.25 3.75];   % SO must fall within this window relative to marker_time (s)

COUPLING.win = 1.2; % spindle peak must be within +/- this many seconds of an SO trough
PHASE.band   = [0.5 2.0];

%% Load data
EEG = pop_loadset('filename', eeg_file);
fs  = EEG.srate;

S = load(spindle_file);
trial_info = S.trial_spindle_info;

labels = {EEG.chanlocs.labels};
n_ch   = numel(labels);
fprintf('Channels: %d | Sampling rate: %.2f Hz\n', n_ch, fs);

trials       = trial_info(strcmp({trial_info.marker_type}, marker_type));
n_trials     = numel(trials);
marker_times = [trials.marker_time];
fprintf('Trials: %d\n', n_trials);

%% Process each channel
result = struct('participant', participant, 'marker_type', marker_type, 'fs', fs, ...
    'params', struct('SO', SO, 'COUPLING', COUPLING, 'PHASE', PHASE), ...
    'channel_results', repmat(struct(), n_ch, 1));

for ch = 1:n_ch
    ch_label = labels{ch};
    fprintf('\n[%d/%d] %s\n', ch, n_ch, ch_label);

    [SO_events, so_meta] = detect_SO_singleChannel(EEG, ch, fs, SO, marker_times);
    fprintf('  SOs detected: %d (threshold: %s, trough >= %.2f uV, t2p >= %.2f uV)\n', ...
        numel(SO_events), so_meta.method, so_meta.thr_trough, so_meta.thr_t2p);

    spindle_events = earliest_spindle_per_trial(trials, ch_label);
    fprintf('  Trials with a spindle on this channel: %d/%d\n', numel(spindle_events), n_trials);

    trial_results = init_trial_results(trials, n_trials, ch_label);
    sp_by_trial = index_by_trial(spindle_events, n_trials);

    if ~isempty(SO_events)
        [~, ord]  = sort([SO_events.trough_time]);
        SO_events = SO_events(ord);
        so_times  = [SO_events.trough_time];
    else
        so_times = [];
    end

    for t = 1:n_trials
        sp = sp_by_trial{t};
        if isempty(sp), continue; end

        trial_results(t).has_spindle       = 1;
        trial_results(t).spindle_peak_time = sp.peak_time;
        trial_results(t).spindle_channel   = sp.channel;

        if isempty(so_times), continue; end

        [min_dist, idx] = min(abs(so_times - sp.peak_time));
        if isempty(idx) || min_dist > COUPLING.win, continue; end
        so = SO_events(idx);

        phi_trough0 = trough_referenced_phase(EEG, so, sp.peak_time, fs, PHASE.band);
        if ~isfinite(phi_trough0), continue; end

        trial_results(t).has_coupling   = 1;
        trial_results(t).relative_time  = sp.peak_time - so.trough_time;
        trial_results(t).so_phase       = wrapToPi(phi_trough0 + pi); % trough=0, upstate=pi
        trial_results(t).so_channel     = so.channel_label;
        trial_results(t).so_trough_time = so.trough_time;
    end

    fprintf('  Coupled: %d/%d\n', sum([trial_results.has_coupling]), n_trials);

    metadata = struct('participant', participant, 'marker_type', marker_type, ...
        'channel', ch_label, 'channel_idx', ch, 'fs', fs, ...
        'SO', SO, 'COUPLING', COUPLING, 'PHASE', PHASE, ...
        'SO_thresholds', so_meta, ...
        'spindle_selection', 'earliest per trial, same channel', ...
        'coupling_rule', sprintf('nearest SO trough within +/-%.2fs', COUPLING.win));

    out_file = fullfile(out_dir, sprintf('%s_%s_%s_LOCAL_SOspindle_TRIALS_PHASE.mat', ...
        participant, marker_type, sanitize_label(ch_label)));
    save(out_file, 'trial_results', 'SO_events', 'spindle_events', 'metadata', '-v7.3');
    fprintf('  Saved: %s\n', out_file);

    result.channel_results(ch) = struct('channel_label', ch_label, 'channel_idx', ch, ...
        'trial_results', trial_results, 'SO_events', SO_events, ...
        'spindle_events', spindle_events, 'metadata', metadata);
end

fprintf('\nDone: %d channel files saved to %s\n', n_ch, out_dir);
end

%% ---------------------------------------------------------------------
function trial_results = init_trial_results(trials, n_trials, ch_label)
trial_results = repmat(struct( ...
    'trial_id', [], 'marker_time', [], 'marker_type', '', 'channel', '', ...
    'has_spindle', 0, 'spindle_peak_time', NaN, 'spindle_channel', '', ...
    'has_coupling', 0, 'relative_time', NaN, ...
    'so_phase', NaN, 'so_channel', '', 'so_trough_time', NaN), n_trials, 1);

for t = 1:n_trials
    trial_results(t).trial_id    = t;
    trial_results(t).marker_time = trials(t).marker_time;
    trial_results(t).marker_type = trials(t).marker_type;
    trial_results(t).channel     = ch_label;
end
end

%% ---------------------------------------------------------------------
function sp_by_trial = index_by_trial(spindle_events, n_trials)
sp_by_trial = cell(n_trials, 1);
for k = 1:numel(spindle_events)
    sp_by_trial{spindle_events(k).trial_idx} = spindle_events(k);
end
end

%% ---------------------------------------------------------------------
function sp_events = earliest_spindle_per_trial(trials, channel_label)
% For each trial, keep only the earliest-peaking spindle on channel_label.
sp_events = struct('trial_idx', {}, 'peak_time', {}, 'channel', {}, 'amplitude', {});
k = 0;

for t = 1:numel(trials)
    tr = trials(t);
    if tr.num_spindles == 0, continue; end

    cand = tr.spindles(strcmp({tr.spindles.channel}, channel_label));
    if isempty(cand), continue; end

    [~, ix] = min([cand.peak]);
    k = k + 1;
    sp_events(k).trial_idx = t;
    sp_events(k).peak_time = cand(ix).peak;
    sp_events(k).channel   = cand(ix).channel;
    if isfield(cand, 'amplitude')
        sp_events(k).amplitude = cand(ix).amplitude;
    else
        sp_events(k).amplitude = NaN;
    end
end
end

%% ---------------------------------------------------------------------
function [SO_events, meta] = detect_SO_singleChannel(EEG, ch, fs, SO, marker_times)
% Detect SOs on one channel: bandpass -> zero-crossing candidate waves ->
% adaptive amplitude thresholds (trough depth, trough-to-peak) from the
% candidate distribution -> keep events exceeding threshold, within
% SO.window of a marker.
labels = {EEG.chanlocs.labels};
[b, a] = butter(SO.order, SO.band / (fs/2), 'bandpass');
x  = filtfilt(b, a, double(EEG.data(ch, :)));
zc = find(diff(sign(x)) ~= 0);

SO_events = struct('channel', {}, 'channel_label', {}, 'trough_time', {}, 'trough_idx', {});
if numel(zc) < 3
    meta = struct('thr_trough', NaN, 'thr_t2p', NaN, 'n_putative', 0, 'method', 'none');
    return
end

% Pass 1: collect candidate events to set adaptive thresholds
all_troughs = [];
all_t2p     = [];
for i = 2:numel(zc)-1
    dur = (zc(i+1) - zc(i-1)) / fs;
    if dur < SO.min_dur || dur > SO.max_dur, continue; end

    seg_idx = zc(i-1):zc(i+1);
    [trough, ix_tr] = min(x(seg_idx));
    trough_idx = seg_idx(1) + ix_tr - 1;
    t_tr = trough_idx / fs;

    if ~any(t_tr >= marker_times + SO.window(1) & t_tr <= marker_times + SO.window(2))
        continue
    end

    peak = max(x(zc(i):zc(i+1)));
    all_troughs(end+1) = abs(trough); %#ok<AGROW>
    all_t2p(end+1)     = peak - trough; %#ok<AGROW>
end

if isempty(all_troughs)
    meta = struct('thr_trough', NaN, 'thr_t2p', NaN, 'n_putative', 0, 'method', 'none');
    return
end

% Threshold: mean * 1.25 (adaptive per channel/session)
thr_trough = mean(all_troughs) * 1.25;
thr_t2p    = mean(all_t2p) * 1.25;
meta = struct('thr_trough', thr_trough, 'thr_t2p', thr_t2p, ...
    'n_putative', numel(all_troughs), 'method', 'mean*1.25');

% Pass 2: keep events exceeding threshold
for i = 2:numel(zc)-1
    dur = (zc(i+1) - zc(i-1)) / fs;
    if dur < SO.min_dur || dur > SO.max_dur, continue; end

    seg_idx = zc(i-1):zc(i+1);
    [trough, ix_tr] = min(x(seg_idx));
    trough_idx = seg_idx(1) + ix_tr - 1;
    t_tr = trough_idx / fs;

    if ~any(t_tr >= marker_times + SO.window(1) & t_tr <= marker_times + SO.window(2))
        continue
    end

    peak = max(x(zc(i):zc(i+1)));
    if abs(trough) < thr_trough || (peak - trough) < thr_t2p, continue; end

    SO_events(end+1) = struct('channel', ch, 'channel_label', labels{ch}, ...
        'trough_time', t_tr, 'trough_idx', trough_idx); %#ok<AGROW>
end
end

%% ---------------------------------------------------------------------
function phi = trough_referenced_phase(EEG, so, spindle_time, fs, band)
% Hilbert phase difference between the spindle peak and the SO trough,
% wrapped to [-pi, pi].
phi = NaN;
[b, a] = butter(2, band / (fs/2), 'bandpass');

pk = round(spindle_time * fs);
L  = max(1, pk - 2*fs);
R  = min(size(EEG.data, 2), pk + 2*fs);

z = hilbert(filtfilt(b, a, double(EEG.data(so.channel, L:R))));
i_pk = pk - L + 1;
i_tr = so.trough_idx - L + 1;
if i_pk < 1 || i_pk > numel(z) || i_tr < 1 || i_tr > numel(z), return; end

phi = angle(exp(1i * (angle(z(i_pk)) - angle(z(i_tr)))));
end

%% ---------------------------------------------------------------------
function s = sanitize_label(label)
% Make a channel label safe for use in a filename.
s = regexprep(label, '\s+', '');
s = regexprep(s, '[^\w\-]', '_');
end
