% D_MYST_ART_CORR.m
% Stage 4: Bad-channel interpolation, average re-referencing, and
% boundary-marked artifact segment removal.
% MYSTI pipeline - Prince Okyere, NSN Lab / KCL
% Last updated: 2026-08-13
%
% Run section-by-section, not all at once - Section "Mark bad
% time segments interactively" pauses for you to mark bad segments in the
% eegplot window before the final section can read them back out as
% boundary events.
%
% Expected inputs:
%   <filepath>/filtered_with_markers.set
%   <badchans_root>/*<sub_tag>_badchans*.xlsx  (columns: chan_idx, chan_label, is_bad)
%   <chanlookup> - a standard channel location file (e.g. standard_1005.elc)

clear all; close all;

%% ==== EDIT FOR YOUR ENVIRONMENT ====
scripts_path  = '/path/to/matlab_scripts';   % folder containing standard_1005.elc
eeglab_path   = '/path/to/eeglab';
chanlookup    = fullfile(scripts_path, 'standard_1005.elc');

data_root     = '/path/to/MYSTI/EEG';
sub_id        = '024';       % subject ID used in folder paths, e.g. '024'
ses_id        = 'NapS001';   % session ID, e.g. 'NapS001'

analysis_sub_id = '03';      % subject ID used in FINAL analysis filenames
                              % (sub-XX_task-nap_eeg.set / _badsamp.mat) -
                              % matches your existing analysis-ready file
                              % naming, which may differ in digit count
                              % from sub_id above. EDIT to match.

badchans_root = data_root;   % folder containing *_badchans*.xlsx
% ====================================

addpath(scripts_path);
addpath(eeglab_path);

sub_tag = ['S', sub_id]; % used for badchans file matching, e.g. *S030_badchans*

filepath     = fullfile(data_root, ['sub-', sub_id], ['sub-P', sub_id], ['ses-', ses_id], 'eeg');
badchans_dir = dir(fullfile(badchans_root, ['*', sub_tag, '_badchans*']));

% Final, analysis-ready output filenames (must match what your
% downstream analysis scripts expect to load):
interp_set_name   = ['sub-', analysis_sub_id, '_task-nap_eeg.set'];              % after interpolation/reref, before bad-segment removal
goodsamps_set_name = ['sub-', analysis_sub_id, '_task-nap_goodsamps_eeg.set'];   % after bad-segment removal
badsamp_mat_name  = ['sub-', analysis_sub_id, '_task-nap_badsamp.mat'];

%% Load EEGLAB and the marker-annotated dataset
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% FIX: originally `EEG2 = EEG; pop_eegplot(EEG2)` ran here BEFORE any
% dataset was loaded, which just plotted an empty EEGLAB structure. The
% dataset is now loaded first, so this plot actually shows the data.
EEG = pop_loadset('filename', 'filtered_with_markers.set', 'filepath', filepath);
EEG = eeg_checkset(EEG);

pop_eegplot(EEG); % inspect and log bad channels in the badchans Excel sheet before continuing

%% Load bad-channel list
if isempty(badchans_dir)
    error('No bad channels file found. Check the file path and name.');
end

% Display all events for inspection (optional)
disp('Displaying all events:');
for i = 1:length(EEG.event)
    fprintf('Event %d: Marker %s, Latency %f\n', i, EEG.event(i).type, EEG.event(i).latency);
end

EEG = pop_saveset(EEG, 'filename', 'filtered_with_markers.set', 'filepath', filepath);

badchans_file = fullfile(badchans_root, badchans_dir(1).name);
if exist(badchans_file, 'file') ~= 2
    error('Bad channels file does not exist. Check the file path and name.');
end

% Using readtable instead of xlsread (deprecated, and less reliable with
% mixed text/number columns) to pull channel index, label, and bad-flag.
badchans_tbl = readtable(badchans_file, 'Sheet', 'Sheet1', 'Range', 'A2:C70', 'ReadVariableNames', false);
badchans_tbl.Properties.VariableNames = {'chan_idx', 'chan_label', 'is_bad'};
badchans_labels = badchans_tbl.chan_label(badchans_tbl.is_bad == 1);

%% Remove auxiliary channel labels, interpolate bad channels, average reference
if length(EEG.chanlocs) >= 18
    EEG.chanlocs(16).labels = 'I1';
    EEG.chanlocs(18).labels = 'I2';
else
    warning('Not enough channels to change labels for indices 16 and 18.');
end

EEG = pop_chanedit(EEG, 'lookup', chanlookup, 'eval', 'chans = pop_chancenter( chans, [],[]);');

originaleeg = EEG; % full channel set, pre-interpolation
for b = 1:length(badchans_labels)
    if ~isempty(badchans_labels{b})
        EEG = pop_select(EEG, 'nochannel', badchans_labels(b));
    else
        warning('Empty bad channel label at index %d.', b);
    end
end

EEG = pop_interp(EEG, originaleeg.chanlocs, 'spherical');
EEG = pop_reref(EEG, []); % average reference

EEG = pop_saveset(EEG, 'filename', interp_set_name, 'filepath', filepath);
disp('EEG preprocessing completed successfully.');

%% Mark bad time segments interactively
% originaleeg here = channel-interpolated, average-referenced continuous
% data (NOT the raw pre-processed data) - kept as the pre-removal
% reference for the final section below.
originaleeg = EEG;

pop_eegplot(EEG, 1, 1, 1);
% --> Mark and REJECT bad segments in the plot window, then register the
%     rejection. EEGLAB inserts a 'boundary' event per rejected segment
%     - the section below reads those back out of EEG.event. Do not run
%     the next section until you've done this.

%% Convert marked boundary events to a bad-sample index and remove
marker    = EEG.event;
samp      = cell2mat({marker.latency});
duration  = {marker.duration};
labels    = {marker.type};
art_marker = find(strcmp(labels, 'boundary') == 1);
art_start  = floor(samp(art_marker));
art_duration = duration(art_marker);

badsamps_all = [];
bad_intervals = []; % Nx2, one row per rejected segment: [start_sample, end_sample]
art_duration2 = 0;

for a = 1:length(art_start)
    if a == 1
        art_start2 = art_start(a);
    else
        art_start2 = art_start(a) + art_duration2;
    end

    art_end2 = art_start2 + art_duration{a};
    badsamps = art_start2:1:art_end2;
    badsamps_all = vertcat(badsamps_all, badsamps');
    bad_intervals = vertcat(bad_intervals, [art_start2, art_end2]);

    art_duration2 = art_duration2 + art_duration{a};
    clear art_start2 art_end2
end

EEG = originaleeg;

goodsamps = 1:size(EEG.data, 2);
EEG.data(:, badsamps_all') = []; % remove bad samples

file2save = goodsamps_set_name;
EEG = pop_saveset(EEG, 'filename', file2save, 'filepath', filepath);

%% Save bad-sample intervals for downstream analysis (epoch-overlap rejection)
% Analysis scripts (MYSTI_WITHIN_SUBJECT_ANALYSIS.m, PSD scripts) load
% this file and check each trial epoch against bad_intervals (Nx2:
% [start_sample, end_sample] per rejected segment) to decide whether to
% skip that trial. Saved under both variable names since different
% analysis scripts load one or the other.
goodsamps(:, badsamps_all') = [];
original_eeg_size = 1:size(originaleeg.data, 2);
goodsamples_eeg = goodsamps;
badsamples_eeg = badsamps_all';
badsamp = bad_intervals; %#ok<NASGU> % alias expected by some analysis scripts

save(fullfile(filepath, badsamp_mat_name), ...
    'original_eeg_size', 'goodsamples_eeg', 'badsamples_eeg', 'bad_intervals', 'badsamp');
