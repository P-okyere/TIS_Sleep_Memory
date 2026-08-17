% A_MYSTI_FILTER.m
% Stage 1: Filter, resample, and export raw EEG (.set for analysis, .edf for scoring).
% MYSTI pipeline - Prince Okyere, NSN Lab / KCL
% Last updated: 2026-08-13
%
% Expected folder layout (edit sub_group derivation below if yours differs):
%   <data_root>/sub-<XX>/sub-P<subID>/ses-<sesID>S001/eeg/

close all; clear all;

%% ==== EDIT FOR YOUR ENVIRONMENT ====
% Paths and subject/session IDs below are placeholders - replace with
% your own before running. This is the only block you should need to
% touch; nothing subject- or machine-specific is hardcoded further down.
eeglab_path = '/path/to/eeglab';
script_path = '/path/to/MYSTI_SCRIPT';
biosig_path = '/path/to/Biosig3.8.3';
data_root   = '/path/to/MYSTI/EEG';

sub_id = '024';       % subject ID, e.g. '024'
ses_id = 'NapS001';   % session ID, e.g. 'NapS001'
% ====================================

addpath(eeglab_path);
addpath(script_path);
addpath(genpath(biosig_path));

%% Build file paths from sub_id/ses_id
% FIX: the outer folder level ('sub-<last two digits>') was previously a
% literal hardcoded value (e.g. 'sub-24') that didn't actually track
% sub_id - changing sub_id to run a different subject silently left this
% folder wrong. It's now derived from sub_id directly; adjust the line
% below if your lab's grouping convention differs.
sub_group = sub_id(end-1:end); % e.g. '24' from '024'

eeg_dir = fullfile(data_root, ['sub-', sub_group], ['sub-P', sub_id], ['ses-', ses_id, 'S001'], 'eeg');

xdf_file        = fullfile(eeg_dir, ['sub-P', sub_id, '_ses-', ses_id, 'S001_task-Default_run-001_eeg.xdf']);
output_set_file = fullfile(eeg_dir, 'filtered.set');
edf_file        = fullfile(eeg_dir, 'filter.edf');

%%
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
EEG = pop_loadxdf(xdf_file, 'streamtype', 'EEG', 'exclude_markerstreams', {});

%% Apply high-pass filter at 0.5 Hz
EEG = pop_eegfiltnew(EEG, 0.5, []);
EEG = eeg_checkset(EEG);

% Apply low-pass filter at 30 Hz
EEG = pop_eegfiltnew(EEG, [], 30);
EEG = eeg_checkset(EEG);

% Apply notch filter at 50 Hz
EEG = pop_eegfiltnew(EEG, 49, 51, [], 1, [], 0);
EEG = eeg_checkset(EEG);

% Downsample the data to 256 Hz
EEG = pop_resample(EEG, 256);
EEG = eeg_checkset(EEG);

%% Save the processed data as a .set file (full channel set, for later EEG analysis)
[set_path, set_name, set_ext] = fileparts(output_set_file);
EEG = pop_saveset(EEG, 'filename', [set_name, set_ext], 'filepath', set_path);

%% Remove auxiliary channels 65-69 (EDF export only - .set above keeps full channel set)
% FIX: comment previously said "Remove channels 32, 33, and 34" but the
% indices removed were [65 66 67 68 69] - comment now matches the code.
channels_to_remove = [65, 66, 67, 68, 69];
EEG_scoring = pop_select(EEG, 'nochannel', channels_to_remove);
EEG_scoring = eeg_checkset(EEG_scoring);

%% Write reduced-channel EDF for sleep scoring
pop_writeeeg(EEG_scoring, edf_file, 'TYPE', 'EDF');
