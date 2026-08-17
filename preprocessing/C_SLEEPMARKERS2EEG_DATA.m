% C_SLEEPMARKERS2EEG_DATA.m
% Stage 3: Import manually scored sleep-profile markers into filtered EEG.
% MYSTI pipeline - Prince Okyere, NSN Lab / KCL
% Last updated: 2026-08-13
%
% Expected folder layout:
%   <deriv_root>/sub-<subID>/sub-P<subID>/ses-<sesID>/eeg/filtered.set
%   <scoring_root>/sub-<subID>/Processed_EEG/Sleep_profile_manual.txt

close all; clear all;

%% ==== EDIT FOR YOUR ENVIRONMENT ====
eeglab_path  = '/path/to/eeglab';
script_path  = '/path/to/MYSTI_SCRIPT';
biosig_path  = '/path/to/Biosig3.8.3';
deriv_root   = '/path/to/MYSTI/EEG';
scoring_root = '/path/to/MYSTI/EEG/MYSTI_SCORING';

sub_id = '05';         % subject ID, e.g. '05'
ses_id = 'test1S001';  % session ID, e.g. 'test1S001'
% ====================================

addpath(eeglab_path);
addpath(script_path);
addpath(genpath(biosig_path));

eeg_filepath = fullfile(deriv_root, ['sub-', sub_id], ['sub-P', sub_id], ['ses-', ses_id], 'eeg');
marker_file  = fullfile(scoring_root, ['sub-', sub_id], 'Processed_EEG', 'Sleep_profile_manual.txt');

%% Initialize EEGLAB and load the EEG dataset
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
EEG = pop_loadset('filename', 'filtered.set', 'filepath', eeg_filepath);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);

%% Read the marker file
fileID = fopen(marker_file, 'r');
if fileID == -1
    error('Failed to open file: %s', marker_file);
end

% Skip the header lines
for i = 1:7
    fgetl(fileID);
end

% Read the data
data = textscan(fileID, '%s %s', 'Delimiter', ';', 'MultipleDelimsAsOne', true);
fclose(fileID);

event_times = data{1};  % timestamps
event_types = data{2};  % event types (sleep stage per 30 s epoch)

disp('Extracted event times and types:');
disp(table(event_times, event_types));

%% Add events at 30 s intervals based on the scored epoch types
interval = 30; % seconds

if ~isfield(EEG, 'event') || isempty(EEG.event)
    EEG.event = struct('type', {}, 'latency', {}, 'duration', {}, 'urevent', {});
end

for i = 1:length(event_types)
    latency = (i - 1) * interval * EEG.srate;
    new_event = struct('type', event_types{i}, 'latency', latency + 1, 'duration', 0, 'urevent', length(EEG.event) + 1);
    EEG.event(end+1) = new_event;
    disp(['Added event: ', event_types{i}, ' at latency (samples): ', num2str(latency)]);
end

EEG = eeg_checkset(EEG, 'eventconsistency');

%% Remove auxiliary channels
channels_to_remove = [65, 66, 67, 68, 69];
EEG = pop_select(EEG, 'nochannel', channels_to_remove);

%% Save the updated .set file
EEG = pop_saveset(EEG, 'filename', 'filtered_with_markers.set', 'filepath', eeg_filepath);

% FIX: the original script reloaded filtered_with_markers.set from a
% DIFFERENT subject/session path right after saving this subject's file -
% a copy-paste leftover that meant the eegplot below showed the wrong
% subject's data. That reload has been removed; EEG already holds what
% was just saved.

%% Plot the EEG data with markers
pop_eegplot(EEG, 1, 1, 1);
