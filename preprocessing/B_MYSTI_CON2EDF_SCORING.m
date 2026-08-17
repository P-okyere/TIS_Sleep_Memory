% B_MYSTI_CON2EDF_SCORING.m
% Stage 2: Re-reference to Cz, derive bipolar scoring channels, export EDF.
% MYSTI pipeline - Prince Okyere, NSN Lab / KCL
% Last updated: 2026-08-13
%
% Expected folder layout:
%   <data_root>/sub-<subID>/sub-P<subID>/ses-<sesID>S001/eeg/
% Loops over every session folder matching ses_pattern within the subject.

clear all; close all;

%% ==== EDIT FOR YOUR ENVIRONMENT ====
eeglab_path = '/path/to/eeglab';
data_root   = '/path/to/MYSTI/EEG';

sub_id      = '020';          % subject ID, e.g. '020'
ses_pattern = 'ses-NapS001*'; % session folder(s) to process for this subject

recording_start = [2024 09 05 16 30 00.000]; % [Y M D H M S] - EDIT per subject/session
% ====================================

addpath(eeglab_path);

Folderpath = fullfile(data_root, ['sub-', sub_id], ['sub-P', sub_id]);
cd(Folderpath);
ListFiles = dir(ses_pattern);

%% EEG processing
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

for ii = 1:length(ListFiles)
    folder = ListFiles(ii).name;

    %% Load EEG data
    file2open = dir(fullfile(folder, 'eeg', '*.edf'));
    if isempty(file2open)
        warning(['No EDF files found in ', folder]);
        continue;
    end
    filename = fullfile(file2open(1).folder, file2open(1).name);
    disp(['Processing file ', filename]);

    EEG = pop_biosig(filename);
    EEG = eeg_checkset(EEG);

    %% Channel indices for this montage
    % A1/A2 are the mastoid channels; their raw label varies depending on
    % signal quality (FT9/TP7 for A1, TP10/TP8 for A2). EDIT A1/A2 below
    % per subject/session based on which channel has the cleanest signal.
    F4 = 29;
    C4 = 24;
    P4 = 19;
    O2 = 18;
    F3 = 3;
    C3 = 8;
    P3 = 14;
    O1 = 16;
    A1 = 9;  % mastoid - FT9 / TP7, make choice depending on signal quality
    A2 = 26; % mastoid - TP10 / TP8, make choice depending on signal quality

    expected_min_chans = max([F4 C4 P4 O2 F3 C3 P3 O1 A1 A2]);
    if EEG.nbchan < expected_min_chans
        error('Channel count mismatch for %s: expected at least %d channels, found %d. Check montage before proceeding.', ...
            folder, expected_min_chans, EEG.nbchan);
    end

    %% Derive bipolar scoring channels
    F4A1 = EEG.data(F4,:) - EEG.data(A1,:);
    C4A1 = EEG.data(C4,:) - EEG.data(A1,:);
    P4A1 = EEG.data(P4,:) - EEG.data(A1,:);
    O2A1 = EEG.data(O2,:) - EEG.data(A1,:);
    F3A2 = EEG.data(F3,:) - EEG.data(A2,:);
    C3A2 = EEG.data(C3,:) - EEG.data(A2,:);
    P3A2 = EEG.data(P3,:) - EEG.data(A2,:);
    O1A2 = EEG.data(O1,:) - EEG.data(A2,:);

    %% Combine into a matrix
    data2 = [F4A1; C4A1; P4A1; O2A1; F3A2; C3A2; P3A2; O1A2];
    data2 = double(data2);

    %% Label the channels
    labels = {'F4A1', 'C4A1', 'P4A1', 'O2A1', 'F3A2', 'C3A2', 'P3A2', 'O1A2'};

    newfolder = fullfile(ListFiles(ii).folder, folder, 'Processed_EEG');
    if ~exist(newfolder, 'dir')
        mkdir(newfolder);
    end

    file2save = fullfile(newfolder, [folder, '_processed.edf']);
    writeeeg_TI(file2save, data2, EEG.srate, 'TYPE', 'EDF', 'Label', labels, 'T0', recording_start);
    disp(['Data saved to: ', file2save]);
end

disp('Processing completed.');
