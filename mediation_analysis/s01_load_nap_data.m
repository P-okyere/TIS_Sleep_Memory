%% s01 — Load nap EEG and behavioural data for all subjects
%
% Reads per-subject EEGLAB .set files, artefact indices and task event
% tables, epochs the stimulation conditions, links each trial to the
% post-nap behavioural response, and stores everything in one table.
%
% Inputs (per subject, from BIDS_ROOT):
%   <BIDS_ROOT>/sub-XX/eeg/sub-XX_task-nap_eeg.set
%   <BIDS_ROOT>/derivatives/artifact_rejection/sub-XX/eeg/sub-XX_task-nap_desc-badsamp_eeg.mat
%   <BIDS_ROOT>/sub-XX/eeg/sub-XX_task-nap_events.tsv
%   <BIDS_ROOT>/derivatives/behaviour/desc-postNapBehaviour_trialdata.xlsx
%
% Output:
%   NapData_AllSubjects.mat   -> tbl_storage, Fs
%
% Conditions:
%   SoundOnly (Cue)  |  TICue (TIS+Cue)  |  TIB4Cue (TIS+Cue_Delayed)
%
% Author: Junheng Li

clear all
clc

%% ============================================================
%  EDIT FOR YOUR ENVIRONMENT
%% ============================================================
BIDS_ROOT = 'E:\BIDS';
% ==============================================================

%% Subject list and constants

% Excluded: s11, s22 (no nap data); s12, s19 (incomplete sleep scoring /
% scoring latency exceeds EEG length)
sbjids = [2,3,4:6,8,12:14,18:25,26,27,29,30,33,36:41];    % Subject indices to use

score_size = 30;    % Sleep scoring epoch length (s)
durstim    = 4;     % Stimulation duration, always 4 s
base_t     = 1;     % Extra pre-stimulus baseline to keep per event (s)

slp_score_labels = {'Wake','N1','N2','N3','REM'};
dict_slpsc       = dictionary(slp_score_labels,[0,1,2,3,5]);

tbl_storage = [];

t_postperf = readtable(fullfile(BIDS_ROOT, 'derivatives', 'behaviour', ...
    'desc-postNapBehaviour_trialdata.xlsx'));    % Post-performance table

basefd = pwd;

%% Loop through subjects

for ns = 1:length(sbjids)

    idxnow = sbjids(ns);
    if idxnow<10
        sbj_str = ['0',num2str(idxnow)];
    else
        sbj_str = [num2str(idxnow)];
    end

    eeg_set_path = fullfile(BIDS_ROOT, ['sub-',sbj_str], 'eeg', ...
        ['sub-',sbj_str,'_task-nap_eeg.set']);
    EEG_nap = load('-mat', eeg_set_path);

    try
        badsamp_path = fullfile(BIDS_ROOT, 'derivatives', 'artifact_rejection', ...
            ['sub-',sbj_str], 'eeg', ['sub-',sbj_str,'_task-nap_desc-badsamp_eeg.mat']);
        load(badsamp_path)      % artefact indices
    catch ME
        bad_intervals = [];
    end

    event_types = {EEG_nap.event.type};

    eeg = double(EEG_nap.data);
    for iii = 1:size(bad_intervals,1)   % Mark bad data as NaN
        eeg(:,bad_intervals(iii,1):bad_intervals(iii,2)) = NaN;
    end

    Fs = EEG_nap.srate;

    eeglen      = length(eeg);
    eeglen_time = floor(eeglen/Fs);            % Length of EEG in time s
    num_scepc   = floor(eeglen_time/score_size);   % Number of sleep stages

    slpsc_idx  = ismember(event_types,slp_score_labels);
    slpsc_this = {event_types{slpsc_idx}};     % Sleep scoring labels
    slpsc_this = {slpsc_this{1:num_scepc}};
    hyp_this   = dict_slpsc(slpsc_this);       % This hypnogram

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read stimulation markers
    % Extract each event individually
    % Record in the meantime, the sleep stage

    eegevents = EEG_nap.event;

    % Sound only (Cue)

    mkrs_sound    = {'MarkerOldSound','MarkerNewSound'};
    idx_soundonly = ismember(event_types,mkrs_sound);
    idx_soundonly = find(idx_soundonly);
    lat_d_soundonly = floor([eegevents(idx_soundonly).latency]);
    % Remove duplicated markers
    diff_latd  = [2*Fs,diff(lat_d_soundonly)];
    redu_makrs = (diff_latd<Fs);
    lat_d_soundonly = lat_d_soundonly(~redu_makrs);
    idx_soundonly   = idx_soundonly(~redu_makrs);

    mkrs_sound_realidx     = {event_types{idx_soundonly}};
    slpstage_soundonly     = floor(lat_d_soundonly/ (score_size*Fs));   % Hypnogram index at which the sound happens
    realslpstage_soundonly = hyp_this(slpstage_soundonly+1);

    % TI With Cue

    mkrs_TICue = {'TICUEE'};
    idx_TICue  = ismember(event_types,mkrs_TICue);
    idx_TICue  = find(idx_TICue);
    lat_d_TICue = floor([eegevents(idx_TICue).latency]);
    % Remove duplicated markers
    diff_latd  = [2*Fs,diff(lat_d_TICue)];
    redu_makrs = (diff_latd<Fs);
    lat_d_TICue = lat_d_TICue(~redu_makrs);
    idx_TICue   = idx_TICue(~redu_makrs);

    slpstage_TICue     = floor(lat_d_TICue/ (score_size*Fs));
    realslpstage_TICue = hyp_this(slpstage_TICue+1);

    % TI Before Cue

    mkrs_TIB4Cue = {'TIB4CUEE'};
    idx_TIB4Cue  = ismember(event_types,mkrs_TIB4Cue);
    idx_TIB4Cue  = find(idx_TIB4Cue);
    lat_d_TIB4Cue = floor([eegevents(idx_TIB4Cue).latency]);
    % Remove duplicated markers
    diff_latd  = [2*Fs,diff(lat_d_TIB4Cue)];
    redu_makrs = (diff_latd<Fs);
    lat_d_TIB4Cue = lat_d_TIB4Cue(~redu_makrs);
    idx_TIB4Cue   = idx_TIB4Cue(~redu_makrs);

    slpstage_TIB4Cue     = floor(lat_d_TIB4Cue/ (score_size*Fs));
    realslpstage_TIB4Cue = hyp_this(slpstage_TIB4Cue+1);

    %%%%%%%%%%%%%%%%%% Extract EEG data for each marker

    durstim_d = durstim*Fs;
    numch     = length(EEG_nap.chanlocs);

    % Add extra EEG baseline to each event
    base_d = floor(base_t*Fs);
    biase  = floor(1*Fs);

    % EEG for sound only
    nume_soundonly   = length(lat_d_soundonly);
    eeg_soundonly    = NaN(nume_soundonly,numch,durstim_d+biase);
    artflag_soundonly = NaN(nume_soundonly,1);
    for ii = 1:nume_soundonly
        eeg_soundonly(ii,:,:) = eeg(:,lat_d_soundonly(ii)+1-base_d:lat_d_soundonly(ii)+durstim_d-base_d + biase);
        artflag_soundonly(ii) = (sum(isnan(eeg_soundonly(ii,1,:)))>0);
        if artflag_soundonly(ii) == 1
            eeg_soundonly(ii,:,:) = NaN;
        end
    end

    % EEG for TICue
    nume_TICue    = length(lat_d_TICue);
    eeg_TICue     = NaN(nume_TICue,numch,durstim_d+biase);
    artflag_TICue = NaN(nume_TICue,1);
    for ii = 1:nume_TICue
        eeg_TICue(ii,:,:) = eeg(:,lat_d_TICue(ii)+1-base_d:lat_d_TICue(ii)+durstim_d-base_d+biase);
        artflag_TICue(ii) = (sum(isnan(eeg_TICue(ii,1,:)))>0);
        if artflag_TICue(ii) == 1
            eeg_TICue(ii,:,:) = NaN;
        end
    end

    % EEG with high-freq before this TI Cue
    nume_TICue      = length(lat_d_TICue);
    eeg_B4TICue     = NaN(nume_TICue,numch,durstim_d);
    artflag_B4TICue = NaN(nume_TICue,1);
    lat_d_B4TICue   = lat_d_TICue - durstim_d+1;
    for ii = 1:nume_TICue
        eeg_B4TICue(ii,:,:) = eeg(:,lat_d_TICue(ii)+1-durstim_d:lat_d_TICue(ii));
        artflag_B4TICue(ii) = (sum(isnan(eeg_B4TICue(ii,1,:)))>0);
        if artflag_B4TICue(ii) == 1
            eeg_B4TICue(ii,:,:) = NaN;
        end
    end
    slpstage_B4TICue     = floor(lat_d_B4TICue/ (score_size*Fs));
    realslpstage_B4TICue = hyp_this(slpstage_B4TICue+1);

    % EEG for TI Before Cue // Taking all 8 s stim
    nume_TIB4Cue    = length(lat_d_TIB4Cue);
    eeg_TIB4Cue     = NaN(nume_TIB4Cue,numch,durstim_d*2);
    artflag_TIB4Cue = NaN(nume_TIB4Cue,1);
    for ii = 1:nume_TIB4Cue
        eeg_TIB4Cue(ii,:,:) = eeg(:,lat_d_TIB4Cue(ii)+1-durstim_d:lat_d_TIB4Cue(ii)+durstim_d);
        artflag_TIB4Cue(ii) = (sum(isnan(eeg_TIB4Cue(ii,1,:)))>0);
        if artflag_TIB4Cue(ii) == 1
            eeg_TIB4Cue(ii,:,:) = NaN;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract periods of sleep where no high-frequency TI was present
    mkrs            = {'empty'};
    idx_empty       = ismember(event_types,mkrs);
    timeafteremtpy  = 16;   % Time to take after empty
    lat_d_empty     = floor([eegevents(idx_empty).latency]);
    slpstage_empty  = floor(lat_d_empty/ (score_size*Fs));
    realslpstage_empty = hyp_this(slpstage_empty+1);

    idxsleep_empty  = find(realslpstage_empty>0);   % Only use those where sleep happens
    idx_empty       = find(idx_empty);
    idx_empty       = idx_empty(idxsleep_empty);
    lat_d_empty     = floor([eegevents(idx_empty).latency]);
    slpstage_empty  = floor(lat_d_empty/ (score_size*Fs));
    realslpstage_empty = hyp_this(slpstage_empty+1);

    % Goes back
    cd(basefd)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract EEG for periods with sleep without TI
    numepcs_empty   = floor(timeafteremtpy/durstim);
    totalepcs_empty = numepcs_empty*length(lat_d_empty);
    eeg_emptysleep  = NaN(totalepcs_empty,numch,durstim_d);
    artflag_emptysleep = NaN(totalepcs_empty,1);
    reallatd_empty  = NaN(totalepcs_empty,1);
    cntee = 0;
    for nemp = 1:length(lat_d_empty)
        for nee_emp = 1:numepcs_empty
            cntee = cntee+1;

            idxstart_nownow = lat_d_empty(nemp)+(nee_emp-1)*durstim_d+1;
            reallatd_empty(cntee)     = idxstart_nownow;
            eeg_emptysleep(cntee,:,:) = eeg(:,idxstart_nownow:idxstart_nownow+durstim_d-1);
            artflag_emptysleep(cntee) = (sum(isnan(eeg_emptysleep(cntee,1,:)))>0);
            if artflag_emptysleep(cntee) == 1
                eeg_emptysleep(cntee,:,:) = NaN;
            end
        end
    end
    % Filter again sleep stages
    idxslpstage_empty_real = floor(reallatd_empty/ (score_size*Fs));
    slpstage_empty_real    = hyp_this(idxslpstage_empty_real+1);
    idx_filt               = (slpstage_empty_real>0);
    slpstage_empty_real    = slpstage_empty_real(idx_filt);
    eeg_emptysleep         = eeg_emptysleep(idx_filt,:,:);
    reallatd_empty         = reallatd_empty(idx_filt);
    artflag_emptysleep     = artflag_emptysleep(idx_filt);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sleep-only periods (4-s epochs, no stimulation)

    mkrs_SleepOnly  = {'NO_STIM'};
    idx_SleepOnly   = ismember(event_types,mkrs_SleepOnly);
    idx_SleepOnly   = find(idx_SleepOnly);
    lat_d_SleepOnly = floor([eegevents(idx_SleepOnly).latency]);

    slpstage_SleepOnly     = floor(lat_d_SleepOnly/ (score_size*Fs));
    realslpstage_SleepOnly = hyp_this(slpstage_SleepOnly+1);

    nume_SleepOnly    = length(lat_d_SleepOnly);
    eeg_SleepOnly     = NaN(nume_SleepOnly,numch,durstim_d);
    artflag_SleepOnly = NaN(nume_SleepOnly,1);
    for ii = 1:nume_SleepOnly
        eeg_SleepOnly(ii,:,:) = eeg(:,lat_d_SleepOnly(ii):lat_d_SleepOnly(ii)+durstim_d-1);
        artflag_SleepOnly(ii) = (sum(isnan(eeg_SleepOnly(ii,1,:)))>0);
        if artflag_SleepOnly(ii) == 1
            eeg_SleepOnly(ii,:,:) = NaN;
        end
    end
    if sum(realslpstage_SleepOnly<2)>0 || sum(realslpstage_SleepOnly>3)>0
        error(['Check Subject ', sbj_str])
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Coarse behavioural performance reading
    allparts       = t_postperf.participant_id;
    allwords       = t_postperf.words;
    all_ori_newold = t_postperf.Corr_resp_Old_New_;
    all_wordacc    = t_postperf.word_selection_Accuracy;
    all_imagacc    = t_postperf.Images_Selection_Accuracy;
    all_blks       = t_postperf.Block;
    ManScore_all   = t_postperf.scaleManual;
    ManScore_all(all_wordacc == 0) = NaN;
    ManScore_all(all_imagacc == 0) = 0;

    idx_thispart = contains(allparts,['sub-',sbj_str]);
    idx_thispart(contains(all_ori_newold,'New')) = 0;  % Muzzle all new trials

    % Taking this participant
    words_thispart      = allwords(idx_thispart);
    ori_newold_thispart = all_ori_newold(idx_thispart);
    wordacc_thispart    = all_wordacc(idx_thispart);
    imagacc_thispart    = all_imagacc(idx_thispart);
    blks_thispart       = all_blks(idx_thispart);
    manscore_thispart   = ManScore_all(idx_thispart);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    Read and extract trial event table
    events_tsv_path = fullfile(BIDS_ROOT, ['sub-',sbj_str], 'eeg', ...
        ['sub-',sbj_str,'_task-nap_events.tsv']);
    t = readtable(events_tsv_path,"FileType","text");

    % The BIDS-enriched events.tsv includes extra "fixation" rows that
    % inherit the Conditions label from the nearest TICUEE/TIB4CUEE
    % marker but have no word/old_new identity (blank). These weren't
    % present in the original standalone file this script was built
    % against, and inflate the per-condition row counts beyond the
    % actual number of stimulus-delivery events in the EEG markers.
    % Remove them here so counts match nume_soundonly/nume_TICue/nume_TIB4Cue.
    t = t(~cellfun(@isempty, t.word), :);

    faceobj       = t.Images_type;
    oldnew        = t.old_new;
    cat_alltrials = t.Conditions;
    allwords      = t.word;

    % Resolves all repetitions
    words_unq = unique(allwords);
    timesrep  = NaN(length(allwords),1);
    for nwod = 1:length(words_unq)
        idxthisword     = find(strcmp(allwords,words_unq{nwod}));
        numrep_thisword = length(idxthisword);
        timesrep(idxthisword) = 1:length(idxthisword);
    end
    t = addvars(t,timesrep,'NewVariableNames','RealWordRepetitions');

    % ---- Sound only (Cue) ----
    idxso_tbl = ismember(cat_alltrials,'Cue+Foil');
    if sum(idxso_tbl) ~= nume_soundonly
        error('Check')
    else
        tinfo_so = t(idxso_tbl,:);
    end
    % Match behaviour performances
    so_words          = tinfo_so.word;
    idxthisword_test  = NaN(length(so_words),1);
    wordacc_post_so   = NaN(length(so_words),1);
    imgacc_post_so    = NaN(length(so_words),1);
    blk_so            = NaN(length(so_words),1);
    manscore_so       = NaN(length(so_words),1);
    oldnew_this       = tinfo_so.old_new;
    for jj = 1:length(so_words)
        thisword = so_words{jj};
        idnow    = find(contains(words_thispart,string(thisword)));
        if ~isempty(idnow)
            idxthisword_test(jj) = find(contains(words_thispart,string(thisword)));
            wordacc_post_so(jj)  = wordacc_thispart(idxthisword_test(jj));
            imgacc_post_so(jj)   = imagacc_thispart(idxthisword_test(jj));
            blk_so(jj)           = blks_thispart(idxthisword_test(jj));
            manscore_so(jj)      = manscore_thispart(idxthisword_test(jj));
        else
            idxthisword_test(jj) = NaN;
            % Check for error
            if strcmp(oldnew_this{jj},{'Old'})
                error('Check code')
            end
        end
    end
    tinfo_so = addvars(tinfo_so,wordacc_post_so,imgacc_post_so,blk_so,manscore_so, ...
        'NewVariableNames',{'WordAcc_Post','ImagAcc_Post','Block_Post','Manual_Score'});

    % ---- TI with Cue ----
    idxticue_tbl = ismember(cat_alltrials,'TI+Cue');
    if sum(idxticue_tbl) ~= nume_TICue
        error('Check')
    else
        tinfo_ticue = t(idxticue_tbl,:);
    end
    % Match behaviour performances
    ticue_words        = tinfo_ticue.word;
    idxthisword_test   = NaN(length(ticue_words),1);
    wordacc_post_ticue = NaN(length(ticue_words),1);
    imgacc_post_ticue  = NaN(length(ticue_words),1);
    blk_ticue          = NaN(length(ticue_words),1);
    manscore_ticue     = NaN(length(ticue_words),1);
    oldnew_this        = tinfo_ticue.old_new;
    for jj = 1:length(ticue_words)
        thisword = ticue_words{jj};
        idnow    = find(contains(words_thispart,string(thisword)));
        if ~isempty(idnow)
            idxthisword_test(jj)  = find(contains(words_thispart,string(thisword)));
            wordacc_post_ticue(jj) = wordacc_thispart(idxthisword_test(jj));
            imgacc_post_ticue(jj)  = imagacc_thispart(idxthisword_test(jj));
            blk_ticue(jj)          = blks_thispart(idxthisword_test(jj));
            manscore_ticue(jj)     = manscore_thispart(idxthisword_test(jj));
        else
            idxthisword_test(jj) = NaN;
            if strcmp(oldnew_this{jj},{'Old'})
                error('Check code')
            end
        end
    end
    tinfo_ticue = addvars(tinfo_ticue,wordacc_post_ticue,imgacc_post_ticue,blk_ticue,manscore_ticue, ...
        'NewVariableNames',{'WordAcc_Post','ImagAcc_Post','Block_Post','Manual_Score'});

    % ---- TI before Cue ----
    idxtib4cue_tbl = ismember(cat_alltrials,'TIB4Cue');
    if sum(idxtib4cue_tbl) ~= nume_TIB4Cue
        error('Check')
    else
        tinfo_tib4cue = t(idxtib4cue_tbl,:);
    end
    % Match behaviour performances
    tib4cue_words        = tinfo_tib4cue.word;
    idxthisword_test     = NaN(length(tib4cue_words),1);
    wordacc_post_tib4cue = NaN(length(tib4cue_words),1);
    imgacc_post_tib4cue  = NaN(length(tib4cue_words),1);
    blk_tib4cue          = NaN(length(tib4cue_words),1);
    oldnew_this          = tinfo_tib4cue.old_new;
    manscore_tib4cue     = NaN(length(tib4cue_words),1);
    for jj = 1:length(tib4cue_words)
        thisword = tib4cue_words{jj};
        idnow    = find(contains(words_thispart,string(thisword)));
        if ~isempty(idnow)
            idxthisword_test(jj)    = find(contains(words_thispart,string(thisword)));
            wordacc_post_tib4cue(jj) = wordacc_thispart(idxthisword_test(jj));
            imgacc_post_tib4cue(jj)  = imagacc_thispart(idxthisword_test(jj));
            blk_tib4cue(jj)          = blks_thispart(idxthisword_test(jj));
            manscore_tib4cue(jj)     = manscore_thispart(idxthisword_test(jj));
        else
            idxthisword_test(jj) = NaN;
            if strcmp(oldnew_this{jj},{'Old'})
                error('Check code')
            end
        end
    end
    tinfo_tib4cue = addvars(tinfo_tib4cue,wordacc_post_tib4cue,imgacc_post_tib4cue,blk_tib4cue,manscore_tib4cue, ...
        'NewVariableNames',{'WordAcc_Post','ImagAcc_Post','Block_Post','Manual_Score'});

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Group table
    tbl_new = table(idxnow,{hyp_this},{lat_d_soundonly},{eeg_soundonly},{slpstage_soundonly},{mkrs_sound_realidx}, ...
        {lat_d_TICue},{eeg_TICue},{slpstage_TICue},{lat_d_TIB4Cue},{eeg_TIB4Cue},{slpstage_TIB4Cue}, ...
        {artflag_soundonly},{artflag_TICue},{artflag_TIB4Cue}, ...
        {idx_soundonly},{idx_TICue},{idx_TIB4Cue},{eegevents},{realslpstage_soundonly},{realslpstage_TICue},{realslpstage_TIB4Cue}, ...
        {eeg_emptysleep},{artflag_emptysleep},{slpstage_empty_real},{reallatd_empty},{eeg_B4TICue},{artflag_B4TICue} ...
        ,{realslpstage_B4TICue},{lat_d_B4TICue},{tinfo_so},{tinfo_ticue},{tinfo_tib4cue}, ...
        {lat_d_SleepOnly},{eeg_SleepOnly},{artflag_SleepOnly},{realslpstage_SleepOnly},'VariableNames',{'Sbjidx','hypnogram' ...
        'Latency_Soundonly','EEG_soundonly','Slpstage_soundonly','Sound_Oldnew', ...
        'Latency_TICue','EEG_TICue','Slpstage_TICue','Latency_TIB4Cue','EEG_TIB4Cue','Slpstage_TIB4Cue', ...
        'Artflag_soundonly','Artflag_TICue','Artflag_TIB4Cue','Idx_soundonly','Idx_TICue','Idx_TIB4Cue','OriginalEvents','RealSlpStage_Soundonly','RealSlpStage_TICue','RealSlpStage_TIB4Cue' ...
        ,'EEG_EmptySleep','Artflag_EmptySleep','SleepStageReal_EmptySleep','EventLats_EmptySleep','EEG_B4TICue','Artflag_B4TICue','SleepStageReal_B4TICue','EventLats_B4TICue', ...
        'Tinfo_SO','Tinfo_TICue','Tinfo_TIB4Cue','Latency_SleepOnly','EEG_SleepOnly','Artflag_SleepOnly','Slpstage_SleepOnly'});

    if isempty(tbl_storage)
        tbl_storage = tbl_new;
    else
        tbl_storage = [tbl_storage; tbl_new];
    end

    clear eeg EEG_nap eegevents

end

%% Save

save('NapData_AllSubjects.mat','tbl_storage','Fs','-v7.3');
fprintf('Saved NapData_AllSubjects.mat (%d subjects).\n', height(tbl_storage));