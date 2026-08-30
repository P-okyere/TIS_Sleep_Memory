clear; clc; close all;

%% ============================
%  EDIT FOR YOUR ENVIRONMENT
%% ============================
bids_root      = 'E:\BIDS';
bids_deriv     = fullfile(bids_root, 'derivatives');
base_path_eeg  = bids_root;
out_root       = fullfile(bids_deriv, 'spindle_amplitude_stats');
batlow_path    = 'C:\Users\po00240\OneDrive - University of Surrey\Desktop\Projects\MYSTI\_matlab_scripts (1)_first\Stats\ScientificColourMaps8\batlow';
% ============================

if ~exist(out_root,'dir'), mkdir(out_root); end

%% ============================
%  SETUP
%% ============================

marker_set = {'MarkerOldSound','TICUEE','TIB4CUEE'};
marker_labels_map = struct( ...
    'MarkerOldSound','Cue', ...
    'TICUEE','TIS+Cue', ...
    'TIB4CUEE','TIS+Cue delayed');

sleep_stages_use = {'N2','N3'};
spindle_types    = {'fast'};
freq_labels      = struct('fast','12–16 Hz');

% Cluster permutation settings
n_perm           = 1000;
cluster_alpha    = 0.05;
cluster_p_thresh = 0.05;
min_cluster_size = 2;

%% ============================
%  LOAD DATA — fast spindle file
%% ============================
load(fullfile(bids_deriv, 'spindle_so_detection', 'desc-fastSpindleMetrics_trialdata.mat'), ...
     'all_trial_data');

%% ============================
%  EEG LAYOUT
%% ============================
example_EEG = pop_loadset('filename', 'sub-41_task-nap_eeg.set', ...
                          'filepath', fullfile(base_path_eeg, 'sub-41', 'eeg'));
ft_data     = eeglab2fieldtrip(example_EEG, 'preprocessing', 'none');

cfg_layout        = [];
cfg_layout.elec   = ft_data.elec;
cfg_layout.layout = 'acticap-64ch-standard2.mat';
cfg_layout.rotate = 0;
layout            = ft_prepare_layout(cfg_layout);

%% === Add missing channels to layout ===
AFz_pos = mean([
    layout.pos(strcmp(layout.label, 'Fz'),  :);
    layout.pos(strcmp(layout.label, 'AF3'), :);
    layout.pos(strcmp(layout.label, 'AF4'), :)
], 1);
FCz_pos = mean([
    layout.pos(strcmp(layout.label, 'Fz'), :);
    layout.pos(strcmp(layout.label, 'Cz'), :)
], 1);
I1_pos = layout.pos(strcmp(layout.label, 'O1'), :);
I2_pos = layout.pos(strcmp(layout.label, 'O2'), :);

new_labels = {'AFz'; 'FCz'; 'I1'; 'I2'};
new_pos    = [AFz_pos; FCz_pos; I1_pos; I2_pos];

layout.label  = [layout.label;  new_labels];
layout.pos    = [layout.pos;    new_pos];
layout.width  = [layout.width;  repmat(mean(layout.width),  4, 1)];
layout.height = [layout.height; repmat(mean(layout.height), 4, 1)];

%% === Rescale electrode positions inward ===
scale_factor = 0.85;
center       = mean(layout.pos, 1);
layout.pos   = center + scale_factor * (layout.pos - center);

fprintf('Layout updated — %d channels\n', numel(layout.label));

channels   = ft_data.label;
n_channels = numel(channels);

%% ============================
%  NEIGHBOUR STRUCTURE
%% ============================
cfg_neigh        = [];
cfg_neigh.method = 'triangulation';
cfg_neigh.layout = 'acticap-64ch-standard2.mat';
neighbours       = ft_prepare_neighbours(cfg_neigh);

% Add missing channels neighbours
idx = numel(neighbours) + 1;
neighbours(idx).label       = 'AFz';
neighbours(idx).neighblabel = {'Fz','AF3','AF4','F1','F2'};
idx = idx + 1;
neighbours(idx).label       = 'FCz';
neighbours(idx).neighblabel = {'Fz','Cz','FC1','FC2','FC3','FC4'};
idx = idx + 1;
neighbours(idx).label       = 'I1';
neighbours(idx).neighblabel = {'O1','PO7','PO3','Oz'};
idx = idx + 1;
neighbours(idx).label       = 'I2';
neighbours(idx).neighblabel = {'O2','PO8','PO4','Oz'};

% Bidirectional links
new_chans = {'AFz','FCz','I1','I2'};
new_neighs = {
    {'Fz','AF3','AF4','F1','F2'}, ...
    {'Fz','Cz','FC1','FC2','FC3','FC4'}, ...
    {'O1','PO7','PO3','Oz'}, ...
    {'O2','PO8','PO4','Oz'} ...
};

for nc = 1:numel(new_chans)
    new_chan = new_chans{nc};
    for nn = 1:numel(new_neighs{nc})
        neigh_name       = new_neighs{nc}{nn};
        neigh_struct_idx = find(strcmp({neighbours.label}, neigh_name));
        if ~isempty(neigh_struct_idx)
            if ~ismember(new_chan, neighbours(neigh_struct_idx).neighblabel)
                neighbours(neigh_struct_idx).neighblabel{end+1} = new_chan;
            end
        end
    end
end

fprintf('Average neighbours per channel: %.1f\n', ...
    mean(arrayfun(@(x) numel(x.neighblabel), neighbours)));
overlap = intersect(channels, {neighbours.label});
fprintf('Channels with neighbours: %d / %d\n', numel(overlap), numel(channels));
fz_idx = find(strcmp({neighbours.label}, 'Fz'));
fprintf('Fz neighbours include AFz: %d\n\n', ...
    ismember('AFz', neighbours(fz_idx).neighblabel));

%% ============================
%  COLORMAP
%% ============================
addpath(batlow_path);
load('batlow.mat');

%% ============================
%  DOT COLOURS
%% ============================
col_sig    = [1.0 1.0 1.0];   % white = significant cluster
col_nonsig = [0.6 0.6 0.6];   % grey  = not significant

%% ============================
%  FIGURE SIZE - small, fixed (matches other reduced topoplot scripts:
%  1.6 x 1.6 inches, 7pt Helvetica, dot sizes 20/4)
%% ============================
font_size  = 4;
font_name  = 'Arial';
fig_width  = 1.2;
fig_height = 1.2;
title_pos  = [0.5, 1.08, 0];   % consistent with all other topoplot scripts

%% ============================
%  LOOP OVER SPINDLE TYPES
%% ============================
for s = 1:numel(spindle_types)

    spindle_type = spindle_types{s};
    freq_range   = freq_labels.(spindle_type);

    fprintf('\n=== Spindle type: %s (%s) ===\n', spindle_type, freq_range);

    plot_dir_corrected = fullfile(out_root, spindle_type, 'topoplot_significant');
    plot_dir_null      = fullfile(out_root, spindle_type, 'topoplot_not_significant');
    if ~exist(plot_dir_corrected,'dir'), mkdir(plot_dir_corrected); end
    if ~exist(plot_dir_null,'dir'),      mkdir(plot_dir_null);      end

    %% ── Filter table ─────────────────────────────────────
    T = all_trial_data( ...
        strcmp(all_trial_data.SpindleType, spindle_type) & ...
        ismember(all_trial_data.SleepStage, sleep_stages_use) & ...
        ismember(all_trial_data.MarkerType, marker_set) & ...
        ~isnan(all_trial_data.MeanAmplitude), :);

    T.MarkerType  = categorical(T.MarkerType, marker_set);
    T.Participant = categorical(T.Participant);

    participants = categories(T.Participant);
    nP           = numel(participants);

    fprintf('Participants: %d | Conditions: %d\n', nP, numel(marker_set));

    %% ====================================================
    %  STEP 1 — Observed F values (channel-wise LME)
    %% ====================================================
    fvals_obs = nan(n_channels, 1);
    pvals_obs = nan(n_channels, 1);

    for c = 1:n_channels
        chan = channels{c};
        tbl  = T(strcmp(T.Channel, chan), :);
        if height(tbl) >= 3 && numel(unique(tbl.Participant)) > 1
            try
                lme           = fitlme(tbl, ...
                    'MeanAmplitude ~ MarkerType + Trial + (1|Participant)');
                a             = anova(lme);
                fvals_obs(c)  = a.FStat(3);
                pvals_obs(c)  = a.pValue(3);
            catch
                warning('LME failed: %s', chan);
            end
        end
    end

    fprintf('Uncorrected significant channels: %d\n', ...
        sum(pvals_obs < cluster_alpha, 'omitnan'));

    %% ====================================================
    %  STEP 2 — Participant-level condition means
    %% ====================================================
    n_conds   = numel(marker_set);
    sub_means = nan(nP, n_channels, n_conds);

    for p = 1:nP
        for c = 1:n_channels
            for k = 1:n_conds
                vals = T.MeanAmplitude( ...
                    T.Participant == participants{p} & ...
                    strcmp(T.Channel, channels{c})  & ...
                    T.MarkerType == marker_set{k});
                if ~isempty(vals)
                    sub_means(p,c,k) = mean(vals, 'omitnan');
                end
            end
        end
    end

    %% ====================================================
    %  STEP 3 — Observed clusters
    %% ====================================================
    sig_mask_obs = pvals_obs < cluster_alpha & ~isnan(fvals_obs);

    [clusters_obs, masses_obs] = find_all_clusters(fvals_obs, sig_mask_obs, ...
        channels, neighbours, min_cluster_size);

    fprintf('Valid clusters (>=%d channels): %d\n', ...
        min_cluster_size, numel(clusters_obs));

    %% ====================================================
    %  STEP 4 — Fast permutation (sign-flip on participant means)
    %% ====================================================
    fprintf('Running %d permutations...\n', n_perm);
    perm_max_mass = nan(n_perm, 1);

    rng(42);

    for perm = 1:n_perm

        if mod(perm, 100) == 0
            fprintf('  Permutation %d / %d\n', perm, n_perm);
        end

        fvals_perm = nan(n_channels, 1);
        pvals_perm = nan(n_channels, 1);

        for c = 1:n_channels
            Y     = squeeze(sub_means(:,c,:));  % nP x n_conds
            valid = all(~isnan(Y), 2);
            Y     = Y(valid,:);
            n     = size(Y,1);
            if n < 5, continue; end

            flip   = (rand(n,1) > 0.5) * 2 - 1;
            Y_flip = Y .* flip;

            grand_mean  = mean(Y_flip(:));
            cond_means  = mean(Y_flip, 1);
            sub_means_c = mean(Y_flip, 2);

            SS_between = n * sum((cond_means - grand_mean).^2);
            SS_subject = n_conds * sum((sub_means_c - grand_mean).^2);
            SS_total   = sum((Y_flip(:) - grand_mean).^2);
            SS_error   = SS_total - SS_between - SS_subject;

            df_between = n_conds - 1;
            df_error   = (n-1) * (n_conds-1);

            if df_error > 0 && SS_error > 0
                F = (SS_between/df_between) / (SS_error/df_error);
                fvals_perm(c) = F;
                pvals_perm(c) = 1 - fcdf(F, df_between, df_error);
            end
        end

        sig_mask_perm       = pvals_perm < cluster_alpha & ~isnan(fvals_perm);
        [~, masses_perm]    = find_all_clusters(fvals_perm, sig_mask_perm, ...
                                                channels, neighbours, min_cluster_size);
        perm_max_mass(perm) = max([0; masses_perm(:)]);

    end

    %% ====================================================
    %  STEP 5 — Test each cluster against null distribution
    %% ====================================================
    all_sig_channels = {};
    n_sig            = 0;

    for cl = 1:numel(clusters_obs)
        cl_pval = mean(perm_max_mass >= masses_obs(cl));
        if cl_pval < cluster_p_thresh
            n_sig = n_sig + 1;
            all_sig_channels = [all_sig_channels, channels(clusters_obs{cl})]; %#ok
            fprintf('  Cluster %d: %d ch, mass=%.3f, p=%.4f ✅\n', ...
                cl, numel(clusters_obs{cl}), masses_obs(cl), cl_pval);
        else
            fprintf('  Cluster %d: %d ch, mass=%.3f, p=%.4f ❌\n', ...
                cl, numel(clusters_obs{cl}), masses_obs(cl), cl_pval);
        end
    end

    all_sig_channels = unique(all_sig_channels);
    any_sig          = ~isempty(all_sig_channels);
    fprintf('Significant cluster channels: %d\n', numel(all_sig_channels));

    %% ====================================================
    %  STEP 6 — FieldTrip data struct
    %% ====================================================
    data        = [];
    data.label  = layout.label;
    data.avg    = nan(numel(layout.label), 1);
    data.time   = 0;
    data.dimord = 'chan_time';

    for j = 1:numel(layout.label)
        idx = strcmp(layout.label{j}, channels);
        if any(idx)
            data.avg(j) = fvals_obs(idx);
        end
    end

    %% ====================================================
    %  STEP 7 — Plot
    %% ====================================================
    fig = figure('Color','w','Units','inches', ...
                 'Position',[1 1 fig_width fig_height]);

    cfg           = [];
    cfg.layout    = layout;
    cfg.parameter = 'avg';
    cfg.colormap  = batlow;
    cfg.comment   = 'no';
    cfg.marker    = 'off';
    cfg.style     = 'both';
    cfg.interplim = 'head';
    cfg.colorbar  = 'yes';

    cfg.zlim = [0 10];
    ft_topoplotER(cfg, data);

    % FIX: thin the head/ears/nose outline (FieldTrip default LineWidth=2),
    % matching the other reduced topoplot scripts
    head_lines = findobj(gca, 'Type', 'Line', 'LineWidth', 2);
    set(head_lines, 'LineWidth', 0.8);

    h = colorbar;
    ylabel(h, 'F value', 'FontSize', font_size);
    h.FontSize = font_size;
    h.FontName = font_name;
    h.Ticks    = 0:2:10;   % explicit intermediate ticks (0,2,4,6,8,10)

    hold on;
    for j = 1:numel(layout.label)
        idx = find(strcmp(layout.label{j}, channels));
        if ~isempty(idx)
            pos      = layout.pos(j,:);
            chan_lbl = layout.label{j};
            if ismember(chan_lbl, all_sig_channels)
                scatter(pos(1), pos(2), 5, col_sig,    'filled');
            else
                scatter(pos(1), pos(2), 3, col_nonsig, 'filled');
            end
        end
    end
    hold off;

    % Title — pushed up, consistent with other topoplot scripts, same
    % font size/name as the other reduced topoplots
    title({'Fast Spindle amplitude', ...
           }, ...
          'FontSize', font_size, 'FontWeight', 'normal', ...
          'FontName', font_name, ...
          'Units', 'normalized', 'Position', title_pos);

    ax = gca; ax.FontSize = font_size; ax.FontName = font_name;

    %% ====================================================
    %  STEP 8 — Save figure
    %% ====================================================
    set(fig,'PaperPositionMode','auto','Renderer','painters');

    fname    = sprintf('%s_absolute_ANOVA_clusterperm', spindle_type);
    save_dir = plot_dir_corrected;
    if ~any_sig
        save_dir = plot_dir_null;
    end

    print(fig, fullfile(save_dir,[fname '.pdf']),  '-dpdf',  '-painters');
    print(fig, fullfile(save_dir,[fname '.eps']),  '-depsc', '-painters');
    print(fig, fullfile(save_dir,[fname '.tiff']), '-dtiff', '-r600');
    print(fig, fullfile(save_dir,[fname '.png']),  '-dpng',  '-r600');
    print(fig, fullfile(save_dir,[fname '.svg']),  '-dsvg',  '-painters');  % ADDED

    %% ====================================================
    %  STEP 9 — Save stats
    %% ====================================================
    if ~isempty(clusters_obs)
        cl_pvals = arrayfun(@(cl) mean(perm_max_mass >= masses_obs(cl)), ...
                            1:numel(clusters_obs));
        cluster_summary = table(...
            (1:numel(clusters_obs))', ...
            masses_obs(:), ...
            cl_pvals(:), ...
            cellfun(@numel, clusters_obs)', ...
            'VariableNames',{'Cluster','Mass','p_value','N_channels'});
        writetable(cluster_summary, ...
                   fullfile(save_dir,[fname '_cluster_summary.csv']));
    end

    in_sig    = ismember(channels, all_sig_channels);
    chan_stats = table(channels(:), fvals_obs(:), pvals_obs(:), in_sig(:), ...
        'VariableNames',{'Channel','F_value','p_uncorrected','In_sig_cluster'});
    writetable(chan_stats, fullfile(save_dir,[fname '_channel_stats.csv']));

    if ~isempty(all_sig_channels)
        writetable(table(all_sig_channels(:), ...
                         'VariableNames',{'SignificantChannel'}), ...
                   fullfile(save_dir,[fname '_sig_channels.csv']));
    end

    % Per-cluster channel export (one CSV per significant cluster, numbered
    % 1, 2, ... in the order found). Downstream plotting scripts should
    % load these files rather than hardcoding channel lists, so a rerun
    % with different permutation results, thresholds, or channel counts
    % can't silently go out of sync with the plots.
    sig_rank = 0;
    for cl = 1:numel(clusters_obs)
        cl_pval = mean(perm_max_mass >= masses_obs(cl));
        if cl_pval >= cluster_p_thresh, continue; end
        sig_rank = sig_rank + 1;

        cl_chans = channels(clusters_obs{cl})';
        cl_table = table(cl_chans(:), ...
            repmat(masses_obs(cl), numel(cl_chans), 1), ...
            repmat(cl_pval,        numel(cl_chans), 1), ...
            'VariableNames', {'Channel','ClusterMass','ClusterPValue'});
        writetable(cl_table, fullfile(save_dir, ...
            sprintf('%s_sig_channels_cluster%d.csv', fname, sig_rank)));
    end

    save(fullfile(save_dir,[fname '_stats.mat']), ...
         'channels','fvals_obs','pvals_obs', ...
         'clusters_obs','masses_obs', ...
         'all_sig_channels','perm_max_mass', ...
         'sub_means','participants','marker_set');

    fprintf('✅ Saved to: %s\n', save_dir);

end

fprintf('\n✅ All done.\n');


%% ============================
%  HELPER FUNCTION
%% ============================
function [clusters, masses] = find_all_clusters(tvals, sig_mask, channels, neighbours, min_size)
    clusters = {};
    masses   = [];
    visited  = false(numel(channels), 1);
    sig_idx  = find(sig_mask);

    for s = 1:numel(sig_idx)
        seed = sig_idx(s);
        if visited(seed), continue; end

        cluster = [];
        queue   = seed;

        while ~isempty(queue)
            node     = queue(1);
            queue(1) = [];
            if visited(node), continue; end
            visited(node)  = true;
            cluster(end+1) = node; %#ok

            chan_name        = channels{node};
            neigh_idx_struct = find(strcmp({neighbours.label}, chan_name));

            if ~isempty(neigh_idx_struct)
                neigh_names = neighbours(neigh_idx_struct).neighblabel;
                for nb = 1:numel(neigh_names)
                    neigh_chan_idx = find(strcmp(channels, neigh_names{nb}));
                    if ~isempty(neigh_chan_idx) && ...
                       sig_mask(neigh_chan_idx)  && ...
                       ~visited(neigh_chan_idx)
                        queue(end+1) = neigh_chan_idx; %#ok
                    end
                end
            end
        end

        if numel(cluster) >= min_size
            clusters{end+1} = cluster;        %#ok
            masses(end+1)   = sum(abs(tvals(cluster))); %#ok
        end
    end
end
