clc; clear; close all;
rng(42);

%% ========================================================================
%% BRAIN-BEHAVIOUR TOPOPLOT: SO-SPINDLE COUPLING RATE -> FORGETTING
%%
%% For each channel: Pearson correlation between coupling rate and forgetting
%% Cluster-based permutation correction across channels
%%
%% Figure size: 1.2 x 1.2 inches (small, fixed - no resizing needed)
%% Font: Arial 4pt
%% White dots = cluster-corrected significant (p<=.05)
%% Grey dots  = not significant
%%
%% LATEST FIXES:
%%   - Head/ears/nose outline thinned (was default LineWidth=2, now 0.8)
%%   - Colorbar reduced to 2 ticks (min/max only) instead of 5
%% ========================================================================

%% -------- EDIT FOR YOUR ENVIRONMENT --------
eeglab_path   = 'C:\Users\po00240\OneDrive - University of Surrey\Desktop\eeglab';
colourmap_path = 'C:\Users\po00240\OneDrive - University of Surrey\Desktop\Projects\MYSTI\ScientificColourMaps8';
bids_root     = 'E:\BIDS';
bids_deriv    = fullfile(bids_root, 'derivatives');

behavior_file = fullfile(bids_deriv, 'behaviour', 'desc-FamRecFidSummary_results.xlsx');
results_file  = fullfile(bids_deriv, 'so_spindle_coupling', 'desc-SOSpindleCouplingRate_metrics.mat');

out_dir = fullfile(bids_deriv, 'coupling_rate_forgetting_stats');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
fprintf('Output: %s\n\n', out_dir);
% --------------------------------------------

%% -------- PARTICIPANTS --------
participants = { ...
    'sub-02','sub-03','sub-04','sub-05','sub-06','sub-08','sub-12', ...
    'sub-13','sub-14','sub-18','sub-19','sub-20','sub-21','sub-22', ...
    'sub-23','sub-24','sub-25','sub-26','sub-27','sub-29','sub-30', ...
    'sub-33','sub-36','sub-37','sub-38','sub-39','sub-40','sub-41'};
n_participants = numel(participants);

%% -------- CONDITIONS --------
cond_labels = {'TIS+Cue','TIS+Cue delayed','Cue'};
n_conds     = numel(cond_labels);

rec_cols_v2 = {'V2_TISCue_recollection_pct', ...
               'V2_TISCue_delayed_recollection_pct', ...
               'V2_Cue_recollection_pct'};

%% -------- SETTINGS --------
n_perms          = 1000;
min_cluster_size = 2;

%% -------- EEG LAYOUT --------
fprintf('Loading EEG layout...\n');
addpath(eeglab_path);
eeglab nogui;

eeg_path = fullfile(bids_root, 'sub-41', 'eeg');
EEG      = pop_loadset('filename','sub-41_task-nap_eeg.set','filepath',eeg_path);
ft_data  = eeglab2fieldtrip(EEG,'preprocessing','none');

cfg_layout        = [];
cfg_layout.elec   = ft_data.elec;
cfg_layout.layout = 'acticap-64ch-standard2.mat';
cfg_layout.rotate = 0;
layout     = ft_prepare_layout(cfg_layout);
layout     = add_missing_channels(layout);
neighbours = build_neighbours();

%% -------- COLORMAP --------
addpath(genpath(colourmap_path));
if exist('batlow.mat','file'), load('batlow.mat'); cmap = batlow;
else, cmap = parula(256); end

%% ========================================================================
%% LOAD COUPLING RATE
%% ========================================================================
fprintf('Loading coupling rate data...\n');
load(results_file, 'p_rate_all', 'channels');
n_channels = size(p_rate_all, 2);
fprintf('Channels: %d | Participants: %d\n\n', n_channels, n_participants);

%% ========================================================================
%% LOAD FORGETTING
%% ========================================================================
fprintf('Loading forgetting data...\n');
behavior_v2 = readtable(behavior_file, 'Sheet', 'EEG_Correlation_Master');

forgetting = nan(n_participants, n_conds);
for p = 1:n_participants
    sub_id = participants{p};
    idx_b  = find(strcmp(behavior_v2.participant_id, sub_id));
    if isempty(idx_b), continue; end
    for cc = 1:n_conds
        rec = behavior_v2.(rec_cols_v2{cc})(idx_b);
        forgetting(p, cc) = 100 - rec;
    end
end

fprintf('Forgetting per condition:\n');
for cc = 1:n_conds
    fprintf('  %-20s: mean=%.1f%%, SD=%.1f%%\n', cond_labels{cc}, ...
        nanmean(forgetting(:,cc)), nanstd(forgetting(:,cc)));
end
fprintf('\n');

%% ========================================================================
%% MAIN ANALYSIS: CHANNEL-WISE PEARSON CORRELATION
%% ========================================================================
fprintf('Running channel-wise correlations...\n\n');

tstat_all       = nan(n_channels, n_conds);
rho_all         = nan(n_channels, n_conds);
p_raw_all       = nan(n_channels, n_conds);
beta_all        = nan(n_channels, n_conds);
h_clust_all     = false(n_channels, n_conds);
h_clust_all_p05 = false(n_channels, n_conds);

df = n_participants - 2;

for cc = 1:n_conds
    forg_vec   = forgetting(:,cc);
    forg_z     = (forg_vec - mean(forg_vec,'omitnan')) ./ std(forg_vec,'omitnan');
    forg_valid = ~isnan(forg_z);
    rate_mat   = squeeze(p_rate_all(:,:,cc));

    [tmap, bmap, rmap] = compute_tmap(rate_mat, forg_z, forg_valid, df);
    pmap = 2*(1 - tcdf(abs(tmap), df));

    fprintf('  [%s] raw p<.05: %d | Permuting...', ...
        cond_labels{cc}, sum(pmap<0.05,'omitnan'));
    t0 = tic;

    [h_clust, h_clust_p05, sig_clusters, cluster_stats] = run_cluster_permutation(...
        tmap, pmap, rate_mat, forg_z, forg_valid, df, ...
        channels, neighbours, n_perms, min_cluster_size);

    fprintf(' %.1fs | Cluster sig: %d\n', toc(t0), sum(h_clust_p05));

    for cl = 1:numel(sig_clusters)
        ch_names = strjoin(channels(sig_clusters{cl}), ', ');
        fprintf('    Cluster %d (%d ch, mass=%.2f, p=%.4f): %s\n', ...
            cl, numel(sig_clusters{cl}), ...
            cluster_stats{cl}(1), cluster_stats{cl}(2), ch_names);
    end

    tstat_all(:,cc)       = tmap;
    rho_all(:,cc)         = rmap;
    p_raw_all(:,cc)       = pmap;
    beta_all(:,cc)        = bmap;
    h_clust_all(:,cc)     = h_clust;
    h_clust_all_p05(:,cc) = h_clust_p05;
end

%% ========================================================================
%% FIGURES (rho, not t-statistic) - single topoplots only
%% ========================================================================

all_r = rho_all(~isnan(rho_all));
clim  = min(1, prctile(abs(all_r),95));
if clim == 0, clim = 0.5; end

for cc = 1:n_conds
    cond_fname = strrep(strrep(cond_labels{cc},'+',' '),' ','_');

    % Uncorrected
    plot_topo(rho_all(:,cc), p_raw_all(:,cc), [], ...
        layout, channels, cmap, cond_labels{cc}, clim, 'uncorrected', out_dir, ...
        sprintf('Topo_CouplingRate_Forgetting_%s_uncorrected', cond_fname));

    % Corrected
    plot_topo(rho_all(:,cc), p_raw_all(:,cc), h_clust_all_p05(:,cc), ...
        layout, channels, cmap, cond_labels{cc}, clim, 'corrected', out_dir, ...
        sprintf('Topo_CouplingRate_Forgetting_%s_corrected', cond_fname));
end

%% -------- SAVE CSV --------
rows = {};
for cc = 1:n_conds
    for ch = 1:n_channels
        rows(end+1,1:8) = {cond_labels{cc}, channels{ch}, ...
            rho_all(ch,cc), beta_all(ch,cc), tstat_all(ch,cc), p_raw_all(ch,cc), ...
            h_clust_all(ch,cc), h_clust_all_p05(ch,cc)}; %#ok
    end
end
writetable(cell2table(rows,'VariableNames', ...
    {'Condition','Channel','Rho','Beta','t_stat','p_raw','sig_Cluster_strict','sig_Cluster_p05'}), ...
    fullfile(out_dir,'PerChannel_CouplingRate_Forgetting.csv'));

save(fullfile(out_dir,'Results_CouplingRate_Forgetting.mat'), ...
    'tstat_all','rho_all','p_raw_all','beta_all','h_clust_all','h_clust_all_p05', ...
    'channels','cond_labels','forgetting');

fprintf('\nDone. Output: %s\n', out_dir);

%% ========================================================================
%% HELPER: COMPUTE T-MAP (also returns rho directly)
%% ========================================================================

function [tmap, bmap, rmap] = compute_tmap(data_mat, forg_z, forg_valid, df)
    n_ch = size(data_mat,2);
    tmap = nan(n_ch,1); bmap = nan(n_ch,1); rmap = nan(n_ch,1);
    forg_v = forg_z(forg_valid); dat_v = data_mat(forg_valid,:);
    forg_c = forg_v - mean(forg_v,'omitnan');
    for ch = 1:n_ch
        y = dat_v(:,ch); v = ~isnan(y);
        if sum(v) < 4; continue; end
        x = forg_c(v); yv = y(v); n = sum(v);
        xc = x - mean(x); yc = yv - mean(yv);
        dx = norm(xc); dy = norm(yc);
        if dx == 0 || dy == 0; continue; end
        r = (xc'*yc)/(dx*dy);
        if isnan(r) || abs(r) >= 1; continue; end
        rmap(ch) = r;
        tmap(ch) = r*sqrt(n-2)/sqrt(1-r^2);
        bmap(ch) = r*(dy/dx);
    end
end

%% ========================================================================
%% HELPER: CLUSTER PERMUTATION (unchanged - still based on t-statistic)
%% ========================================================================

function [h_clust, h_clust_p05, sig_clusters, cluster_stats] = run_cluster_permutation(...
        tmap, pmap, data_mat, forg_z, forg_valid, df, ...
        channels, neighbours, n_perms, min_cluster_size)

    sig_clusters  = {};
    cluster_stats = {};

    sig_pos_obs = pmap < 0.05 & tmap > 0 & ~isnan(tmap);
    sig_neg_obs = pmap < 0.05 & tmap < 0 & ~isnan(tmap);

    [obs_cp, obs_mp] = find_all_clusters(tmap,      sig_pos_obs, channels, neighbours, min_cluster_size);
    [obs_cn, obs_mn] = find_all_clusters(abs(tmap), sig_neg_obs, channels, neighbours, min_cluster_size);

    null_pos = zeros(n_perms,1);
    null_neg = zeros(n_perms,1);

    for perm = 1:n_perms
        forg_perm = forg_z(randperm(numel(forg_z)));
        [tp, ~, ~] = compute_tmap(data_mat, forg_perm, forg_valid, df);
        pp = 2*(1 - tcdf(abs(tp), df));
        sig_pos_p = pp < 0.05 & tp > 0 & ~isnan(tp);
        sig_neg_p = pp < 0.05 & tp < 0 & ~isnan(tp);
        [~,mp] = find_all_clusters(tp,      sig_pos_p, channels, neighbours, min_cluster_size);
        [~,mn] = find_all_clusters(abs(tp), sig_neg_p, channels, neighbours, min_cluster_size);
        null_pos(perm) = max([0, mp]);
        null_neg(perm) = max([0, mn]);
    end

    h_clust     = false(numel(channels),1);
    h_clust_p05 = false(numel(channels),1);

    for cl = 1:numel(obs_cp)
        p_clust  = mean(null_pos >= obs_mp(cl));
        ch_names = strjoin(channels(obs_cp{cl}), ', ');
        if p_clust <= 0.05
            h_clust(obs_cp{cl})     = true;
            h_clust_p05(obs_cp{cl}) = true;
            sig_clusters{end+1}     = obs_cp{cl}; %#ok
            cluster_stats{end+1}    = [obs_mp(cl), p_clust]; %#ok
            fprintf('\n      Pos cluster %d: %d ch, mass=%.2f, p=%.4f SIGNIFICANT — %s', ...
                cl, numel(obs_cp{cl}), obs_mp(cl), p_clust, ch_names);
        elseif p_clust < 0.10
            fprintf('\n      Pos cluster %d: %d ch, mass=%.2f, p=%.4f (trend) — %s', ...
                cl, numel(obs_cp{cl}), obs_mp(cl), p_clust, ch_names);
        else
            fprintf('\n      Pos cluster %d: %d ch, mass=%.2f, p=%.4f — %s', ...
                cl, numel(obs_cp{cl}), obs_mp(cl), p_clust, ch_names);
        end
    end
    for cl = 1:numel(obs_cn)
        p_clust  = mean(null_neg >= obs_mn(cl));
        ch_names = strjoin(channels(obs_cn{cl}), ', ');
        if p_clust <= 0.05
            h_clust(obs_cn{cl})     = true;
            h_clust_p05(obs_cn{cl}) = true;
            sig_clusters{end+1}     = obs_cn{cl}; %#ok
            cluster_stats{end+1}    = [obs_mn(cl), p_clust]; %#ok
            fprintf('\n      Neg cluster %d: %d ch, mass=%.2f, p=%.4f SIGNIFICANT — %s', ...
                cl, numel(obs_cn{cl}), obs_mn(cl), p_clust, ch_names);
        elseif p_clust < 0.10
            fprintf('\n      Neg cluster %d: %d ch, mass=%.2f, p=%.4f (trend) — %s', ...
                cl, numel(obs_cn{cl}), obs_mn(cl), p_clust, ch_names);
        else
            fprintf('\n      Neg cluster %d: %d ch, mass=%.2f, p=%.4f — %s', ...
                cl, numel(obs_cn{cl}), obs_mn(cl), p_clust, ch_names);
        end
    end
    fprintf('\n');
end

%% ========================================================================
%% HELPER: FIND CLUSTERS
%% ========================================================================

function [clusters, masses] = find_all_clusters(tvals, sig_mask, channels, neighbours, min_size)
    clusters = {}; masses = [];
    visited  = false(numel(channels),1);
    sig_idx  = find(sig_mask);
    for s = 1:numel(sig_idx)
        seed = sig_idx(s);
        if visited(seed), continue; end
        cluster = []; queue = seed;
        while ~isempty(queue)
            node = queue(1); queue(1) = [];
            if visited(node), continue; end
            visited(node) = true; cluster(end+1) = node; %#ok
            si = find(strcmp({neighbours.label}, channels{node}));
            if ~isempty(si)
                for n = 1:numel(neighbours(si).neighblabel)
                    ni = find(strcmp(channels, neighbours(si).neighblabel{n}));
                    if ~isempty(ni) && sig_mask(ni) && ~visited(ni)
                        queue(end+1) = ni; %#ok
                    end
                end
            end
        end
        if numel(cluster) >= min_size
            clusters{end+1} = cluster; %#ok
            masses(end+1)   = sum(abs(tvals(cluster))); %#ok
        end
    end
end

%% ========================================================================
%% HELPER: PLOT SINGLE TOPOPLOT (small fixed size, rho, thin head outline,
%% 2-tick colorbar)
%% ========================================================================

function plot_topo(rmap, pmap, h_clust, layout, channels, cmap, ...
        con_name, clim, mode, out_dir, fname)

    font_size = 4;
    font_name = 'Arial';
    fig_w     = 1.2;
    fig_h     = 1.2;

    fig = figure('Color','w','Units','inches','Position',[1 1 fig_w fig_h]);

    ft_plot        = [];
    ft_plot.label  = layout.label;
    ft_plot.avg    = nan(numel(layout.label),1);
    ft_plot.time   = 0;
    ft_plot.dimord = 'chan_time';
    for j = 1:numel(layout.label)
        idx = strcmp(layout.label{j}, channels);
        if any(idx), ft_plot.avg(j) = rmap(idx); end
    end

    cfg           = [];
    cfg.layout    = layout;
    cfg.parameter = 'avg';
    cfg.colormap  = cmap;
    cfg.comment   = 'no';
    cfg.marker    = 'off';
    cfg.style     = 'both';
    cfg.interplim = 'head';
    cfg.zlim      = [-clim clim];
    cfg.figure    = 'gca';
    cfg.colorbar  = 'yes';
    ft_topoplotER(cfg, ft_plot);

    % FIX: thin the head/ears/nose outline. FieldTrip draws these as Line
    % objects with LineWidth=2 by default (thicker than the iso-contour
    % lines, which are ~0.5) - find and thin just those.
    head_lines = findobj(gca, 'Type', 'Line', 'LineWidth', 2);
    set(head_lines, 'LineWidth', 0.8);

    cb = findobj(fig,'Type','ColorBar');
    if ~isempty(cb)
        cb.Label.String            = 'Pearson r';
        cb.Label.FontSize          = font_size;
        cb.Label.FontName          = font_name;
        cb.Label.Rotation          = 270;
        cb.Label.VerticalAlignment = 'bottom';
        cb.FontSize                = font_size;
        cb.FontName                = font_name;
        % Colorbar ticks left at MATLAB default (auto, shows intermediate
        % values like -0.4/-0.2/0/0.2/0.4) - per your request to keep the
        % original scale style
    end

    title({'Coupling Rate ~ Forgetting', con_name}, ...
        'FontSize', font_size+1, 'FontName', font_name, 'FontWeight', 'normal', ...
        'Units', 'normalized', 'Position', [0.5, 1.08, 0]);

    hold on;
    for j = 1:numel(layout.label)
        idx = find(strcmp(layout.label{j}, channels),1);
        if isempty(idx), continue; end
        pos = layout.pos(j,:);
        if strcmp(mode,'uncorrected')
            if ~isnan(pmap(idx)) && pmap(idx) < 0.05
                scatter(pos(1),pos(2),5,'w','filled','MarkerEdgeColor','k','LineWidth',0.4);
            else
                scatter(pos(1),pos(2),3,[0.45 0.45 0.45],'filled');
            end
        else
            if h_clust(idx)
                scatter(pos(1),pos(2),5,'w','filled','MarkerEdgeColor','k','LineWidth',0.4);
            else
                scatter(pos(1),pos(2),3,[0.45 0.45 0.45],'filled');
            end
        end
    end
    hold off;

    ax = gca; ax.FontSize = font_size; ax.FontName = font_name;
    fig.PaperUnits    = 'inches';
    fig.PaperSize     = [fig_w fig_h];
    fig.PaperPosition = [0 0 fig_w fig_h];
    exportgraphics(fig, fullfile(out_dir,[fname '.png']), 'Resolution',300);
    exportgraphics(fig, fullfile(out_dir,[fname '.pdf']), 'ContentType','vector');
    print(fig, fullfile(out_dir,[fname '.svg']), '-dsvg', '-painters');
    savefig(fig, fullfile(out_dir,[fname '.fig']));
    fprintf('  Saved: %s\n', fname);
    close(fig);
end

%% ========================================================================
%% HELPER: ADD MISSING CHANNELS
%% ========================================================================

function layout = add_missing_channels(layout)
    AFz_pos = mean([layout.pos(strcmp(layout.label,'Fz'),:); ...
                    layout.pos(strcmp(layout.label,'AF3'),:); ...
                    layout.pos(strcmp(layout.label,'AF4'),:)],1);
    FCz_pos = mean([layout.pos(strcmp(layout.label,'Fz'),:); ...
                    layout.pos(strcmp(layout.label,'Cz'),:)],1);
    I1_pos  = layout.pos(strcmp(layout.label,'O1'),:);
    I2_pos  = layout.pos(strcmp(layout.label,'O2'),:);
    layout.label  = [layout.label;  {'AFz';'FCz';'I1';'I2'}];
    layout.pos    = [layout.pos;    AFz_pos; FCz_pos; I1_pos; I2_pos];
    layout.width  = [layout.width;  repmat(mean(layout.width),4,1)];
    layout.height = [layout.height; repmat(mean(layout.height),4,1)];
    center     = mean(layout.pos,1);
    layout.pos = center + 0.85*(layout.pos - center);
end

%% ========================================================================
%% HELPER: BUILD NEIGHBOURS
%% ========================================================================

function neighbours = build_neighbours()
    cfg_neigh        = [];
    cfg_neigh.method = 'triangulation';
    cfg_neigh.layout = 'acticap-64ch-standard2.mat';
    neighbours       = ft_prepare_neighbours(cfg_neigh);
    extras = {'AFz',{'Fz','AF3','AF4','F1','F2'};
              'FCz',{'Fz','Cz','FC1','FC2','FC3','FC4'};
              'I1', {'O1','PO7','PO3','Oz'};
              'I2', {'O2','PO8','PO4','Oz'}};
    for e = 1:size(extras,1)
        idx                         = numel(neighbours)+1;
        neighbours(idx).label       = extras{e,1};
        neighbours(idx).neighblabel = extras{e,2};
    end
    new_chans  = {'AFz','FCz','I1','I2'};
    new_neighs = {{'Fz','AF3','AF4','F1','F2'}, ...
                  {'Fz','Cz','FC1','FC2','FC3','FC4'}, ...
                  {'O1','PO7','PO3','Oz'}, ...
                  {'O2','PO8','PO4','Oz'}};
    for nc = 1:numel(new_chans)
        for nn = 1:numel(new_neighs{nc})
            si = find(strcmp({neighbours.label}, new_neighs{nc}{nn}));
            if ~isempty(si) && ~ismember(new_chans{nc}, neighbours(si).neighblabel)
                neighbours(si).neighblabel{end+1} = new_chans{nc};
            end
        end
    end
end
