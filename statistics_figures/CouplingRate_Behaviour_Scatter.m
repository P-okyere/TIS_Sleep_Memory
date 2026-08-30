%% CouplingRate_Behaviour_Scatter.m
%
% SO-SPINDLE COUPLING RATE x BEHAVIOUR SCATTER PLOTS
% Coupling rate (averaged across significant cluster channels) ->
% Associative recall forgetting, per condition
%
% MATLAB conversion of CouplingRate_Behaviour_Scatter.R - same logic,
% same significant-channel lists, same output structure.
%
% Only conditions with a significant cluster get a scatter panel
% (Cue and TIS+Cue; TIS+Cue delayed had no significant cluster).

close all; clear all;

%% ============================================================
%  EDIT FOR YOUR ENVIRONMENT
%% ============================================================

bids_deriv   = 'E:\BIDS\derivatives';
behav_v2     = fullfile(bids_deriv, 'behaviour', 'desc-FamRecFidSummary_results.xlsx');
rate_results = fullfile(bids_deriv, 'so_spindle_coupling', 'desc-SOSpindleCouplingRate_metrics.mat');

out_dir = fullfile(bids_deriv, 'CouplingRate_Behaviour_Scatter');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
fprintf('Output folder: %s\n', out_dir);

font_size = 10;

% Significant channels, confirmed directly from the topoplot's
% cluster-permutation output (both clusters combined per condition).
% TIS+Cue delayed intentionally omitted - no significant cluster
% (both clusters were only trend-level: p=0.148, p=0.165).
sig_channels_map = containers.Map();
sig_channels_map('Cue')     = {'Pz','CP2','P2','CP4','P4','C6','T8','Fp2','AF8','AF4','AFz'};
sig_channels_map('TIS+Cue') = {'T7','C5','TP9','TP7','FC3','P7','P5','PO3','C4','FC4','C6','FT8','T8','TP8'};

conditions_with_clusters = {'Cue', 'TIS+Cue'};

% Colours - matches COL_COND from the Chapter 3 behaviour script
condition_colors = containers.Map();
condition_colors('TIS+Cue')         = [0.608 0.561 0.831];  % #9B8FD4
condition_colors('TIS+Cue delayed') = [0.290 0.565 0.851];  % #4A90D9
condition_colors('Cue')             = [0.106 0.769 0.769];  % #1BC4C4

%% ============================================================
%  LOAD COUPLING RATE (participants x channels x conditions)
%% ============================================================

rate_data = load(rate_results, 'p_rate_all', 'channels');
p_rate_all = rate_data.p_rate_all;
channels   = rate_data.channels;

participants = { ...
    'sub-02','sub-03','sub-04','sub-05','sub-06','sub-08','sub-12', ...
    'sub-13','sub-14','sub-18','sub-19','sub-20','sub-21','sub-22', ...
    'sub-23','sub-24','sub-25','sub-26','sub-27','sub-29','sub-30', ...
    'sub-33','sub-36','sub-37','sub-38','sub-39','sub-40','sub-41' ...
};
cond_labels = {'TIS+Cue', 'TIS+Cue delayed', 'Cue'};

fprintf('Coupling rate loaded: %d participants x %d channels x %d conditions\n', ...
    size(p_rate_all, 1), size(p_rate_all, 2), size(p_rate_all, 3));

%% ============================================================
%  LOAD BEHAVIOUR - associative recall forgetting (V2)
%  forgetting_pct = 100 - recollection_pct
%% ============================================================

behav = readtable(behav_v2, 'Sheet', 'EEG_Correlation_Master');

%% ============================================================
%  BUILD PER-CONDITION SCATTER DATA + STATS + PLOTS
%% ============================================================

cor_table = table();
plots_data = struct();

for c = 1:length(conditions_with_clusters)
    cond_name = conditions_with_clusters{c};
    cond_idx  = find(strcmp(cond_labels, cond_name));

    sig_ch = sig_channels_map(cond_name);
    ch_idx = find(ismember(channels, sig_ch));

    if isempty(ch_idx)
        error('None of the specified channels matched the loaded channel list for %s.', cond_name);
    end

    fprintf('  %s: averaging over %d significant channels (%s)\n', ...
        cond_name, length(ch_idx), strjoin(channels(ch_idx), ', '));

    % Average coupling rate across sig channels, per participant
    rate_roi = squeeze(mean(p_rate_all(:, ch_idx, cond_idx), 2, 'omitnan'));

    % Match behaviour column to this condition
    switch cond_name
        case 'Cue'
            forgetting = 100 - behav.V2_Cue_recollection_pct;
        case 'TIS+Cue'
            forgetting = 100 - behav.V2_TISCue_recollection_pct;
        case 'TIS+Cue delayed'
            forgetting = 100 - behav.V2_TISCue_delayed_recollection_pct;
    end

    % Match participant order (behav table may not be in the same
    % order as `participants`)
    [~, match_idx] = ismember(participants, behav.participant_id);
    forgetting_matched = forgetting(match_idx);

    keep = ~isnan(rate_roi) & ~isnan(forgetting_matched);
    x = rate_roi(keep);
    y = forgetting_matched(keep);
    n = sum(keep);

    [r, p] = corr(x, y, 'Type', 'Pearson');

    cor_table = [cor_table; table({cond_name}, r, p, n, ...
        'VariableNames', {'Condition', 'r_pearson', 'p_pearson', 'n'})]; %#ok

    plots_data.(matlab.lang.makeValidName(cond_name)) = struct('x', x, 'y', y);
end

disp(cor_table);
writetable(cor_table, fullfile(out_dir, 'CouplingRate_Forgetting_correlations.csv'));

%% ============================================================
%  SCATTER PLOTS - one per condition, with regression line + CI band
%  (manually computed prediction interval, since MATLAB has no direct
%  geom_smooth equivalent)
%% ============================================================

fig_combined = figure('Color', 'w', 'Units', 'inches', ...
    'Position', [1 1 3.2*length(conditions_with_clusters) 3.2]);

for c = 1:length(conditions_with_clusters)
    cond_name  = conditions_with_clusters{c};
    field_name = matlab.lang.makeValidName(cond_name);
    x = plots_data.(field_name).x;
    y = plots_data.(field_name).y;
    col = condition_colors(cond_name);

    % Individual figure
    fig_ind = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 3.2 3.2]);
    plot_scatter_with_ci(x, y, col, cond_name, font_size);
    exportgraphics(fig_ind, fullfile(out_dir, sprintf('scatter_%s.png', strrep(cond_name, '+', '_'))), 'Resolution', 600);
    exportgraphics(fig_ind, fullfile(out_dir, sprintf('scatter_%s.pdf', strrep(cond_name, '+', '_'))), 'ContentType', 'vector');
    print(fig_ind, fullfile(out_dir, sprintf('scatter_%s.svg', strrep(cond_name, '+', '_'))), '-dsvg', '-painters');

    % Add to combined figure
    figure(fig_combined);
    subplot(1, length(conditions_with_clusters), c);
    plot_scatter_with_ci(x, y, col, cond_name, font_size);

    close(fig_ind);
end

exportgraphics(fig_combined, fullfile(out_dir, 'fig_couplingrate_behaviour.png'), 'Resolution', 600);
exportgraphics(fig_combined, fullfile(out_dir, 'fig_couplingrate_behaviour.pdf'), 'ContentType', 'vector');
print(fig_combined, fullfile(out_dir, 'fig_couplingrate_behaviour.svg'), '-dsvg', '-painters');

fprintf('\nDone!\n');
fprintf('Output: %s\n', out_dir);
fprintf('  CouplingRate_Forgetting_correlations.csv - r/p/n per condition\n');
fprintf('  fig_couplingrate_behaviour.png/.pdf/.svg - combined panel\n');
fprintf('  scatter_<condition>.png/.pdf/.svg - individual panels\n');

%% ============================================================
%  HELPER: scatter + linear fit + shaded 95%% CI band, matching the
%  ggplot geom_smooth(method="lm", se=TRUE) look
%% ============================================================

function plot_scatter_with_ci(x, y, col, title_str, font_size)
    hold on;

    % Linear fit
    p_fit = polyfit(x, y, 1);
    x_line = linspace(min(x), max(x), 100)';
    y_line = polyval(p_fit, x_line);

    % Prediction interval (95% CI for the mean prediction, matching
    % ggplot's default confidence band, not a full prediction interval)
    n = length(x);
    y_fit_orig = polyval(p_fit, x);
    residuals = y - y_fit_orig;
    mse = sum(residuals.^2) / (n - 2);
    x_mean = mean(x);
    sxx = sum((x - x_mean).^2);
    se_line = sqrt(mse * (1/n + (x_line - x_mean).^2 / sxx));
    t_crit = tinv(0.975, n - 2);
    ci_upper = y_line + t_crit * se_line;
    ci_lower = y_line - t_crit * se_line;

    % Shaded CI band
    fill([x_line; flipud(x_line)], [ci_upper; flipud(ci_lower)], col, ...
        'FaceAlpha', 0.10, 'EdgeColor', 'none');

    % Regression line
    plot(x_line, y_line, '-', 'Color', col, 'LineWidth', 0.8);

    % Scatter points
    scatter(x, y, 60, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.5, 'MarkerFaceAlpha', 0.9);

    xlabel({'SO-spindle coupling rate (%)', '(mean across significant cluster channels)'}, ...
        'FontSize', font_size);
    ylabel('Associative recall forgetting (%)', 'FontSize', font_size);
    title(title_str, 'FontSize', font_size + 1, 'FontWeight', 'bold');
    set(gca, 'FontSize', font_size, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 0.3);
    hold off;
end
