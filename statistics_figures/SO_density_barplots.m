clear; clc; close all;

%% ============================
%  EDIT FOR YOUR ENVIRONMENT
%% ============================
bids_root  = 'E:\BIDS';
bids_deriv = fullfile(bids_root, 'derivatives');

% Must match out_root/topoplot_significant from
% SO_density_cluster_permutation_stats.m, and fname must match its
% 'fname' variable - this is where that script saved
% <fname>_sig_channels_clusterN.csv
stats_dir = fullfile(bids_deriv, 'SO_density_stats', 'topoplot_significant');
stats_fname_prefix = 'SO_density_Gamma_posonly_clusterperm';

output_folder = fullfile(bids_deriv, 'SO_density_stats', 'barplots');

if ~exist(output_folder,'dir'), mkdir(output_folder); end
% ============================

%% ============================
%  LOAD SIGNIFICANT CHANNELS FROM STATS SCRIPT OUTPUT
%% ============================
% Loads <stats_fname_prefix>_sig_channels_cluster1.csv, _cluster2.csv, ...
% saved by SO_density_cluster_permutation_stats.m — no channel list is
% hardcoded here. As of the last run this covered all 64 channels (a
% whole-scalp effect), so all cluster files found are pooled together.
cluster_files = dir(fullfile(stats_dir, [stats_fname_prefix '_sig_channels_cluster*.csv']));
if isempty(cluster_files)
    error('No significant cluster files found in %s. Run SO_density_cluster_permutation_stats.m first.', stats_dir);
end

sig_channels = {};
for i = 1:numel(cluster_files)
    tbl = readtable(fullfile(cluster_files(i).folder, cluster_files(i).name));
    sig_channels = [sig_channels, tbl.Channel']; %#ok
end
sig_channels = unique(sig_channels);
fprintf('Loaded %d significant channels across %d cluster(s)\n', ...
    numel(sig_channels), numel(cluster_files));

%% ============================
%  LOAD DATA
%% ============================
load(fullfile(bids_deriv, 'spindle_so_detection', 'desc-SOMetrics_trialdata.mat'), ...
     'all_so_trial_data');

%% ============================
%  PARAMETERS
%% ============================
marker_set    = {'MarkerOldSound','TICUEE','TIB4CUEE'};
marker_labels = {'Cue','TIS+Cue','TIS+Cue delayed'};

% Bar colours
bar_colours = [ 0.0941 0.7804 0.8235;   % Cue             — teal
                0.7804 0.5451 0.8980;   % TIS+Cue         — purple
                0.1176 0.5647 1.0000];  % TIS+Cue delayed — blue

pairs = [1 2; 1 3; 2 3];

%% ============================
%  FILTER DATA — positive density only
%% ============================
T = all_so_trial_data( ...
    ismember(all_so_trial_data.SleepStage,  {'N2','N3'}) & ...
    ismember(all_so_trial_data.MarkerType,  marker_set)  & ...
    ismember(all_so_trial_data.Channel,     sig_channels) & ...
    all_so_trial_data.Density > 0, :);

T.MarkerType  = categorical(T.MarkerType, marker_set);
T.Participant = categorical(T.Participant);

fprintf('Positive density trials used: %d\n', height(T));

%% ============================
%  PARTICIPANT-LEVEL MEANS
%% ============================
avg_tbl = groupsummary(T, {'Participant','MarkerType'}, 'mean', 'Density');
avg_tbl.Properties.VariableNames{end} = 'AvgDensity';

participants = categories(avg_tbl.Participant);
nP  = numel(participants);
nC  = numel(marker_set);

plot_data = nan(nP, nC);
for p = 1:nP
    for c = 1:nC
        idx = avg_tbl.Participant == participants{p} & ...
              avg_tbl.MarkerType  == marker_set{c};
        if any(idx)
            plot_data(p,c) = avg_tbl.AvgDensity(idx);
        end
    end
end

group_means = mean(plot_data, 'omitnan');
group_sems  = std(plot_data, 'omitnan') ./ sqrt(sum(~isnan(plot_data)));

%% ============================
%  PAIRWISE PAIRED T-TESTS
%  On participant-level means — df = n-1 = 27
%% ============================
pvals  = nan(3,1);
tstats = nan(3,1);
dfs    = nan(3,1);

pair_combos = {[1 2], [1 3], [2 3]};
pair_names  = {'Cue vs TIS+Cue', 'Cue vs TIS+Cue delayed', 'TIS+Cue vs TIS+Cue delayed'};

fprintf('\n=== Pairwise paired t-tests (participant-level means) ===\n');
for i = 1:3
    c1 = pair_combos{i}(1);
    c2 = pair_combos{i}(2);

    x1 = plot_data(:,c1);
    x2 = plot_data(:,c2);

    % Remove participants missing either condition
    valid = ~isnan(x1) & ~isnan(x2);
    x1 = x1(valid);
    x2 = x2(valid);

    [~, p, ~, stats] = ttest(x1, x2);
    pvals(i)  = p;
    tstats(i) = stats.tstat;
    dfs(i)    = stats.df;

    fprintf('%s: t(%d) = %.3f, p = %.4f\n', ...
        pair_names{i}, stats.df, stats.tstat, p);
end

%% ============================
%  PLOT
%% ============================
fig = figure('Color','w','Units','inches','Position',[1 1 1.5 2.0]);
b = bar(group_means, 'FaceColor','flat'); hold on;
b.CData = bar_colours;

% Connecting lines for each participant across conditions - matches
% the behavioural plots (Roi/Derk-Jan noted the inconsistency: lines
% were present for behaviour but not EEG panels).
for p = 1:nP
    y_line = plot_data(p,:);
    valid_line = ~isnan(y_line);
    if sum(valid_line) > 1
        plot(find(valid_line), y_line(valid_line), '-', ...
            'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 0.45);
    end
end

errorbar(1:nC, group_means, group_sems, ...
    'k.', 'LineWidth',1.6, 'MarkerSize',1);

for c = 1:nC
    y = plot_data(:,c);
    y = y(~isnan(y));
    x = c + (rand(size(y)) - 0.5) * 0.25;
    scatter(x, y, 10, bar_colours(c,:), ...
        'filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.65);
end

% Significance brackets
y_data_max = max(plot_data(:), [], 'omitnan');
y_bar_max  = max(group_means + group_sems);
y_base     = max([y_data_max, y_bar_max]);
step       = 0.07 * y_base;
h_tick     = 0.02 * y_base;

for i = 1:3
    if ~isnan(pvals(i)) && pvals(i) < 0.05
        x1 = pairs(i,1); x2 = pairs(i,2);
        y  = y_base + i * step;
        plot([x1 x1 x2 x2],[y y+h_tick y+h_tick y],'k','LineWidth',1.0);
        if     pvals(i) < 0.001, s = '***';
        elseif pvals(i) < 0.01,  s = '**';
        elseif pvals(i) < 0.05,  s = '*';
        else,                    s = '†';
        end
        text(mean([x1 x2]), y+h_tick*1.4, s, ...
            'HorizontalAlignment','center','FontSize',7);
    end
end

ylim([0, y_base + 4*step])
xlim([0.5, nC + 0.5])
xticks(1:nC); xticklabels(marker_labels)
ylabel('SO density (SOs/s)', 'FontSize',7)

ax = gca;
ax.Box       = 'off';
ax.LineWidth = 1.4;
ax.XAxis.TickLength = [0.02 0.02];
ax.FontSize  = 7;

%% ============================
%  SAVE
%% ============================
if ~isvalid(fig)
    error('fig is not a valid figure handle - make sure the whole figure-creation block above ran without any intervening close/clear.');
end

set(fig,'PaperPositionMode','auto','Renderer','painters')
print(fig, fullfile(output_folder,'SO_density_pairedttest_sigchannels.png'), '-dpng','-r600');
print(fig, fullfile(output_folder,'SO_density_pairedttest_sigchannels.pdf'), '-dpdf','-painters');
print(fig, fullfile(output_folder,'SO_density_pairedttest_sigchannels.svg'), '-dsvg','-painters');
savefig(fig, fullfile(output_folder,'SO_density_pairedttest_sigchannels.fig'));
fprintf('\n✅ Saved to: %s\n', output_folder);
