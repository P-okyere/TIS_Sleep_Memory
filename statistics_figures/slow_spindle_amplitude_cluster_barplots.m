clear; clc; close all;

%% ============================
%  EDIT FOR YOUR ENVIRONMENT
%% ============================
base_path = '/path/to/MYSTI';

% Must match out_root/spindle_type/topoplot_significant from
% slow_spindle_cluster_permutation_stats.m, and fname must match its
% 'fname' variable - this is where that script saved
% <fname>_sig_channels_clusterN.csv
stats_dir = fullfile(base_path, 'derivatives', 'spindle_amplitude_stats', 'slow', 'topoplot_significant');
stats_fname_prefix = 'slow_absolute_ANOVA_clusterperm';

output_folder = fullfile(base_path, 'derivatives', 'spindle_amplitude_stats', 'slow', 'barplots');

if ~exist(output_folder,'dir'), mkdir(output_folder); end
% ============================

%% ============================
%  LOAD SIGNIFICANT CLUSTERS FROM STATS SCRIPT OUTPUT
%% ============================
% Loads <stats_fname_prefix>_sig_channels_cluster1.csv, _cluster2.csv, ...
% saved by slow_spindle_cluster_permutation_stats.m — no channel lists
% or p-values are hardcoded here.
cluster_files = dir(fullfile(stats_dir, [stats_fname_prefix '_sig_channels_cluster*.csv']));
if isempty(cluster_files)
    error('No significant cluster files found in %s. Run slow_spindle_cluster_permutation_stats.m first.', stats_dir);
end

cluster_nums = zeros(numel(cluster_files),1);
for i = 1:numel(cluster_files)
    tok = regexp(cluster_files(i).name, 'cluster(\d+)\.csv$', 'tokens');
    cluster_nums(i) = str2double(tok{1}{1});
end
[~, ord] = sort(cluster_nums);
cluster_files = cluster_files(ord);

n_clusters = numel(cluster_files);
all_chan_sets = cell(n_clusters, 1);
cluster_pvals = nan(n_clusters, 1);

for i = 1:n_clusters
    tbl = readtable(fullfile(cluster_files(i).folder, cluster_files(i).name));
    all_chan_sets{i} = tbl.Channel';
    cluster_pvals(i) = tbl.ClusterPValue(1);
end

combined_chans = unique(horzcat(all_chan_sets{:}));
all_chan_sets{end+1} = combined_chans; % combined panel = last entry
cluster_pvals(end+1) = NaN;

% Auto-generated labels from cluster number + p-value. Edit by hand for
% anatomical descriptions (e.g. "Central/Centroparietal") once you know
% which cluster is which from the topoplot.
cluster_labels = cell(n_clusters+1, 1);
for i = 1:n_clusters
    cluster_labels{i} = sprintf('Cluster %d (p = %.3f)', i, cluster_pvals(i));
end
cluster_labels{end} = 'All Significant Channels';

n_panels = n_clusters + 1;

%% ============================
%  LOAD DATA
%% ============================
load(fullfile(base_path, ...
    'trial_density_slowspindle_conditions_clean3.5s.mat'), ...
    'all_trial_data');

%% ============================
%  PARAMETERS
%% ============================
spindle_type  = 'slow';
measure_field = 'MeanAmplitude';

marker_set    = {'MarkerOldSound','TICUEE','TIB4CUEE'};
marker_labels = {'Cue','TIS+Cue','TIS+Cue delayed'};

bar_colours = [ 0.0941 0.7804 0.8235;
                0.7804 0.5451 0.8980;
                0.1176 0.5647 1.0000];

pair_list = {
    {'MarkerOldSound','TICUEE'}
    {'MarkerOldSound','TIB4CUEE'}
    {'TICUEE','TIB4CUEE'}
};
pairs = [1 2; 1 3; 2 3];

all_stats     = [];

%% ============================
%  LOOP OVER ALL 3 PLOTS
%% ============================
for cl = 1:n_panels

    sig_channels = all_chan_sets{cl};
    cl_label     = cluster_labels{cl};

    fprintf('\n=== %s ===\n', cl_label);

    T = all_trial_data( ...
        strcmp(all_trial_data.SpindleType, spindle_type) & ...
        ismember(all_trial_data.SleepStage, {'N2','N3'}) & ...
        ismember(all_trial_data.MarkerType, marker_set) & ...
        ismember(all_trial_data.Channel, sig_channels) & ...
        ~isnan(all_trial_data.(measure_field)), :);

    T.MarkerType  = categorical(T.MarkerType, marker_set);
    T.Participant = categorical(T.Participant);

    avg_tbl = groupsummary(T, {'Participant','MarkerType'}, ...
                           'mean', measure_field);
    avg_tbl.Properties.VariableNames{end} = 'AvgValue';

    participants = categories(avg_tbl.Participant);
    nP = numel(participants);
    nC = numel(marker_set);

    plot_data = nan(nP, nC);
    for p = 1:nP
        for c = 1:nC
            idx = avg_tbl.Participant == participants{p} & ...
                  avg_tbl.MarkerType  == marker_set{c};
            if any(idx), plot_data(p,c) = avg_tbl.AvgValue(idx); end
        end
    end

    group_means = mean(plot_data, 'omitnan');
    group_sems  = std(plot_data,  'omitnan') ./ sqrt(sum(~isnan(plot_data)));

    full_lme  = fitlme(avg_tbl, 'AvgValue ~ MarkerType + (1|Participant)');
    anova_res = anova(full_lme);
    fprintf('Omnibus F(%d,%d) = %.3f, p = %.4f\n', ...
        anova_res.DF1(2), anova_res.DF2(2), ...
        anova_res.FStat(2), anova_res.pValue(2));

    pvals = nan(3,1);
    for i = 1:3
        Tpair = avg_tbl(ismember(avg_tbl.MarkerType, pair_list{i}), :);
        Tpair.MarkerType = removecats(Tpair.MarkerType);
        lme   = fitlme(Tpair, 'AvgValue ~ MarkerType + (1|Participant)');
        pvals(i) = lme.Coefficients.pValue(2);
        fprintf('%s vs %s: t(%d)=%.3f, p=%.4f\n', ...
            marker_labels{pairs(i,1)}, marker_labels{pairs(i,2)}, ...
            lme.Coefficients.DF(2), lme.Coefficients.tStat(2), pvals(i));

        all_stats = [all_stats; struct( ...
            'Cluster',    cl_label, ...
            'Comparison', sprintf('%s vs %s', marker_labels{pairs(i,1)}, ...
                                              marker_labels{pairs(i,2)}), ...
            'Estimate',   lme.Coefficients.Estimate(2), ...
            'SE',         lme.Coefficients.SE(2), ...
            'tStat',      lme.Coefficients.tStat(2), ...
            'DF',         lme.Coefficients.DF(2), ...
            'pValue',     lme.Coefficients.pValue(2))]; %#ok
    end

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
        y = plot_data(:,c); y = y(~isnan(y));
        x = c + (rand(size(y)) - 0.5) * 0.25;
        scatter(x, y, 10, bar_colours(c,:), ...
            'filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.65);
    end

    y_data_max = max(plot_data(:), [], 'omitnan');
    y_bar_max  = max(group_means + group_sems);
    y_base     = max([y_data_max, y_bar_max]);
    step       = 0.07 * y_base;
    h_tick     = 0.02 * y_base;

    for i = 1:3
        if pvals(i) < 0.10
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

    ylim([0, y_base + 4.5*step])
    xlim([0.5, nC + 0.5])
    xticks(1:nC); xticklabels(marker_labels)
    ylabel('Slow spindle amplitude (µV)', 'FontSize',7)
    % NO TITLE (single-panel figures are untitled by design)

    ax = gca;
    ax.Box       = 'off';
    ax.LineWidth = 1.4;
    ax.XAxis.TickLength = [0.02 0.02];
    ax.FontSize  = 7;

    cl_safe = strrep(strrep(cl_label,'/','_'),' ','_');
    set(fig,'PaperPositionMode','auto','Renderer','painters')
    print(fig, fullfile(output_folder, ...
        sprintf('slow_spindle_cluster%d_%s.png', cl, cl_safe)), '-dpng','-r600');
    print(fig, fullfile(output_folder, ...
        sprintf('slow_spindle_cluster%d_%s.pdf', cl, cl_safe)), '-dpdf','-painters');
    print(fig, fullfile(output_folder, ...
        sprintf('slow_spindle_cluster%d_%s.svg', cl, cl_safe)), '-dsvg','-painters');
    savefig(fig, fullfile(output_folder, ...
        sprintf('slow_spindle_cluster%d_%s.fig', cl, cl_safe)));
    fprintf('✅ Plot %d saved.\n', cl);

end

%% ============================
%  COMBINED PANEL (3 side by side) — NO TITLES
%% ============================
fprintf('\nGenerating 3-panel combined figure...\n');

fig_panel = figure('Color','w','Units','inches','Position',[1 1 3.5*n_panels 3.5]);

for cl = 1:n_panels

    sig_channels = all_chan_sets{cl};

    T = all_trial_data( ...
        strcmp(all_trial_data.SpindleType, spindle_type) & ...
        ismember(all_trial_data.SleepStage, {'N2','N3'}) & ...
        ismember(all_trial_data.MarkerType, marker_set) & ...
        ismember(all_trial_data.Channel, sig_channels) & ...
        ~isnan(all_trial_data.(measure_field)), :);

    T.MarkerType  = categorical(T.MarkerType, marker_set);
    T.Participant = categorical(T.Participant);

    avg_tbl = groupsummary(T, {'Participant','MarkerType'}, ...
                           'mean', measure_field);
    avg_tbl.Properties.VariableNames{end} = 'AvgValue';

    participants = categories(avg_tbl.Participant);
    nP = numel(participants);
    nC = numel(marker_set);

    plot_data = nan(nP, nC);
    for p = 1:nP
        for c = 1:nC
            idx = avg_tbl.Participant == participants{p} & ...
                  avg_tbl.MarkerType  == marker_set{c};
            if any(idx), plot_data(p,c) = avg_tbl.AvgValue(idx); end
        end
    end

    group_means = mean(plot_data, 'omitnan');
    group_sems  = std(plot_data,  'omitnan') ./ sqrt(sum(~isnan(plot_data)));

    pvals = nan(3,1);
    for i = 1:3
        Tpair = avg_tbl(ismember(avg_tbl.MarkerType, pair_list{i}), :);
        Tpair.MarkerType = removecats(Tpair.MarkerType);
        lme   = fitlme(Tpair, 'AvgValue ~ MarkerType + (1|Participant)');
        pvals(i) = lme.Coefficients.pValue(2);
    end

    subplot(1,n_panels,cl);
    b = bar(group_means,'FaceColor','flat'); hold on;
    b.CData = bar_colours;

    % Connecting lines for each participant - matches behavioural plots
    for p = 1:nP
        y_line = plot_data(p,:);
        valid_line = ~isnan(y_line);
        if sum(valid_line) > 1
            plot(find(valid_line), y_line(valid_line), '-', ...
                'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 0.45);
        end
    end

    errorbar(1:nC, group_means, group_sems, ...
        'k.','LineWidth',1.6,'MarkerSize',1);

    for c = 1:nC
        y = plot_data(:,c); y = y(~isnan(y));
        x = c + (rand(size(y))-0.5)*0.25;
        scatter(x, y, 20, bar_colours(c,:), ...
            'filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.65);
    end

    y_data_max = max(plot_data(:),[],'omitnan');
    y_bar_max  = max(group_means + group_sems);
    y_base     = max([y_data_max, y_bar_max]);
    step       = 0.07 * y_base;
    h_tick     = 0.02 * y_base;

    for i = 1:3
        if pvals(i) < 0.10
            x1 = pairs(i,1); x2 = pairs(i,2);
            y  = y_base + i * step;
            plot([x1 x1 x2 x2],[y y+h_tick y+h_tick y],'k','LineWidth',1.4);
            if     pvals(i) < 0.001, s = '***';
            elseif pvals(i) < 0.01,  s = '**';
            elseif pvals(i) < 0.05,  s = '*';
            else,                    s = '†';
            end
            text(mean([x1 x2]), y+h_tick*1.4, s, ...
                'HorizontalAlignment','center','FontSize',11);
        end
    end

    ylim([0, y_base + 4.5*step])
    xlim([0.5, nC + 0.5])
    xticks(1:nC); xticklabels(marker_labels)
    ylabel('Slow spindle amplitude (µV)','FontSize',8)

    title(cluster_labels{cl}, 'FontSize', 10, 'FontWeight', 'normal');
    ax = gca; ax.Box = 'off'; ax.LineWidth = 1.4;
    ax.XAxis.TickLength = [0.02 0.02];

end

if ~isvalid(fig_panel)
    error('fig_panel is not a valid figure handle - make sure you run the whole panel-creation block (figure(...) through the cl loop) in one go, without closing figures in between.');
end

set(fig_panel,'PaperPositionMode','auto','Renderer','painters')
print(fig_panel, fullfile(output_folder,'slow_spindle_ALL3_panels.png'), '-dpng','-r600');
print(fig_panel, fullfile(output_folder,'slow_spindle_ALL3_panels.pdf'), '-dpdf','-painters');
print(fig_panel, fullfile(output_folder,'slow_spindle_ALL3_panels.svg'), '-dsvg','-painters');
savefig(fig_panel, fullfile(output_folder,'slow_spindle_ALL3_panels.fig'));
fprintf('✅ 3-panel combined figure saved.\n');

%% ============================
%  SAVE ALL STATS
%% ============================
writetable(struct2table(all_stats), ...
    fullfile(output_folder,'pairwise_statistics_all_clusters.csv'));

fprintf('\n✅ All done. Saved to:\n%s\n', output_folder);