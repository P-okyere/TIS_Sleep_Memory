clear; clc; close all;

%% ============================
%  EDIT FOR YOUR ENVIRONMENT
%% ============================
base_path = '/path/to/MYSTI';

% Must match out_root/spindle_type/topoplot_significant from
% fast_spindle_cluster_permutation_stats.m, and fname must match its
% 'fname' variable - this is where that script saved
% <fname>_sig_channels_clusterN.csv
stats_dir = fullfile(base_path, 'derivatives', 'spindle_amplitude_stats', 'fast', 'topoplot_significant');
stats_fname_prefix = 'fast_absolute_ANOVA_clusterperm';

output_folder = fullfile(base_path, 'derivatives', 'spindle_amplitude_stats', 'fast', 'barplots');

if ~exist(output_folder,'dir'), mkdir(output_folder); end
% ============================

%% ============================
%  LOAD SIGNIFICANT CLUSTERS FROM f2.m OUTPUT
%% ============================
% Loads <stats_fname_prefix>_sig_channels_cluster1.csv, _cluster2.csv, ...
% saved by fast_spindle_cluster_permutation_stats.m — no channel lists are hardcoded here.
cluster_files = dir(fullfile(stats_dir, [stats_fname_prefix '_sig_channels_cluster*.csv']));
if isempty(cluster_files)
    error('No significant cluster files found in %s. Run fast_spindle_cluster_permutation_stats.m first.', stats_dir);
end

% Sort numerically by cluster number in the filename
cluster_nums = zeros(numel(cluster_files),1);
for i = 1:numel(cluster_files)
    tok = regexp(cluster_files(i).name, 'cluster(\d+)\.csv$', 'tokens');
    cluster_nums(i) = str2double(tok{1}{1});
end
[~, ord] = sort(cluster_nums);
cluster_files = cluster_files(ord);

n_clusters   = numel(cluster_files);
cluster_list = cell(n_clusters, 1);
cluster_pval = nan(n_clusters, 1);

for i = 1:n_clusters
    tbl = readtable(fullfile(cluster_files(i).folder, cluster_files(i).name));
    cluster_list{i} = tbl.Channel';
    cluster_pval(i) = tbl.ClusterPValue(1);
end

combined_chans = unique(horzcat(cluster_list{:}));
cluster_list{end+1} = combined_chans; % combined panel = last entry

% Auto-generated titles from the loaded p-values. Edit these labels by
% hand if you want anatomical descriptions (e.g. "Left Temporal/Central")
% once you know which cluster is which from the topoplot.
cluster_titles = cell(n_clusters+1, 1);
for i = 1:n_clusters
    cluster_titles{i} = {sprintf('Fast spindle amplitude'), ...
                          sprintf('Cluster %d (cluster p = %.3f)', i, cluster_pval(i))};
end
cluster_titles{end} = {'Fast spindle amplitude','All significant channels'};

cluster_fnames = cell(n_clusters+1, 1);
for i = 1:n_clusters
    cluster_fnames{i} = sprintf('fast_spindle_cluster%d', i);
end
cluster_fnames{end} = 'fast_spindle_all_clusters_combined';

n_panels = n_clusters + 1;

%% ============================
%  LOAD DATA
%% ============================
load(fullfile(base_path,'trial_density_all_conditions_clean3.5s.mat'), ...
     'all_trial_data');

%% ============================
%  PARAMETERS
%% ============================
spindle_type  = 'fast';
measure_field = 'MeanAmplitude';

marker_set    = {'MarkerOldSound','TICUEE','TIB4CUEE'};
marker_labels = {'Cue','TIS+Cue','TIS+Cue delayed'};

% Bar colours
bar_colours = [ 0.0941 0.7804 0.8235;
                0.7804 0.5451 0.8980;
                0.1176 0.5647 1.0000];

pairs      = [1 2; 1 3; 2 3];
pair_names = {'Cue vs TIS+Cue','Cue vs TIS+Cue delayed','TIS+Cue vs TIS+Cue delayed'};

%% ============================
%  LOOP OVER CLUSTERS — individual plots
%% ============================
all_plot_data  = cell(n_panels,1);
all_pvals      = cell(n_panels,1);
all_tstats     = cell(n_panels,1);
all_dfs        = cell(n_panels,1);
all_means      = cell(n_panels,1);
all_sems       = cell(n_panels,1);

for cl = 1:n_panels

    sig_channels = cluster_list{cl};

    %% Filter
    T = all_trial_data( ...
        strcmp(all_trial_data.SpindleType, spindle_type) & ...
        ismember(all_trial_data.SleepStage, {'N2','N3'}) & ...
        ismember(all_trial_data.MarkerType, marker_set) & ...
        ismember(all_trial_data.Channel,    sig_channels) & ...
        ~isnan(all_trial_data.(measure_field)), :);

    T.MarkerType  = categorical(T.MarkerType, marker_set);
    T.Participant = categorical(T.Participant);

    %% Participant-level means
    avg_tbl = groupsummary(T, {'Participant','MarkerType'}, 'mean', measure_field);
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

    %% Paired t-tests on participant-level means
    pvals  = nan(3,1);
    tstats = nan(3,1);
    dfs    = nan(3,1);

    fprintf('\n--- Cluster %d ---\n', cl);
    for i = 1:3
        x1 = plot_data(:, pairs(i,1));
        x2 = plot_data(:, pairs(i,2));
        valid = ~isnan(x1) & ~isnan(x2);
        [~, pvals(i), ~, stats] = ttest(x1(valid), x2(valid));
        tstats(i) = stats.tstat;
        dfs(i)    = stats.df;
        fprintf('  %s: t(%d)=%.3f, p=%.4f\n', pair_names{i}, dfs(i), tstats(i), pvals(i));
    end

    % Store for panel plot
    all_plot_data{cl} = plot_data;
    all_pvals{cl}     = pvals;
    all_tstats{cl}    = tstats;
    all_dfs{cl}       = dfs;
    all_means{cl}     = group_means;
    all_sems{cl}      = group_sems;

    %% Plot
    fig = figure('Color','w','Units','inches','Position',[1 1 2.0 2.5]);
    b = bar(group_means, 'FaceColor','flat'); hold on;
    b.CData = bar_colours;

    errorbar(1:nC, group_means, group_sems, ...
        'k.', 'LineWidth',1.6, 'MarkerSize',1);

    for c = 1:nC
        y = plot_data(:,c);
        y = y(~isnan(y));
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
        if pvals(i) < 0.05
            x1 = pairs(i,1); x2 = pairs(i,2);
            yy = y_base + i * step;
            plot([x1 x1 x2 x2],[yy yy+h_tick yy+h_tick yy],'k','LineWidth',1.0);
            if     pvals(i) < 0.001, s = '***';
            elseif pvals(i) < 0.01,  s = '**';
            elseif pvals(i) < 0.05,  s = '*';
            else,                    s = '†';
            end
            text(mean([x1 x2]), yy+h_tick*1.4, s, ...
                'HorizontalAlignment','center','FontSize',10);
        end
    end

    ylim([0, y_base + 4*step])
    xlim([0.5, nC + 0.5])
    xticks(1:nC); xticklabels(marker_labels)
    ylabel('Fast spindle amplitude (µV)', 'FontSize',10)


    ax = gca;
    ax.Box       = 'off';
    ax.LineWidth = 1.4;
    ax.XAxis.TickLength = [0.02 0.02];
    ax.FontSize  = 10;

    set(fig,'PaperPositionMode','auto','Renderer','painters')
    print(fig, fullfile(output_folder,[cluster_fnames{cl} '.png']), '-dpng','-r600');
    print(fig, fullfile(output_folder,[cluster_fnames{cl} '.pdf']), '-dpdf','-painters');
    fprintf('✅ Saved: %s\n', cluster_fnames{cl});

end

%% ============================
%  PLOT 4 — Side-by-side panel
%% ============================
panel_titles = {
    {'Left Temporal/Central','(cluster p < .001)'}, ...
    {'Right Temporal/Central','(cluster p < .001)'}, ...
    {'All Significant Channels','(both clusters p < .001)'} ...
};

fig_panel = figure('Color','w','Units','inches','Position',[1 1 3.5*n_panels 3.5]);

for cl = 1:n_panels

    plot_data   = all_plot_data{cl};
    pvals       = all_pvals{cl};
    group_means = all_means{cl};
    group_sems  = all_sems{cl};
    nC          = numel(marker_set);

    subplot(1,n_panels,cl);
    b = bar(group_means,'FaceColor','flat'); hold on;
    b.CData = bar_colours;

    errorbar(1:nC, group_means, group_sems, ...
        'k.','LineWidth',1.6,'MarkerSize',1);

    for c = 1:nC
        y = plot_data(:,c); y = y(~isnan(y));
        x = c + (rand(size(y))-0.5)*0.25;
        scatter(x, y, 8, bar_colours(c,:), ...
            'filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.65);
    end

    y_data_max = max(plot_data(:),[],'omitnan');
    y_bar_max  = max(group_means + group_sems);
    y_base     = max([y_data_max, y_bar_max]);
    step       = 0.07 * y_base;
    h_tick     = 0.02 * y_base;

    for i = 1:3
        if pvals(i) < 0.05
            x1 = pairs(i,1); x2 = pairs(i,2);
            yy = y_base + i * step;
            plot([x1 x1 x2 x2],[yy yy+h_tick yy+h_tick yy],'k','LineWidth',1.4);
            if     pvals(i) < 0.001, s = '***';
            elseif pvals(i) < 0.01,  s = '**';
            elseif pvals(i) < 0.05,  s = '*';
            else,                    s = '†';
            end
            text(mean([x1 x2]), yy+h_tick*1.4, s, ...
                'HorizontalAlignment','center','FontSize',11);
        end
    end

    ylim([0, y_base + 4*step])
    xlim([0.5, nC + 0.5])
    xticks(1:nC); xticklabels(marker_labels)
    ylabel('Fast spindle amplitude (µV)','FontSize',11)

    ax = gca; ax.Box = 'off'; ax.LineWidth = 1.4; ax.FontSize = 10;

end

set(fig_panel,'PaperPositionMode','auto','Renderer','painters')
print(fig_panel, fullfile(output_folder,'fast_spindle_ALL_panels.png'), '-dpng','-r600');
print(fig_panel, fullfile(output_folder,'fast_spindle_ALL_panels.pdf'), '-dpdf','-painters');
fprintf('✅ All panels combined saved.\n');

fprintf('\n✅ All figures saved to:\n%s\n', output_folder);