%% s03 — Feature collinearity reduction (hierarchical clustering)
%
% Computes a Spearman correlation matrix across all 51 EEG features,
% performs hierarchical agglomerative clustering on the distance matrix
% (d = 1 - |r|), cuts the dendrogram at CLUSTER_THRESHOLD, then selects
% the one feature per cluster whose correlation profile is closest to the
% cluster centroid.
%
% Input : EEG_Features_Behaviour_RelBase.csv  (from s02)
% Output: Results_SelectedFeatures_RelBase.csv  (metadata + representatives)
%         Feature_Clustering_Results.mat
%         Feature_Clustering_RelBase.svg        (Supplementary figure)
%
% Author: Junheng Li

clear all; clc; close all;

%% ========================================================================
%  USER SETTINGS
%  ========================================================================

CLUSTER_THRESHOLD = 0.5;        % distance threshold for cluster cutting
                                % distance = 1 - |r|, so 0.5 <-> |r| > 0.5

CORR_METHOD       = 'Spearman'; % 'Pearson' or 'Spearman'

LINKAGE_METHOD    = 'average';  % 'average' (UPGMA) | 'complete' | 'ward'
                                % 'average' is robust for correlation-based clustering

%  Columns to EXCLUDE from feature analysis (metadata / behavioural)
META_COLS = {'global_trial_idx','trial_idx','subject','condition', ...
             'word_acc_post','imag_acc_post','old_new','word', ...
             'block_post','NewManscore1'};

INPUT_FILE    = 'EEG_Features_Behaviour_RelBase.csv';
OUT_STEM      = 'Results_SelectedFeatures_RelBase';
FIGURE_OUTPUT = 'Feature_Clustering_RelBase.svg';

%% ========================================================================
%  LOAD DATA
%  ========================================================================

results_tbl = readtable(INPUT_FILE);
table_label = 'Relative (post-stim - baseline)';
fprintf('Loaded %s: %d rows.\n', INPUT_FILE, height(results_tbl));

%% ========================================================================
%  EXTRACT NUMERIC EEG FEATURE COLUMNS
%  ========================================================================

all_cols     = results_tbl.Properties.VariableNames;
feature_mask = ~ismember(all_cols, META_COLS);

% Keep only numeric columns that are not metadata
feature_cols = {};
for c = 1:length(all_cols)
    if feature_mask(c) && isnumeric(results_tbl.(all_cols{c}))
        feature_cols{end+1} = all_cols{c}; %#ok<SAGROW>
    end
end

fprintf('Found %d EEG feature columns.\n', length(feature_cols));

% Build feature matrix [trials x features], remove rows with any NaN
feat_mat   = table2array(results_tbl(:, feature_cols));
valid_rows = all(isfinite(feat_mat), 2);
feat_mat   = feat_mat(valid_rows, :);
fprintf('Using %d / %d trials (rows with complete data).\n', ...
    sum(valid_rows), length(valid_rows));

% Remove features with zero variance (constant columns)
feat_std     = std(feat_mat, 0, 1);
nonconst     = feat_std > 0;
feat_mat     = feat_mat(:, nonconst);
feature_cols = feature_cols(nonconst);
fprintf('%d features after removing zero-variance columns.\n', length(feature_cols));

n_features = length(feature_cols);

%% ========================================================================
%  CORRELATION MATRIX
%  ========================================================================

fprintf('Computing %s correlation matrix (%d features)...\n', ...
    CORR_METHOD, n_features);

if strcmpi(CORR_METHOD, 'Spearman')
    [R, P] = corr(feat_mat, 'Type', 'Spearman', 'Rows', 'complete');
else
    [R, P] = corr(feat_mat, 'Type', 'Pearson',  'Rows', 'complete');
end

% Distance matrix for clustering: distance = 1 - |r|
D = 1 - abs(R);
% Ensure symmetry and zero diagonal (guard against floating point drift)
D = (D + D') / 2;
D(1:n_features+1:end) = 0;

%% ========================================================================
%  HIERARCHICAL CLUSTERING
%  ========================================================================

fprintf('Running hierarchical clustering (linkage: %s)...\n', LINKAGE_METHOD);

% squareform converts square distance matrix -> condensed vector form
Z          = linkage(squareform(D, 'tovector'), LINKAGE_METHOD);
cluster_id = cluster(Z, 'Cutoff', CLUSTER_THRESHOLD, 'Criterion', 'distance');

n_clusters = max(cluster_id);
fprintf('Cut at distance %.2f -> %d clusters.\n', CLUSTER_THRESHOLD, n_clusters);

%% ========================================================================
%  SELECT REPRESENTATIVE FEATURE PER CLUSTER (closest to centroid)
%  ========================================================================
% Centroid = mean correlation profile across all members of the cluster.
% Proximity = Euclidean distance in correlation-profile space.

selected_features = cell(n_clusters, 1);
selected_indices  = NaN(n_clusters, 1);
cluster_sizes     = NaN(n_clusters, 1);

for k = 1:n_clusters
    members = find(cluster_id == k);
    cluster_sizes(k) = length(members);

    if length(members) == 1
        % Singleton cluster - the only member is automatically representative
        selected_indices(k)  = members(1);
        selected_features{k} = feature_cols{members(1)};
        continue;
    end

    % Centroid: mean absolute-correlation profile restricted to cluster members
    centroid = mean(abs(R(members, :)), 1);  % 1 x n_features

    % Distance of each member from centroid
    dists_to_centroid = NaN(length(members), 1);
    for mi = 1:length(members)
        dists_to_centroid(mi) = norm(abs(R(members(mi), :)) - centroid);
    end

    [~, closest] = min(dists_to_centroid);
    selected_indices(k)  = members(closest);
    selected_features{k} = feature_cols{members(closest)};
end

%% ========================================================================
%  CONSOLE SUMMARY
%  ========================================================================

fprintf('\n========================================\n');
fprintf('SELECTED FEATURES  (%d clusters -> %d representatives)\n', ...
    n_clusters, n_clusters);
fprintf('========================================\n');
for k = 1:n_clusters
    members     = find(cluster_id == k);
    member_list = strjoin(feature_cols(members), ', ');
    fprintf('Cluster %2d  [n=%d]  -> %s\n       members: %s\n\n', ...
        k, cluster_sizes(k), selected_features{k}, member_list);
end

%% ========================================================================
%  REORDER FEATURES ACCORDING TO DENDROGRAM LEAF ORDER
%  ========================================================================

leaf_order     = optimalleaforder(Z, squareform(D, 'tovector'));
R_ordered      = R(leaf_order, leaf_order);
labels_ordered = feature_cols(leaf_order);
cid_ordered    = cluster_id(leaf_order);

% Boolean mask: is this (reordered) feature the selected representative?
is_selected = ismember(leaf_order, selected_indices);

%% ========================================================================
%  FIGURE
%  ========================================================================
% ---- Colour palette (black-and-white theme) ----
bg_color      = [1.00 1.00 1.00];   % white background
ax_text_color = [0.10 0.10 0.10];   % near-black text
cmap_heat     = flipud(gray(256));  % white = low |r|, black = high |r|
cluster_box_c = [0.00 0.00 0.00];   % black cluster boundary squares
centroid_box_c= [0.15 0.40 0.80];   % blue diagonal box for representatives
dendro_line_c = [0.30 0.30 0.30];   % dark-grey dendrogram lines
thresh_line_c = [0.80 0.10 0.10];   % red threshold line

% ---- Prepare feature labels: replace underscores with spaces ----
labels_dendro = cellfun(@(s) strrep(s, '_', ' '), feature_cols, 'UniformOutput', false);

% ---- Figure geometry ----
fig = figure('Color', bg_color, 'Units', 'normalized', ...
    'Position', [0.02 0.02 0.96 0.94]);
set(fig, 'DefaultAxesFontName', 'Helvetica', ...
         'DefaultTextFontName', 'Helvetica', ...
         'DefaultColorbarFontName', 'Helvetica');

% Axes positions: [left bottom width height] in normalised units
dendro_ax = axes('Parent', fig, 'Position', [0.01 0.08 0.13 0.82]);
heat_ax   = axes('Parent', fig, 'Position', [0.15 0.08 0.62 0.82]);
cbar_ax   = axes('Parent', fig, 'Position', [0.79 0.08 0.020 0.82]);

% ---- Dendrogram ----
axes(dendro_ax);
dendrogram(Z, 0, ...
    'Labels',         labels_dendro, ...
    'Orientation',    'left', ...
    'ColorThreshold', 0, ...          % single colour, we style manually
    'Reorder',        leaf_order);

% Uniform dark-grey lines, no coloured clusters in dendrogram
hl = findobj(dendro_ax, 'Type', 'Line');
set(hl, 'Color', dendro_line_c, 'LineWidth', 0.9);

set(dendro_ax, 'Color', bg_color, ...
    'XColor', ax_text_color, 'YColor', 'none', ...
    'TickDir', 'out', 'FontSize', 7, 'XDir', 'reverse', ...
    'Box', 'off');
hold(dendro_ax, 'on');

% Threshold dashed line
plot(dendro_ax, [CLUSTER_THRESHOLD CLUSTER_THRESHOLD], ...
    [0.5 n_features+0.5], '--', 'Color', thresh_line_c, 'LineWidth', 1.2);

xlabel(dendro_ax, 'Distance', 'Color', ax_text_color, 'FontSize', 8);
title(dendro_ax, sprintf('%.1f', CLUSTER_THRESHOLD), ...
    'Color', thresh_line_c, 'FontSize', 8, 'FontWeight', 'normal');

% ---- Heatmap (absolute correlation values for B&W legibility) ----
axes(heat_ax);
imagesc(heat_ax, abs(R_ordered), [0 1]);
colormap(heat_ax, cmap_heat);
axis(heat_ax, 'square');
set(heat_ax, 'YDir', 'normal', 'TickDir', 'out', ...
    'Color', bg_color, 'XColor', 'none', 'YColor', 'none', ...
    'Box', 'off', 'FontSize', 7);
hold(heat_ax, 'on');

% ---- Solid black squares around each cluster block ----
% Identify contiguous runs of the same cluster label in leaf order
k_boundaries = [0; find(diff(cid_ordered) ~= 0); n_features];
for kb = 1:length(k_boundaries)-1
    r_start = k_boundaries(kb)   + 1;   % first row/col of cluster
    r_end   = k_boundaries(kb+1);       % last  row/col of cluster
    bw = r_end - r_start + 1;           % block width/height

    rectangle(heat_ax, ...
        'Position',   [r_start-0.5, r_start-0.5, bw, bw], ...
        'EdgeColor',  cluster_box_c, ...
        'FaceColor',  'none', ...
        'LineWidth',  1.6);
end

% ---- Blue box on the diagonal cell for each cluster representative ----
for k = 1:n_clusters
    orig_idx  = selected_indices(k);
    reord_pos = find(leaf_order == orig_idx);   % position in heatmap

    rectangle(heat_ax, ...
        'Position',   [reord_pos-0.5, reord_pos-0.5, 1, 1], ...
        'EdgeColor',  centroid_box_c, ...
        'FaceColor',  'none', ...
        'LineWidth',  2.0);
end

% ---- Feature name labels on right side of heatmap ----
labels_ordered_display = cellfun(@(s) strrep(s, '_', ' '), labels_ordered, 'UniformOutput', false);

for fi = 1:n_features
    label_str = labels_ordered_display{fi};

    if is_selected(fi)
        text(heat_ax, n_features + 0.7, fi, label_str, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
            'FontSize', 6.5, 'FontWeight', 'bold', ...
            'Color', centroid_box_c);          % blue = representative
    else
        text(heat_ax, n_features + 0.7, fi, label_str, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
            'FontSize', 6, 'FontWeight', 'normal', ...
            'Color', ax_text_color);
    end
end

xlim(heat_ax, [0.5 n_features + 0.5]);
ylim(heat_ax, [0.5 n_features + 0.5]);

% ---- Colourbar ----
axes(cbar_ax);
imagesc(cbar_ax, 1, linspace(0, 1, 256), linspace(0, 1, 256)');
colormap(cbar_ax, cmap_heat);
set(cbar_ax, 'YDir', 'normal', 'XTick', [], ...
    'YAxisLocation', 'right', 'TickDir', 'out', ...
    'Color', bg_color, 'YColor', ax_text_color, ...
    'FontSize', 7, 'Box', 'off');
yticks(cbar_ax,      [0 0.25 0.5 0.75 1]);
yticklabels(cbar_ax, {'0','0.25','0.5','0.75','1'});
ylabel(cbar_ax, sprintf('%s r', CORR_METHOD), ...
    'Color', ax_text_color, 'FontSize', 8);

% ---- Super-title ----
annotation(fig, 'textbox', [0 0.96 1 0.04], ...
    'String', 'EEG Feature Collinearity: Hierarchical Clustering and Representative Selection', ...
    'Color', ax_text_color, 'BackgroundColor', 'none', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');

%% ========================================================================
%  SAVE FIGURE
%  ========================================================================

if ~isempty(FIGURE_OUTPUT)
    print(fig, FIGURE_OUTPUT, '-dsvg', '-vector');
    fprintf('\nFigure saved to: %s\n', FIGURE_OUTPUT);
end

%% ========================================================================
%  SAVE RESULTS TO .mat
%  ========================================================================

clustering_results = struct();
clustering_results.selected_features = selected_features;
clustering_results.selected_indices  = selected_indices;
clustering_results.cluster_id        = cluster_id;
clustering_results.cluster_sizes     = cluster_sizes;
clustering_results.feature_cols      = feature_cols;
clustering_results.n_clusters        = n_clusters;
clustering_results.R                 = R;
clustering_results.Z                 = Z;
clustering_results.leaf_order        = leaf_order;
clustering_results.threshold         = CLUSTER_THRESHOLD;

save('Feature_Clustering_Results.mat', 'clustering_results', 'selected_features');
fprintf('Clustering results saved to Feature_Clustering_Results.mat\n');

%% ========================================================================
%  SAVE SELECTED-FEATURES TABLE FOR DOWNSTREAM ANALYSIS
%  ========================================================================
% Combine metadata columns (if present) with the selected feature columns.

% Only keep meta cols that actually exist in the original table
existing_meta = META_COLS(ismember(META_COLS, results_tbl.Properties.VariableNames));

% Build the reduced table: metadata + one representative per cluster
selected_col_names = [existing_meta, selected_features'];
results_selected   = results_tbl(:, selected_col_names);

writetable(results_selected, [OUT_STEM '.csv']);

fprintf('Selected-features table (%d features + %d meta cols, %d rows) saved to:\n', ...
    length(selected_features), length(existing_meta), height(results_selected));
fprintf('  %s.csv\n', OUT_STEM);
