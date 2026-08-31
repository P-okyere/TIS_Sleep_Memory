%% s06 - Path-a t value heatmap (standalone path-a analysis)
% Reads : Data for reproducing final plots/path_a_results.csv (from s04)
% Writes: Plots/Output/path_a_heatmap.svg and .png
% Shows path-a t values (treatment to EEG feature) as a features x contrasts
% heatmap. Features grouped by region (frontal to occipital, cross channel last).
%
% Significance is judged on the FDR-adjusted p value (path_a_p_adj: BH within
% contrast), which is the p value reported in the manuscript. Cells with
% adjusted p < 0.05 are boxed and carry asterisks (* < .05, ** < .01,
% *** < .001). Uncorrected p values are read for reference but are NOT used
% for marking. Three cells survive FDR in the released results table.
%
% Outcome variable does NOT enter path-a, so no outcome dimension here.
%
% Junheng Li

clear; clc; close all

%% Settings

% Paths are resolved against the repository root (the parent of this script's
% folder) so the script runs correctly from any working directory.
INPUT_FILE = 'Data for reproducing final plots/path_a_results.csv';
OUTPUT_DIR = 'Results/Plots';

CONTRAST_LABELS  = {'TICue_vs_SoundOnly', 'TIB4Cue_vs_SoundOnly', 'TICue_vs_TIB4Cue'};
CONTRAST_DISPLAY = {'TI Cue vs SO', 'TI Before Cue vs SO', 'TI Cue vs TI Before Cue'};

SIG_THRESH_FDR = 0.05;   % FDR-adjusted threshold, drives all markers
SIG_THRESH     = 0.05;   % fallback only, if no path_a_p_adj column is present

% Region order: frontal to occipital, cross channel last
REGION_ORDER   = {'frontal', 'central', 'left_temporal', 'right_temporal', ...
                  'posterior', 'occipital', 'frontal_central'};
REGION_DISPLAY = {'Frontal', 'Central', 'L Temporal', 'R Temporal', ...
                  'Posterior', 'Occipital', 'Frontal Central'};

% Feature type order and display labels
FTYPE_ORDER   = {'spindle_power', 'spindle_peak_freq', 'theta_power', ...
                 'sigma_burst_amp', 'so_spindle', 'delta_spindle', ...
                 'so_delta_diff', 'aperiodic'};
FTYPE_DISPLAY = {'Spindle power', 'Spindle peak freq', 'Theta power', ...
                 'Spindle band amplitude', 'SO spindle coupling', ...
                 'Delta spindle coupling', ...
                 'SO delta diff', 'Aperiodic slope'};

%% Load data

T = readtable(INPUT_FILE);
all_features = unique(T.Feature, 'stable');

% Check if FDR column exists
has_fdr = ismember('path_a_p_adj', T.Properties.VariableNames);
if has_fdr
    fprintf('FDR-adjusted p-values found; markers use path_a_p_adj.\n');
else
    warning(['No path_a_p_adj column found. Falling back to uncorrected ' ...
             'p-values, which does NOT match the manuscript.']);
end

%% Parse each feature into region + type

[~, len_idx] = sort(cellfun(@length, REGION_ORDER), 'descend');
region_order_sorted = REGION_ORDER(len_idx);

feat_region = cell(size(all_features));
feat_type   = cell(size(all_features));

for fi = 1:length(all_features)
    fname = all_features{fi};
    matched = false;
    for ri = 1:length(region_order_sorted)
        prefix = [region_order_sorted{ri}, '_'];
        if startsWith(fname, prefix)
            feat_region{fi} = region_order_sorted{ri};
            feat_type{fi}   = fname(length(prefix)+1 : end);
            matched = true;
            break
        end
    end
    if ~matched
        feat_region{fi} = 'unknown';
        feat_type{fi}   = fname;
    end
end

%% Build canonical sorted order (region then feature type)

region_rank = NaN(length(all_features), 1);
type_rank   = NaN(length(all_features), 1);

for fi = 1:length(all_features)
    idx_r = find(strcmp(REGION_ORDER, feat_region{fi}));
    if ~isempty(idx_r); region_rank(fi) = idx_r;
    else;               region_rank(fi) = length(REGION_ORDER) + 1; end

    idx_t = find(strcmp(FTYPE_ORDER, feat_type{fi}));
    if ~isempty(idx_t); type_rank(fi) = idx_t;
    else;               type_rank(fi) = length(FTYPE_ORDER) + 1; end
end

[~, canon_idx] = sortrows([region_rank, type_rank], [1 2]);
features_sorted = all_features(canon_idx);
regions_sorted  = feat_region(canon_idx);
types_sorted    = feat_type(canon_idx);

n_feat = length(features_sorted);
n_con  = length(CONTRAST_LABELS);

%% Build display labels from feature type names

type_display_map = containers.Map(FTYPE_ORDER, FTYPE_DISPLAY);

feat_display = cell(n_feat, 1);
for fi = 1:n_feat
    if type_display_map.isKey(types_sorted{fi})
        feat_display{fi} = type_display_map(types_sorted{fi});
    else
        feat_display{fi} = strrep(types_sorted{fi}, '_', ' ');
    end
end

%% Region group boundaries

region_breaks    = [];
region_label_pos = [];
region_label_str = {};

prev_region = '';
group_start = 1;
for fi = 1:n_feat
    if ~strcmp(regions_sorted{fi}, prev_region)
        if fi > 1
            region_breaks(end+1) = fi - 0.5;
            region_label_pos(end+1) = (group_start + fi - 1) / 2;
            idx_r = find(strcmp(REGION_ORDER, prev_region));
            if ~isempty(idx_r); region_label_str{end+1} = REGION_DISPLAY{idx_r};
            else;               region_label_str{end+1} = prev_region; end
            group_start = fi;
        end
        prev_region = regions_sorted{fi};
    end
end
region_label_pos(end+1) = (group_start + n_feat) / 2;
idx_r = find(strcmp(REGION_ORDER, prev_region));
if ~isempty(idx_r); region_label_str{end+1} = REGION_DISPLAY{idx_r};
else;               region_label_str{end+1} = prev_region; end

%% Diverging colormap: blue white red

n_colors = 256;
half = n_colors / 2;
blue_to_white = [linspace(0.2, 1, half)', linspace(0.4, 1, half)', linspace(0.8, 1, half)'];
white_to_red  = [linspace(1, 0.8, half)', linspace(1, 0.2, half)', linspace(1, 0.2, half)'];
cmap = [blue_to_white; white_to_red];

%% Build t-value, p-value, and FDR-p matrices

t_mat   = NaN(n_feat, n_con);
p_mat   = NaN(n_feat, n_con);
p_mat_a = NaN(n_feat, n_con);   % FDR-adjusted

for fi = 1:n_feat
    for ci = 1:n_con
        row = T(strcmp(T.Feature, features_sorted{fi}) & ...
                strcmp(T.Contrast, CONTRAST_LABELS{ci}), :);
        if ~isempty(row)
            t_mat(fi, ci) = row.path_a_t;
            p_mat(fi, ci) = row.path_a_p;
            if has_fdr
                p_mat_a(fi, ci) = row.path_a_p_adj;
            end
        end
    end
end

%% Generate heatmap

fprintf('Generating path-a heatmap ...\n');

fig_h = max(450, n_feat * 22 + 150);
fig = figure('Position', [50 50 860 fig_h], 'Color', 'w');

% Left margin holds two label columns: italic region names, then feature names
ax_main = axes('Position', [0.40 0.12 0.46 0.78]);
imagesc(t_mat, 'AlphaData', ~isnan(t_mat));
hold on

colormap(ax_main, cmap);

c_lim = max(abs(t_mat(:)));
if isnan(c_lim) || c_lim == 0; c_lim = 1; end
clim([-c_lim, c_lim]);

cb = colorbar;
cb.Label.String = 'Path-a t value';
cb.Label.FontSize = 11;

set(ax_main, 'Color', [0.93 0.93 0.93]);

% Significance markers. The manuscript reports FDR-adjusted p values, so the
% adjusted p matrix drives both the asterisks and the box.
if has_fdr
    p_sig      = p_mat_a;
    sig_thresh = SIG_THRESH_FDR;
else
    p_sig      = p_mat;
    sig_thresh = SIG_THRESH;
end

n_sig = 0;
for fi = 1:n_feat
    for ci = 1:n_con
        if isnan(p_sig(fi, ci)); continue; end
        if p_sig(fi, ci) >= sig_thresh; continue; end

        if     p_sig(fi, ci) < 0.0005; marker = '***';
        elseif p_sig(fi, ci) < 0.005;  marker = '**';
        else;                          marker = '*'; end

        if abs(t_mat(fi, ci)) > c_lim * 0.6
            txt_color = 'w';
        else
            txt_color = 'k';
        end

        text(ci, fi, marker, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment',   'middle', ...
            'FontSize', 10, 'FontWeight', 'bold', 'Color', txt_color);

        % Box every cell that passes the reported threshold
        rectangle('Position', [ci-0.5, fi-0.5, 1, 1], ...
            'EdgeColor', [0 0 0], 'LineWidth', 1.8, 'LineStyle', '-');

        n_sig = n_sig + 1;
        fprintf('  sig: %-24s %-32s t = %+.3f, p_adj = %.3g\n', ...
            CONTRAST_LABELS{ci}, features_sorted{fi}, ...
            t_mat(fi, ci), p_sig(fi, ci));
    end
end
fprintf('Significant cells marked: %d\n', n_sig);

% Region separators
for bi = 1:length(region_breaks)
    plot([0.5 n_con+0.5], [region_breaks(bi) region_breaks(bi)], ...
        'k-', 'LineWidth', 1.5);
end

% Light grid
for fi = 0.5:1:(n_feat+0.5)
    plot([0.5 n_con+0.5], [fi fi], '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.3);
end
for ci = 0.5:1:(n_con+0.5)
    plot([ci ci], [0.5 n_feat+0.5], '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.3);
end

% Axes
set(ax_main, 'XTick', 1:n_con, 'XTickLabel', CONTRAST_DISPLAY, ...
    'XTickLabelRotation', 25, ...
    'YTick', 1:n_feat, 'YTickLabel', feat_display, ...
    'FontSize', 9, 'TickLength', [0 0], 'YDir', 'reverse');

if has_fdr
    title_str = 'Path-a: treatment effect on EEG features (boxed = FDR p_{adj} < 0.05)';
else
    title_str = 'Path-a: treatment effect on EEG features (boxed = uncorrected p < 0.05)';
end
title(title_str, 'FontSize', 13);

hold off

% Region labels (overlay on transparent twin axes)
ax_region = axes('Position', ax_main.Position, 'Color', 'none');
set(ax_region, 'XLim', [0 1], 'YLim', [0.5 n_feat+0.5], ...
    'YDir', 'reverse', 'Visible', 'off');
hold on
for ri = 1:length(region_label_str)
    text(-0.42, region_label_pos(ri), region_label_str{ri}, ...
        'HorizontalAlignment', 'right', 'FontSize', 10, ...
        'FontWeight', 'bold', 'FontAngle', 'italic', 'Parent', ax_region);
end
hold off

%% Save

if ~exist(OUTPUT_DIR, 'dir'); mkdir(OUTPUT_DIR); end

fname_svg = fullfile(OUTPUT_DIR, 'path_a_heatmap.svg');
fname_png = fullfile(OUTPUT_DIR, 'path_a_heatmap.png');
saveas(fig, fname_svg, 'svg');
saveas(fig, fname_png, 'png');
fprintf('Saved: %s\n', fname_svg);
fprintf('Saved: %s\n', fname_png);

fprintf('Done.\n');
