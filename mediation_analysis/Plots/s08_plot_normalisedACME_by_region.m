%% s08 - Per-feature normalised ACME bar chart, grouped by brain region
% Reads: Mediation_JointModel_RelBase/multivariate_mediation_results.csv (from s05)
% Shows ACME_CoV per feature, grouped by region, coloured by direction.
% Markers: * / ** / *** = QB_p uncorrected;  † = Perm_p uncorrected
%
% Junheng Li

clear; clc; close all

%% Settings

INPUT_FILE = 'Data for reproducing final plots/multivariate_mediation_results.csv';
OUTPUT_DIR = 'Results/Plots';

OUTCOMES = {'Image Accuracy'};
CONTRAST = 'TICue_vs_SoundOnly';

COLOR_POS = [0.91 0.39 0.27];
COLOR_NEG = [0.13 0.70 0.67];

SIG_QB_THRESH   = 0.05;
SIG_PERM_THRESH = 0.05;

REGION_ORDER   = {'frontal', 'central', 'left_temporal', 'right_temporal', ...
                  'posterior', 'occipital', 'frontal_central'};
REGION_DISPLAY = {'Frontal', 'Central', 'L Temporal', 'R Temporal', ...
                  'Posterior', 'Occipital', 'Frontal Central'};

FTYPE_ORDER   = {'spindle_power', 'spindle_peak_freq', 'theta_power', ...
                 'sigma_burst_amp', 'so_spindle', 'delta_spindle', ...
                 'so_delta_diff', 'aperiodic'};
FTYPE_DISPLAY = {'Spindle power', 'Spindle peak freq', 'Theta power', ...
                 'Spindle band amplitude', 'SO spindle coupling', ...
                 'Delta spindle coupling', 'SO delta diff', 'Aperiodic slope'};

%% Load

if ~exist(OUTPUT_DIR, 'dir'); mkdir(OUTPUT_DIR); end
T = readtable(INPUT_FILE);
has_perm = ismember('Perm_p', T.Properties.VariableNames);

%% Helper: parse feature name → region + type

[~, len_idx] = sort(cellfun(@length, REGION_ORDER), 'descend');
region_sorted = REGION_ORDER(len_idx);

%% Main loop

for oi = 1:length(OUTCOMES)
    outcome = OUTCOMES{oi};

    mask = strcmp(T.Outcome, outcome) & strcmp(T.Contrast, CONTRAST);
    Tsub = T(mask, :);
    if height(Tsub) == 0
        fprintf('No data for %s | %s\n', outcome, CONTRAST); continue
    end

    n_feat = height(Tsub);
    feat_region = cell(n_feat, 1);
    feat_type   = cell(n_feat, 1);
    for fi = 1:n_feat
        [feat_region{fi}, feat_type{fi}] = parse_feature(Tsub.Feature{fi}, region_sorted);
    end

    % Sort by region order, then feature-type order
    region_rank = NaN(n_feat, 1); type_rank = NaN(n_feat, 1);
    for fi = 1:n_feat
        idx_r = find(strcmp(REGION_ORDER, feat_region{fi}));
        idx_t = find(strcmp(FTYPE_ORDER,   feat_type{fi}));
        if ~isempty(idx_r); region_rank(fi) = idx_r; else; region_rank(fi) = 99; end
        if ~isempty(idx_t); type_rank(fi)   = idx_t; else; type_rank(fi)   = 99; end
    end
    [~, sort_idx] = sortrows([region_rank, type_rank]);

    cov_vals  = Tsub.ACME_CoV(sort_idx);
    qb_p      = Tsub.QB_p(sort_idx);
    if has_perm; perm_p = Tsub.Perm_p(sort_idx);
    else;        perm_p = NaN(n_feat, 1); end
    regions_s = feat_region(sort_idx);
    types_s   = feat_type(sort_idx);

    % Bar colours by sign
    bar_colors = zeros(n_feat, 3);
    for fi = 1:n_feat
        if cov_vals(fi) >= 0; bar_colors(fi,:) = COLOR_POS;
        else;                  bar_colors(fi,:) = COLOR_NEG; end
    end

    % Region group boundaries
    group_starts = []; group_ends = []; group_names = {};
    prev = '';
    for fi = 1:n_feat
        r = regions_s{fi};
        if ~strcmp(r, prev)
            group_starts(end+1) = fi;       %#ok
            group_names{end+1}  = r;        %#ok
            if fi > 1; group_ends(end+1) = fi - 1; end  %#ok
            prev = r;
        end
    end
    group_ends(end+1) = n_feat;

    % Feature type display labels
    type_display_map = containers.Map(FTYPE_ORDER, FTYPE_DISPLAY);
    y_labels = cell(n_feat, 1);
    for fi = 1:n_feat
        if type_display_map.isKey(types_s{fi})
            y_labels{fi} = type_display_map(types_s{fi});
        else
            y_labels{fi} = strrep(types_s{fi}, '_', ' ');
        end
    end

    % Figure
    fig_h = max(380, n_feat * 28 + 120);
    fig = figure('Position', [100 100 560 fig_h], 'Color', 'w');
    ax  = axes('Position', [0.36 0.10 0.54 0.82]);
    hold on

    xline(0, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0);

    for fi = 1:n_feat
        cv = cov_vals(fi);
        if isnan(cv); continue; end

        fc = bar_colors(fi,:);
        barh(fi, cv, 0.62, 'FaceColor', fc, 'EdgeColor', 'none', 'LineWidth', 0.8);

        % CoV value label at bar tip
        x_text = cv + sign(cv) * 0.02;
        ha = 'left'; if cv < 0; ha = 'right'; end
        text(x_text, fi, sprintf('%.2f', cv), ...
            'HorizontalAlignment', ha, 'VerticalAlignment', 'middle', ...
            'FontSize', 7.5, 'Color', [0.2 0.2 0.2]);

        % Significance markers
        mark_str = '';
        if has_perm && ~isnan(perm_p(fi)) && perm_p(fi) < SIG_PERM_THRESH
            if     qb_p(fi) < 0.001; mark_str = '***';
            elseif qb_p(fi) < 0.01;  mark_str = '**';
            else;                     mark_str = '*';
            end
        end
        if ~isempty(mark_str)
            x_mark = cv + sign(cv) * 0.10;
            ha_m = 'left'; if cv < 0; ha_m = 'right'; end
            text(x_mark, fi, mark_str, ...
                'HorizontalAlignment', ha_m, 'VerticalAlignment', 'middle', ...
                'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.1 0.1 0.1]);
        end
    end

    % Region separators + labels
    for gi = 1:length(group_starts)
        gs = group_starts(gi); ge = group_ends(gi);
        if gi > 1
            plot(xlim, [gs-0.5, gs-0.5], '-', 'Color', [0.80 0.80 0.80], 'LineWidth', 0.8);
        end
        y_mid = (gs + ge) / 2;
        idx_r = find(strcmp(REGION_ORDER, group_names{gi}));
        if ~isempty(idx_r); rlab = REGION_DISPLAY{idx_r};
        else;                rlab = strrep(group_names{gi}, '_', ' '); end
        text(ax.XLim(1) - 0.05, y_mid, rlab, ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
            'FontSize', 9, 'FontWeight', 'bold', 'FontAngle', 'italic', 'Units', 'data');
    end

    yticks(1:n_feat); yticklabels(y_labels);
    set(ax, 'YDir', 'reverse', 'FontSize', 8, 'Box', 'off', 'TickDir', 'out');
    ylim([0.5, n_feat + 0.5]);

    xlabel('ACME ', 'FontSize', 11);
    title(sprintf('%s    %s ', outcome, strrep(CONTRAST, '_', ' ')), ...
        'FontSize', 11, 'FontWeight', 'normal');

    finite_cv = cov_vals(isfinite(cov_vals));
    if ~isempty(finite_cv)
        x_pad = range(finite_cv) * 0.25 + 0.1;
        xlim([min(finite_cv) - x_pad, max(finite_cv) + x_pad]);
    end

    hold off

    % Legend
    ax_leg = axes('Position', [0.02 0.04 0.30 0.18], 'Visible', 'off');
    hold(ax_leg, 'on')
    h1 = patch(ax_leg, NaN, NaN, COLOR_POS, 'EdgeColor', COLOR_POS*0.6);
    h2 = patch(ax_leg, NaN, NaN, COLOR_NEG, 'EdgeColor', COLOR_NEG*0.6);
    legend(ax_leg, [h1,h2], {'Positive CoV','Negative CoV'}, ...
        'Location', 'northwest', 'FontSize', 8, 'Box', 'off');
    leg_note = '* QB p<0.05';
    text(0, -0.35, leg_note, 'FontSize', 7.5, 'Color', [0.25 0.25 0.25], ...
        'Units', 'normalized', 'Parent', ax_leg);
    hold(ax_leg, 'off')

    fname = fullfile(OUTPUT_DIR, sprintf('CoV_grouped_%s_ensemble.svg', ...
        lower(strrep(outcome, ' ', '_'))));
    saveas(fig, fname, 'svg');
    fprintf('Saved: %s\n', fname);
end

fprintf('Done.\n');

%% ── Local helper ──────────────────────────────────────────────────────────

function [reg, ftype] = parse_feature(fname, region_sorted)
    reg = ''; ftype = fname;
    for ri = 1:length(region_sorted)
        prefix = [region_sorted{ri}, '_'];
        if startsWith(fname, prefix)
            reg   = region_sorted{ri};
            ftype = fname(length(prefix)+1:end);
            return
        end
    end
end
