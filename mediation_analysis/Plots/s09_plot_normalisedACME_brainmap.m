%% s09 - Brain bubble map: mean normalised ACME per region, joint model
% Fixed-size bubbles, colour = mean ACME_mean across features in that region.
% Diverging colormap: blue (negative) — white (0) — red (positive).
% Dot in bubble centre if any feature in that region has Perm_p_FDR < 0.05.
%
% Junheng Li

clear; clc; close all

%% ── Settings ──────────────────────────────────────────────────────────────

INPUT_FILE = 'Data for reproducing final plots/multivariate_mediation_results.csv';
OUTPUT_DIR = 'Results/Plots';

CONTRAST = 'TICue_vs_SoundOnly';
OUTCOMES = {'Image Accuracy'};

BUBBLE_RADIUS = 0.14;
BUBBLE_ALPHA  = 0.90;
DOT_RADIUS    = 0.028;   % inner dot if any feature in region is perm-sig

SIG_PERM_THRESH = 0.05;  % Perm_p_FDR threshold

COLOR_HEAD = [0.93 0.93 0.93];
COLOR_EDGE = [0.50 0.50 0.50];

REGION_ORDER   = {'frontal','central','left_temporal','right_temporal', ...
                  'posterior','occipital','frontal_central'};
REGION_DISPLAY = {'Frontal','Central','L Temporal','R Temporal', ...
                  'Posterior','Occipital','FC'};

REGION_XY = [
     0.00,  0.58;   % frontal
     0.00,  0.15;   % central
    -0.62,  0.05;   % left_temporal
     0.62,  0.05;   % right_temporal
     0.00, -0.32;   % posterior
     0.00, -0.62;   % occipital
     0.22,  0.37;   % frontal_central (nudged right)
];

%% ── Diverging colormap (blue — white — red) ───────────────────────────────

n_colors = 256;
half = n_colors / 2;
blue_to_white = [linspace(0.17, 1, half)', linspace(0.35, 1, half)', linspace(0.80, 1, half)'];
white_to_red  = [linspace(1, 0.80, half)', linspace(1, 0.15, half)', linspace(1, 0.15, half)'];
CMAP = [blue_to_white; white_to_red];

%% ── Helper: parse feature name into region ────────────────────────────────

[~, len_idx] = sort(cellfun(@length, REGION_ORDER), 'descend');
region_by_len = REGION_ORDER(len_idx);

%% ── Load & compute mean ACME per region ──────────────────────────────────

if ~exist(OUTPUT_DIR, 'dir'); mkdir(OUTPUT_DIR); end
T = readtable(INPUT_FILE);

has_perm = ismember('Perm_p_FDR', T.Properties.VariableNames);

n_regions  = length(REGION_ORDER);
n_outcomes = length(OUTCOMES);
mean_acme  = NaN(n_regions, n_outcomes);
any_perm_sig = false(n_regions, n_outcomes);

for oi = 1:n_outcomes
    mask = strcmp(T.Outcome, OUTCOMES{oi}) & strcmp(T.Contrast, CONTRAST);
    Tsub = T(mask, :);

    for ri = 1:n_regions
        feat_mask = false(height(Tsub), 1);
        for fi = 1:height(Tsub)
            feat_mask(fi) = strcmp(get_region(Tsub.Feature{fi}, region_by_len), ...
                                   REGION_ORDER{ri});
        end
        vals = Tsub.ACME_CoV(feat_mask & isfinite(Tsub.ACME_mean));
        if ~isempty(vals)
            mean_acme(ri, oi) = mean(vals);
        end
        if has_perm
            pv = Tsub.Perm_p(feat_mask);
            if any(~isnan(pv) & pv < SIG_PERM_THRESH)
                any_perm_sig(ri, oi) = true;
            end
        end
    end
end

fprintf('Mean ACME per region (TICue vs SO):\n');
for ri = 1:n_regions
    fprintf('  %-22s  %+.3f%s\n', REGION_ORDER{ri}, mean_acme(ri,1), ...
        if_else(any_perm_sig(ri,1), '  [perm sig]', ''));
end

%% ── Colour mapping ────────────────────────────────────────────────────────

clim_abs = max(abs(mean_acme(:)), [], 'omitnan');
clim_abs = ceil(clim_abs * 10) / 10 + 0.02;

cov_to_color = @(cv) interp1(linspace(-clim_abs, clim_abs, n_colors), ...
                              CMAP, max(-clim_abs, min(clim_abs, cv)), 'linear');

%% ── Draw ──────────────────────────────────────────────────────────────────

fig = figure('Position', [80 80 520 500], 'Color', 'w');
theta_b = linspace(0, 2*pi, 120);

for oi = 1:n_outcomes
    ax = subplot(1, n_outcomes, oi);
    hold on; axis equal; axis off;

    % Head outline
    theta_h = linspace(0, 2*pi, 300);
    head_rx = 0.82; head_ry = 0.88;
    patch(head_rx*cos(theta_h), head_ry*sin(theta_h), ...
        COLOR_HEAD, 'EdgeColor', COLOR_EDGE, 'LineWidth', 1.8);

    % Nose
    patch([-0.09, 0.00, 0.09], [0.87, 1.00, 0.87], ...
        COLOR_HEAD, 'EdgeColor', COLOR_EDGE, 'LineWidth', 1.8);

    % Ears
    theta_e = linspace(pi*0.55, pi*1.45, 60);
    ear_rx = 0.07; ear_ry = 0.13; ear_cy = 0.02;
    patch(-head_rx + ear_rx*cos(theta_e), ear_cy + ear_ry*sin(theta_e), ...
        COLOR_HEAD, 'EdgeColor', COLOR_EDGE, 'LineWidth', 1.8);
    theta_e = linspace(-pi*0.45, pi*0.45, 60);
    patch( head_rx + ear_rx*cos(theta_e), ear_cy + ear_ry*sin(theta_e), ...
        COLOR_HEAD, 'EdgeColor', COLOR_EDGE, 'LineWidth', 1.8);

    % Bubbles
    for ri = 1:n_regions
        cv = mean_acme(ri, oi);
        if isnan(cv); continue; end

        xc = REGION_XY(ri, 1);
        yc = REGION_XY(ri, 2);
        fc = cov_to_color(cv);

        % Main bubble
        patch(xc + BUBBLE_RADIUS*cos(theta_b), yc + BUBBLE_RADIUS*sin(theta_b), ...
            fc, 'EdgeColor', fc*0.55, 'LineWidth', 1.3, 'FaceAlpha', BUBBLE_ALPHA);

        % % Inner dot if any feature in region is perm-significant
        % if any_perm_sig(ri, oi)
        %     dot_col = if_else(cv >= 0, [0.6 0.1 0.1], [0.1 0.1 0.6]);
        %     patch(xc + DOT_RADIUS*cos(theta_b), yc + DOT_RADIUS*sin(theta_b), ...
        %         dot_col, 'EdgeColor', dot_col, 'LineWidth', 0.5, 'FaceAlpha', 1.0);
        % end

    end

    title(OUTCOMES{oi}, 'FontSize', 11, 'FontWeight', 'bold');
    xlim([-1.15  1.25]); ylim([-1.10  1.15]);
    hold off
end

% Colorbar
ax_cb = axes('Position', [0.88 0.15 0.020 0.70]);
imagesc(ax_cb, 1, linspace(-clim_abs, clim_abs, n_colors), reshape(CMAP, n_colors, 1, 3));
set(ax_cb, 'YDir', 'normal', 'XTick', [], 'YAxisLocation', 'right', ...
    'TickDir', 'out', 'FontSize', 8, 'Box', 'off');
n_ticks   = 5;
tick_vals = linspace(-clim_abs, clim_abs, n_ticks);
yticks(ax_cb, tick_vals);
yticklabels(ax_cb, arrayfun(@(v) sprintf('%+.3f', v), tick_vals, 'UniformOutput', false));
ylabel(ax_cb, 'Mean ACME', 'FontSize', 9);

% Super-title
annotation(fig, 'textbox', [0 0.93 1 0.06], ...
    'String', sprintf('Regional ACME (joint model) — %s', strrep(CONTRAST,'_',' ')), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', 11, 'FontWeight', 'normal', ...
    'EdgeColor', 'none', 'BackgroundColor', 'none');

% % Legend note
% if has_perm
%     annotation(fig, 'textbox', [0.01 0.01 0.70 0.06], ...
%         'String', sprintf('Dot = any feature in region Perm_p_{FDR} < %.2f', SIG_PERM_THRESH), ...
%         'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
%         'FontSize', 7.5, 'EdgeColor', 'none', 'BackgroundColor', 'none');
% end

%% ── Save ──────────────────────────────────────────────────────────────────

fname = fullfile(OUTPUT_DIR, 'BrainBubble_ACME_Ensemble.svg');
saveas(fig, fname, 'svg');
fprintf('\nSaved: %s\n', fname);

%% ── Local helper ──────────────────────────────────────────────────────────

function out = if_else(cond, a, b)
    if cond; out = a; else; out = b; end
end


function reg = get_region(fname, region_by_len)
    reg = '';
    for ri = 1:length(region_by_len)
        if startsWith(fname, [region_by_len{ri}, '_'])
            reg = region_by_len{ri};
            return
        end
    end
end

