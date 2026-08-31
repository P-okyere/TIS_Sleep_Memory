%% s07 - Normalised ACME (ACME mean/SD) violin plots, joint model
% Reads: Mediation_JointModel_RelBase/multivariate_mediation_results.csv (from s05)
% Produces one figure: 3 violins (one per contrast), CoV distribution
%
% Junheng Li

clear; clc; close all

%% ── Settings ──────────────────────────────────────────────────────────────

INPUT_FILE = 'Data for reproducing final plots/multivariate_mediation_results.csv';
OUTPUT_DIR = 'Results/Plots';

OUTCOME_LABELS   = {'Image Accuracy'};
CONTRAST_LABELS  = {'TICue_vs_SoundOnly', 'TIB4Cue_vs_SoundOnly', 'TICue_vs_TIB4Cue'};
CONTRAST_DISPLAY = {'TICue vs SO', 'TIB4Cue vs SO', 'TICue vs TIB4Cue'};

VIOLIN_COLORS = [0.40 0.76 0.65;   % teal
                 0.99 0.55 0.38;   % coral
                 0.55 0.63 0.80];  % steel blue

%% ── Load ──────────────────────────────────────────────────────────────────

if ~exist(OUTPUT_DIR, 'dir'); mkdir(OUTPUT_DIR); end
T = readtable(INPUT_FILE);

%% ── Plot ──────────────────────────────────────────────────────────────────

for oi = 1:length(OUTCOME_LABELS)
    outcome = OUTCOME_LABELS{oi};

    % Collect CoV per contrast
    cov_cells = cell(1, length(CONTRAST_LABELS));
    max_n = 0;
    for ci = 1:length(CONTRAST_LABELS)
        mask = strcmp(T.Outcome, outcome) & strcmp(T.Contrast, CONTRAST_LABELS{ci});
        vals = T.ACME_CoV(mask);
        vals = vals(~isnan(vals));
        cov_cells{ci} = vals;
        max_n = max(max_n, length(vals));
    end

    % NaN-padded matrix
    data_mat = NaN(max_n, length(CONTRAST_LABELS));
    for ci = 1:length(CONTRAST_LABELS)
        n = length(cov_cells{ci});
        data_mat(1:n, ci) = cov_cells{ci};
    end

    fig = figure('Position', [100 100 600 450], 'Color', 'w');
    hold on

    violins = violinplot(data_mat, CONTRAST_DISPLAY, ...
        'ViolinColor',  {VIOLIN_COLORS}, ...
        'ViolinAlpha',  {0.35}, ...
        'ShowMean',     true, ...
        'ShowMedian',   true, ...
        'ShowBox',      true, ...
        'ShowWhiskers', true, ...
        'ShowData',     true, ...
        'MarkerSize',   20, ...
        'Width',        0.35);

    yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);

    % One-sample t-test annotations (two-sided, then mark if significant)
    yl = ylim;
    y_top = yl(2) * 0.93;
    for ci = 1:length(CONTRAST_LABELS)
        vals = cov_cells{ci};
        if length(vals) < 3; continue; end
        [~, p_test] = ttest(vals, 0);
        if p_test < 0.001
            sig_str = '***';
        elseif p_test < 0.01
            sig_str = '**';
        elseif p_test < 0.05
            sig_str = '*';
        else
            continue
        end
        text(ci, y_top, sig_str, ...
            'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
    end

    ylabel('CoV  (ACME mean / SD)', 'FontSize', 12);
    title(sprintf('ACME CoV across features:  %s  [joint model]', outcome), 'FontSize', 13);
    set(gca, 'FontSize', 11, 'Box', 'off', 'TickDir', 'out', ...
        'TickLength', 2*get(gca,'TickLength'), 'LineWidth', 2);
    xlim([0.3  length(CONTRAST_LABELS) + 0.7]);
    ylim([-0.6, 1]);

    hold off

    fname = fullfile(OUTPUT_DIR, sprintf('cov_violin_%s_ensemble.svg', ...
        lower(strrep(outcome, ' ', '_'))));
    saveas(fig, fname, 'svg');
    fprintf('Saved: %s\n', fname);
end

fprintf('Done.\n');


%% ── Summary stats for manuscript reporting ────────────────────────────────

fprintf('\n========== CoV Summary Stats ==========\n');
fprintf('%-22s  %8s  %8s  %8s\n', 'Contrast', 'Mean', 'SD', 't-test p');
fprintf('%s\n', repmat('-', 1, 52));

for oi = 1:length(OUTCOME_LABELS)
    outcome = OUTCOME_LABELS{oi};
    fprintf('Outcome: %s\n', outcome);
    for ci = 1:length(CONTRAST_LABELS)
        mask = strcmp(T.Outcome, outcome) & strcmp(T.Contrast, CONTRAST_LABELS{ci});
        vals = T.ACME_CoV(mask);
        vals = vals(~isnan(vals));
        if length(vals) < 3; continue; end
        [~, p_val, ~, stats] = ttest(vals, 0);
        fprintf('  %-20s  %+8.4f  %8.4f  %8.4f  (t=%.3f, df=%d, n=%d)\n', ...
            CONTRAST_DISPLAY{ci}, mean(vals), std(vals), p_val, ...
            stats.tstat, stats.df, length(vals));
    end
end
fprintf('========================================\n');