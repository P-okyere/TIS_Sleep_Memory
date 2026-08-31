# ============================================================================
# s04 - PATH-A ANALYSIS: Treatment -> EEG Feature
#
# Tests whether each EEG feature differs between conditions, independent of
# any behavioural outcome. This is the X -> M arm of the mediation model,
# fitted as a standalone univariate mixed model per feature.
#
# Model: feature_z ~ treat + (1 | subject)
# Fitted with lmerTest::lmer(); Satterthwaite df for the treat coefficient.
#
# Filtering: only mediator (feature) non-NaN trials. Outcome is NOT used,
# so n per cell will typically exceed that of the full mediation pipeline.
# MIN_TRIALS_PER_CELL is applied on mediator-valid trial counts.
#
# Outputs:
#   path_a_results.csv : input to s06_plot_path_a_heatmap.m
#   path_a_summary.txt : paper-ready stats summary
#
# Author: Junheng Li
# ============================================================================

library(lme4)
library(lmerTest)
library(dplyr)

# ============================================================================
# USER SETTINGS
# ============================================================================

DATA_FILE  <- "Results_SelectedFeatures_RelBase.csv"
OUTPUT_DIR <- "PathA_RelBase"

# Metadata columns (must match MATLAB META_COLS)
META_COLS <- c("global_trial_idx", "trial_idx", "subject", "condition",
               "word_acc_post", "imag_acc_post", "old_new", "word",
               "block_post", "Manscore", "NewManscore", "NewManscore1")

# Pairwise contrasts
CONTRASTS <- list(
  list(treat = "TICue",   control = "SoundOnly", label = "TICue_vs_SoundOnly"),
  list(treat = "TIB4Cue", control = "SoundOnly", label = "TIB4Cue_vs_SoundOnly"),
  list(treat = "TICue",   control = "TIB4Cue",   label = "TICue_vs_TIB4Cue")
)

MIN_TRIALS_PER_CELL <- 5
set.seed(2026)

# ============================================================================
# SETUP
# ============================================================================

dir.create(OUTPUT_DIR, showWarnings = FALSE)

cat("════════════════════════════════════════════════════════════════════════\n")
cat("   PATH-A ANALYSIS: Treatment -> EEG Feature\n")
cat("════════════════════════════════════════════════════════════════════════\n\n")
cat("Outcome variable is NOT involved. Filtering on mediator non-NaN only.\n\n")

# ============================================================================
# LOAD DATA & AUTO-DETECT FEATURES
# ============================================================================

data_raw <- read.csv(DATA_FILE)

existing_meta <- intersect(META_COLS, names(data_raw))
all_numeric   <- names(data_raw)[sapply(data_raw, is.numeric)]
feature_cols  <- setdiff(all_numeric, existing_meta)

cat(sprintf("Data: %d trials, %d subjects, %d features\n\n",
            nrow(data_raw), length(unique(data_raw$subject)), length(feature_cols)))

data_raw$subject <- factor(data_raw$subject)

# ============================================================================
# HELPER: fit single path-a model
# ============================================================================

fit_path_a <- function(data, feat_name, treat_val, control_val, min_trials_per_cell) {
  
  d <- data %>%
    filter(condition %in% c(treat_val, control_val)) %>%
    mutate(
      treat    = as.integer(condition == treat_val),
      mediator = as.numeric(scale(!!sym(feat_name)))
    ) %>%
    filter(!is.na(mediator))
  
  # Trial-count filter on mediator-valid trials only
  trial_counts <- d %>% count(subject, treat, name = "n_trials")
  valid_subjects <- trial_counts %>%
    group_by(subject) %>%
    filter(n() == 2, all(n_trials >= min_trials_per_cell)) %>%
    pull(subject) %>%
    unique()
  
  n_excluded <- length(setdiff(unique(d$subject), valid_subjects))
  d <- d %>% filter(subject %in% valid_subjects) %>% droplevels()
  
  if (nrow(d) < 30 || length(unique(d$treat)) < 2 || length(unique(d$subject)) < 3) {
    return(NULL)
  }
  
  mod <- tryCatch(
    lmerTest::lmer(mediator ~ treat + (1|subject),
                   data = d,
                   control = lmerControl(optimizer = "bobyqa",
                                         optCtrl = list(maxfun = 100000))),
    error = function(e) NULL
  )
  if (is.null(mod)) return(NULL)
  
  sing <- isSingular(mod)
  
  cs <- summary(mod)$coefficients
  beta   <- cs["treat", "Estimate"]
  se     <- cs["treat", "Std. Error"]
  t_val  <- cs["treat", "t value"]
  df_val <- cs["treat", "df"]
  p_val  <- cs["treat", "Pr(>|t|)"]
  
  # 95% CI from Satterthwaite df
  t_crit <- qt(0.975, df_val)
  ci_lo  <- beta - t_crit * se
  ci_hi  <- beta + t_crit * se
  
  # Cohen's d analogue: since mediator is z-scored globally, beta is approximately
  # the standardised mean difference between conditions on the feature scale.
  cohens_d <- beta
  
  data.frame(
    N            = nrow(d),
    N_subjects   = length(unique(d$subject)),
    N_excluded   = n_excluded,
    singular     = sing,
    path_a_coef  = beta,
    path_a_se    = se,
    path_a_t     = t_val,
    path_a_df    = df_val,
    path_a_p     = p_val,
    path_a_CI_lo = ci_lo,
    path_a_CI_hi = ci_hi,
    path_a_d     = cohens_d,
    stringsAsFactors = FALSE
  )
}

# ============================================================================
# MAIN LOOP
# ============================================================================

all_results  <- data.frame()
model_count  <- 0
total_models <- length(CONTRASTS) * length(feature_cols)

t_start <- Sys.time()
cat(sprintf("Fitting %d models (%d features x %d contrasts) ...\n\n",
            total_models, length(feature_cols), length(CONTRASTS)))

for (con in CONTRASTS) {
  cat(sprintf("Contrast: %s\n", con$label))
  
  for (feat in feature_cols) {
    model_count <- model_count + 1
    cat(sprintf("  [%3d/%3d] %-45s ", model_count, total_models, feat))
    
    res <- tryCatch(
      fit_path_a(data_raw, feat, con$treat, con$control, MIN_TRIALS_PER_CELL),
      error = function(e) {
        cat(sprintf("FAIL (%s)\n", substr(e$message, 1, 40)))
        NULL
      }
    )
    
    if (is.null(res)) { cat("SKIP\n"); next }
    
    res$Feature  <- feat
    res$Contrast <- con$label
    all_results  <- rbind(all_results, res)
    
    if (res$singular) cat("[sing] ")
    sig_flag <- if (res$path_a_p < 0.001) "***" else if (res$path_a_p < 0.01) "**" else if (res$path_a_p < 0.05) "*" else ""
    cat(sprintf("t(%.0f) = %+6.2f, p = %.4f %s\n",
                res$path_a_df, res$path_a_t, res$path_a_p, sig_flag))
  }
  cat("\n")
}

t_total <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
cat(sprintf("Complete: %d models in %.1f s\n\n", nrow(all_results), t_total))

# ============================================================================
# FDR CORRECTION (BH, within contrast)
# ============================================================================

all_results <- all_results %>%
  group_by(Contrast) %>%
  mutate(path_a_p_adj = p.adjust(path_a_p, method = "fdr")) %>%
  ungroup()

# Singularity summary
n_sing <- sum(all_results$singular, na.rm = TRUE)
cat(sprintf("Singularity: %d / %d models flagged isSingular()\n\n",
            n_sing, nrow(all_results)))

# ============================================================================
# SAVE STATS TABLE (compatible with PathA_Heatmap_PathAOnly.m)
# ============================================================================

col_order <- c("Contrast", "Feature", "N", "N_subjects", "N_excluded", "singular",
               "path_a_coef", "path_a_se", "path_a_t", "path_a_df",
               "path_a_p", "path_a_p_adj",
               "path_a_CI_lo", "path_a_CI_hi", "path_a_d")
all_results <- all_results %>% dplyr::select(all_of(col_order))

out_csv <- file.path(OUTPUT_DIR, "path_a_results.csv")
write.csv(all_results, out_csv, row.names = FALSE)
cat(sprintf("Saved: %s\n", out_csv))

# ============================================================================
# PAPER-READY SUMMARY
# ============================================================================

cat("\n════════════════════════════════════════════════════════════════════════\n")
cat("PATH-A SUMMARY (for paper)\n")
cat("════════════════════════════════════════════════════════════════════════\n")

summary_lines <- c(
  sprintf("Path-a analysis: treatment -> EEG feature"),
  sprintf("Model: feature_z ~ treat + (1|subject), fitted via lmerTest::lmer"),
  sprintf("Satterthwaite df. FDR (BH) within contrast."),
  sprintf("MIN_TRIALS_PER_CELL = %d", MIN_TRIALS_PER_CELL),
  ""
)

for (con in CONTRASTS) {
  d_sub <- all_results %>% filter(Contrast == con$label) %>% arrange(path_a_p)
  
  n_sig_unc <- sum(d_sub$path_a_p < 0.05, na.rm = TRUE)
  n_sig_adj <- sum(d_sub$path_a_p_adj < 0.05, na.rm = TRUE)
  n_total   <- nrow(d_sub)
  
  hdr <- sprintf("--- %s (n = %d features) ---", con$label, n_total)
  cat("\n", hdr, "\n", sep = "")
  summary_lines <- c(summary_lines, hdr)
  
  l1 <- sprintf("  Uncorrected p < 0.05:  %d / %d (%.1f%%)",
                n_sig_unc, n_total, 100 * n_sig_unc / n_total)
  l2 <- sprintf("  FDR-adjusted p < 0.05: %d / %d (%.1f%%)",
                n_sig_adj, n_total, 100 * n_sig_adj / n_total)
  cat(l1, "\n", l2, "\n", sep = "")
  summary_lines <- c(summary_lines, l1, l2)
  
  if (n_sig_unc > 0) {
    cat("\n  Significant features (sorted by p):\n")
    summary_lines <- c(summary_lines, "", "  Significant features (sorted by p):")
    sig_features <- d_sub %>% filter(path_a_p < 0.05)
    
    for (i in seq_len(nrow(sig_features))) {
      r <- sig_features[i, ]
      fdr_mark <- if (r$path_a_p_adj < 0.05) " [FDR<.05]" else ""
      line <- sprintf("    %-45s  beta = %+6.3f  (95%% CI: %+6.3f to %+6.3f),  t(%5.1f) = %+6.2f,  p = %.4f%s",
                      r$Feature, r$path_a_coef, r$path_a_CI_lo, r$path_a_CI_hi,
                      r$path_a_df, r$path_a_t, r$path_a_p, fdr_mark)
      cat(line, "\n")
      summary_lines <- c(summary_lines, line)
    }
  }
}

out_txt <- file.path(OUTPUT_DIR, "path_a_summary.txt")
writeLines(summary_lines, out_txt)
cat(sprintf("\nSaved: %s\n\n", out_txt))

# ============================================================================
# DONE
# ============================================================================

cat("════════════════════════════════════════════════════════════════════════\n")
cat("PATH-A ANALYSIS COMPLETE\n")
cat("════════════════════════════════════════════════════════════════════════\n")
cat(sprintf("Output directory: %s/\n", OUTPUT_DIR))
cat("Files:\n")
cat("  * path_a_results.csv   Stats table (input to s06_plot_path_a_heatmap.m)\n")
cat("  * path_a_summary.txt   Paper-ready stats summary\n")
