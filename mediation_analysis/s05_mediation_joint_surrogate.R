# ============================================================================
# s05 - MULTIVARIATE MEDIATION MODEL — ALL FEATURES SIMULTANEOUSLY + SURROGATE TEST
#
# Joint path-b: outcome ~ s_feat1 + ... + s_featK + treat + (1|subject)
# Path-a:       per-feature (feature_i ~ treat + (1|subject))
# ACME per feature = unique mediated effect, controlling for all others.
#
# Surrogate: within-subject treat-label permutation builds empirical null.
# Perm p-value: P(null_ACME >= observed_ACME), one-sided right-tail.
#
# Output table: ~18 manuscript-ready columns (path-a excluded; CIs removed).
#
# Runtime note: 1000 perms × 3 contrasts × K features is CPU-intensive.
# Parallelise the permutation for-loop with parallel::mclapply for speedup.
#
# Author: Junheng Li
# ============================================================================

library(lme4)
library(lmerTest)
library(mediation)
library(dplyr)

# ============================================================================
# USER SETTINGS
# ============================================================================

DATA_FILE  <- "Results_SelectedFeatures_RelBase.csv"
OUTPUT_DIR <- "Mediation_JointModel_RelBase"

OUTCOMES       <- c("imag_acc_post")
OUTCOME_LABELS <- c("Image Accuracy")

META_COLS <- c("global_trial_idx", "trial_idx", "subject", "condition",
               "word_acc_post", "imag_acc_post", "old_new", "word",
               "block_post", "NewManscore1")

CONTRASTS <- list(
  list(treat = "TICue",   control = "SoundOnly", label = "TICue_vs_SoundOnly"),
  list(treat = "TIB4Cue", control = "SoundOnly", label = "TIB4Cue_vs_SoundOnly"),
  list(treat = "TICue",   control = "TIB4Cue",   label = "TICue_vs_TIB4Cue")
)

N_SIMS <- 1000       # QB sims for observed models
set.seed(2026)

# ── Surrogate ────────────────────────────────────────────────────────────────
RUN_SURROGATE    <- TRUE
N_PERMS          <- 1000    # permutations per contrast
N_SIMS_PERM      <- 200     # reduced QB sims per permutation refit
PERM_SEED        <- 2026    # per-perm seed = PERM_SEED + p

# ── Misc ─────────────────────────────────────────────────────────────────────
SAVE_SIM_DRAWS      <- TRUE
MIN_TRIALS_PER_CELL <- 5

# ============================================================================
# SETUP
# ============================================================================

dir.create(OUTPUT_DIR, showWarnings = FALSE)

cat("════════════════════════════════════════════════════════════════════════\n")
cat("   MULTIVARIATE MEDIATION + SURROGATE TESTING\n")
cat("════════════════════════════════════════════════════════════════════════\n\n")
cat("Joint path-b: outcome ~ s_feat1 + ... + s_featK + treat + (1|subject)\n")
cat("Each ACME = unique mediated effect (all other features as covariates)\n\n")
if (RUN_SURROGATE) {
  cat(sprintf("Surrogate: %d within-subject permutations per contrast (%d QB sims each)\n",
              N_PERMS, N_SIMS_PERM))
  cat("           Parallelise the perm for-loop for large speedup.\n\n")
}

# ============================================================================
# LOAD DATA & AUTO-DETECT FEATURES
# ============================================================================

data_raw      <- read.csv(DATA_FILE)
existing_meta <- intersect(META_COLS, names(data_raw))
all_numeric   <- names(data_raw)[sapply(data_raw, is.numeric)]
feature_cols  <- setdiff(all_numeric, existing_meta)

cat(sprintf("Data: %d trials, %d subjects\n",
            nrow(data_raw), length(unique(data_raw$subject))))
cat(sprintf("Features (%d): %s\n\n",
            length(feature_cols), paste(feature_cols, collapse = ", ")))

data_raw$subject    <- factor(data_raw$subject)
data_raw$block_post <- factor(data_raw$block_post)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

apply_trial_filter <- function(data, treat_val, control_val, min_trials) {
  d <- data %>%
    filter(condition %in% c(treat_val, control_val)) %>%
    mutate(treat = as.integer(condition == treat_val))
  tc <- d %>% count(subject, treat, name = "n_trials")
  valid_subj <- tc %>%
    group_by(subject) %>%
    filter(n() == 2, all(n_trials >= min_trials)) %>%
    pull(subject) %>% unique()
  list(
    data      = filter(d, subject %in% valid_subj) %>% droplevels(),
    n_excluded = length(setdiff(unique(d$subject), valid_subj))
  )
}

fit_path_a <- function(d, feat_name) {
  d$mediator <- as.numeric(scale(d[[feat_name]]))
  m <- tryCatch(
    lme4::lmer(mediator ~ treat + (1|subject), data = d,
               control = lmerControl(optimizer = "bobyqa",
                                     optCtrl = list(maxfun = 100000))),
    error = function(e) NULL)
  if (is.null(m)) return(NULL)
  # lmerTest for Satterthwaite p-value
  mt <- tryCatch(
    lmerTest::lmer(mediator ~ treat + (1|subject), data = d,
                   control = lmerControl(optimizer = "bobyqa",
                                         optCtrl = list(maxfun = 100000))),
    error = function(e) NULL)
  cf  <- if (!is.null(mt)) summary(mt)$coefficients else summary(m)$coefficients
  pa_p <- tryCatch(cf["treat", "Pr(>|t|)"],
                   error = function(e) 2 * pnorm(-abs(cf["treat", "t value"])))
  list(model = m, coef = cf["treat","Estimate"], se = cf["treat","Std. Error"],
       t_val = cf["treat","t value"], p = pa_p, singular = isSingular(m))
}

fit_joint_path_b <- function(d, feature_names, outcome_var) {
  for (f in feature_names)
    d[[paste0("s_", f)]] <- as.numeric(scale(d[[f]]))
  rhs <- paste(c(paste0("s_", feature_names), "treat"), collapse = " + ")
  formula_str <- paste0(outcome_var, " ~ ", rhs, " + (1|subject)")
  m <- tryCatch(
    glmer(as.formula(formula_str), data = d, family = binomial(link = "logit"),
          control = glmerControl(optimizer = "bobyqa",
                                 optCtrl = list(maxfun = 200000))),
    error = function(e) NULL)
  if (is.null(m)) return(NULL)
  list(model = m, singular = isSingular(m))
}

run_feature_mediation <- function(d, feat_name, path_a_result, joint_model,
                                  scaled_feat_name, n_sims) {
  med <- tryCatch(
    mediate(path_a_result$model, joint_model$model,
            treat = "treat", mediator = scaled_feat_name,
            sims = n_sims, boot = FALSE),
    error = function(e) NULL)
  if (is.null(med)) return(NULL)

  sims     <- med$d0.sims
  acme_m   <- mean(sims)
  acme_sd  <- sd(sims)
  acme_cov <- ifelse(acme_sd > .Machine$double.eps, acme_m / acme_sd, NA)

  cf_y    <- summary(joint_model$model)$coefficients
  pb_est  <- tryCatch(cf_y[scaled_feat_name, "Estimate"],  error = function(e) NA)
  pb_p    <- tryCatch(cf_y[scaled_feat_name, "Pr(>|z|)"],  error = function(e) NA)

  list(
    summary = data.frame(
      PathB_estimate = pb_est,  PathB_p    = pb_p,
      ACME_mean      = acme_m,  ACME_sd    = acme_sd,  ACME_CoV = acme_cov,
      QB_p           = med$d0.p,
      Direct_effect  = med$z0,  ADE_p      = med$z0.p,
      Total_effect   = med$tau.coef, Total_p = med$tau.p,
      Prop_mediated  = med$n0,
      stringsAsFactors = FALSE
    ),
    acme_sims = sims
  )
}

# ============================================================================
# MAIN LOOP
# ============================================================================

cat("════════════════════════════════════════════════════════════════════════\n")
cat("RUNNING\n")
cat("════════════════════════════════════════════════════════════════════════\n\n")

all_results    <- data.frame()
all_sim_draws  <- list()
perm_null_acme <- list()   # keyed by "outcome|contrast|feature"

for (out_idx in seq_along(OUTCOMES)) {
  outcome_var   <- OUTCOMES[out_idx]
  outcome_label <- OUTCOME_LABELS[out_idx]

  cat(sprintf("OUTCOME: %s\n%s\n\n", outcome_label, strrep("─", 70)))

  for (con in CONTRASTS) {
    cat(sprintf("  Contrast: %s\n", con$label))

    # Step 1: Filter
    filt <- apply_trial_filter(data_raw, con$treat, con$control, MIN_TRIALS_PER_CELL)
    d    <- filt$data
    if (filt$n_excluded > 0)
      cat(sprintf("    [-%d subj filtered]\n", filt$n_excluded))
    if (nrow(d) < 50 || length(unique(d$treat)) < 2 || length(unique(d$subject)) < 3) {
      cat("    SKIP: insufficient data\n\n"); next
    }
    cat(sprintf("    N = %d trials, %d subjects\n", nrow(d), length(unique(d$subject))))

    # Step 2: Path-a per feature
    path_a_results <- list(); valid_features <- c()
    for (feat in feature_cols) {
      pa <- fit_path_a(d, feat)
      if (!is.null(pa)) {
        path_a_results[[feat]] <- pa
        valid_features <- c(valid_features, feat)
        if (pa$singular) cat(sprintf("      %s [path-a singular]\n", feat))
      }
    }
    if (length(valid_features) == 0) { cat("    SKIP: no path-a models\n\n"); next }
    cat(sprintf("    %d/%d path-a OK\n", length(valid_features), length(feature_cols)))

    # Step 3: Joint path-b
    cat("    Fitting joint path-b...\n")
    joint <- fit_joint_path_b(d, valid_features, outcome_var)
    if (is.null(joint)) { cat("    SKIP: joint path-b failed\n\n"); next }
    if (joint$singular) cat("    ⚠ Joint path-b singular\n") else cat("    Joint path-b OK\n")

    # Step 4: Mediation per feature
    cat("    ACME estimation:\n")
    for (feat in valid_features) {
      cat(sprintf("      %s...", feat))
      res <- tryCatch(
        run_feature_mediation(d, feat, path_a_results[[feat]], joint,
                              paste0("s_", feat), N_SIMS),
        error = function(e) NULL)
      if (is.null(res)) { cat(" SKIP\n"); next }

      row <- data.frame(
        Outcome = outcome_label, Contrast = con$label, Feature = feat,
        N_trials = nrow(d), N_subjects = length(unique(d$subject)),
        Joint_singular = joint$singular,
        stringsAsFactors = FALSE
      )
      row <- cbind(row, res$summary)
      all_results <- rbind(all_results, row)
      all_sim_draws[[paste(outcome_label, con$label, feat, sep = "|")]] <- res$acme_sims
      cat(" OK\n")
    }

    # Step 5: Surrogate permutation testing
    if (RUN_SURROGATE && nrow(all_results) > 0) {
      cat(sprintf("\n    Surrogate: %d permutations (N_SIMS_PERM = %d)...\n",
                  N_PERMS, N_SIMS_PERM))

      perm_null_tmp <- setNames(
        lapply(valid_features, function(f) rep(NA_real_, N_PERMS)),
        valid_features
      )
      n_joint_ok <- 0L
      t0 <- Sys.time()

      for (p in seq_len(N_PERMS)) {
        set.seed(PERM_SEED + p)

        # Within-subject treat shuffle
        d_perm <- d %>%
          group_by(subject) %>%
          mutate(treat = sample(treat)) %>%
          ungroup()

        # Path-a
        perm_pa <- list(); perm_vf <- c()
        for (feat in valid_features) {
          pa_p <- suppressWarnings(suppressMessages(fit_path_a(d_perm, feat)))
          if (!is.null(pa_p)) { perm_pa[[feat]] <- pa_p; perm_vf <- c(perm_vf, feat) }
        }
        if (length(perm_vf) == 0) next

        # Joint path-b
        perm_joint <- suppressWarnings(suppressMessages(
          fit_joint_path_b(d_perm, perm_vf, outcome_var)
        ))
        if (is.null(perm_joint)) next
        n_joint_ok <- n_joint_ok + 1L

        # Mediation per feature
        for (feat in perm_vf) {
          res_p <- tryCatch(
            suppressWarnings(suppressMessages(
              run_feature_mediation(d_perm, feat, perm_pa[[feat]],
                                    perm_joint, paste0("s_", feat), N_SIMS_PERM)
            )),
            error = function(e) NULL)
          if (!is.null(res_p))
            perm_null_tmp[[feat]][p] <- res_p$summary$ACME_mean
        }

        if (p %% 200 == 0 || p == N_PERMS) {
          dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
          cat(sprintf("      perm %4d/%d  %.1f min  %d joint-b OK\n",
                      p, N_PERMS, dt, n_joint_ok))
        }
      }

      for (feat in valid_features) {
        key <- paste(outcome_label, con$label, feat, sep = "|")
        perm_null_acme[[key]] <- perm_null_tmp[[feat]]
      }
      valid_n <- sapply(perm_null_tmp, function(x) sum(!is.na(x)))
      cat(sprintf("    Done: %d joint-b converged. Valid draws/feature: %d–%d\n\n",
                  n_joint_ok, min(valid_n), max(valid_n)))
    }

    cat("\n")
  }
}

cat(sprintf("Total feature-contrast models: %d\n\n", nrow(all_results)))

# ============================================================================
# PERMUTATION P-VALUES
# ============================================================================

if (RUN_SURROGATE && length(perm_null_acme) > 0) {
  cat("Computing permutation p-values (one-sided: observed >= null)...\n")
  all_results$Perm_p       <- NA_real_
  all_results$N_perm_valid <- NA_integer_

  for (i in seq_len(nrow(all_results))) {
    key <- paste(all_results$Outcome[i], all_results$Contrast[i],
                 all_results$Feature[i], sep = "|")
    nv <- perm_null_acme[[key]]
    if (!is.null(nv)) {
      nv_valid <- sum(!is.na(nv))
      if (nv_valid > 0) {
        all_results$Perm_p[i]       <- (1 + sum(nv >= all_results$ACME_mean[i], na.rm = TRUE)) /
                                        (1 + nv_valid)
        all_results$N_perm_valid[i] <- nv_valid
      }
    }
  }
  cat("Done.\n\n")
}

# ============================================================================
# FDR CORRECTION (within each Outcome × Contrast)
# ============================================================================

all_results <- all_results %>%
  group_by(Outcome, Contrast) %>%
  mutate(
    PathB_p_FDR = p.adjust(PathB_p, method = "fdr"),
    QB_p_FDR    = p.adjust(QB_p,    method = "fdr")
  ) %>%
  ungroup()

if (RUN_SURROGATE && "Perm_p" %in% names(all_results)) {
  all_results <- all_results %>%
    group_by(Outcome, Contrast) %>%
    mutate(Perm_p_FDR = p.adjust(Perm_p, method = "fdr")) %>%
    ungroup()
}

# ============================================================================
# STREAMLINED OUTPUT TABLE  (~18 manuscript-ready columns)
# ============================================================================
# Removed: path-a stats, CIs (ACME/ADE/Total/PropMed), ACME_median,
#          ADE_p, Total_p, PropMed CIs & p, path_cp, N_excluded,
#          singular_path_a.

results_clean <- all_results %>%
  dplyr::select(
    Outcome, Contrast, Feature,
    N_trials, N_subjects,
    Joint_singular,
    PathB_estimate, PathB_p, PathB_p_FDR,
    ACME_mean, ACME_sd, ACME_CoV,
    QB_p, QB_p_FDR,
    any_of(c("Perm_p", "Perm_p_FDR")),
    Direct_effect, Prop_mediated
  )

write.csv(results_clean,
          file.path(OUTPUT_DIR, "multivariate_mediation_results.csv"),
          row.names = FALSE)
cat(sprintf("Main table saved: %s/multivariate_mediation_results.csv\n",
            OUTPUT_DIR))
cat(sprintf("  %d rows × %d columns\n\n", nrow(results_clean), ncol(results_clean)))

# Per-outcome split
for (ol in OUTCOME_LABELS) {
  d_out <- filter(results_clean, Outcome == ol)
  write.csv(d_out,
            file.path(OUTPUT_DIR, paste0("multivar_med_",
                                          gsub(" ", "_", tolower(ol)), ".csv")),
            row.names = FALSE)
}

# ============================================================================
# SAVE SIMULATION DRAWS
# ============================================================================

if (SAVE_SIM_DRAWS && length(all_sim_draws) > 0) {
  saveRDS(all_sim_draws, file.path(OUTPUT_DIR, "acme_simulation_draws.rds"))
  draws_mat <- do.call(cbind, lapply(all_sim_draws, function(x) {
    v <- as.numeric(x); if (length(v) != N_SIMS) v <- rep(NA, N_SIMS); v
  }))
  colnames(draws_mat) <- names(all_sim_draws)
  write.csv(draws_mat, file.path(OUTPUT_DIR, "acme_simulation_draws_matrix.csv"),
            row.names = FALSE)
}

if (RUN_SURROGATE && length(perm_null_acme) > 0) {
  saveRDS(perm_null_acme, file.path(OUTPUT_DIR, "perm_null_acme.rds"))
  null_mat <- do.call(cbind, lapply(perm_null_acme, function(x) {
    v <- as.numeric(x); if (length(v) != N_PERMS) v <- rep(NA, N_PERMS); v
  }))
  colnames(null_mat) <- names(perm_null_acme)
  write.csv(null_mat, file.path(OUTPUT_DIR, "perm_null_acme_matrix.csv"),
            row.names = FALSE)
  cat("Perm null draws saved.\n\n")
}

# ============================================================================
# CoV ANALYSIS (omnibus: is CoV > 0 across all features per contrast?)
# ============================================================================

cat("════════════════════════════════════════════════════════════════════════\n")
cat("CoV ANALYSIS\n")
cat("════════════════════════════════════════════════════════════════════════\n\n")

cov_tests <- data.frame()
for (ol in OUTCOME_LABELS) {
  for (cl in sapply(CONTRASTS, `[[`, "label")) {
    sub <- results_clean %>% filter(Outcome == ol, Contrast == cl, !is.na(ACME_CoV))
    if (nrow(sub) < 3) next
    cv   <- sub$ACME_CoV
    t_r  <- tryCatch(t.test(cv, mu = 0, alternative = "greater"),     error = function(e) NULL)
    w_r  <- tryCatch(wilcox.test(cv, mu = 0, alternative = "greater",
                                  exact = FALSE),                       error = function(e) NULL)
    sw_r <- tryCatch(shapiro.test(cv),                                  error = function(e) NULL)
    cov_tests <- rbind(cov_tests, data.frame(
      Outcome = ol, Contrast = cl, N_features = length(cv),
      CoV_mean = mean(cv), CoV_sd = sd(cv), CoV_median = median(cv),
      t_stat   = if (!is.null(t_r))  t_r$statistic  else NA,
      t_p      = if (!is.null(t_r))  t_r$p.value    else NA,
      wilcox_V = if (!is.null(w_r))  w_r$statistic  else NA,
      wilcox_p = if (!is.null(w_r))  w_r$p.value    else NA,
      shapiro_p = if (!is.null(sw_r)) sw_r$p.value  else NA,
      stringsAsFactors = FALSE
    ))
    cat(sprintf("%s | %s\n  N=%d  CoV mean=%+.3f  t-p=%.4f  Wilcoxon-p=%.4f\n\n",
                ol, cl, length(cv), mean(cv),
                if (!is.null(t_r)) t_r$p.value else NA,
                if (!is.null(w_r)) w_r$p.value else NA))
  }
}
write.csv(cov_tests, file.path(OUTPUT_DIR, "cov_test_results.csv"), row.names = FALSE)

# ============================================================================
# CONSOLE SUMMARY TABLE
# ============================================================================

cat("════════════════════════════════════════════════════════════════════════\n")
cat("ACME SUMMARY (sorted by QB p-value)\n")
cat("════════════════════════════════════════════════════════════════════════\n\n")

for (ol in OUTCOME_LABELS) {
  for (cl in sapply(CONTRASTS, `[[`, "label")) {
    sub <- results_clean %>% filter(Outcome == ol, Contrast == cl) %>% arrange(QB_p)
    if (nrow(sub) == 0) next
    cat(sprintf("%s | %s\n", ol, cl))
    disp_cols <- intersect(c("Feature","ACME_mean","ACME_sd","ACME_CoV",
                              "QB_p","QB_p_FDR","Perm_p","Perm_p_FDR"), names(sub))
    print(sub[, disp_cols], digits = 4, row.names = FALSE)
    cat("\n")
  }
}

# ============================================================================
# FINAL
# ============================================================================

cat("\n════════════════════════════════════════════════════════════════════════\n")
cat("COMPLETE\n")
cat(sprintf("Output: %s/\n\n", OUTPUT_DIR))
cat("Output table columns:\n")
cat(paste0("  ", names(results_clean), collapse = "\n"), "\n\n")
cat("Key:\n")
cat("  PathB_*         Unique feature coefficient in joint outcome model\n")
cat("  ACME_mean/sd    Average Causal Mediation Effect (unique, +/- QB sim SD)\n")
cat("  ACME_CoV        Consistency index (ACME_mean / ACME_sd)\n")
cat("  QB_p/FDR        Quasi-Bayesian mediation p-value (Imai et al. 2010)\n")
cat("  Perm_p/FDR      Within-subject permutation p-value for ACME\n")
cat("  Direct_effect   Average Direct Effect (ADE)\n")
cat("  Prop_mediated   Proportion of total effect mediated\n")
cat("════════════════════════════════════════════════════════════════════════\n")
