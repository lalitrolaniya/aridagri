# =============================================================================
# BASIC DESIGNS: CRD, RCBD AND COMBINED RCBD OVER ENVIRONMENTS
# Package: aridagri
# =============================================================================

# Build the treatment means table shared by CRD and RCBD results.
.means_table <- function(d, response, treatment, lsm, ph, methods) {
  y <- d[[response]]
  g <- d[[treatment]]
  sdv <- tapply(y, g, stats::sd)
  tab <- data.frame(Treatment = names(lsm$means),
                    Mean = as.numeric(lsm$means),
                    SD = as.numeric(sdv[names(lsm$means)]),
                    N = as.numeric(lsm$n[names(lsm$means)]),
                    SE = sqrt(pmax(diag(lsm$vcov), 0)),
                    stringsAsFactors = FALSE)
  first <- setdiff(methods, "dunnett")
  if (length(first) > 0 && !is.null(ph[[first[1]]])) {
    tab$Group <- ph[[first[1]]]$letters[tab$Treatment]
  }
  tab <- tab[order(-tab$Mean), , drop = FALSE]
  rownames(tab) <- NULL
  tab
}

# Name of the method whose letters appear in the means table.
.letters_method <- function(methods) {
  m <- setdiff(methods, "dunnett")
  if (length(m) == 0) return(NA_character_)
  switch(m[1], lsd = "Fisher's LSD", duncan = "Duncan's MRT",
         tukey = "Tukey HSD", snk = "Student-Newman-Keuls",
         scheffe = "Scheffe", bonferroni = "Bonferroni")
}

#' Completely Randomized Design (CRD) ANOVA
#'
#' @description
#' Analysis of variance for a completely randomized design with standard
#' errors, critical difference, coefficient of variation, post-hoc letter
#' groupings and assumption checks. Unequal replication is handled
#' exactly.
#'
#' @param data Data frame containing the experimental data.
#' @param response Name of the response variable (character).
#' @param treatment Name of the treatment column.
#' @param posthoc Post-hoc test(s): \code{"lsd"}, \code{"duncan"},
#'   \code{"tukey"}, \code{"snk"}, \code{"scheffe"}, \code{"bonferroni"},
#'   \code{"dunnett"}, a vector of these, \code{"all"} or \code{"none"}.
#'   The first method supplies the letters in the means table.
#' @param alpha Significance level for the CD and the post-hoc tests
#'   (default 0.05).
#' @param verbose Logical. If TRUE (default), prints the report.
#' @param control Control level for Dunnett's test (default: first level).
#'
#' @details
#' SEm = sqrt(MSE / r), SEd = sqrt(2 MSE / r) and CD = t(alpha/2, error
#' df) x SEd, where r is the number of replications (the harmonic mean
#' when replication is unequal; the letters then use the exact SE of each
#' pair). CV = 100 sqrt(MSE) / grand mean. Rows with a missing response
#' are dropped with a message.
#'
#' @return An object of class \code{aridagri_anova} with
#'   \code{anova_table}, \code{model}, \code{grand_mean}, \code{cv},
#'   \code{ms_error}, \code{df_error}, \code{se_mean} (SEm),
#'   \code{se_diff} (SEd), \code{cd} (at \code{alpha}), \code{cd_5},
#'   \code{cd_1}, \code{treatment_means} (mean, SD, N, SE and letters),
#'   \code{posthoc}, \code{diagnostics}, \code{significance} and
#'   \code{balanced}.
#'
#' @references
#' Gomez, K.A. and Gomez, A.A. (1984). \emph{Statistical Procedures for
#' Agricultural Research}, 2nd ed. John Wiley and Sons, New York.
#'
#' @examples
#' data <- data.frame(
#'   treatment = rep(c("T1", "T2", "T3", "T4"), each = 5),
#'   yield = c(1210, 1185, 1242, 1198, 1225, 1340, 1362, 1318, 1375, 1350,
#'             1105, 1090, 1122, 1131, 1098, 1445, 1462, 1430, 1471, 1452)
#' )
#' res <- anova_crd(data, response = "yield", treatment = "treatment",
#'                  posthoc = "duncan")
#' res$treatment_means
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
anova_crd <- function(data, response, treatment, posthoc = "lsd", alpha = 0.05,
                      verbose = TRUE, control = NULL) {
  .check_alpha(alpha)
  methods <- .resolve_methods(posthoc)
  d <- .prep_data(data, response, treatment, verbose)
  trt <- d[[treatment]]
  k <- nlevels(trt)
  if (k < 2) stop("Need at least 2 treatments.", call. = FALSE)
  n_i <- as.numeric(table(trt))
  if (nrow(d) - k < 1) {
    stop("No error degrees of freedom: the treatments need replication.",
         call. = FALSE)
  }
  model <- stats::aov(.formula(response, .bq(treatment)), data = d)
  anova_table <- .clean_rownames(stats::anova(model))
  ms_e <- anova_table[["Mean Sq"]][2]
  df_e <- anova_table[["Df"]][2]
  p_trt <- anova_table[["Pr(>F)"]][1]
  balanced <- length(unique(n_i)) == 1
  r_h <- 1 / mean(1 / n_i)
  gm <- mean(d[[response]])
  sem <- sqrt(ms_e / r_h)
  sed <- sqrt(2) * sem
  lsm <- .ls_means_model(model, treatment)
  ph <- .posthoc_set(model, treatment, methods, alpha, control = control,
                     response = response, lsm = lsm)
  means_df <- .means_table(d, response, treatment, lsm, ph, methods)
  diagnostics <- check_assumptions(model, factor = treatment, verbose = FALSE)
  result <- list(
    design = "CRD",
    response = response,
    factors = c(treatment = treatment),
    anova_table = anova_table,
    model = model,
    grand_mean = gm,
    cv = 100 * sqrt(ms_e) / gm,
    ms_error = ms_e,
    df_error = df_e,
    se_mean = sem,
    se_diff = sed,
    cd = stats::qt(1 - alpha / 2, df_e) * sed,
    cd_5 = stats::qt(0.975, df_e) * sed,
    cd_1 = stats::qt(0.995, df_e) * sed,
    alpha = alpha,
    treatment_means = means_df,
    letters_method = .letters_method(methods),
    posthoc = ph,
    diagnostics = diagnostics,
    significance = .sig_code(p_trt),
    balanced = balanced,
    replications = stats::setNames(n_i, levels(trt))
  )
  class(result) <- c("aridagri_anova", "list")
  if (verbose) print(result)
  invisible(result)
}


#' Randomized Complete Block Design (RCBD) ANOVA
#'
#' @description
#' Analysis of variance for a randomized complete block design with
#' standard errors, critical difference, CV, relative efficiency over a
#' CRD, post-hoc letter groupings and assumption checks. Missing plots are
#' handled by least squares: treatments are tested after blocks, means are
#' least-squares (adjusted) means and every pair is compared with its own
#' standard error.
#'
#' @inheritParams anova_crd
#' @param block Name of the block (replication) column.
#'
#' @details
#' Relative efficiency of the RCBD over a CRD follows Gomez and Gomez
#' (1984): \deqn{RE = \frac{(r - 1) E_b + r (t - 1) E_e}{(r t - 1) E_e}}
#' where \eqn{E_b} and \eqn{E_e} are the block and error mean squares,
#' multiplied by the precision adjustment factor
#' \eqn{k = \frac{[(r-1)(t-1)+1][t(r-1)+3]}{[(r-1)(t-1)+3][t(r-1)+1]}}
#' when the error degrees of freedom are below 20. Values above 100\% mean
#' that blocking reduced the error.
#'
#' @return An object of class \code{aridagri_anova}. In addition to the
#'   components described in \code{\link{anova_crd}} it contains
#'   \code{relative_efficiency} (percent, adjusted),
#'   \code{relative_efficiency_unadjusted} and \code{missing_plots}.
#'
#' @references
#' Gomez, K.A. and Gomez, A.A. (1984). \emph{Statistical Procedures for
#' Agricultural Research}, 2nd ed. John Wiley and Sons, New York.
#'
#' @examples
#' data(mothbean_rbd)
#' res <- anova_rbd(mothbean_rbd, response = "seed_yield",
#'                  treatment = "treatment", block = "block",
#'                  posthoc = "duncan")
#' res$treatment_means
#' res$relative_efficiency
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
anova_rbd <- function(data, response, treatment, block, posthoc = "lsd",
                      alpha = 0.05, verbose = TRUE, control = NULL) {
  .check_alpha(alpha)
  methods <- .resolve_methods(posthoc)
  d <- .prep_data(data, response, c(treatment, block), verbose)
  trt <- d[[treatment]]
  blk <- d[[block]]
  k <- nlevels(trt)
  r <- nlevels(blk)
  if (k < 2) stop("Need at least 2 treatments.", call. = FALSE)
  if (r < 2) stop("Need at least 2 blocks.", call. = FALSE)
  cell <- table(blk, trt)
  if (any(colSums(cell) == 0)) {
    stop("Some treatments have no observations.", call. = FALSE)
  }
  balanced <- all(cell == 1)
  missing_plots <- sum(cell == 0)
  if (!balanced && verbose) {
    if (missing_plots > 0) {
      message(missing_plots, " missing plot(s): treatments are tested after ",
              "blocks and means are least-squares means.")
    }
    if (any(cell > 1)) {
      message("Some treatments occur more than once in a block; the error ",
              "term includes block x treatment variation.")
    }
  }
  model <- stats::aov(.formula(response, paste(.bq(block), "+",
                                               .bq(treatment))), data = d)
  anova_table <- .clean_rownames(stats::anova(model))
  if (anova_table[["Df"]][3] < 1) {
    stop("No error degrees of freedom left.", call. = FALSE)
  }
  ms_b <- anova_table[["Mean Sq"]][1]
  ms_e <- anova_table[["Mean Sq"]][3]
  df_e <- anova_table[["Df"]][3]
  p_trt <- anova_table[["Pr(>F)"]][2]
  n_i <- colSums(cell)
  r_h <- 1 / mean(1 / n_i)
  gm <- mean(d[[response]])
  sem <- sqrt(ms_e / r_h)
  sed <- sqrt(2) * sem
  re_raw <- NA_real_
  re_adj <- NA_real_
  if (balanced) {
    re_raw <- 100 * ((r - 1) * ms_b + r * (k - 1) * ms_e) /
      ((r * k - 1) * ms_e)
    adj <- if (df_e < 20) {
      (((r - 1) * (k - 1) + 1) * (k * (r - 1) + 3)) /
        (((r - 1) * (k - 1) + 3) * (k * (r - 1) + 1))
    } else 1
    re_adj <- re_raw * adj
  }
  lsm <- .ls_means_model(model, treatment)
  ph <- .posthoc_set(model, treatment, methods, alpha, control = control,
                     response = response, lsm = lsm)
  means_df <- .means_table(d, response, treatment, lsm, ph, methods)
  diagnostics <- check_assumptions(model, factor = treatment, verbose = FALSE)
  result <- list(
    design = "RBD",
    response = response,
    factors = c(treatment = treatment, block = block),
    anova_table = anova_table,
    model = model,
    grand_mean = gm,
    cv = 100 * sqrt(ms_e) / gm,
    ms_error = ms_e,
    df_error = df_e,
    se_mean = sem,
    se_diff = sed,
    cd = stats::qt(1 - alpha / 2, df_e) * sed,
    cd_5 = stats::qt(0.975, df_e) * sed,
    cd_1 = stats::qt(0.995, df_e) * sed,
    alpha = alpha,
    relative_efficiency = re_adj,
    relative_efficiency_unadjusted = re_raw,
    treatment_means = means_df,
    letters_method = .letters_method(methods),
    posthoc = ph,
    diagnostics = diagnostics,
    significance = .sig_code(p_trt),
    balanced = balanced,
    missing_plots = missing_plots,
    replications = n_i
  )
  class(result) <- c("aridagri_anova", "list")
  if (verbose) print(result)
  invisible(result)
}


#' Print Method for CRD and RCBD Results
#'
#' @param x An object of class \code{aridagri_anova}.
#' @param digits Number of decimal places for means. Default 2.
#' @param ... Unused.
#'
#' @return Invisibly returns \code{x}. Called for its side effect of
#'   printing the ANOVA report.
#'
#' @export
print.aridagri_anova <- function(x, digits = 2, ...) {
  if (!is.null(x$design) && x$design %in% c("CRD", "RBD")) {
    .print_crd_rbd(x, digits)
  } else {
    .print_generic_anova(x, digits)
  }
  invisible(x)
}

.print_crd_rbd <- function(x, digits) {
  f <- x$factors
  title <- if (x$design == "CRD") "COMPLETELY RANDOMIZED DESIGN (CRD)" else
    "RANDOMIZED COMPLETE BLOCK DESIGN (RCBD)"
  info <- c(Response = x$response, Treatment = f[["treatment"]],
            Treatments = length(x$replications))
  if (x$design == "RBD") {
    info <- c(info, Blocks = nlevels(stats::model.frame(x$model)[[f[["block"]]]]))
  }
  info <- c(info, Observations = sum(x$replications))
  .print_info(title, info)
  tab <- x$anova_table
  labs <- if (x$design == "CRD") c("Treatment", "Error") else
    c("Block", "Treatment", "Error")
  .print_anova_df(tab, labs)
  st <- c(`Grand mean` = x$grand_mean, `CV (%)` = x$cv,
          SEm = x$se_mean, SEd = x$se_diff)
  st[[sprintf("CD (%g%%)", 100 * x$alpha)]] <- x$cd
  .print_stats_line(st, 3)
  if (!isTRUE(x$balanced)) {
    cat(" (SEm, SEd and CD use the harmonic mean replication; letters use",
        "the exact SE of each pair.)\n")
  }
  if (x$design == "RBD" && !is.na(x$relative_efficiency)) {
    cat(sprintf(" Relative efficiency of RCBD over CRD = %.1f%%\n",
                x$relative_efficiency))
  }
  ttl <- if (is.na(x$letters_method)) "Treatment means" else
    sprintf("Treatment means (%s, alpha = %g)", x$letters_method, x$alpha)
  .print_means(x$treatment_means, ttl, digits)
  p_trt <- x$anova_table[["Pr(>F)"]][if (x$design == "CRD") 1 else 2]
  if (!is.na(x$letters_method)) {
    cat(" Means sharing a letter do not differ significantly.\n")
    if (!is.na(p_trt) && p_trt >= x$alpha) {
      cat(" Note: the treatment F-test is not significant, so pairwise",
          "differences\n should not be interpreted.\n")
    }
  }
  if (!is.null(x$posthoc$dunnett)) .print_posthoc(x$posthoc["dunnett"], digits)
  invisible(NULL)
}

# Fallback printer for result lists that carry an ANOVA table.
.print_generic_anova <- function(x, digits) {
  .heading(if (!is.null(x$design)) toupper(x$design) else "ANOVA")
  if (!is.null(x$anova_table)) {
    print(.add_sig(.round_anova(x$anova_table)))
  }
  if (!is.null(x$grand_mean)) {
    cat(sprintf("\n Grand mean = %.3f\n", x$grand_mean))
  }
  if (!is.null(x$cv)) {
    cv <- x$cv
    cat(" CV (%):", paste(names(cv), round(cv, 2), sep = " ", collapse = "; "),
        "\n")
  }
  invisible(NULL)
}


#' Combined (Pooled) Analysis of RCBD Experiments over Environments
#'
#' @description
#' Combined analysis of variance for the same set of treatments tested in a
#' randomized complete block design in several environments (years,
#' seasons or locations). Gives the analysis of each environment,
#' Bartlett's test of homogeneity of error variances, the combined ANOVA
#' with the correct error term for every source, SE and CD for treatment,
#' environment and interaction means, and letter groupings for the
#' treatment means.
#'
#' @inheritParams anova_crd
#' @param environment Name of the environment (year, season or location)
#'   column.
#' @param block Name of the block column (blocks are nested within
#'   environments; labels may repeat across environments).
#' @param env_effect Either \code{"fixed"} (default; years or seasons of
#'   interest) or \code{"random"} (locations or years treated as a sample).
#'
#' @details
#' Sources and their F-tests (McIntosh, 1983; Gomez and Gomez, 1984):
#' environment is tested against blocks within environments; the
#' treatment x environment interaction against the pooled error; the
#' treatment effect against the pooled error when environments are fixed
#' and against the treatment x environment mean square when they are
#' random. With random environments the usual restricted mixed model is
#' assumed, so environment is still tested against blocks within
#' environments. The SE and CD of treatment means use the same denominator
#' as the treatment F-test. Bartlett's chi-square tests whether the error
#' variances of the environments can be pooled; when it is significant,
#' interpret the combined analysis with care.
#'
#' @return An object of class \code{aridagri_pooled} with
#'   \code{pooled_anova} (combined ANOVA with correct F-tests),
#'   \code{individual_anovas}, \code{environment_summary} (mean, MSE, CV,
#'   F, p and CD of each environment), \code{bartlett_test},
#'   \code{treatment_means} (with letters), \code{interaction_means},
#'   \code{se_treatment}, \code{cd_treatment}, \code{se_environment},
#'   \code{cd_environment}, \code{se_interaction}, \code{cd_interaction},
#'   \code{posthoc}, \code{cv}, \code{grand_mean}, \code{mse_values} and
#'   \code{model}.
#'
#' @references
#' McIntosh, M.S. (1983). Analysis of combined experiments.
#' \emph{Agronomy Journal}, 75(1), 153-155.
#' \doi{10.2134/agronj1983.00021962007500010041x}
#'
#' Gomez, K.A. and Gomez, A.A. (1984). \emph{Statistical Procedures for
#' Agricultural Research}, 2nd ed. John Wiley and Sons, New York.
#'
#' @examples
#' data(pearlmillet_met)
#' res <- anova_rbd_pooled(pearlmillet_met, response = "grain_yield",
#'                         treatment = "genotype", environment = "environment",
#'                         block = "rep", verbose = FALSE)
#' res$pooled_anova
#' res$bartlett_test
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
anova_rbd_pooled <- function(data, response, treatment, environment, block,
                             posthoc = "lsd", alpha = 0.05, verbose = TRUE,
                             env_effect = c("fixed", "random"),
                             control = NULL) {
  .check_alpha(alpha)
  env_effect <- match.arg(env_effect)
  methods <- setdiff(.resolve_methods(posthoc), "dunnett")
  d <- .prep_data(data, response, c(treatment, environment, block), verbose)
  trt <- d[[treatment]]
  env <- d[[environment]]
  k <- nlevels(trt)
  e <- nlevels(env)
  if (k < 2) stop("Need at least 2 treatments.", call. = FALSE)
  if (e < 2) stop("Need at least 2 environments.", call. = FALSE)
  d$.blk <- interaction(env, d[[block]], drop = TRUE, sep = ":")
  r_env <- tapply(d$.blk, env, function(b) length(unique(b)))
  if (any(r_env < 2)) {
    stop("Every environment needs at least 2 blocks.", call. = FALSE)
  }
  balanced <- all(table(d$.blk, trt) == 1) && length(unique(r_env)) == 1
  if (!balanced) {
    warning("The data are not balanced (every treatment once in every block ",
            "of every environment). The combined ANOVA uses sequential sums ",
            "of squares and is approximate.", call. = FALSE)
  }
  r <- as.numeric(stats::median(r_env))
  y <- d[[response]]

  # Analysis of each environment
  env_levels <- levels(env)
  individual <- list()
  summ <- vector("list", e)
  for (i in seq_len(e)) {
    di <- d[env == env_levels[i], , drop = FALSE]
    di$.blk <- droplevels(di$.blk)
    di[[treatment]] <- droplevels(di[[treatment]])
    fi <- stats::aov(.formula(response, paste(".blk +", .bq(treatment))),
                     data = di)
    ai <- stats::anova(fi)
    rownames(ai) <- c(block, treatment, "Residuals")
    individual[[env_levels[i]]] <- ai
    mi <- mean(di[[response]])
    sed_i <- sqrt(2 * ai[["Mean Sq"]][3] / r_env[[i]])
    summ[[i]] <- data.frame(
      Environment = env_levels[i], Mean = mi, MSE = ai[["Mean Sq"]][3],
      df_error = ai[["Df"]][3], CV = 100 * sqrt(ai[["Mean Sq"]][3]) / mi,
      F_value = ai[["F value"]][2], p_value = ai[["Pr(>F)"]][2],
      Sig = .sig_code(ai[["Pr(>F)"]][2]),
      CD = stats::qt(1 - alpha / 2, ai[["Df"]][3]) * sed_i,
      stringsAsFactors = FALSE)
  }
  env_summary <- do.call(rbind, summ)
  mse_values <- env_summary$MSE
  df_values <- env_summary$df_error

  # Bartlett's test of homogeneity of error variances
  s2p <- sum(df_values * mse_values) / sum(df_values)
  chi <- (sum(df_values) * log(s2p) - sum(df_values * log(mse_values))) /
    (1 + (sum(1 / df_values) - 1 / sum(df_values)) / (3 * (e - 1)))
  p_bart <- stats::pchisq(chi, e - 1, lower.tail = FALSE)
  bartlett <- list(chi_sq = chi, df = e - 1, p_value = p_bart,
                   variance_ratio = max(mse_values) / min(mse_values),
                   homogeneous = p_bart >= 0.05)

  # Combined analysis
  rhs <- paste(.bq(environment), "+ .blk +", .bq(treatment), "+",
               paste0(.bq(environment), ":", .bq(treatment)))
  model <- stats::aov(.formula(response, rhs), data = d)
  at <- stats::anova(model)
  ss <- at[["Sum Sq"]]
  dfv <- at[["Df"]]
  ms <- at[["Mean Sq"]]
  # rows: env, blocks(env), trt, env:trt, residuals
  den <- if (env_effect == "fixed") c(2, 5, 5, 5) else c(2, 5, 4, 5)
  num <- c(1, 2, 3, 4)
  fv <- ms[num] / ms[den]
  pv <- stats::pf(fv, dfv[num], dfv[den], lower.tail = FALSE)
  pooled <- data.frame(Df = dfv, `Sum Sq` = ss, `Mean Sq` = ms,
                       `F value` = c(fv, NA), `Pr(>F)` = c(pv, NA),
                       check.names = FALSE)
  rownames(pooled) <- c(environment, paste0(block, "(", environment, ")"),
                        treatment, paste0(environment, ":", treatment),
                        "Residuals")
  ms_e <- ms[5]
  df_e <- dfv[5]
  ms_te <- ms[4]
  df_te <- dfv[4]
  ms_re <- ms[2]
  df_re <- dfv[2]
  gm <- mean(y)

  ms_t <- if (env_effect == "fixed") ms_e else ms_te
  df_t <- if (env_effect == "fixed") df_e else df_te
  se_t <- sqrt(ms_t / (e * r))
  cd_t <- stats::qt(1 - alpha / 2, df_t) * sqrt(2) * se_t
  se_env <- sqrt(ms_re / (k * r))
  cd_env <- stats::qt(1 - alpha / 2, df_re) * sqrt(2) * se_env
  se_int <- sqrt(ms_e / r)
  cd_int <- stats::qt(1 - alpha / 2, df_e) * sqrt(2) * se_int

  tm <- tapply(y, trt, mean)
  ph <- list()
  if (length(methods) > 0) {
    V <- diag(se_t^2, k)
    dimnames(V) <- list(names(tm), names(tm))
    for (m in methods) {
      ph[[m]] <- .compare(tm, V, df_t, ms_t, m, alpha,
                          n = table(trt), response = response)
    }
  }
  treatment_means <- data.frame(Treatment = names(tm), Mean = as.numeric(tm),
                                stringsAsFactors = FALSE)
  if (length(methods) > 0) {
    treatment_means$Group <- ph[[methods[1]]]$letters[treatment_means$Treatment]
  }
  treatment_means <- treatment_means[order(-treatment_means$Mean), ,
                                     drop = FALSE]
  rownames(treatment_means) <- NULL
  im <- tapply(y, list(trt, env), mean)
  interaction_table <- data.frame(Treatment = rownames(im),
                                  as.data.frame.matrix(im),
                                  Mean = as.numeric(tm[rownames(im)]),
                                  check.names = FALSE,
                                  stringsAsFactors = FALSE)
  rownames(interaction_table) <- NULL
  env_means <- data.frame(Environment = env_levels,
                          Mean = as.numeric(tapply(y, env, mean)[env_levels]),
                          stringsAsFactors = FALSE)

  result <- list(
    design = "Pooled RBD",
    response = response,
    factors = c(treatment = treatment, environment = environment,
                block = block),
    env_effect = env_effect,
    pooled_anova = pooled,
    anova_table = pooled,
    individual_anovas = individual,
    environment_summary = env_summary,
    model = model,
    grand_mean = gm,
    cv = 100 * sqrt(ms_e) / gm,
    bartlett_test = bartlett,
    ms_error = ms_e,
    df_error = df_e,
    se_treatment = se_t,
    cd_treatment = cd_t,
    se_environment = se_env,
    cd_environment = cd_env,
    se_interaction = se_int,
    cd_interaction = cd_int,
    alpha = alpha,
    treatment_means = treatment_means,
    env_means = env_means,
    interaction_means = interaction_table,
    letters_method = .letters_method(methods),
    posthoc = ph,
    mse_values = stats::setNames(mse_values, env_levels),
    balanced = balanced
  )
  class(result) <- c("aridagri_pooled", "list")
  if (verbose) print(result)
  invisible(result)
}


#' Print Method for Combined RCBD Analysis
#'
#' @param x An object of class \code{aridagri_pooled}.
#' @param digits Number of decimal places for means. Default 2.
#' @param ... Unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.aridagri_pooled <- function(x, digits = 2, ...) {
  f <- x$factors
  .print_info("COMBINED ANALYSIS OF RCBD OVER ENVIRONMENTS",
              c(Response = x$response, Treatment = f[["treatment"]],
                Environment = paste0(f[["environment"]], " (", x$env_effect,
                                     ")"),
                Treatments = nrow(x$treatment_means),
                Environments = nrow(x$environment_summary)))
  cat("\n Analysis of each environment\n")
  es <- x$environment_summary
  es$Mean <- round(es$Mean, digits)
  es$MSE <- signif(es$MSE, 5)
  es$CV <- round(es$CV, 2)
  es$F_value <- round(es$F_value, 3)
  es$p_value <- signif(es$p_value, 3)
  es$CD <- round(es$CD, 3)
  print(es, row.names = FALSE)
  b <- x$bartlett_test
  cat(sprintf(paste0("\n Bartlett's test of error variances: chi-square = ",
                     "%.3f, df = %d, p = %.4f (%s)\n"),
              b$chi_sq, b$df, b$p_value,
              if (b$homogeneous) "homogeneous, pooling is valid" else
                "heterogeneous, interpret pooling with care"))
  den_t <- if (x$env_effect == "fixed") "pooled error" else "T x E"
  labs <- c("Environment", "Blocks within env.", "Treatment",
            "Treatment x Environment", "Pooled error")
  .print_anova_df(x$pooled_anova, labs,
                  title = paste0("Combined analysis of variance (treatment ",
                                 "tested against ", den_t, ")"))
  cat(sprintf("\n Grand mean = %.3f   Pooled CV = %.2f%%\n", x$grand_mean,
              x$cv))
  a <- 100 * x$alpha
  cat(sprintf(" Treatment means   : SEm = %.3f   CD (%g%%) = %.3f\n",
              x$se_treatment, a, x$cd_treatment))
  cat(sprintf(" Environment means : SEm = %.3f   CD (%g%%) = %.3f\n",
              x$se_environment, a, x$cd_environment))
  cat(sprintf(" T within E        : SEm = %.3f   CD (%g%%) = %.3f\n",
              x$se_interaction, a, x$cd_interaction))
  ttl <- if (is.na(x$letters_method)) "Treatment means over environments" else
    sprintf("Treatment means over environments (%s)", x$letters_method)
  .print_means(x$treatment_means, ttl, digits)
  invisible(x)
}
