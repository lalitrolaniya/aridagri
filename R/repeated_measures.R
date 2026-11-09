# =============================================================================
# REPEATED MEASURES (SPLIT-PLOT IN TIME)
# Package: aridagri
# =============================================================================

#' Repeated Measures ANOVA (Split-Plot in Time)
#'
#' @description
#' Analysis of experiments in which the same experimental units (plots,
#' animals, trees) are measured on several occasions: plant height or LAI
#' at successive dates, soil moisture through the season, or body weight,
#' heart rate or blood metabolites over time. Gives the univariate
#' split-plot-in-time ANOVA, Mauchly's test of sphericity with
#' Greenhouse-Geisser and Huynh-Feldt corrected p-values, a comparison of
#' covariance structures by mixed models (compound symmetry, first-order
#' autoregressive and unstructured), SE and CD for every comparison, letter
#' groupings, and a separate analysis at each time.
#'
#' @param data Data frame in long format (one row per unit per time).
#' @param response Name of the response variable.
#' @param treatment Name of the treatment column.
#' @param time Name of the time (occasion) column. Times are taken in the
#'   order of their levels (numeric times are sorted).
#' @param subject Name of the column identifying the experimental unit that
#'   is measured repeatedly. If \code{NULL} and \code{block} is given, each
#'   block x treatment combination is one unit.
#' @param block Name of the block column for units laid out in randomized
#'   complete blocks, or \code{NULL} (completely randomized units).
#' @param posthoc Method for the letters (default \code{"lsd"}).
#' @param alpha Significance level. Default 0.05.
#' @param covariance Logical. If TRUE (default), fits mixed models with
#'   different covariance structures using \pkg{nlme}.
#' @param verbose Logical. If TRUE (default), prints the report.
#'
#' @details
#' \strong{Univariate analysis.} Treatment (and block) are tested against
#' the between-unit error, Error(a). Time and the treatment x time
#' interaction are tested against the within-unit error, Error(b). With
#' blocks, the block x time interaction is fitted as a separate
#' within-unit term, as in the univariate tests of SAS PROC GLM with a
#' REPEATED statement, so that the sphericity corrections apply to the
#' same error term. This is
#' valid when the covariance of the repeated measurements is spherical.
#' Mauchly's test checks this, and the Greenhouse-Geisser (1959) and
#' Huynh-Feldt (1976) epsilons, the latter with the correction of Lecoutre
#' (1991) as in \code{anova.mlm()}, scale the degrees of freedom of the
#' within-unit tests when it fails. A common rule is to use the
#' Greenhouse-Geisser p-value when epsilon is below 0.75 and the
#' Huynh-Feldt p-value otherwise.
#'
#' \strong{Mixed models.} Generalized least squares models with the same
#' fixed effects are fitted by REML with compound symmetry (CS),
#' first-order autoregressive (AR1, on the time order) and unstructured
#' (UN, with separate variances per time) covariance, and compared by AIC
#' and BIC. The marginal F-tests of the best structure are reported (with
#' the residual degrees of freedom used by \pkg{nlme}, so tests of
#' between-unit effects are approximate; the univariate test of treatment
#' is exact for balanced data). Structures that fail to converge are
#' skipped.
#'
#' \strong{Comparisons.} With \eqn{r} units per treatment, \eqn{a}
#' treatments and \eqn{t} times, the SEd are as in a split-plot design
#' with treatments in main plots and times in sub-plots: treatment means
#' \eqn{\sqrt{2E_a/(rt)}}, time means \eqn{\sqrt{2E_b/(ra)}}, times within
#' a treatment \eqn{\sqrt{2E_b/r}}, and treatments at the same time
#' \eqn{\sqrt{2[(t-1)E_b + E_a]/(rt)}} with the weighted t.
#'
#' @return An object of class \code{aridagri_repeated} with
#'   \code{anova_table} (univariate), \code{sphericity} (Mauchly's W and
#'   p, GG and HF epsilons), \code{corrected} (within-unit tests with GG
#'   and HF p-values), \code{covariance_models} (AIC/BIC table),
#'   \code{best_structure}, \code{mixed_anova}, \code{se_cd},
#'   \code{means} (treatment x time table), \code{treatment_means},
#'   \code{time_means}, \code{letters}, \code{by_time} (ANOVA at each time)
#'   and \code{models}.
#'
#' @references
#' Greenhouse, S.W. and Geisser, S. (1959). On methods in the analysis of
#' profile data. \emph{Psychometrika}, 24(2), 95-112.
#' \doi{10.1007/BF02289823}
#'
#' Lecoutre, B. (1991). A correction for the epsilon approximate test in
#' repeated measures designs with two or more independent groups.
#' \emph{Journal of Educational Statistics}, 16(4), 371-372.
#' \doi{10.3102/10769986016004371}
#'
#' Huynh, H. and Feldt, L.S. (1976). Estimation of the Box correction for
#' degrees of freedom from sample data in randomized block and split-plot
#' designs. \emph{Journal of Educational Statistics}, 1(1), 69-82.
#' \doi{10.3102/10769986001001069}
#'
#' Littell, R.C., Henry, P.R. and Ammerman, C.B. (1998). Statistical
#' analysis of repeated measures data using SAS procedures. \emph{Journal
#' of Animal Science}, 76(4), 1216-1231. \doi{10.2527/1998.7641216x}
#'
#' @examples
#' data(mungbean_repeated)
#' rm <- anova_repeated(mungbean_repeated, response = "spad",
#'                      treatment = "treatment", time = "das", block = "block")
#' rm$sphericity
#' rm$covariance_models
#'
#' @export
anova_repeated <- function(data, response, treatment, time, subject = NULL,
                           block = NULL, posthoc = "lsd", alpha = 0.05,
                           covariance = TRUE, verbose = TRUE) {
  .check_alpha(alpha)
  methods <- setdiff(.resolve_methods(posthoc), "dunnett")
  if (is.null(subject) && is.null(block)) {
    stop("Give 'subject' (the unit measured repeatedly) or 'block'.",
         call. = FALSE)
  }
  .check_cols(data, c(response, treatment, time, subject, block))
  d0 <- as.data.frame(data, stringsAsFactors = FALSE)
  if (is.numeric(d0[[time]])) {
    d0[[time]] <- factor(d0[[time]], levels = sort(unique(d0[[time]])))
  }
  if (is.null(subject)) {
    d0$.unit <- interaction(d0[[block]], d0[[treatment]], drop = TRUE)
    subject <- ".unit"
  }
  d <- .prep_data(d0, response, unique(c(treatment, time, subject, block)),
                  verbose)
  trt <- d[[treatment]]
  tm <- d[[time]]
  su <- d[[subject]]
  a <- nlevels(trt)
  k <- nlevels(tm)
  if (a < 2) stop("Need at least 2 treatments.", call. = FALSE)
  if (k < 2) stop("Need at least 2 times.", call. = FALSE)
  # each unit belongs to one treatment (and block)
  ut <- tapply(as.character(trt), su, function(v) length(unique(v)))
  if (any(ut > 1)) {
    stop("Each subject must receive a single treatment.", call. = FALSE)
  }
  cnt <- table(su, tm)
  complete <- all(cnt == 1)
  n_unit <- nlevels(su)
  r_t <- table(tapply(as.character(trt), su, `[`, 1))
  balanced <- complete && length(unique(as.numeric(r_t))) == 1
  if (!complete) {
    warning("Not every unit is measured exactly once at every time. The ",
            "univariate analysis and sphericity tests need complete data ",
            "and are skipped; only the mixed models are fitted.",
            call. = FALSE)
  }
  res <- list(design = "Repeated measures (split-plot in time)",
              response = response, letters = list(),
              factors = c(treatment = treatment, time = time,
                          subject = subject, block = block %||% NA),
              alpha = alpha, letters_method = .letters_method(methods),
              n_units = n_unit, n_times = k, complete = complete)
  y <- d[[response]]

  if (complete) {
    btw <- if (is.null(block)) .bq(treatment) else
      paste(.bq(block), "+", .bq(treatment))
    bxt <- if (is.null(block)) "" else
      paste0(" + ", .bq(block), ":", .bq(time))
    f <- .formula(response, paste0(btw, " + ", .bq(time), bxt, " + ",
                                   .bq(treatment), ":", .bq(time),
                                   " + Error(", .bq(subject), ")"))
    fit <- suppressWarnings(stats::aov(f, data = d))
    s <- summary(fit)
    st_names <- names(s)
    t_btw <- s[[which(st_names == paste0("Error: ", subject))]][[1]]
    t_win <- s[[which(st_names == "Error: Within")]][[1]]
    rn_b <- trimws(gsub("`", "", rownames(t_btw)))
    rn_w <- trimws(gsub("`", "", rownames(t_win)))
    rn_b[rn_b == "Residuals"] <- "Error(a)"
    rn_w[rn_w == "Residuals"] <- "Error(b)"
    tab <- rbind(as.data.frame(t_btw), as.data.frame(t_win))
    rownames(tab) <- make.unique(c(rn_b, rn_w))
    names(tab) <- c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
    Ea <- tab["Error(a)", "Mean Sq"]
    dfa <- tab["Error(a)", "Df"]
    Eb <- tab["Error(b)", "Mean Sq"]
    dfb <- tab["Error(b)", "Df"]
    res$anova_table <- tab
    res$model <- fit
    res$cv <- c(a = 100 * sqrt(Ea) / mean(y), b = 100 * sqrt(Eb) / mean(y))
    res$grand_mean <- mean(y)

    # Sphericity: mlm on the wide data
    wide <- stats::reshape(d[, c(subject, treatment, block, time, response)],
                           idvar = subject, timevar = time,
                           v.names = response, direction = "wide")
    ycols <- paste0(response, ".", levels(tm))
    Y <- as.matrix(wide[, ycols])
    wide[[treatment]] <- factor(wide[[treatment]])
    rhs <- if (is.null(block)) .bq(treatment) else
      paste(.bq(block), "+", .bq(treatment))
    mlm <- stats::lm(stats::as.formula(paste("Y ~", rhs)), data = wide)
    n_res <- mlm$df.residual
    C <- stats::contr.poly(k)
    R <- stats::residuals(mlm)
    B <- t(C) %*% crossprod(R) %*% C
    lam <- Re(eigen(B, only.values = TRUE)$values)
    p <- k - 1
    gg <- sum(lam)^2 / (p * sum(lam^2))
    hf <- min(1, ((n_res + 1) * p * gg - 2) / (p * (n_res - p * gg)))
    mt <- tryCatch(stats::mauchly.test(mlm, X = ~1), error = function(e) NULL)
    res$sphericity <- data.frame(
      Mauchly_W = if (is.null(mt)) NA_real_ else unname(mt$statistic),
      p_value = if (is.null(mt)) NA_real_ else mt$p.value,
      GG_epsilon = gg, HF_epsilon = hf, lower_bound = 1 / p)
    win <- tab[rownames(tab) %in% c(time, paste0(treatment, ":", time)), ,
               drop = FALSE]
    win <- win[!is.na(win[["F value"]]), , drop = FALSE]
    corr <- data.frame(Source = rownames(win), Df = win$Df,
                       F_value = win[["F value"]],
                       p_uncorrected = win[["Pr(>F)"]],
                       p_GG = stats::pf(win[["F value"]], win$Df * gg,
                                        dfb * gg, lower.tail = FALSE),
                       p_HF = stats::pf(win[["F value"]], win$Df * hf,
                                        dfb * hf, lower.tail = FALSE),
                       stringsAsFactors = FALSE)
    corr$recommended <- if (gg < 0.75) corr$p_GG else corr$p_HF
    res$corrected <- corr

    # SE / CD and letters (balanced)
    if (balanced) {
      r <- as.numeric(r_t[1])
      facs <- c(treatment, time)
      lev_of <- stats::setNames(c(1L, 2L), facs)
      nl <- stats::setNames(c(a, k), facs)
      sc <- .split_comparisons(facs, lev_of, nl, r, c(Ea, Eb), c(dfa, dfb),
                               alpha)
      sc$Error <- sub("Eb", "E(b)", sub("Ea", "E(a)", sc$Error))
      res$se_cd <- sc
      lets <- list()
      if (length(methods) > 0) {
        c1 <- sc[sc$Effect == treatment, ]
        c2 <- sc[sc$Effect == time, ]
        c3 <- sc[sc$Effect == paste(treatment, "x", time) & sc$stratum == 1, ]
        c4 <- sc[sc$Effect == paste(treatment, "x", time) & sc$stratum == 2, ]
        lets$treatment <- .split_letters(d, response, treatment, NULL,
                                         c1$SEd, c1$df, c1$t_value, methods,
                                         alpha, Ea)
        lets$time <- .split_letters(d, response, time, NULL, c2$SEd, c2$df,
                                    c2$t_value, methods, alpha, Eb)
        lets$treatment_within_time <- .split_letters(
          d, response, treatment, time, c3$SEd, c3$df, c3$t_value, methods,
          alpha, Eb)
        lets$time_within_treatment <- .split_letters(
          d, response, time, treatment, c4$SEd, c4$df, c4$t_value, methods,
          alpha, Eb)
      }
      res$letters <- lets
    }
  }
  tmn <- tapply(y, trt, mean)
  res$treatment_means <- data.frame(Treatment = names(tmn),
                                    Mean = as.numeric(tmn),
                                    stringsAsFactors = FALSE)
  if (!is.null(res[["letters"]][["treatment"]])) {
    res$treatment_means$Group <- res[["letters"]][["treatment"]]$Group
  }
  tim <- tapply(y, tm, mean)
  res$time_means <- data.frame(Time = names(tim), Mean = as.numeric(tim),
                               stringsAsFactors = FALSE)
  if (!is.null(res[["letters"]][["time"]])) {
    res$time_means$Group <- res[["letters"]][["time"]]$Group
  }
  res$means <- tapply(y, list(trt, tm), mean)
  res$se_table <- tapply(y, list(trt, tm), function(v) {
    if (length(v) > 1) stats::sd(v) / sqrt(length(v)) else NA_real_
  })

  # Mixed models with different covariance structures
  if (covariance && requireNamespace("nlme", quietly = TRUE)) {
    dm <- d
    dm$.y <- y
    dm$.trt <- trt
    dm$.time <- tm
    dm$.tix <- as.integer(tm)
    dm$.sub <- su
    # sum-to-zero contrasts so that the marginal (type III) tests of main
    # effects are averaged over the levels of the other factor
    stats::contrasts(dm$.trt) <- stats::contr.sum(nlevels(dm$.trt))
    stats::contrasts(dm$.time) <- stats::contr.sum(nlevels(dm$.time))
    if (!is.null(block)) {
      dm[[block]] <- factor(dm[[block]])
      stats::contrasts(dm[[block]]) <- stats::contr.sum(nlevels(dm[[block]]))
    }
    fx <- if (is.null(block)) ".y ~ .trt * .time" else
      paste(".y ~", .bq(block), "+ .trt * .time")
    fx <- stats::as.formula(fx)
    fits <- list(
      CS = tryCatch(nlme::gls(fx, data = dm, method = "REML",
                              correlation = nlme::corCompSymm(
                                form = ~ .tix | .sub)),
                    error = function(e) NULL),
      AR1 = tryCatch(nlme::gls(fx, data = dm, method = "REML",
                               correlation = nlme::corAR1(
                                 form = ~ .tix | .sub)),
                     error = function(e) NULL),
      UN = if (k <= 8) tryCatch(nlme::gls(fx, data = dm, method = "REML",
                                          correlation = nlme::corSymm(
                                            form = ~ .tix | .sub),
                                          weights = nlme::varIdent(
                                            form = ~ 1 | .time)),
                                error = function(e) NULL) else NULL)
    fits <- fits[!vapply(fits, is.null, logical(1))]
    if (length(fits) > 0) {
      cm <- data.frame(Structure = names(fits),
                       Parameters = vapply(fits, function(m) attr(stats::logLik(m), "df"), numeric(1)),
                       AIC = vapply(fits, stats::AIC, numeric(1)),
                       BIC = vapply(fits, stats::BIC, numeric(1)),
                       logLik = vapply(fits, function(m) as.numeric(stats::logLik(m)), numeric(1)),
                       stringsAsFactors = FALSE)
      rownames(cm) <- NULL
      best <- cm$Structure[which.min(cm$AIC)]
      an <- as.data.frame(stats::anova(fits[[best]], type = "marginal"))
      an <- an[rownames(an) != "(Intercept)", , drop = FALSE]
      rn <- rownames(an)
      rn <- sub("^\\.trt:\\.time$", paste0(treatment, ":", time), rn)
      rn <- sub("^\\.trt$", treatment, rn)
      rn <- sub("^\\.time$", time, rn)
      rownames(an) <- rn
      res$covariance_models <- cm
      res$best_structure <- best
      res$mixed_anova <- an
      res$models <- fits
    }
  }

  # Separate analysis at each time
  bt <- list()
  for (tl in levels(tm)) {
    di <- d[tm == tl, , drop = FALSE]
    ai <- tryCatch(
      if (is.null(block)) {
        anova_crd(di, response, treatment, posthoc = posthoc, alpha = alpha,
                  verbose = FALSE)
      } else {
        anova_rbd(di, response, treatment, block, posthoc = posthoc,
                  alpha = alpha, verbose = FALSE)
      }, error = function(e) NULL)
    if (is.null(ai)) next
    p_i <- ai$anova_table[["Pr(>F)"]][if (is.null(block)) 1 else 2]
    bt[[tl]] <- list(summary = data.frame(
      Time = tl, Mean = ai$grand_mean, F_value =
        ai$anova_table[["F value"]][if (is.null(block)) 1 else 2],
      p_value = p_i, Sig = .sig_code(p_i), CV = ai$cv, SEm = ai$se_mean,
      CD = ai$cd, stringsAsFactors = FALSE),
      means = ai$treatment_means)
  }
  if (length(bt) > 0) {
    res$by_time <- do.call(rbind, lapply(bt, `[[`, "summary"))
    rownames(res$by_time) <- NULL
    res$by_time_means <- lapply(bt, `[[`, "means")
  }
  class(res) <- c("aridagri_repeated", "list")
  if (verbose) print(res)
  invisible(res)
}


#' Print Method for Repeated Measures Results
#'
#' @param x An object of class \code{aridagri_repeated}.
#' @param digits Number of decimal places for means. Default 2.
#' @param ... Unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.aridagri_repeated <- function(x, digits = 2, ...) {
  f <- x$factors
  .print_info("REPEATED MEASURES ANOVA (SPLIT-PLOT IN TIME)",
              c(Response = x$response, Treatment = f[["treatment"]],
                Time = f[["time"]], Units = x$n_units, Times = x$n_times))
  if (!is.null(x$anova_table)) {
    .print_anova_df(x$anova_table, rownames(x$anova_table),
                    title = "Univariate analysis of variance")
    cat(sprintf("\n Grand mean = %.3f   CV(a) = %.2f%%   CV(b) = %.2f%%\n",
                x$grand_mean, x$cv[1], x$cv[2]))
    sp <- x$sphericity
    cat(sprintf(paste0("\n Mauchly's test of sphericity: W = %.4f, p = %.4f\n",
                       " Greenhouse-Geisser epsilon = %.4f   Huynh-Feldt ",
                       "epsilon = %.4f\n"),
                sp$Mauchly_W, sp$p_value, sp$GG_epsilon, sp$HF_epsilon))
    co <- x$corrected
    for (cl in c("F_value")) co[[cl]] <- round(co[[cl]], 3)
    for (cl in c("p_uncorrected", "p_GG", "p_HF", "recommended")) {
      co[[cl]] <- signif(co[[cl]], 4)
    }
    cat("\n Within-unit tests with corrected p-values\n")
    print(co, row.names = FALSE)
  }
  if (!is.null(x$covariance_models)) {
    cat("\n Mixed models: covariance structures (lower AIC is better)\n")
    cm <- x$covariance_models
    cm[, 3:5] <- round(cm[, 3:5], 2)
    print(cm, row.names = FALSE)
    cat("\n F-tests with the", x$best_structure, "structure (residual df;",
        "for treatment, the univariate test above is exact)\n")
    ma <- x$mixed_anova
    ma[["F-value"]] <- round(ma[["F-value"]], 3)
    ma[["p-value"]] <- signif(ma[["p-value"]], 4)
    print(ma)
  }
  if (!is.null(x$se_cd)) {
    cat("\n SEm, SEd and CD (alpha = ", x$alpha, ")\n", sep = "")
    sc <- x$se_cd[, c("Comparison", "SEm", "SEd", "CD", "Error")]
    sc[, 2:4] <- round(sc[, 2:4], 3)
    print(sc, row.names = FALSE, right = FALSE)
  }
  cat("\n Treatment x time means\n")
  print(round(x$means, digits))
  if (!is.null(x$by_time)) {
    cat("\n Analysis at each time\n")
    bt <- x$by_time
    bt$Mean <- round(bt$Mean, digits)
    bt$F_value <- round(bt$F_value, 3)
    bt$p_value <- signif(bt$p_value, 3)
    bt$CV <- round(bt$CV, 2)
    bt$SEm <- round(bt$SEm, 3)
    bt$CD <- round(bt$CD, 3)
    print(bt, row.names = FALSE)
  }
  invisible(x)
}
