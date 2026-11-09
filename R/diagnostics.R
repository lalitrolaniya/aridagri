# =============================================================================
# ANOVA ASSUMPTION CHECKS
# =============================================================================

#' Check ANOVA Assumptions
#'
#' @description
#' Checks the assumptions behind an analysis of variance: normality of
#' residuals (Shapiro-Wilk), homogeneity of variances across treatments
#' (Bartlett and Levene/Brown-Forsythe), additivity of block and treatment
#' effects (Tukey's one-degree-of-freedom test, for two-way models),
#' serial independence (Durbin-Watson, in data order) and outliers
#' (studentized residuals).
#'
#' @param model A fitted \code{aov} or \code{lm} model.
#' @param verbose Logical. If TRUE (default), prints the results.
#' @param factor Name(s) of the factor(s) defining the groups for the
#'   variance tests. Default \code{NULL} uses the last main-effect factor
#'   in the model (the treatment in \code{y ~ block + treatment}). Give
#'   several names to use their combinations.
#' @param plot Logical. If TRUE, draws the four standard residual plots
#'   with base graphics. Default FALSE.
#'
#' @details
#' The variance tests use the model residuals, so block effects do not
#' inflate the group variances. Levene's test is the Brown-Forsythe form
#' (absolute deviations from the group median). Tukey's test for
#' non-additivity adds the squared fitted values as a one-degree-of-freedom
#' covariate to a model with two additive factors; a significant result
#' suggests a transformation (Gomez and Gomez, 1984, ch. 7). The
#' Durbin-Watson statistic is only meaningful when the rows are in field
#' or time order.
#'
#' @return A list with \code{shapiro}, \code{bartlett}, \code{levene},
#'   \code{nonadditivity} (where applicable), \code{durbin_watson},
#'   \code{outliers} (row indices with |studentized residual| > 3) and
#'   \code{groups_factor}.
#'
#' @references
#' Tukey, J.W. (1949). One degree of freedom for non-additivity.
#' \emph{Biometrics}, 5(3), 232-242. \doi{10.2307/3001938}
#'
#' Brown, M.B. and Forsythe, A.B. (1974). Robust tests for the equality of
#' variances. \emph{Journal of the American Statistical Association},
#' 69(346), 364-367. \doi{10.1080/01621459.1974.10482955}
#'
#' @examples
#' data(mothbean_rbd)
#' fit <- aov(seed_yield ~ block + treatment, data = mothbean_rbd)
#' check_assumptions(fit)
#'
#' @export
check_assumptions <- function(model, verbose = TRUE, factor = NULL,
                              plot = FALSE) {
  res <- stats::residuals(model)
  fit_v <- stats::fitted(model)
  n <- length(res)
  mf <- stats::model.frame(model)
  stud <- tryCatch(stats::rstandard(model), error = function(e) {
    s <- stats::sd(res)
    if (is.finite(s) && s > 0) res / s else res * 0
  })
  stud[!is.finite(stud)] <- 0
  out <- list()

  # Normality
  if (n >= 3 && n <= 5000 && stats::var(res) > 0) {
    out$shapiro <- stats::shapiro.test(res)
  } else {
    out$shapiro <- NULL
  }

  # Groups for the variance tests
  facs <- names(mf)[vapply(mf, function(v) is.factor(v) || is.character(v),
                           logical(1))]
  if (is.null(factor)) {
    labs <- attr(stats::terms(model), "term.labels")
    main <- labs[labs %in% facs]
    factor <- if (length(main)) main[length(main)] else
      if (length(facs)) facs[length(facs)] else NULL
  }
  if (!is.null(factor)) {
    miss <- setdiff(factor, names(mf))
    if (length(miss) > 0) {
      stop("Factor not in the model: ", paste(miss, collapse = ", "),
           call. = FALSE)
    }
    grp <- interaction(mf[factor], drop = TRUE)
    ok <- table(grp)
    if (nlevels(grp) >= 2 && all(ok >= 2)) {
      out$bartlett <- tryCatch(stats::bartlett.test(res, grp),
                               error = function(e) NULL)
      dev <- abs(res - stats::ave(res, grp, FUN = stats::median))
      lv <- stats::anova(stats::lm(dev ~ grp))
      out$levene <- list(statistic = lv[["F value"]][1],
                         df = c(lv[["Df"]][1], lv[["Df"]][2]),
                         p.value = lv[["Pr(>F)"]][1],
                         method = "Levene's test (Brown-Forsythe, median)")
    }
  }
  out$groups_factor <- factor

  # Tukey's one-degree-of-freedom test for non-additivity
  labs <- attr(stats::terms(model), "term.labels")
  if (length(labs) == 2 && all(labs %in% facs)) {
    yname <- names(mf)[1]
    d2 <- mf
    d2$.q <- fit_v^2
    f2 <- stats::as.formula(paste(.bq(yname), "~", .bq(labs[1]), "+",
                                  .bq(labs[2]), "+ .q"))
    a2 <- tryCatch(stats::anova(stats::lm(f2, data = d2)),
                   error = function(e) NULL)
    if (!is.null(a2) && ".q" %in% rownames(a2) &&
        a2[["Df"]][nrow(a2)] > 0) {
      out$nonadditivity <- list(F_value = a2[".q", "F value"],
                                df = c(1, a2[["Df"]][nrow(a2)]),
                                p.value = a2[".q", "Pr(>F)"])
    }
  }

  # Independence and outliers
  out$durbin_watson <- if (n > 1 && sum(res^2) > 0) {
    sum(diff(res)^2) / sum(res^2)
  } else NA_real_
  out$outliers <- which(abs(stud) > 3)

  if (verbose) {
    .heading("ANOVA ASSUMPTION CHECKS")
    if (!is.null(out$shapiro)) {
      cat(sprintf(" Normality (Shapiro-Wilk)       : W = %.4f, p = %.4f  %s\n",
                  out$shapiro$statistic, out$shapiro$p.value,
                  if (out$shapiro$p.value >= 0.05) "OK" else
                    "non-normal; consider a transformation"))
    } else {
      cat(" Normality (Shapiro-Wilk)       : not run (needs 3 to 5000",
          "residuals)\n")
    }
    if (!is.null(out$bartlett)) {
      cat(sprintf(" Equal variances (Bartlett)     : K2 = %.3f, p = %.4f  %s\n",
                  out$bartlett$statistic, out$bartlett$p.value,
                  if (out$bartlett$p.value >= 0.05) "OK" else "unequal"))
    }
    if (!is.null(out$levene)) {
      cat(sprintf(" Equal variances (Levene)       : F = %.3f, p = %.4f  %s\n",
                  out$levene$statistic, out$levene$p.value,
                  if (out$levene$p.value >= 0.05) "OK" else "unequal"))
    }
    if (!is.null(out$nonadditivity)) {
      cat(sprintf(" Additivity (Tukey 1 df)        : F = %.3f, p = %.4f  %s\n",
                  out$nonadditivity$F_value, out$nonadditivity$p.value,
                  if (out$nonadditivity$p.value >= 0.05) "OK" else
                    "non-additive; consider a transformation"))
    }
    cat(sprintf(" Durbin-Watson (data order)     : %.3f\n", out$durbin_watson))
    if (length(out$outliers) > 0) {
      cat(" Possible outliers (|studentized residual| > 3): rows",
          paste(out$outliers, collapse = ", "), "\n")
    } else {
      cat(" No studentized residual exceeds 3 in absolute value.\n")
    }
    if (!is.null(factor)) {
      cat(" Variance tests compare groups of:", paste(factor, collapse = " x "),
          "\n")
    }
  }

  if (plot) {
    op <- graphics::par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1))
    on.exit(graphics::par(op), add = TRUE)
    graphics::plot(fit_v, res, xlab = "Fitted values", ylab = "Residuals",
                   main = "Residuals vs fitted", pch = 19,
                   col = grDevices::adjustcolor("#2a78d6", 0.7))
    graphics::abline(h = 0, lty = 2, col = "grey40")
    stats::qqnorm(stud, main = "Normal Q-Q", pch = 19,
                  col = grDevices::adjustcolor("#2a78d6", 0.7))
    stats::qqline(stud, col = "grey40", lty = 2)
    graphics::hist(res, main = "Histogram of residuals", xlab = "Residuals",
                   col = "#9ec5f4", border = "white")
    graphics::plot(fit_v, sqrt(abs(stud)), xlab = "Fitted values",
                   ylab = "sqrt(|studentized residual|)",
                   main = "Scale-location", pch = 19,
                   col = grDevices::adjustcolor("#2a78d6", 0.7))
  }
  invisible(out)
}
