# =============================================================================
# MEAN COMPARISON ENGINE
# Exact pairwise tests and compact letter displays for LSD, Tukey HSD,
# Duncan (DMRT), Student-Newman-Keuls, Scheffe and Bonferroni. Works with
# equal or unequal replication and with least-squares means from models
# with missing plots. No extra packages are needed; Dunnett's test uses
# 'multcomp' when it is installed.
# =============================================================================

.mc_methods <- c("lsd", "duncan", "tukey", "snk", "scheffe", "bonferroni")

# Resolve the 'posthoc' argument to a vector of method names.
.resolve_methods <- function(posthoc) {
  if (is.null(posthoc) || identical(posthoc, "none") || length(posthoc) == 0) {
    return(character(0))
  }
  posthoc <- tolower(posthoc)
  if ("all" %in% posthoc) return(c(.mc_methods, "dunnett"))
  bad <- setdiff(posthoc, c(.mc_methods, "dunnett", "hsd", "dmrt"))
  if (length(bad) > 0) {
    stop("Unknown post-hoc method: ", paste(bad, collapse = ", "),
         ". Use one or more of: ", paste(c(.mc_methods, "dunnett"),
                                         collapse = ", "), ", or \"all\".",
         call. = FALSE)
  }
  posthoc[posthoc == "hsd"] <- "tukey"
  posthoc[posthoc == "dmrt"] <- "duncan"
  unique(posthoc)
}

# Pairwise significance for one method.
#   means : named numeric vector
#   vcov  : covariance matrix of the means
#   df    : error degrees of freedom
#   t_crit: optional critical t replacing qt(1 - alpha/2, df) for LSD
# Returns a list with the logical matrix 'sig', a long table of pairwise
# comparisons and the critical difference(s).
.mc_pairs <- function(means, vcov, df, method, alpha, t_crit = NULL) {
  k <- length(means)
  lev <- names(means)
  v <- as.matrix(vcov)
  sed <- sqrt(pmax(outer(diag(v), diag(v), "+") - 2 * v, 0))
  dif <- outer(means, means, "-")
  m <- k * (k - 1) / 2
  ord <- order(means, decreasing = TRUE)
  pos <- integer(k)
  pos[ord] <- seq_len(k)
  span <- abs(outer(pos, pos, "-")) + 1
  se_q <- sed / sqrt(2)
  q <- abs(dif) / se_q
  tval <- abs(dif) / sed
  tc <- if (is.null(t_crit)) stats::qt(1 - alpha / 2, df) else t_crit
  p <- switch(method,
    lsd = 2 * stats::pt(-tval, df),
    bonferroni = pmin(1, 2 * stats::pt(-tval, df) * m),
    tukey = stats::ptukey(q, k, df, lower.tail = FALSE),
    scheffe = stats::pf(dif^2 / ((k - 1) * sed^2), k - 1, df,
                        lower.tail = FALSE),
    snk = stats::ptukey(q, pmax(span, 2), df, lower.tail = FALSE),
    duncan = {
      pr <- stats::ptukey(q, pmax(span, 2), df, lower.tail = FALSE)
      1 - (1 - pr)^(1 / pmax(span - 1, 1))
    })
  p <- matrix(p, k, k, dimnames = list(lev, lev))
  diag(p) <- 1
  p[!is.finite(p)] <- 1
  sig <- p < alpha
  if (method == "lsd" && !is.null(t_crit)) sig <- tval > t_crit
  diag(sig) <- FALSE
  sig[is.na(sig)] <- FALSE
  if (method %in% c("snk", "duncan") && k > 2) {
    # Closure: two means differ only if every range of means containing them
    # differs. The reported p-value is the largest p of those ranges, so
    # p < alpha agrees with the letters.
    p_ord <- p[ord, ord]
    adj <- matrix(1, k, k)
    for (w in (k - 1):1) {
      for (a in seq_len(k - w)) {
        b <- a + w
        v <- p_ord[a, b]
        if (a > 1) v <- max(v, adj[a - 1, b])
        if (b < k) v <- max(v, adj[a, b + 1])
        adj[a, b] <- adj[b, a] <- v
      }
    }
    p[ord, ord] <- adj
    diag(p) <- 1
    sig <- p < alpha
    diag(sig) <- FALSE
  }
  # Critical differences (exact when every pair has the same SEd).
  sed_u <- sed[upper.tri(sed)]
  equal <- length(sed_u) > 0 && diff(range(sed_u)) <= 1e-8 * max(sed_u)
  s0 <- if (length(sed_u)) mean(sed_u) else NA_real_
  crit <- switch(method,
    lsd = tc * s0,
    bonferroni = stats::qt(1 - alpha / (2 * m), df) * s0,
    tukey = stats::qtukey(1 - alpha, k, df) * s0 / sqrt(2),
    scheffe = sqrt((k - 1) * stats::qf(1 - alpha, k - 1, df)) * s0,
    snk = stats::qtukey(1 - alpha, 2:k, df) * s0 / sqrt(2),
    duncan = stats::qtukey((1 - alpha)^(1:(k - 1)), 2:k, df) * s0 / sqrt(2))
  if (method %in% c("snk", "duncan")) names(crit) <- paste0("p=", 2:k)
  ij <- which(upper.tri(sig), arr.ind = TRUE)
  ij <- ij[order(ij[, 1], ij[, 2]), , drop = FALSE]
  comp <- data.frame(
    comparison = paste(lev[ij[, 1]], "-", lev[ij[, 2]]),
    difference = dif[ij],
    SE_diff = sed[ij],
    p_value = p[ij],
    significant = sig[ij],
    stringsAsFactors = FALSE)
  list(sig = sig, comparisons = comp, critical = crit, equal_se = equal,
       sed = s0)
}

# Compact letter display from a logical matrix of significant pairs
# (insert-and-absorb, Piepho 2004). The group holding the highest mean
# gets the first letter.
.cld <- function(means, sig) {
  k <- length(means)
  cols <- matrix(TRUE, nrow = k, ncol = 1)
  if (k > 1) {
    for (i in seq_len(k - 1)) {
      for (j in (i + 1):k) {
        if (!isTRUE(sig[i, j]) && !isTRUE(sig[j, i])) next
        hit <- which(cols[i, ] & cols[j, ])
        if (length(hit) == 0) next
        add <- cols[, hit, drop = FALSE]
        cols[i, hit] <- FALSE
        add[j, ] <- FALSE
        cols <- .absorb(cbind(cols, add))
      }
    }
  }
  top <- apply(cols, 2, function(cl) max(means[cl]))
  avg <- apply(cols, 2, function(cl) mean(means[cl]))
  cols <- cols[, order(-top, -avg), drop = FALSE]
  labs <- .cld_labels(ncol(cols))
  out <- apply(cols, 1, function(rw) paste(labs[rw], collapse = ""))
  names(out) <- names(means)
  out
}

# Drop columns that are subsets of (or duplicate) another column.
.absorb <- function(cols) {
  n <- ncol(cols)
  if (n < 2) return(cols)
  m <- cols * 1
  out_of <- crossprod(m, 1 - m)
  sub <- out_of == 0
  diag(sub) <- FALSE
  drop <- rep(FALSE, n)
  for (a in seq_len(n)) {
    for (b in which(sub[a, ] & !drop)) {
      if (!sub[b, a] || b < a) {
        drop[a] <- TRUE
        break
      }
    }
  }
  cols[, !drop, drop = FALSE]
}

.cld_labels <- function(n) {
  base <- c(letters, LETTERS)
  if (n <= length(base)) return(base[seq_len(n)])
  extra <- paste0(letters, rep(seq_len(ceiling((n - 52) / 26)), each = 26))
  c(base, extra)[seq_len(n)]
}

# Full comparison for one method: letters, agricolae-style group table,
# pairwise table and statistics.
.compare <- function(means, vcov, df, mse, method, alpha, n = NULL,
                     response = "mean", t_crit = NULL) {
  pr <- .mc_pairs(means, vcov, df, method, alpha, t_crit = t_crit)
  let <- .cld(means, pr$sig)
  ord <- order(means, decreasing = TRUE)
  groups <- data.frame(means[ord], groups = let[ord],
                       row.names = names(means)[ord],
                       stringsAsFactors = FALSE)
  names(groups)[1] <- response
  se_m <- sqrt(pmax(diag(as.matrix(vcov)), 0))
  tab <- data.frame(level = names(means)[ord], mean = as.numeric(means[ord]),
                    SE = se_m[ord], group = let[ord],
                    stringsAsFactors = FALSE)
  if (!is.null(n)) tab$n <- as.numeric(n[names(means)[ord]])
  rownames(tab) <- NULL
  label <- switch(method, lsd = "Fisher's LSD", duncan = "Duncan's MRT",
                  tukey = "Tukey HSD", snk = "Student-Newman-Keuls",
                  scheffe = "Scheffe", bonferroni = "Bonferroni")
  stats_df <- data.frame(test = label, alpha = alpha, df_error = df,
                         MSE = mse, SE_diff = pr$sed,
                         equal_SE = pr$equal_se, stringsAsFactors = FALSE)
  if (length(pr$critical) == 1) stats_df$critical_difference <- pr$critical
  list(method = method, groups = groups, table = tab,
       letters = let, comparisons = pr$comparisons,
       critical_range = pr$critical, statistics = stats_df)
}

# Least-squares means of one factor, or of a combination of factors, from a
# fitted lm/aov, averaging over the levels of the other factors (numeric
# covariates held at their mean). Means that are not estimable (for example
# because of empty cells) are returned as NA and flagged in 'estimable'.
.ls_means_model <- function(fit, factor_name) {
  mf <- stats::model.frame(fit)
  lmfit <- if (inherits(fit, "lm") && !inherits(fit, "aov")) fit else
    stats::lm(stats::formula(fit), data = mf)
  miss <- setdiff(factor_name, names(mf))
  if (length(miss) > 0) {
    stop("Factor '", paste(miss, collapse = ", "), "' is not in the model.",
         call. = FALSE)
  }
  g <- interaction(mf[factor_name], drop = TRUE, sep = " : ")
  lev <- levels(g)
  n <- as.numeric(table(g))
  names(n) <- lev
  df <- lmfit$df.residual
  if (df < 1) {
    stop("No error degrees of freedom left in the model.", call. = FALSE)
  }
  mse <- sum(stats::residuals(lmfit)^2) / df
  b <- stats::coef(lmfit)
  keep <- !is.na(b)
  preds <- setdiff(names(mf), names(mf)[1])
  grid <- lapply(preds, function(p) {
    x <- mf[[p]]
    if (is.factor(x) || is.character(x)) levels(factor(x)) else mean(x)
  })
  names(grid) <- preds
  nd <- expand.grid(grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = TRUE)
  ndg <- as.character(interaction(nd[factor_name], drop = FALSE, sep = " : "))
  tt <- stats::delete.response(stats::terms(lmfit))
  X <- stats::model.matrix(tt, nd, xlev = lmfit$xlevels,
                           contrasts.arg = lmfit$contrasts)
  L <- t(vapply(lev, function(l) colMeans(X[ndg == l, , drop = FALSE]),
                numeric(ncol(X))))
  if (ncol(X) == 1) L <- matrix(L, ncol = 1)
  est <- rep(TRUE, length(lev))
  if (!all(keep)) {
    Xm <- stats::model.matrix(lmfit)
    sv <- svd(Xm)
    N <- sv$v[, sv$d < max(sv$d) * 1e-8, drop = FALSE]
    if (ncol(N) > 0) {
      est <- apply(abs(L %*% N), 1, max) < 1e-6 * max(1, max(abs(L)))
    }
  }
  Lk <- L[, keep, drop = FALSE]
  means <- as.numeric(Lk %*% b[keep])
  names(means) <- lev
  V <- Lk %*% stats::vcov(lmfit, complete = FALSE) %*% t(Lk)
  dimnames(V) <- list(lev, lev)
  means[!est] <- NA_real_
  list(means = means, vcov = V, n = n, mse = mse, df = df,
       exact = all(est), estimable = stats::setNames(est, lev))
}

# Keep only the estimable least-squares means (with a warning).
.drop_nonestimable <- function(lsm, what = "means") {
  if (all(lsm$estimable)) return(lsm)
  bad <- names(lsm$estimable)[!lsm$estimable]
  warning("Least-squares ", what, " are not estimable (empty cells) for: ",
          paste(bad, collapse = ", "), ". They are left out of the ",
          "comparisons.", call. = FALSE)
  ok <- lsm$estimable
  lsm$means <- lsm$means[ok]
  lsm$vcov <- lsm$vcov[ok, ok, drop = FALSE]
  lsm$n <- lsm$n[ok]
  lsm$estimable <- lsm$estimable[ok]
  lsm
}

# Dunnett's test against a control level (needs 'multcomp'). Differences
# are averaged over interacting factors; 'mse' and 'df' replace the model's
# error term when given.
.dunnett <- function(fit, factor_name, control = NULL, alpha = 0.05,
                     mse = NULL, df = NULL) {
  if (!requireNamespace("multcomp", quietly = TRUE)) {
    return(NULL)
  }
  mf <- stats::model.frame(fit)
  if (!is.null(control)) {
    if (!control %in% levels(factor(mf[[factor_name]]))) {
      stop("'control' level '", control, "' not found in '", factor_name,
           "'.", call. = FALSE)
    }
    mf[[factor_name]] <- stats::relevel(factor(mf[[factor_name]]),
                                        ref = control)
  }
  refit <- stats::lm(stats::formula(fit), data = mf)
  mse_fit <- sum(stats::residuals(refit)^2) / refit$df.residual
  scale <- if (is.null(mse)) 1 else mse / mse_fit
  df_use <- if (is.null(df)) refit$df.residual else df
  mcp_arg <- c(stats::setNames(list("Dunnett"), factor_name),
               list(interaction_average = TRUE, covariate_average = TRUE))
  g <- multcomp::glht(refit, linfct = do.call(multcomp::mcp, mcp_arg),
                      vcov. = function(m) stats::vcov(m) * scale,
                      df = df_use)
  s <- summary(g)
  data.frame(comparison = names(s$test$coefficients),
             estimate = as.numeric(s$test$coefficients),
             SE = as.numeric(s$test$sigma),
             t_value = as.numeric(s$test$tstat),
             p_adjusted = as.numeric(s$test$pvalues),
             significant = as.numeric(s$test$pvalues) < alpha,
             stringsAsFactors = FALSE)
}

# Run the requested comparisons for one factor of a fitted model.
# Returns a named list (one element per method) for the result object.
.posthoc_set <- function(fit, factor_name, methods, alpha, mse = NULL,
                         df = NULL, control = NULL, response = "mean",
                         lsm = NULL) {
  if (length(methods) == 0) return(list())
  if (is.null(lsm)) lsm <- .ls_means_model(fit, factor_name)
  lsm <- .drop_nonestimable(lsm)
  V <- lsm$vcov
  if (!is.null(mse) && is.finite(mse) && lsm$mse > 0) V <- V * mse / lsm$mse
  mse_use <- if (is.null(mse)) lsm$mse else mse
  df_use <- if (is.null(df)) lsm$df else df
  out <- list()
  for (m in methods) {
    if (m == "dunnett") {
      dn <- .dunnett(fit, factor_name, control = control, alpha = alpha,
                     mse = mse, df = df)
      if (!is.null(dn)) out$dunnett <- dn
      next
    }
    out[[m]] <- .compare(lsm$means, V, df_use, mse_use, m, alpha,
                         n = lsm$n, response = response)
  }
  out
}

# Print the result of .posthoc_set().
.print_posthoc <- function(ph, digits = 2) {
  for (m in names(ph)) {
    x <- ph[[m]]
    if (m == "dunnett") {
      cat("\nDunnett's test (each treatment vs control)\n")
      d <- x
      d$estimate <- round(d$estimate, digits)
      d$SE <- round(d$SE, 3)
      d$t_value <- round(d$t_value, 3)
      d$p_adjusted <- signif(d$p_adjusted, 4)
      print(d, row.names = FALSE)
      next
    }
    st <- x$statistics
    cat("\n", st$test, " (alpha = ", st$alpha, ", error df = ",
        round(st$df_error, 2), ")\n", sep = "")
    if (!is.null(st$critical_difference) && isTRUE(st$equal_SE)) {
      cat("Critical difference: ", round(st$critical_difference, 3), "\n",
          sep = "")
    } else if (length(x$critical_range) > 1 && isTRUE(st$equal_SE)) {
      cat("Critical ranges: ",
          paste(names(x$critical_range), round(x$critical_range, 3),
                sep = " ", collapse = ", "), "\n", sep = "")
    } else if (!isTRUE(st$equal_SE)) {
      cat("Unequal replication: each pair tested with its own SE.\n")
    }
    tb <- x$table
    tb$mean <- round(tb$mean, digits)
    tb$SE <- round(tb$SE, 3)
    print(tb, row.names = FALSE)
  }
  if (any(names(ph) != "dunnett")) {
    cat("\nMeans sharing a letter do not differ significantly.\n")
  }
  invisible(NULL)
}


#' Post-Hoc Mean Comparisons After ANOVA
#'
#' @description
#' Compares treatment means after an analysis of variance and returns
#' compact letter displays (means sharing a letter do not differ). Fisher's
#' LSD, Duncan's multiple range test (DMRT), Tukey's HSD,
#' Student-Newman-Keuls, Scheffe and Bonferroni are computed inside the
#' package, so no extra package is needed. Dunnett's test against a
#' control uses \pkg{multcomp} when it is installed.
#'
#' @param model A fitted \code{aov} or \code{lm} model containing the
#'   treatment factor.
#' @param data Data frame used to fit the model. Kept for compatibility;
#'   the model frame is used for the calculations.
#' @param response Name of the response variable (used to label output).
#' @param treatment Name of the treatment factor in the model.
#' @param mse Error mean square. Default \code{NULL} takes it from the
#'   model. Supply it to compare means against a different error term
#'   (for example the main-plot error of a split-plot design).
#' @param df_error Error degrees of freedom matching \code{mse}. Default
#'   \code{NULL} takes it from the model.
#' @param posthoc Method(s): \code{"lsd"}, \code{"duncan"},
#'   \code{"tukey"}, \code{"snk"}, \code{"scheffe"}, \code{"bonferroni"},
#'   \code{"dunnett"}, a vector of these, or \code{"all"}.
#' @param alpha Significance level. Default 0.05.
#' @param control Control level for Dunnett's test. Default \code{NULL}
#'   uses the first level of the treatment factor.
#' @param verbose Logical. If TRUE (default), prints the results.
#'
#' @details
#' Means are least-squares means of the treatment factor averaged over the
#' other factors in the model. With complete, balanced data they equal the
#' arithmetic means. With missing plots or unequal replication every pair
#' is tested with its own standard error of difference, so the letters stay
#' exact (Tukey-Kramer form for the range tests).
#'
#' Duncan's and the Student-Newman-Keuls tests use the studentized range
#' for the number of means spanned by each pair, with the protection level
#' \eqn{1 - (1 - \alpha)^{p - 1}} for Duncan's test, and apply the closure
#' rule: two means are declared different only if every range of means
#' containing them is also significant (Duncan, 1955). The p-value reported
#' for a pair is the largest p-value of the ranges containing it, so it
#' agrees with the letters.
#'
#' Least-squares means that cannot be estimated (for example a level whose
#' cells are all empty in a model with interactions) are left out with a
#' warning. Dunnett's test averages differences over interacting factors
#' and uses \code{mse} and \code{df_error} when they are supplied.
#'
#' Letters are assigned with the insert-and-absorb algorithm of Piepho
#' (2004), which gives a valid display for any pattern of significant
#' pairs. The group containing the highest mean gets the letter "a".
#'
#' @return A list with one element per method. Each element contains
#'   \code{groups} (means and letters, ordered from highest to lowest, with
#'   treatments as row names), \code{table} (level, mean, SE, letter,
#'   replications), \code{comparisons} (all pairwise differences with
#'   their SE, p-value and significance), \code{critical_range} and
#'   \code{statistics}. The Dunnett element is a data frame of comparisons
#'   with the control.
#'
#' @references
#' Duncan, D.B. (1955). Multiple range and multiple F tests.
#' \emph{Biometrics}, 11(1), 1-42. \doi{10.2307/3001478}
#'
#' Piepho, H.-P. (2004). An algorithm for a letter-based representation of
#' all-pairwise comparisons. \emph{Journal of Computational and Graphical
#' Statistics}, 13(2), 456-466. \doi{10.1198/1061860043515}
#'
#' Kramer, C.Y. (1956). Extension of multiple range tests to group means
#' with unequal numbers of replications. \emph{Biometrics}, 12(3),
#' 307-310. \doi{10.2307/3001469}
#'
#' @examples
#' data(mothbean_rbd)
#' fit <- aov(seed_yield ~ block + treatment, data = mothbean_rbd)
#' ph <- perform_posthoc(fit, mothbean_rbd, response = "seed_yield",
#'                       treatment = "treatment", posthoc = c("lsd", "duncan"))
#' ph$duncan$groups
#'
#' @seealso \code{\link{compare_means}} for summary statistics,
#'   \code{\link{anova_rbd}}.
#' @export
perform_posthoc <- function(model, data, response, treatment, mse = NULL,
                            df_error = NULL, posthoc = "lsd", alpha = 0.05,
                            control = NULL, verbose = TRUE) {
  .check_alpha(alpha)
  methods <- .resolve_methods(posthoc)
  if (inherits(model, "aovlist")) {
    stop("Models with Error() strata are not supported here. Use the ",
         "split-plot functions (for example anova_spd()), which give ",
         "letters with the correct error term.", call. = FALSE)
  }
  if (!inherits(model, "lm")) {
    stop("'model' must be a fitted aov() or lm() model.", call. = FALSE)
  }
  mf <- stats::model.frame(model)
  if (!treatment %in% names(mf)) {
    stop("Treatment factor '", treatment, "' is not in the model.",
         call. = FALSE)
  }
  if (!is.factor(mf[[treatment]]) && !is.character(mf[[treatment]])) {
    stop("'", treatment, "' must be a factor in the model.", call. = FALSE)
  }
  if (!is.null(mse) && is.null(df_error)) {
    warning("'mse' given without 'df_error'; the residual degrees of freedom ",
            "of the model are used.", call. = FALSE)
  }
  res <- .posthoc_set(model, treatment, methods, alpha, mse = mse,
                      df = df_error, control = control, response = response)
  if ("dunnett" %in% methods && is.null(res$dunnett) && verbose) {
    message("Install the 'multcomp' package for Dunnett's test.")
  }
  if (verbose && length(res) > 0) {
    .heading("POST-HOC COMPARISONS")
    .print_posthoc(res)
  }
  invisible(res)
}


#' Mean Separation from Summary Statistics
#'
#' @description
#' Assigns significance letters to treatment means when only the means,
#' the error mean square and its degrees of freedom are available (for
#' example from a published table or from another program). Uses the same
#' tests as \code{\link{perform_posthoc}}.
#'
#' @param means Named numeric vector of treatment means.
#' @param mse Error mean square from the ANOVA.
#' @param df_error Error degrees of freedom.
#' @param reps Number of replications per mean: a single number or a
#'   vector with one value per mean.
#' @param method One of \code{"lsd"}, \code{"duncan"}, \code{"tukey"},
#'   \code{"snk"}, \code{"scheffe"} or \code{"bonferroni"}.
#' @param alpha Significance level. Default 0.05.
#' @param verbose Logical. If TRUE (default), prints the table.
#'
#' @return The comparison object described in
#'   \code{\link{perform_posthoc}}: \code{groups}, \code{table},
#'   \code{comparisons}, \code{critical_range} and \code{statistics}.
#'
#' @examples
#' m <- c(T1 = 1180, T2 = 1265, T3 = 1342, T4 = 1410, T5 = 1238)
#' compare_means(m, mse = 3600, df_error = 12, reps = 4, method = "duncan")
#'
#' @export
compare_means <- function(means, mse, df_error, reps, method = "lsd",
                          alpha = 0.05, verbose = TRUE) {
  .check_alpha(alpha)
  if (!is.numeric(means) || length(means) < 2) {
    stop("'means' must be a numeric vector of at least two means.",
         call. = FALSE)
  }
  if (anyNA(means)) stop("'means' must not contain NA.", call. = FALSE)
  if (is.null(names(means))) names(means) <- paste0("T", seq_along(means))
  if (anyDuplicated(names(means))) {
    stop("Names of 'means' must be unique.", call. = FALSE)
  }
  if (!is.numeric(mse) || length(mse) != 1 || mse <= 0) {
    stop("'mse' must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(df_error) || length(df_error) != 1 || df_error < 1) {
    stop("'df_error' must be a single number of at least 1.", call. = FALSE)
  }
  if (!is.numeric(reps) || any(reps <= 0) ||
      !(length(reps) %in% c(1, length(means)))) {
    stop("'reps' must be positive: one value or one per mean.",
         call. = FALSE)
  }
  method <- .resolve_methods(method)
  if (length(method) != 1 || method == "dunnett") {
    stop("Choose a single method other than Dunnett.", call. = FALSE)
  }
  n <- rep_len(reps, length(means))
  names(n) <- names(means)
  V <- diag(mse / n, length(means))
  dimnames(V) <- list(names(means), names(means))
  out <- .compare(means, V, df_error, mse, method, alpha, n = n)
  if (verbose) .print_posthoc(stats::setNames(list(out), method))
  invisible(out)
}
