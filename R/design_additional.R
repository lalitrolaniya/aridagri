# =============================================================================
# STRIP-PLOT, AUGMENTED BLOCK AND ALPHA LATTICE DESIGNS
# Package: aridagri
# =============================================================================

#' Strip-Plot (Strip-Block) Design ANOVA
#'
#' @description
#' Analysis of variance for a strip-plot design, where factor A is applied
#' to horizontal strips and factor B to vertical strips in each
#' replication. A is tested against Error(a) = Rep x A, B against
#' Error(b) = Rep x B, and A x B against Error(c). Gives CV(a), CV(b),
#' CV(c), SEm, SEd and CD for the four standard comparisons, and letter
#' groupings.
#'
#' @param data Data frame containing the experimental data.
#' @param response Name of the response variable.
#' @param horizontal_factor Factor applied in horizontal strips (A).
#' @param vertical_factor Factor applied in vertical strips (B).
#' @param replication Name of the replication column.
#' @param alpha Significance level. Default 0.05.
#' @param verbose Logical. If TRUE (default), prints the report.
#' @param posthoc Method for the letters (default \code{"lsd"}).
#'
#' @details
#' Standard errors of difference (Gomez and Gomez, 1984, table 3.15), with
#' \eqn{a} and \eqn{b} levels of A and B and \eqn{r} replications:
#' A means \eqn{\sqrt{2E_a/(rb)}}; B means \eqn{\sqrt{2E_b/(ra)}}; A means
#' at the same level of B \eqn{\sqrt{2[(b-1)E_c + E_a]/(rb)}}; B means at
#' the same level of A \eqn{\sqrt{2[(a-1)E_c + E_b]/(ra)}}. The last two
#' use the weighted t of Cochran and Cox.
#'
#' @return An object of class \code{aridagri_strip} with the corrected
#'   \code{anova_table}, \code{cv} (a, b, c), error mean squares and df,
#'   \code{se_cd}, \code{horizontal_means} and \code{vertical_means} (with
#'   letters), \code{interaction_means}, \code{letters} and \code{model}.
#'
#' @references
#' Gomez, K.A. and Gomez, A.A. (1984). \emph{Statistical Procedures for
#' Agricultural Research}, 2nd ed. John Wiley and Sons, New York.
#'
#' @examples
#' data <- expand.grid(rep = 1:4, irrigation = c("I1", "I2", "I3"),
#'                     tillage = c("CT", "MT", "ZT"))
#' set.seed(9)
#' data$yield <- rnorm(nrow(data), 1200, 70)
#' anova_strip(data, response = "yield", horizontal_factor = "irrigation",
#'             vertical_factor = "tillage", replication = "rep")
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram
#'   Choudhary, Hemant Poonia
anova_strip <- function(data, response, horizontal_factor, vertical_factor,
                        replication, alpha = 0.05, verbose = TRUE,
                        posthoc = "lsd") {
  .check_alpha(alpha)
  methods <- setdiff(.resolve_methods(posthoc), "dunnett")
  A <- horizontal_factor
  B <- vertical_factor
  R <- replication
  d <- .prep_data(data, response, c(R, A, B), verbose)
  .check_split_balance(d, R, c(A, B))
  a <- nlevels(d[[A]])
  b <- nlevels(d[[B]])
  r <- nlevels(d[[R]])
  if (a < 2 || b < 2) stop("Each factor needs at least 2 levels.",
                           call. = FALSE)
  if (r < 2) stop("Need at least 2 replications.", call. = FALSE)
  rhs <- paste(.bq(R), "+", .bq(A), "+", paste0(.bq(R), ":", .bq(A)), "+",
               .bq(B), "+", paste0(.bq(R), ":", .bq(B)), "+",
               paste0(.bq(A), ":", .bq(B)))
  model <- stats::aov(.formula(response, rhs), data = d)
  at <- stats::anova(model)
  rn <- trimws(gsub("`", "", rownames(at)))
  get <- function(parts) {
    key <- vapply(strsplit(rn, ":", fixed = TRUE),
                  function(z) paste(sort(z), collapse = "|"), character(1))
    which(key == paste(sort(parts), collapse = "|"))
  }
  i_r <- get(R)
  i_a <- get(A)
  i_ea <- get(c(R, A))
  i_b <- get(B)
  i_eb <- get(c(R, B))
  i_ab <- get(c(A, B))
  i_ec <- which(rn == "Residuals")
  ss <- at[["Sum Sq"]]
  dfv <- at[["Df"]]
  ms <- at[["Mean Sq"]]
  Ea <- ms[i_ea]
  Eb <- ms[i_eb]
  Ec <- ms[i_ec]
  fA <- ms[i_a] / Ea
  fB <- ms[i_b] / Eb
  fAB <- ms[i_ab] / Ec
  tab <- data.frame(
    Df = dfv[c(i_r, i_a, i_ea, i_b, i_eb, i_ab, i_ec)],
    `Sum Sq` = ss[c(i_r, i_a, i_ea, i_b, i_eb, i_ab, i_ec)],
    `Mean Sq` = ms[c(i_r, i_a, i_ea, i_b, i_eb, i_ab, i_ec)],
    `F value` = c(NA, fA, NA, fB, NA, fAB, NA),
    `Pr(>F)` = c(NA, stats::pf(fA, dfv[i_a], dfv[i_ea], lower.tail = FALSE),
                 NA, stats::pf(fB, dfv[i_b], dfv[i_eb], lower.tail = FALSE),
                 NA, stats::pf(fAB, dfv[i_ab], dfv[i_ec], lower.tail = FALSE),
                 NA),
    check.names = FALSE)
  rownames(tab) <- c(R, A, "Error(a)", B, "Error(b)", paste0(A, ":", B),
                     "Error(c)")
  gm <- mean(d[[response]])
  E <- c(Ea, Eb, Ec)
  DFE <- dfv[c(i_ea, i_eb, i_ec)]
  mk <- function(lab, co) {
    use <- co > 0
    sed <- sqrt(sum(co * E))
    tw <- .weighted_t(co[use], E[use], DFE[use], alpha)
    data.frame(Comparison = lab, SEm = sed / sqrt(2), SEd = sed,
               t_value = tw, df = .satterthwaite(co[use], E[use], DFE[use]),
               CD = tw * sed,
               Error = paste(c("Ea", "Eb", "Ec")[use], collapse = "+"),
               stringsAsFactors = FALSE)
  }
  se_cd <- rbind(
    mk(paste(A, "means"), c(2 / (r * b), 0, 0)),
    mk(paste(B, "means"), c(0, 2 / (r * a), 0)),
    mk(paste(A, "levels at the same", B), c(2 / (r * b), 0,
                                             2 * (b - 1) / (r * b))),
    mk(paste(B, "levels at the same", A), c(0, 2 / (r * a),
                                             2 * (a - 1) / (r * a))))
  hm <- .split_letters(d, response, A, NULL, se_cd$SEd[1], se_cd$df[1],
                       se_cd$t_value[1], methods, alpha, Ea)
  vm <- .split_letters(d, response, B, NULL, se_cd$SEd[2], se_cd$df[2],
                       se_cd$t_value[2], methods, alpha, Eb)
  if (is.null(hm)) {
    m <- tapply(d[[response]], d[[A]], mean)
    hm <- data.frame(Level = names(m), Mean = as.numeric(m))
    m <- tapply(d[[response]], d[[B]], mean)
    vm <- data.frame(Level = names(m), Mean = as.numeric(m))
  }
  lets <- list()
  if (length(methods) > 0) {
    lets[[paste(A, "within", B)]] <- .split_letters(
      d, response, A, B, se_cd$SEd[3], se_cd$df[3], se_cd$t_value[3],
      methods, alpha, Ec)
    lets[[paste(B, "within", A)]] <- .split_letters(
      d, response, B, A, se_cd$SEd[4], se_cd$df[4], se_cd$t_value[4],
      methods, alpha, Ec)
  }
  result <- list(
    design = "Strip Plot Design",
    response = response,
    factors = c(horizontal = A, vertical = B, replication = R),
    anova_table = tab,
    model = model,
    grand_mean = gm,
    cv = c(a = 100 * sqrt(Ea) / gm, b = 100 * sqrt(Eb) / gm,
           c = 100 * sqrt(Ec) / gm),
    ms_error_a = Ea, ms_error_b = Eb, ms_error_c = Ec,
    df_error_a = DFE[1], df_error_b = DFE[2], df_error_c = DFE[3],
    alpha = alpha,
    se_cd = se_cd,
    horizontal_means = hm,
    vertical_means = vm,
    interaction_means = tapply(d[[response]], list(d[[A]], d[[B]]), mean),
    letters = lets,
    letters_method = .letters_method(methods),
    r = r
  )
  class(result) <- c("aridagri_strip", "list")
  if (verbose) print(result)
  invisible(result)
}

#' @rdname print.aridagri_split
#' @export
print.aridagri_strip <- function(x, digits = 2, ...) {
  f <- x$factors
  .print_info("STRIP-PLOT DESIGN",
              c(Response = x$response, `Horizontal factor (A)` = f[[1]],
                `Vertical factor (B)` = f[[2]], Replications = x$r))
  .print_anova_df(x$anova_table, rownames(x$anova_table))
  cat(sprintf("\n Grand mean = %.3f   CV(a) = %.2f%%   CV(b) = %.2f%%   CV(c) = %.2f%%\n",
              x$grand_mean, x$cv[1], x$cv[2], x$cv[3]))
  cat("\n SEm, SEd and CD (alpha = ", x$alpha, ")\n", sep = "")
  sc <- x$se_cd[, c("Comparison", "SEm", "SEd", "CD", "Error")]
  sc[, 2:4] <- round(sc[, 2:4], 3)
  print(sc, row.names = FALSE, right = FALSE)
  lm_lab <- if (is.na(x$letters_method)) "" else
    paste0(" (", x$letters_method, ")")
  .print_means(x$horizontal_means, paste0("Means of ", f[[1]], lm_lab), digits)
  .print_means(x$vertical_means, paste0("Means of ", f[[2]], lm_lab), digits)
  cat("\n Interaction means\n")
  print(round(x$interaction_means, digits))
  invisible(x)
}


#' Augmented Randomized Complete Block Design
#'
#' @description
#' Analysis of an augmented randomized complete block design (Federer,
#' 1956), where a few check varieties are replicated in every block and
#' the new (test) entries appear once. Gives the ANOVA with the treatment
#' sum of squares split into checks, tests and checks vs tests, adjusted
#' means of all entries, the four standard errors of difference, a
#' comparison of every test entry with the checks, and genetic parameters
#' of the test entries.
#'
#' @param data Data frame containing the data.
#' @param response Name of the response variable.
#' @param genotype Name of the entry (genotype) column.
#' @param block Name of the block column.
#' @param check_names Character vector with the names of the check
#'   entries.
#' @param verbose Logical. If TRUE (default), prints the report.
#' @param alpha Significance level for CD. Default 0.05.
#' @param k Standardized selection differential for genetic advance.
#'   Default 2.063 (5\% selection intensity).
#'
#' @details
#' Block effects are estimated from the checks and removed from the test
#' entries: adjusted value = observed value - (mean of checks in its block -
#' overall check mean). This equals the least-squares mean from the model
#' \code{response ~ block + genotype}, which is used here, so the standard
#' errors are exact. With \eqn{b} blocks and \eqn{c} checks the four
#' standard errors of difference are: two check means
#' \eqn{\sqrt{2 MSE / b}}; two test entries in the same block
#' \eqn{\sqrt{2 MSE}}; two test entries in different blocks
#' \eqn{\sqrt{2 MSE (1 + 1/c)}}; a test entry and a check mean
#' \eqn{\sqrt{MSE (1 + 1/b + 1/c - 1/(bc))}}.
#'
#' Two ANOVA tables are given. In the first, treatments are fitted after
#' blocks and their sum of squares is split into tests vs checks, checks
#' and tests, each adjusted for blocks. In the second, blocks are fitted
#' after treatments and the treatment sum of squares (ignoring blocks) is
#' split into checks, tests and checks vs tests (Federer, 1961).
#'
#' Genetic parameters of the test entries use the block-adjusted mean
#' square among tests from the first table: phenotypic variance
#' \eqn{\sigma^2_p} = MS(tests), error variance \eqn{\sigma^2_e} = MSE,
#' genotypic variance \eqn{\sigma^2_g = \sigma^2_p - \sigma^2_e},
#' broad-sense heritability \eqn{h^2 = \sigma^2_g / \sigma^2_p}, GCV and
#' PCV in percent of the mean of the adjusted test values, and genetic
#' advance \eqn{GA = k h^2 \sigma_p}.
#'
#' @return An object of class \code{aridagri_augmented} with
#'   \code{anova_adjusted} (treatments after blocks),
#'   \code{anova_partitioned} (blocks after treatments, treatments split),
#'   \code{adjusted_means} (all entries with type, raw and adjusted mean,
#'   SE and block), \code{se_diff} and \code{cd} (the four comparisons),
#'   \code{test_vs_checks} (difference of each test entry from the check
#'   mean and the best check, with significance), \code{block_effects},
#'   \code{genetic_parameters}, \code{mse}, \code{df_error},
#'   \code{anova_checks} and \code{model}.
#'
#' @references
#' Federer, W.T. (1956). Augmented (or hoonuiaku) designs.
#' \emph{Hawaiian Planters' Record}, 55, 191-208.
#'
#' Federer, W.T. (1961). Augmented designs with one-way elimination of
#' heterogeneity. \emph{Biometrics}, 17(3), 447-473.
#' \doi{10.2307/2527837}
#'
#' @examples
#' data(chickpea_augmented)
#' res <- anova_augmented(chickpea_augmented, response = "seed_yield",
#'                        genotype = "entry", block = "block",
#'                        check_names = c("C1", "C2", "C3", "C4"))
#' head(res$test_vs_checks)
#' res$genetic_parameters
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram
#'   Choudhary, Hemant Poonia
anova_augmented <- function(data, response, genotype, block, check_names,
                            verbose = TRUE, alpha = 0.05, k = 2.063) {
  .check_alpha(alpha)
  d <- .prep_data(data, response, c(genotype, block), verbose)
  g <- d[[genotype]]
  blk <- d[[block]]
  check_names <- as.character(check_names)
  absent <- setdiff(check_names, levels(g))
  if (length(absent) > 0) {
    stop("Check(s) not found in '", genotype, "': ",
         paste(absent, collapse = ", "), call. = FALSE)
  }
  is_chk <- as.character(g) %in% check_names
  b <- nlevels(blk)
  cc <- length(check_names)
  if (b < 2) stop("Need at least 2 blocks.", call. = FALSE)
  if (cc < 2) stop("Need at least 2 checks.", call. = FALSE)
  ct <- table(droplevels(blk[is_chk]), droplevels(g[is_chk]))
  if (nrow(ct) < b || any(ct == 0)) {
    stop("Every check must appear in every block.", call. = FALSE)
  }
  if (any(ct > 1)) {
    warning("Some checks appear more than once in a block.", call. = FALSE)
  }
  tests <- setdiff(levels(g), check_names)
  if (length(tests) < 1) stop("No test entries found.", call. = FALSE)
  rep_t <- table(droplevels(g[!is_chk]))
  if (any(rep_t > 1)) {
    message("Some test entries appear more than once; they are analysed ",
            "with their replication.")
  }
  y <- d[[response]]
  model <- stats::lm(.formula(response, paste(.bq(block), "+",
                                              .bq(genotype))), data = d)
  a1 <- stats::anova(model)
  df_e <- a1[["Df"]][3]
  if (df_e < 1) stop("No error degrees of freedom.", call. = FALSE)
  mse <- a1[["Mean Sq"]][3]

  # Treatments ignoring blocks, split into checks, tests, checks vs tests.
  d$.type <- factor(ifelse(is_chk, "Check", "Test"), levels = c("Check",
                                                                "Test"))
  d$.chk <- factor(ifelse(is_chk, as.character(g), ".test"))
  m_g <- stats::lm(.formula(response, .bq(genotype)), data = d)
  ss_trt <- sum((stats::fitted(m_g) - mean(y))^2)
  yc <- y[is_chk]
  gc <- droplevels(g[is_chk])
  ss_chk <- sum(tapply(yc, gc, function(v) length(v) * (mean(v) - mean(yc))^2))
  yt <- y[!is_chk]
  gt <- droplevels(g[!is_chk])
  ss_tst <- sum(tapply(yt, gt, function(v) length(v) * (mean(v) - mean(yt))^2))
  ss_cvt <- ss_trt - ss_chk - ss_tst
  ss_blk_adj <- sum(stats::residuals(m_g)^2) - sum(stats::residuals(model)^2)
  nt <- length(tests)
  df_p <- c(b - 1, cc + nt - 1, cc - 1, nt - 1, 1, df_e)
  ss_p <- c(ss_blk_adj, ss_trt, ss_chk, ss_tst, ss_cvt, mse * df_e)
  ms_p <- ss_p / df_p
  f_p <- c(ms_p[1:5] / mse, NA)
  part <- data.frame(Df = df_p, `Sum Sq` = ss_p, `Mean Sq` = ms_p,
                     `F value` = f_p,
                     `Pr(>F)` = stats::pf(f_p, df_p, df_e, lower.tail = FALSE),
                     check.names = FALSE)
  rownames(part) <- c("Blocks (eliminating treatments)",
                      "Treatments (ignoring blocks)", "  Checks", "  Tests",
                      "  Checks vs Tests", "Residuals")
  # Treatments eliminating blocks, split into tests vs checks, checks and
  # tests (each adjusted for blocks and for the sources before it).
  rss <- function(f) sum(stats::residuals(stats::lm(f, data = d))^2)
  rss0 <- rss(.formula(response, .bq(block)))
  rss_b <- rss(.formula(response, paste(.bq(block), "+ .type")))
  rss_a <- rss(.formula(response, paste(.bq(block), "+ .chk")))
  rss_f <- sum(stats::residuals(model)^2)
  df_a <- c(b - 1, cc + nt - 1, 1, cc - 1, nt - 1, df_e)
  ss_a <- c(a1[["Sum Sq"]][1], rss0 - rss_f, rss0 - rss_b, rss_b - rss_a,
            rss_a - rss_f, rss_f)
  ms_a <- ss_a / df_a
  f_a <- c(ms_a[1:5] / mse, NA)
  adj_tab <- data.frame(Df = df_a, `Sum Sq` = ss_a, `Mean Sq` = ms_a,
                        `F value` = f_a,
                        `Pr(>F)` = stats::pf(f_a, df_a, df_e,
                                             lower.tail = FALSE),
                        check.names = FALSE)
  rownames(adj_tab) <- c("Blocks (ignoring treatments)",
                         "Treatments (eliminating blocks)",
                         "  Tests vs Checks", "  Checks", "  Tests",
                         "Residuals")
  ms_tests <- ms_a[5]

  # Adjusted means
  lsm <- .ls_means_model(model, genotype)
  raw <- tapply(y, g, mean)
  blk_of <- tapply(as.character(blk), g, function(v) paste(unique(v),
                                                           collapse = ","))
  adj <- data.frame(Entry = names(lsm$means),
                    Type = ifelse(names(lsm$means) %in% check_names, "Check",
                                  "Test"),
                    Block = as.character(blk_of[names(lsm$means)]),
                    Raw_mean = as.numeric(raw[names(lsm$means)]),
                    Adjusted_mean = as.numeric(lsm$means),
                    SE = sqrt(pmax(diag(lsm$vcov), 0)),
                    stringsAsFactors = FALSE)
  adj$Block[adj$Type == "Check"] <- "all"
  adj <- adj[order(adj$Type, -adj$Adjusted_mean), ]
  rownames(adj) <- NULL

  cm <- tapply(yc, droplevels(blk[is_chk]), mean)
  block_effects <- cm - mean(yc)

  se <- c(`Two check means` = sqrt(2 * mse / b),
          `Two test entries (same block)` = sqrt(2 * mse),
          `Two test entries (different blocks)` = sqrt(2 * mse * (1 + 1 / cc)),
          `Test entry vs check mean` = sqrt(mse * (1 + 1 / b + 1 / cc -
                                                    1 / (b * cc))))
  tcrit <- stats::qt(1 - alpha / 2, df_e)
  cd <- tcrit * se

  # Each test entry against the checks
  V <- lsm$vcov
  chk_means <- lsm$means[check_names]
  best <- names(which.max(chk_means))
  tst <- names(lsm$means)[!names(lsm$means) %in% check_names]
  w <- rep(1 / cc, cc)
  tv <- lapply(tst, function(tn) {
    d_mean <- lsm$means[[tn]] - mean(chk_means)
    v_mean <- V[tn, tn] + as.numeric(t(w) %*% V[check_names, check_names] %*%
                                       w) - 2 * sum(w * V[tn, check_names])
    d_best <- lsm$means[[tn]] - chk_means[[best]]
    v_best <- V[tn, tn] + V[best, best] - 2 * V[tn, best]
    data.frame(Entry = tn, Adjusted_mean = lsm$means[[tn]],
               Diff_vs_check_mean = d_mean,
               Pct_over_check_mean = 100 * d_mean / mean(chk_means),
               p_vs_check_mean = 2 * stats::pt(-abs(d_mean / sqrt(v_mean)),
                                               df_e),
               Diff_vs_best_check = d_best,
               Pct_over_best_check = 100 * d_best / chk_means[[best]],
               p_vs_best_check = 2 * stats::pt(-abs(d_best / sqrt(v_best)),
                                               df_e),
               stringsAsFactors = FALSE)
  })
  tvc <- do.call(rbind, tv)
  tvc$Sig_vs_best <- ifelse(tvc$p_vs_best_check < alpha,
                            ifelse(tvc$Diff_vs_best_check > 0, "higher",
                                   "lower"), "ns")
  tvc <- tvc[order(-tvc$Adjusted_mean), ]
  rownames(tvc) <- NULL

  # Genetic parameters of the test entries, from the block-adjusted mean
  # square among tests: sigma2_p = MS(tests), sigma2_g = MS(tests) - MSE.
  at <- lsm$means[tst]
  gp <- NULL
  if (length(at) >= 3) {
    vp <- ms_tests
    vg <- vp - mse
    mt <- mean(at)
    h2 <- if (vp > 0) vg / vp else NA_real_
    gp <- data.frame(
      Parameter = c("Mean of tests", "Range (min)", "Range (max)",
                    "Phenotypic variance", "Genotypic variance",
                    "Error variance", "PCV (%)", "GCV (%)",
                    "Heritability (broad sense, %)", "Genetic advance",
                    "GA as % of mean"),
      Value = c(mt, min(at), max(at), vp, vg, mse, 100 * sqrt(vp) / mt,
                if (vg > 0) 100 * sqrt(vg) / mt else NA_real_, 100 * h2,
                k * h2 * sqrt(vp), 100 * k * h2 * sqrt(vp) / mt),
      stringsAsFactors = FALSE)
  }
  chk_d <- d[is_chk, , drop = FALSE]
  chk_d[[genotype]] <- droplevels(chk_d[[genotype]])
  chk_d[[block]] <- droplevels(chk_d[[block]])
  anova_checks <- stats::anova(stats::lm(.formula(response, paste(
    .bq(block), "+", .bq(genotype))), data = chk_d))

  result <- list(
    design = "Augmented Block Design",
    response = response,
    factors = c(genotype = genotype, block = block),
    anova_adjusted = adj_tab,
    anova_partitioned = part,
    anova_checks = anova_checks,
    adjusted_means = adj,
    se_diff = se,
    cd = cd,
    alpha = alpha,
    test_vs_checks = tvc,
    best_check = best,
    block_effects = block_effects,
    genetic_parameters = gp,
    mse = mse,
    df_error = df_e,
    cd_check = cd[[1]],
    grand_mean = mean(y),
    cv = 100 * sqrt(mse) / mean(y),
    n_checks = cc, n_tests = nt, n_blocks = b,
    model = model
  )
  class(result) <- c("aridagri_augmented", "list")
  if (verbose) print(result)
  invisible(result)
}

#' Print Method for Augmented Design Results
#'
#' @param x An object of class \code{aridagri_augmented}.
#' @param digits Number of decimal places for means. Default 2.
#' @param n Number of top test entries to show. Default 10.
#' @param ... Unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.aridagri_augmented <- function(x, digits = 2, n = 10, ...) {
  .print_info("AUGMENTED RANDOMIZED COMPLETE BLOCK DESIGN",
              c(Response = x$response, Blocks = x$n_blocks,
                Checks = x$n_checks, `Test entries` = x$n_tests))
  .print_anova_df(x$anova_adjusted, rownames(x$anova_adjusted),
                  total = FALSE,
                  title = "ANOVA (treatments adjusted for blocks)")
  .print_anova_df(x$anova_partitioned, rownames(x$anova_partitioned),
                  total = FALSE,
                  title = "ANOVA (blocks adjusted; treatments partitioned)")
  cat(sprintf("\n Grand mean = %.3f   CV = %.2f%%   MSE = %.4f (df = %d)\n",
              x$grand_mean, x$cv, x$mse, as.integer(x$df_error)))
  cat("\n Standard errors of difference and CD (alpha = ", x$alpha, ")\n",
      sep = "")
  sx <- data.frame(Comparison = names(x$se_diff), SEd = round(x$se_diff, 3),
                   CD = round(x$cd, 3), row.names = NULL)
  print(sx, row.names = FALSE, right = FALSE)
  am <- x$adjusted_means
  .print_means(am[am$Type == "Check", c("Entry", "Raw_mean", "SE")],
               "Check means", digits)
  top <- utils::head(x$test_vs_checks[, c("Entry", "Adjusted_mean",
                                          "Pct_over_check_mean",
                                          "Pct_over_best_check",
                                          "Sig_vs_best")], n)
  .print_means(top, paste0("Top ", nrow(top), " test entries (adjusted; best ",
                           "check = ", x$best_check, ")"), digits)
  if (!is.null(x$genetic_parameters)) {
    gp <- x$genetic_parameters
    gp$Value <- round(gp$Value, 3)
    cat("\n Genetic parameters of test entries\n")
    print(gp, row.names = FALSE, right = FALSE)
  }
  invisible(x)
}


#' Alpha Lattice Design ANOVA
#'
#' @description
#' Analysis of an alpha lattice (resolvable incomplete block) design.
#' Gives the intra-block ANOVA, genotype means adjusted for incomplete
#' blocks with exact standard errors and letter groupings, and the
#' efficiency of the design relative to a randomized complete block design.
#' By default replications and blocks are treated as random (REML), so
#' that inter-block information is recovered (Patterson and Williams,
#' 1976); \code{analysis = "intra"} gives the classical intra-block
#' estimates.
#'
#' @param data Data frame containing the data.
#' @param response Name of the response variable.
#' @param genotype Name of the genotype column.
#' @param replication Name of the replication column.
#' @param block Name of the incomplete block column (nested within
#'   replication; labels may repeat across replications).
#' @param verbose Logical. If TRUE (default), prints the report.
#' @param analysis \code{"reml"} (default; replications and incomplete
#'   blocks random, fitted with \pkg{nlme}, so that inter-block information
#'   is recovered) or \code{"intra"} (blocks fixed; the classical
#'   intra-block analysis). If the REML fit fails, the intra-block analysis
#'   is used with a message.
#' @param posthoc Method for the letters (default \code{"lsd"}).
#' @param alpha Significance level. Default 0.05.
#'
#' @details
#' The intra-block ANOVA is the sequential analysis of the model
#' replication + block(replication) + genotype: genotypes are tested after
#' blocks. Blocks within replications are also shown adjusted for
#' genotypes. Adjusted means are least-squares means. The efficiency
#' relative to an RCBD is the average variance of a difference under the
#' RCBD analysis (\eqn{2 MSE_{RCBD} / r}) divided by the average variance
#' of a difference of adjusted means. With REML, letters and CD use the
#' intra-block error degrees of freedom, an approximation.
#'
#' @return An object of class \code{aridagri_alpha} with
#'   \code{anova_table}, \code{genotype_means} (adjusted means with SE and
#'   letters, highest first), \code{avg_sed}, \code{cd} (from the average
#'   SEd), \code{efficiency} (percent, relative to RCBD), \code{cv},
#'   \code{grand_mean}, \code{posthoc}, \code{variance_components} (REML)
#'   and \code{model}.
#'
#' @references
#' Patterson, H.D. and Williams, E.R. (1976). A new class of resolvable
#' incomplete block designs. \emph{Biometrika}, 63(1), 83-92.
#' \doi{10.1093/biomet/63.1.83}
#'
#' @examples
#' data(cowpea_alpha)
#' res <- anova_alpha_lattice(cowpea_alpha, response = "grain_yield",
#'                            genotype = "genotype", replication = "rep",
#'                            block = "block")
#' head(res$genotype_means)
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram
#'   Choudhary, Hemant Poonia
anova_alpha_lattice <- function(data, response, genotype, replication, block,
                                verbose = TRUE,
                                analysis = c("reml", "intra"),
                                posthoc = "lsd", alpha = 0.05) {
  .check_alpha(alpha)
  analysis <- match.arg(analysis)
  methods <- setdiff(.resolve_methods(posthoc), "dunnett")
  d <- .prep_data(data, response, c(genotype, replication, block), verbose)
  d$.ib <- interaction(d[[replication]], d[[block]], drop = TRUE, sep = ":")
  g <- nlevels(d[[genotype]])
  r <- nlevels(d[[replication]])
  if (r < 2) stop("Need at least 2 replications.", call. = FALSE)
  if (g < 3) stop("Need at least 3 genotypes.", call. = FALSE)
  y <- d[[response]]
  gm <- mean(y)
  f1 <- .formula(response, paste(.bq(replication), "+ .ib +", .bq(genotype)))
  f2 <- .formula(response, paste(.bq(replication), "+", .bq(genotype),
                                 "+ .ib"))
  m1 <- stats::lm(f1, data = d)
  m2 <- stats::lm(f2, data = d)
  a1 <- stats::anova(m1)
  a2 <- stats::anova(m2)
  df_e <- a1[["Df"]][4]
  if (df_e < 1) stop("No error degrees of freedom.", call. = FALSE)
  mse <- a1[["Mean Sq"]][4]
  tab <- rbind(a1[1, ], a1[2, ], a2[3, ], a1[3, ], a1[4, ])
  tab[["F value"]] <- c(NA, NA, tab[["Mean Sq"]][3:4] / mse, NA)
  tab[["Pr(>F)"]] <- c(NA, NA, stats::pf(tab[["F value"]][3:4],
                                         tab[["Df"]][3:4], df_e,
                                         lower.tail = FALSE), NA)
  rownames(tab) <- c(replication, "Blocks within rep (unadjusted)",
                     "Blocks within rep (adjusted)",
                     paste(genotype, "(adjusted)"), "Intra-block error")
  # RCBD analysis for efficiency
  m0 <- stats::lm(.formula(response, paste(.bq(replication), "+",
                                           .bq(genotype))), data = d)
  mse_rcbd <- sum(stats::residuals(m0)^2) / m0$df.residual
  vc <- NULL
  fit_reml <- NULL
  if (analysis == "reml") {
    d$.rep <- d[[replication]]
    d$.gen <- d[[genotype]]
    d$.y <- y
    fit_reml <- tryCatch(
      nlme::lme(.y ~ 0 + .gen, random = ~ 1 | .rep / .ib, data = d,
                method = "REML"),
      error = function(e) NULL)
    if (is.null(fit_reml)) {
      if (verbose) message("REML fit failed; using the intra-block analysis.")
      analysis <- "intra"
    }
  }
  if (analysis == "intra") {
    mi <- stats::lm(.formula(response, paste(".ib +", .bq(genotype))),
                    data = d)
    lsm <- .ls_means_model(mi, genotype)
    if (!isTRUE(lsm$exact)) {
      warning("The design is disconnected: some genotype contrasts are not ",
              "estimable in the intra-block analysis. Use analysis = ",
              "\"reml\".", call. = FALSE)
    }
    V <- lsm$vcov
    df_cmp <- df_e
    mse_cmp <- mse
  } else {
    fit <- fit_reml
    b <- nlme::fixef(fit)
    V <- as.matrix(stats::vcov(fit))
    lev <- sub("^\\.gen", "", names(b))
    names(b) <- lev
    dimnames(V) <- list(lev, lev)
    lsm <- list(means = b, vcov = V, n = table(d[[genotype]])[lev])
    vcomp <- nlme::VarCorr(fit)
    vv <- suppressWarnings(as.numeric(vcomp[, "Variance"]))
    vc <- data.frame(Component = c(replication, "Block within rep",
                                   "Residual"),
                     Variance = vv[!is.na(vv)], stringsAsFactors = FALSE)
    df_cmp <- df_e
    mse_cmp <- vc$Variance[3]
  }
  sed <- sqrt(pmax(outer(diag(V), diag(V), "+") - 2 * V, 0))
  avg_v <- mean(sed[upper.tri(sed)]^2)
  avg_sed <- sqrt(avg_v)
  eff <- 100 * (2 * mse_rcbd / r) / avg_v
  ph <- list()
  for (m in methods) {
    ph[[m]] <- .compare(lsm$means, V, df_cmp, mse_cmp, m, alpha, n = lsm$n,
                        response = response)
  }
  raw <- tapply(y, d[[genotype]], mean)
  gmeans <- data.frame(Genotype = names(lsm$means),
                       Mean = as.numeric(raw[names(lsm$means)]),
                       Adjusted_mean = as.numeric(lsm$means),
                       SE = sqrt(pmax(diag(V), 0)),
                       stringsAsFactors = FALSE)
  if (length(methods) > 0) {
    gmeans$Group <- ph[[methods[1]]]$letters[gmeans$Genotype]
  }
  gmeans <- gmeans[order(-gmeans$Adjusted_mean), ]
  rownames(gmeans) <- NULL
  result <- list(
    design = "Alpha Lattice Design",
    response = response,
    analysis = analysis,
    factors = c(genotype = genotype, replication = replication,
                block = block),
    anova_table = tab,
    model = m1,
    grand_mean = gm,
    cv = 100 * sqrt(mse) / gm,
    ms_error = mse,
    df_error = df_e,
    mse_rcbd = mse_rcbd,
    efficiency = eff,
    avg_sed = avg_sed,
    cd = stats::qt(1 - alpha / 2, df_cmp) * avg_sed,
    alpha = alpha,
    genotype_means = gmeans,
    variance_components = vc,
    letters_method = .letters_method(methods),
    posthoc = ph,
    n_blocks = nlevels(d$.ib), r = r
  )
  class(result) <- c("aridagri_alpha", "list")
  if (verbose) print(result)
  invisible(result)
}

#' Print Method for Alpha Lattice Results
#'
#' @param x An object of class \code{aridagri_alpha}.
#' @param digits Number of decimal places for means. Default 2.
#' @param n Number of genotypes to show. Default 15.
#' @param ... Unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.aridagri_alpha <- function(x, digits = 2, n = 15, ...) {
  .print_info("ALPHA LATTICE DESIGN",
              c(Response = x$response,
                Analysis = if (x$analysis == "intra") "intra-block (blocks fixed)" else
                  "REML (blocks random)",
                Genotypes = nrow(x$genotype_means), Replications = x$r,
                `Incomplete blocks` = x$n_blocks))
  .print_anova_df(x$anova_table, rownames(x$anova_table), total = FALSE,
                  title = "Intra-block analysis of variance")
  cat(sprintf(paste0("\n Grand mean = %.3f   CV = %.2f%%   Average SEd = ",
                     "%.3f   CD (%g%%) = %.3f\n"),
              x$grand_mean, x$cv, x$avg_sed, 100 * x$alpha, x$cd))
  cat(sprintf(" Efficiency relative to RCBD = %.1f%%\n", x$efficiency))
  if (!is.null(x$variance_components)) {
    cat("\n REML variance components\n")
    vc <- x$variance_components
    vc$Variance <- signif(vc$Variance, 4)
    print(vc, row.names = FALSE)
  }
  gm <- utils::head(x$genotype_means, n)
  ttl <- paste0("Adjusted genotype means (top ", nrow(gm), if (!is.na(
    x$letters_method)) paste0(", ", x$letters_method) else "", ")")
  .print_means(gm, ttl, digits)
  invisible(x)
}
