# =============================================================================
# FACTORIAL AND LATIN SQUARE DESIGNS
# Package: aridagri
# =============================================================================

# Precision adjustment factor for relative efficiency when error df < 20
# (Fisher's information correction, as used by Gomez and Gomez, 1984).
.re_adjust <- function(df1, df2) {
  if (df1 >= 20) return(1)
  ((df1 + 1) * (df2 + 3)) / ((df1 + 3) * (df2 + 1))
}

# Means of one effect (single factor or combination) with letters. These
# are least-squares means from the full factorial model, so they are
# adjusted for the other factors when data are unbalanced.
.effect_means <- function(model, response, effect, methods, alpha, ms_e,
                          df_e) {
  lsm <- .ls_means_model(model, effect)
  lsm <- .drop_nonestimable(lsm, paste(effect, collapse = " x "))
  V <- lsm$vcov * ms_e / lsm$mse
  ph <- list()
  for (m in setdiff(methods, "dunnett")) {
    ph[[m]] <- .compare(lsm$means, V, df_e, ms_e, m, alpha, n = lsm$n,
                        response = response)
  }
  lev <- names(lsm$means)
  parts <- do.call(rbind, strsplit(lev, " : ", fixed = TRUE))
  if (is.null(dim(parts))) parts <- matrix(parts, ncol = length(effect))
  tab <- as.data.frame(parts, stringsAsFactors = FALSE)
  names(tab) <- effect
  tab$Mean <- as.numeric(lsm$means)
  tab$SE <- sqrt(pmax(diag(V), 0))
  tab$Obs <- as.numeric(lsm$n)
  if (length(ph) > 0) tab$Group <- ph[[1]]$letters[lev]
  rownames(tab) <- NULL
  list(table = tab, posthoc = ph)
}

# SEm / SEd / CD table for the effects of a single-stratum factorial.
.factorial_comparisons <- function(effects, n_per_mean, ms_e, df_e, alpha) {
  sem <- sqrt(ms_e / n_per_mean)
  data.frame(Effect = effects, Obs_per_mean = n_per_mean, SEm = sem,
             SEd = sqrt(2) * sem,
             CD = stats::qt(1 - alpha / 2, df_e) * sqrt(2) * sem,
             stringsAsFactors = FALSE)
}

#' Two-Factor Factorial ANOVA
#'
#' @description
#' Analysis of variance for a two-factor factorial experiment in a
#' randomized complete block design (when \code{replication} is given) or
#' a completely randomized design, with SEm, SEd and CD for both main
#' effects and the interaction, and letter groupings for factor means and
#' treatment combinations.
#'
#' @param data Data frame containing the data.
#' @param response Name of the response variable.
#' @param factor1 Name of the first factor (A).
#' @param factor2 Name of the second factor (B).
#' @param replication Name of the replication (block) column, or
#'   \code{NULL} for a completely randomized layout.
#' @param posthoc Post-hoc method for the letters (see
#'   \code{\link{perform_posthoc}}). Default \code{"lsd"}.
#' @param alpha Significance level for CD and letters. Default 0.05.
#' @param verbose Logical. If TRUE (default), prints the report.
#'
#' @return An object of class \code{aridagri_factorial} with
#'   \code{anova_table}, \code{model}, \code{grand_mean}, \code{cv},
#'   \code{ms_error}, \code{df_error}, \code{comparisons} (SEm, SEd and CD
#'   for A, B and A x B), \code{se_a}, \code{cd_a}, \code{se_b},
#'   \code{cd_b}, \code{se_ab}, \code{cd_ab}, \code{mean_a},
#'   \code{mean_b} (with letters), \code{interaction_means} (all
#'   combinations with letters), \code{mean_ab} (two-way table) and
#'   \code{posthoc}.
#'
#' @examples
#' data <- expand.grid(rep = 1:4, nitrogen = c("N0", "N40", "N80"),
#'                     phosphorus = c("P0", "P30"))
#' set.seed(1)
#' data$yield <- 900 + 3 * as.numeric(sub("N", "", data$nitrogen)) +
#'   2 * as.numeric(sub("P", "", data$phosphorus)) + rnorm(nrow(data), 0, 40)
#' res <- anova_factorial(data, response = "yield", factor1 = "nitrogen",
#'                        factor2 = "phosphorus", replication = "rep")
#' res$interaction_means
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
anova_factorial <- function(data, response, factor1, factor2,
                            replication = NULL, posthoc = "lsd", alpha = 0.05,
                            verbose = TRUE) {
  .check_alpha(alpha)
  methods <- setdiff(.resolve_methods(posthoc), "dunnett")
  facs <- c(factor1, factor2)
  d <- .prep_data(data, response, c(facs, replication), verbose)
  a <- nlevels(d[[factor1]])
  b <- nlevels(d[[factor2]])
  if (a < 2 || b < 2) stop("Each factor needs at least 2 levels.",
                           call. = FALSE)
  if (!is.null(replication)) {
    r <- nlevels(d[[replication]])
    balanced <- .is_balanced(d, replication, facs)
    rhs <- paste(.bq(replication), "+", .bq(factor1), "*", .bq(factor2))
  } else {
    rc <- .reps_per_cell(d, facs)
    balanced <- !is.na(rc)
    r <- if (balanced) rc else nrow(d) / (a * b)
    rhs <- paste(.bq(factor1), "*", .bq(factor2))
  }
  if (!balanced) {
    warning("Unbalanced data: sums of squares are sequential (A, then B, ",
            "then A x B) and SE/CD use the average replication. Means are ",
            "least-squares means and letters use their exact standard ",
            "errors.", call. = FALSE)
  }
  model <- stats::aov(.formula(response, rhs), data = d)
  anova_table <- .clean_rownames(stats::anova(model))
  nr <- nrow(anova_table)
  ms_e <- anova_table[["Mean Sq"]][nr]
  df_e <- anova_table[["Df"]][nr]
  if (df_e < 1) stop("No error degrees of freedom.", call. = FALSE)
  gm <- mean(d[[response]])
  comp <- .factorial_comparisons(c(factor1, factor2,
                                   paste(factor1, "x", factor2)),
                                 c(b * r, a * r, r), ms_e, df_e, alpha)
  ea <- .effect_means(model, response, factor1, methods, alpha, ms_e, df_e)
  eb <- .effect_means(model, response, factor2, methods, alpha, ms_e, df_e)
  eab <- .effect_means(model, response, facs, methods, alpha, ms_e, df_e)
  mean_a <- ea$table
  names(mean_a)[1] <- "Level"
  mean_b <- eb$table
  names(mean_b)[1] <- "Level"
  im <- eab$table
  wide <- stats::reshape(im[, c(factor1, factor2, "Mean")], idvar = factor1,
                         timevar = factor2, direction = "wide")
  names(wide) <- sub("^Mean\\.", "", names(wide))
  names(wide)[1] <- "A"
  rownames(wide) <- NULL
  result <- list(
    design = "Two-Factor Factorial",
    layout = if (is.null(replication)) "CRD" else "RCBD",
    response = response,
    factors = c(A = factor1, B = factor2),
    anova_table = anova_table,
    model = model,
    grand_mean = gm,
    cv = 100 * sqrt(ms_e) / gm,
    ms_error = ms_e,
    df_error = df_e,
    alpha = alpha,
    comparisons = comp,
    se_a = comp$SEm[1], cd_a = comp$CD[1],
    se_b = comp$SEm[2], cd_b = comp$CD[2],
    se_ab = comp$SEm[3], cd_ab = comp$CD[3],
    mean_a = mean_a,
    mean_b = mean_b,
    interaction_means = im,
    mean_ab = wide,
    letters_method = .letters_method(methods),
    posthoc = list(A = ea$posthoc, B = eb$posthoc, AB = eab$posthoc),
    balanced = balanced
  )
  class(result) <- c("aridagri_factorial", "list")
  if (verbose) print(result)
  invisible(result)
}


#' Three-Factor Factorial ANOVA
#'
#' @description
#' Analysis of variance for a three-factor factorial experiment in a
#' randomized complete block design (when \code{replication} is given) or
#' a completely randomized design, with SEm, SEd and CD for every main
#' effect and interaction and letter groupings for main-effect means and
#' treatment combinations.
#'
#' @inheritParams anova_factorial
#' @param factor3 Name of the third factor (C).
#'
#' @return An object of class \code{aridagri_factorial3} with
#'   \code{anova_table}, \code{model}, \code{grand_mean}, \code{cv},
#'   \code{ms_error}, \code{df_error}, \code{comparisons} (SEm, SEd and CD
#'   for all seven effects), \code{means} (a list of tables for A, B, C,
#'   A x B, A x C, B x C and A x B x C, with letters) and \code{posthoc}.
#'
#' @examples
#' data <- expand.grid(rep = 1:3, A = c("A1", "A2"), B = c("B1", "B2"),
#'                     C = c("C1", "C2", "C3"))
#' set.seed(2)
#' data$yield <- rnorm(nrow(data), 1200, 60)
#' res <- anova_factorial_3way(data, "yield", "A", "B", "C", "rep")
#' res$comparisons
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram
#'   Choudhary, Hemant Poonia
anova_factorial_3way <- function(data, response, factor1, factor2, factor3,
                                 replication = NULL, alpha = 0.05,
                                 verbose = TRUE, posthoc = "lsd") {
  .check_alpha(alpha)
  methods <- setdiff(.resolve_methods(posthoc), "dunnett")
  facs <- c(factor1, factor2, factor3)
  d <- .prep_data(data, response, c(facs, replication), verbose)
  lv <- vapply(facs, function(f) nlevels(d[[f]]), numeric(1))
  if (any(lv < 2)) stop("Each factor needs at least 2 levels.", call. = FALSE)
  if (!is.null(replication)) {
    r <- nlevels(d[[replication]])
    balanced <- .is_balanced(d, replication, facs)
    rhs <- paste(.bq(replication), "+",
                 paste(.bq(facs), collapse = " * "))
  } else {
    rc <- .reps_per_cell(d, facs)
    balanced <- !is.na(rc)
    r <- if (balanced) rc else nrow(d) / prod(lv)
    rhs <- paste(.bq(facs), collapse = " * ")
  }
  if (!balanced) {
    warning("Unbalanced data: sums of squares are sequential and SE/CD use ",
            "the average replication. Means are least-squares means.",
            call. = FALSE)
  }
  model <- stats::aov(.formula(response, rhs), data = d)
  anova_table <- .clean_rownames(stats::anova(model))
  nr <- nrow(anova_table)
  ms_e <- anova_table[["Mean Sq"]][nr]
  df_e <- anova_table[["Df"]][nr]
  if (df_e < 1) stop("No error degrees of freedom.", call. = FALSE)
  gm <- mean(d[[response]])
  sets <- list(1, 2, 3, c(1, 2), c(1, 3), c(2, 3), c(1, 2, 3))
  labs <- vapply(sets, function(s) paste(facs[s], collapse = " x "),
                 character(1))
  n_mean <- vapply(sets, function(s) r * prod(lv[-s]), numeric(1))
  comp <- .factorial_comparisons(labs, n_mean, ms_e, df_e, alpha)
  means <- list()
  ph <- list()
  for (i in seq_along(sets)) {
    em <- .effect_means(model, response, facs[sets[[i]]], methods, alpha,
                        ms_e, df_e)
    means[[labs[i]]] <- em$table
    ph[[labs[i]]] <- em$posthoc
  }
  result <- list(
    design = "Three-Factor Factorial",
    layout = if (is.null(replication)) "CRD" else "RCBD",
    response = response,
    factors = c(A = factor1, B = factor2, C = factor3),
    anova_table = anova_table,
    model = model,
    grand_mean = gm,
    cv = 100 * sqrt(ms_e) / gm,
    ms_error = ms_e,
    df_error = df_e,
    alpha = alpha,
    comparisons = comp,
    means = means,
    mean_a = means[[1]], mean_b = means[[2]], mean_c = means[[3]],
    letters_method = .letters_method(methods),
    posthoc = ph,
    balanced = balanced
  )
  class(result) <- c("aridagri_factorial3", "list")
  if (verbose) print(result)
  invisible(result)
}


#' Print Methods for Factorial Results
#'
#' @param x An object of class \code{aridagri_factorial} or
#'   \code{aridagri_factorial3}.
#' @param digits Number of decimal places for means. Default 2.
#' @param ... Unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @name print.aridagri_factorial
#' @export
print.aridagri_factorial <- function(x, digits = 2, ...) {
  f <- x$factors
  .print_info(paste0("TWO-FACTOR FACTORIAL (", x$layout, ")"),
              c(Response = x$response, `Factor A` = f[["A"]],
                `Factor B` = f[["B"]],
                Observations = nrow(stats::model.frame(x$model))))
  tab <- x$anova_table
  labs <- trimws(gsub("`", "", rownames(tab)))
  labs[labs == "Residuals"] <- "Error"
  .print_anova_df(tab, labs)
  cat(sprintf("\n Grand mean = %.3f   CV = %.2f%%\n", x$grand_mean, x$cv))
  cmp <- x$comparisons
  cmp[, 3:5] <- round(cmp[, 3:5], 3)
  cat("\n SE and CD (alpha = ", x$alpha, ")\n", sep = "")
  print(cmp, row.names = FALSE)
  lm_lab <- if (is.na(x$letters_method)) "" else
    paste0(" (", x$letters_method, ")")
  .print_means(x$mean_a, paste0("Means of ", f[["A"]], lm_lab), digits)
  .print_means(x$mean_b, paste0("Means of ", f[["B"]], lm_lab), digits)
  .print_means(x$interaction_means,
               paste0("Treatment combinations", lm_lab), digits)
  invisible(x)
}

#' @rdname print.aridagri_factorial
#' @export
print.aridagri_factorial3 <- function(x, digits = 2, ...) {
  f <- x$factors
  .print_info(paste0("THREE-FACTOR FACTORIAL (", x$layout, ")"),
              c(Response = x$response, `Factor A` = f[["A"]],
                `Factor B` = f[["B"]], `Factor C` = f[["C"]],
                Observations = nrow(stats::model.frame(x$model))))
  tab <- x$anova_table
  labs <- trimws(gsub("`", "", rownames(tab)))
  labs[labs == "Residuals"] <- "Error"
  .print_anova_df(tab, labs)
  cat(sprintf("\n Grand mean = %.3f   CV = %.2f%%\n", x$grand_mean, x$cv))
  cmp <- x$comparisons
  cmp[, 3:5] <- round(cmp[, 3:5], 3)
  cat("\n SE and CD (alpha = ", x$alpha, ")\n", sep = "")
  print(cmp, row.names = FALSE)
  for (i in 1:3) {
    .print_means(x$means[[i]], paste0("Means of ", names(x$means)[i]), digits)
  }
  cat("\n Two- and three-way tables are in $means.\n")
  invisible(x)
}


#' Latin Square Design ANOVA
#'
#' @description
#' Analysis of variance for a Latin square design with row and column
#' blocking, SEm, SEd, CD, CV, the efficiency of the Latin square relative
#' to a CRD and to an RCBD, and post-hoc letter groupings.
#'
#' @param data Data frame containing the data.
#' @param response Name of the response variable.
#' @param treatment Name of the treatment column.
#' @param row Name of the row column.
#' @param column Name of the column column.
#' @param posthoc Post-hoc method(s); see \code{\link{perform_posthoc}}.
#' @param alpha Significance level. Default 0.05.
#' @param verbose Logical. If TRUE (default), prints the report.
#' @param control Control level for Dunnett's test.
#'
#' @details
#' Relative efficiencies (Cochran and Cox, 1957; Gomez and Gomez, 1984)
#' use the row (\eqn{E_r}), column (\eqn{E_c}) and error (\eqn{E_e}) mean
#' squares: versus a CRD, \eqn{[E_r + E_c + (t-1)E_e] / [(t+1)E_e]};
#' versus an RCBD with rows as blocks, \eqn{[E_c + (t-1)E_e] / (t E_e)};
#' versus an RCBD with columns as blocks, \eqn{[E_r + (t-1)E_e] / (t E_e)},
#' each multiplied by the precision adjustment factor when the error
#' degrees of freedom are below 20.
#'
#' @return An object of class \code{aridagri_anova} with the components
#'   of \code{\link{anova_crd}} and \code{relative_efficiency} (a named
#'   vector, percent).
#'
#' @examples
#' data <- data.frame(
#'   row = rep(1:5, each = 5),
#'   col = rep(1:5, 5),
#'   treatment = c("A","B","C","D","E", "B","C","D","E","A",
#'                 "C","D","E","A","B", "D","E","A","B","C",
#'                 "E","A","B","C","D"),
#'   yield = c(1180, 1265, 1342, 1410, 1238, 1270, 1335, 1402, 1251, 1176,
#'             1340, 1415, 1245, 1190, 1258, 1405, 1232, 1178, 1262, 1339,
#'             1236, 1183, 1271, 1348, 1398)
#' )
#' anova_latin(data, response = "yield", treatment = "treatment",
#'             row = "row", column = "col")
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
anova_latin <- function(data, response, treatment, row, column,
                        posthoc = "lsd", alpha = 0.05, verbose = TRUE,
                        control = NULL) {
  .check_alpha(alpha)
  methods <- .resolve_methods(posthoc)
  d <- .prep_data(data, response, c(treatment, row, column), verbose)
  k <- nlevels(d[[treatment]])
  latin <- nlevels(d[[row]]) == k && nlevels(d[[column]]) == k &&
    all(table(d[[row]], d[[treatment]]) == 1) &&
    all(table(d[[column]], d[[treatment]]) == 1) && nrow(d) == k * k
  if (!latin) {
    warning("The layout is not a complete Latin square (each treatment once ",
            "in every row and column). Treatments are tested after rows ",
            "and columns and means are least-squares means.", call. = FALSE)
  }
  rhs <- paste(.bq(row), "+", .bq(column), "+", .bq(treatment))
  model <- stats::aov(.formula(response, rhs), data = d)
  anova_table <- .clean_rownames(stats::anova(model))
  if (anova_table[["Df"]][4] < 1) {
    stop("No error degrees of freedom.", call. = FALSE)
  }
  ms <- anova_table[["Mean Sq"]]
  ms_e <- ms[4]
  df_e <- anova_table[["Df"]][4]
  gm <- mean(d[[response]])
  n_i <- as.numeric(table(d[[treatment]]))
  r_h <- 1 / mean(1 / n_i)
  sem <- sqrt(ms_e / r_h)
  sed <- sqrt(2) * sem
  re <- c(CRD = NA_real_, RCBD_rows_as_blocks = NA_real_,
          RCBD_columns_as_blocks = NA_real_)
  if (latin) {
    er <- ms[1]
    ec <- ms[2]
    re[["CRD"]] <- 100 * (er + ec + (k - 1) * ms_e) / ((k + 1) * ms_e) *
      .re_adjust(df_e, k * (k - 1))
    re[["RCBD_rows_as_blocks"]] <- 100 * (ec + (k - 1) * ms_e) / (k * ms_e) *
      .re_adjust(df_e, (k - 1)^2)
    re[["RCBD_columns_as_blocks"]] <- 100 * (er + (k - 1) * ms_e) /
      (k * ms_e) * .re_adjust(df_e, (k - 1)^2)
  }
  lsm <- .ls_means_model(model, treatment)
  ph <- .posthoc_set(model, treatment, methods, alpha, control = control,
                     response = response, lsm = lsm)
  means_df <- .means_table(d, response, treatment, lsm, ph, methods)
  result <- list(
    design = "Latin Square Design",
    response = response,
    factors = c(treatment = treatment, row = row, column = column),
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
    relative_efficiency = re,
    treatment_means = means_df,
    letters_method = .letters_method(methods),
    posthoc = ph,
    diagnostics = check_assumptions(model, factor = treatment, verbose = FALSE),
    significance = .sig_code(anova_table[["Pr(>F)"]][3]),
    balanced = latin
  )
  class(result) <- c("aridagri_latin", "aridagri_anova", "list")
  if (verbose) print(result)
  invisible(result)
}

#' @rdname print.aridagri_anova
#' @export
print.aridagri_latin <- function(x, digits = 2, ...) {
  f <- x$factors
  .print_info("LATIN SQUARE DESIGN",
              c(Response = x$response, Treatment = f[["treatment"]],
                Rows = f[["row"]], Columns = f[["column"]],
                Treatments = nrow(x$treatment_means)))
  .print_anova_df(x$anova_table, c("Row", "Column", "Treatment", "Error"))
  st <- c(`Grand mean` = x$grand_mean, `CV (%)` = x$cv, SEm = x$se_mean,
          SEd = x$se_diff)
  st[[sprintf("CD (%g%%)", 100 * x$alpha)]] <- x$cd
  .print_stats_line(st, 3)
  if (!anyNA(x$relative_efficiency)) {
    cat(sprintf(paste0(" Relative efficiency (%%): vs CRD %.1f; vs RCBD ",
                       "(rows as blocks) %.1f; vs RCBD (columns as blocks) ",
                       "%.1f\n"),
                x$relative_efficiency[1], x$relative_efficiency[2],
                x$relative_efficiency[3]))
  }
  ttl <- if (is.na(x$letters_method)) "Treatment means" else
    sprintf("Treatment means (%s, alpha = %g)", x$letters_method, x$alpha)
  .print_means(x$treatment_means, ttl, digits)
  if (!is.null(x$posthoc$dunnett)) .print_posthoc(x$posthoc["dunnett"], digits)
  invisible(x)
}
