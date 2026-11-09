# =============================================================================
# SPLIT-PLOT FAMILY OF DESIGNS
# Package: aridagri
# =============================================================================

# Assemble the result list shared by the split-plot family.
.spd_result <- function(core, design, response, alpha, methods, class_name) {
  res <- list(
    design = design,
    response = response,
    anova_table = core$anova_table,
    grand_mean = core$grand_mean,
    cv = core$cv,
    error_ms = core$error_ms,
    error_df = core$error_df,
    alpha = alpha,
    se_cd = core$se_cd,
    comparisons = core$comparisons,
    factor_means = core$factor_means,
    interaction_means = core$interaction_means,
    interaction_tables = core$interaction_means,
    letters = core$letters,
    letters_method = .letters_method(methods),
    model = core$model,
    core = core[c("r", "levels", "lev_of")]
  )
  K <- length(core$error_ms)
  for (k in seq_len(K)) {
    res[[paste0("ms_error_", letters[k])]] <- core$error_ms[[k]]
    res[[paste0("df_error_", letters[k])]] <- core$error_df[[k]]
  }
  class(res) <- c(class_name, "aridagri_split", "list")
  res
}

#' Split-Plot Design ANOVA
#'
#' @description
#' Analysis of variance for a split-plot design in randomized complete
#' blocks, with the main-plot factor tested against Error(a) and the
#' sub-plot factor and interaction against Error(b). Gives CV(a), CV(b),
#' SEm, SEd and CD for all four standard comparisons and letter groupings
#' for main-plot means, sub-plot means and each factor within the levels
#' of the other.
#'
#' @param data Data frame containing the experimental data.
#' @param response Name of the response variable.
#' @param main_plot Name of the main-plot factor (A).
#' @param sub_plot Name of the sub-plot factor (B).
#' @param replication Name of the replication (block) column.
#' @param posthoc Method for the letters (see
#'   \code{\link{perform_posthoc}}), or \code{"none"}. Default
#'   \code{"lsd"}.
#' @param alpha Significance level for CD and letters. Default 0.05.
#' @param verbose Logical. If TRUE (default), prints the report.
#'
#' @details
#' Standard errors of difference (Gomez and Gomez, 1984, table 3.10):
#' main-plot means \eqn{\sqrt{2E_a/(rb)}}, sub-plot means
#' \eqn{\sqrt{2E_b/(ra)}}, sub-plot means at the same main-plot level
#' \eqn{\sqrt{2E_b/r}}, and main-plot means at the same or different
#' sub-plot level \eqn{\sqrt{2[(b-1)E_b + E_a]/(rb)}}. For the last
#' comparison the CD uses the weighted t value
#' \eqn{t' = [(b-1)E_b t_b + E_a t_a] / [(b-1)E_b + E_a]}.
#' The data must be complete (every combination once in each
#' replication).
#'
#' @return An object of class \code{aridagri_spd} with the corrected
#'   \code{anova_table}, \code{cv_a}, \code{cv_b}, error mean squares and
#'   df, \code{se_main}, \code{cd_main}, \code{se_sub}, \code{cd_sub},
#'   \code{se_sub_same}, \code{cd_sub_same}, \code{se_main_same_sub},
#'   \code{cd_main_same_sub}, \code{se_cd} (all comparisons),
#'   \code{main_means} and \code{sub_means} (with letters),
#'   \code{interaction_means}, \code{letters} (each factor within the
#'   levels of the other) and \code{model}.
#'
#' @references
#' Gomez, K.A. and Gomez, A.A. (1984). \emph{Statistical Procedures for
#' Agricultural Research}, 2nd ed. John Wiley and Sons, New York.
#'
#' @examples
#' data(clusterbean_spd)
#' res <- anova_spd(clusterbean_spd, response = "seed_yield",
#'                  main_plot = "irrigation", sub_plot = "nitrogen",
#'                  replication = "rep")
#' res$se_cd
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
anova_spd <- function(data, response, main_plot, sub_plot, replication,
                      posthoc = "lsd", alpha = 0.05, verbose = TRUE) {
  .check_alpha(alpha)
  methods <- setdiff(.resolve_methods(posthoc), "dunnett")
  groups <- list(main = main_plot, sub = sub_plot)
  core <- .aridagri_spd_anova(data, response, replication, groups, alpha,
                              methods, verbose)
  res <- .spd_result(core, "Split Plot Design", response, alpha, methods,
                     "aridagri_spd")
  sc <- core$se_cd
  pick <- function(eff, s) sc[sc$Effect == eff & sc$stratum == s, ]
  ra <- pick(main_plot, 1)
  rb <- pick(sub_plot, 2)
  rab_b <- pick(paste(main_plot, "x", sub_plot), 2)
  rab_a <- pick(paste(main_plot, "x", sub_plot), 1)
  d <- core$data
  im <- stats::aggregate(d[[response]], by = list(d[[main_plot]],
                                                  d[[sub_plot]]), FUN = mean)
  names(im) <- c("MainPlot", "SubPlot", "Mean")
  res$cv_a <- core$cv[["main"]]
  res$cv_b <- core$cv[["sub"]]
  res$se_main <- ra$SEm
  res$cd_main <- ra$CD
  res$se_sub <- rb$SEm
  res$cd_sub <- rb$CD
  res$se_sub_same <- rab_b$SEm
  res$cd_sub_same <- rab_b$CD
  res$se_main_same_sub <- rab_a$SEm
  res$cd_main_same_sub <- rab_a$CD
  res$main_means <- core$factor_means[[main_plot]]
  res$sub_means <- core$factor_means[[sub_plot]]
  res$interaction_means <- im
  p <- core$anova_table[["Pr(>F)"]]
  rn <- rownames(core$anova_table)
  res$significance <- list(
    main = .sig_code(p[rn == main_plot]),
    sub = .sig_code(p[rn == sub_plot]),
    interaction = .sig_code(p[rn %in% c(paste0(main_plot, ":", sub_plot),
                                        paste0(sub_plot, ":", main_plot))]))
  res$factors <- c(main_plot = main_plot, sub_plot = sub_plot)
  if (verbose) print(res)
  invisible(res)
}


#' Split-Plot Design with Two Factors in the Main Plot
#'
#' @description
#' Split-plot design where the main plots carry the factorial combination
#' of two factors (A x B) and the sub-plots carry factor C. A, B and A x B
#' are tested against Error(a); C and every interaction with C against
#' Error(b).
#'
#' @param data Data frame containing the experimental data.
#' @param response Name of the response variable.
#' @param main_factor1 First main-plot factor (A).
#' @param main_factor2 Second main-plot factor (B).
#' @param sub_plot Sub-plot factor (C).
#' @param replication Name of the replication column.
#' @param posthoc Method for the letters (default \code{"lsd"}).
#' @param alpha Significance level. Default 0.05.
#' @param verbose Logical. If TRUE (default), prints the report.
#'
#' @return An object of class \code{aridagri_spd_ab} with the corrected
#'   \code{anova_table}, \code{cv}, error mean squares and df,
#'   \code{se_cd} (SEm, SEd and CD for every main effect and interaction
#'   comparison), \code{factor_means} (with letters), \code{letters},
#'   \code{interaction_means} and \code{model}.
#'
#' @examples
#' data <- expand.grid(rep = 1:3, irrigation = c("I1", "I2"),
#'                     variety = c("V1", "V2", "V3"),
#'                     nitrogen = c("N0", "N1", "N2", "N3"))
#' set.seed(3)
#' data$yield <- rnorm(nrow(data), 1200, 80)
#' res <- anova_spd_ab_main(data, response = "yield",
#'                          main_factor1 = "irrigation",
#'                          main_factor2 = "variety",
#'                          sub_plot = "nitrogen", replication = "rep",
#'                          verbose = FALSE)
#' res$anova_table
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
anova_spd_ab_main <- function(data, response, main_factor1, main_factor2,
                              sub_plot, replication, posthoc = "lsd",
                              alpha = 0.05, verbose = TRUE) {
  .check_alpha(alpha)
  methods <- setdiff(.resolve_methods(posthoc), "dunnett")
  groups <- list(main = c(main_factor1, main_factor2), sub = sub_plot)
  core <- .aridagri_spd_anova(data, response, replication, groups, alpha,
                              methods, verbose)
  res <- .spd_result(core, "Split-Plot (A x B main, C sub)", response, alpha,
                     methods, "aridagri_spd_ab")
  if (verbose) print(res)
  invisible(res)
}


#' Split-Plot Design with Two Factors in the Sub-Plot
#'
#' @description
#' Split-plot design where the main plots carry factor C and the sub-plots
#' carry the factorial combination of two factors (A x B). C is tested
#' against Error(a); A, B, A x B and every interaction with C against
#' Error(b).
#'
#' @param data Data frame containing the experimental data.
#' @param response Name of the response variable.
#' @param main_plot Main-plot factor (C).
#' @param sub_factor1 First sub-plot factor (A).
#' @param sub_factor2 Second sub-plot factor (B).
#' @param replication Name of the replication column.
#' @param posthoc Method for the letters (default \code{"lsd"}).
#' @param alpha Significance level. Default 0.05.
#' @param verbose Logical. If TRUE (default), prints the report.
#'
#' @return An object of class \code{aridagri_spd_cab}; see
#'   \code{\link{anova_spd_ab_main}} for the components.
#'
#' @examples
#' data <- expand.grid(rep = 1:3, irrigation = c("I1", "I2", "I3"),
#'                     variety = c("V1", "V2"),
#'                     nitrogen = c("N1", "N2", "N3", "N4"))
#' set.seed(4)
#' data$yield <- rnorm(nrow(data), 1200, 80)
#' res <- anova_spd_c_main_ab_sub(data, response = "yield",
#'                                main_plot = "irrigation",
#'                                sub_factor1 = "variety",
#'                                sub_factor2 = "nitrogen",
#'                                replication = "rep", verbose = FALSE)
#' res$anova_table
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
anova_spd_c_main_ab_sub <- function(data, response, main_plot, sub_factor1,
                                    sub_factor2, replication, posthoc = "lsd",
                                    alpha = 0.05, verbose = TRUE) {
  .check_alpha(alpha)
  methods <- setdiff(.resolve_methods(posthoc), "dunnett")
  groups <- list(main = main_plot, sub = c(sub_factor1, sub_factor2))
  core <- .aridagri_spd_anova(data, response, replication, groups, alpha,
                              methods, verbose)
  res <- .spd_result(core, "Split-Plot (C main, A x B sub)", response, alpha,
                     methods, "aridagri_spd_cab")
  if (verbose) print(res)
  invisible(res)
}


#' Split-Plot Design with Factorial Main and Sub-Plots
#'
#' @description
#' Split-plot design where the main plots carry A x B and the sub-plots
#' carry C x D. Main-plot sources are tested against Error(a), all other
#' sources against Error(b).
#'
#' @param data Data frame containing the experimental data.
#' @param response Name of the response variable.
#' @param main_factor1 First main-plot factor (A).
#' @param main_factor2 Second main-plot factor (B).
#' @param sub_factor1 First sub-plot factor (C).
#' @param sub_factor2 Second sub-plot factor (D).
#' @param replication Name of the replication column.
#' @param verbose Logical. If TRUE (default), prints the report.
#' @param posthoc Method for the letters (default \code{"lsd"}).
#' @param alpha Significance level. Default 0.05.
#'
#' @return An object of class \code{aridagri_spd_abcd}; see
#'   \code{\link{anova_spd_ab_main}} for the components.
#'
#' @examples
#' data <- expand.grid(rep = 1:3, irrigation = c("I1", "I2"),
#'                     tillage = c("CT", "ZT"), variety = c("V1", "V2"),
#'                     nitrogen = c("N1", "N2", "N3"))
#' set.seed(5)
#' data$yield <- rnorm(nrow(data), 1200, 80)
#' res <- anova_spd_ab_cd(data, response = "yield",
#'                        main_factor1 = "irrigation", main_factor2 = "tillage",
#'                        sub_factor1 = "variety", sub_factor2 = "nitrogen",
#'                        replication = "rep", verbose = FALSE)
#' res$anova_table
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
anova_spd_ab_cd <- function(data, response, main_factor1, main_factor2,
                            sub_factor1, sub_factor2, replication,
                            verbose = TRUE, posthoc = "lsd", alpha = 0.05) {
  .check_alpha(alpha)
  methods <- setdiff(.resolve_methods(posthoc), "dunnett")
  groups <- list(main = c(main_factor1, main_factor2),
                 sub = c(sub_factor1, sub_factor2))
  core <- .aridagri_spd_anova(data, response, replication, groups, alpha,
                              methods, verbose)
  res <- .spd_result(core, "Split-Plot (A x B main, C x D sub)", response,
                     alpha, methods, "aridagri_spd_abcd")
  if (verbose) print(res)
  invisible(res)
}


#' Split-Split-Plot Design ANOVA
#'
#' @description
#' Analysis of variance for a split-split-plot design with main-plot
#' sources tested against Error(a), sub-plot sources against Error(b) and
#' sub-sub-plot sources against Error(c). Gives CVs for the three strata,
#' SEm, SEd and CD for every type of comparison (including those that mix
#' error terms) and letter groupings.
#'
#' @param data Data frame containing the experimental data.
#' @param response Name of the response variable.
#' @param main_plot Name of the main-plot factor (A).
#' @param sub_plot Name of the sub-plot factor (B).
#' @param sub_sub_plot Name of the sub-sub-plot factor (C).
#' @param replication Name of the replication column.
#' @param verbose Logical. If TRUE (default), prints the report.
#' @param posthoc Method for the letters (default \code{"lsd"}).
#' @param alpha Significance level. Default 0.05.
#'
#' @return An object of class \code{aridagri_sspd} with the corrected
#'   \code{anova_table}, \code{cv}, \code{ms_error_a/b/c},
#'   \code{df_error_a/b/c}, \code{se_cd}, \code{factor_means},
#'   \code{letters}, \code{interaction_means} and \code{model}.
#'
#' @references
#' Gomez, K.A. and Gomez, A.A. (1984). \emph{Statistical Procedures for
#' Agricultural Research}, 2nd ed. John Wiley and Sons, New York.
#'
#' @examples
#' data <- expand.grid(rep = 1:3, A = c("A1", "A2"), B = c("B1", "B2"),
#'                     C = c("C1", "C2", "C3"))
#' set.seed(6)
#' data$yield <- rnorm(nrow(data), 1200, 80)
#' res <- anova_sspd(data, response = "yield", main_plot = "A",
#'                   sub_plot = "B", sub_sub_plot = "C", replication = "rep")
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
anova_sspd <- function(data, response, main_plot, sub_plot, sub_sub_plot,
                       replication, verbose = TRUE, posthoc = "lsd",
                       alpha = 0.05) {
  .check_alpha(alpha)
  methods <- setdiff(.resolve_methods(posthoc), "dunnett")
  groups <- list(main = main_plot, sub = sub_plot, subsub = sub_sub_plot)
  core <- .aridagri_spd_anova(data, response, replication, groups, alpha,
                              methods, verbose)
  res <- .spd_result(core, "Split-Split-Plot Design", response, alpha,
                     methods, "aridagri_sspd")
  if (verbose) print(res)
  invisible(res)
}


#' Print Method for Split-Plot Family Results
#'
#' @param x A result from \code{\link{anova_spd}} or a related split-plot
#'   function.
#' @param digits Number of decimal places for means. Default 2.
#' @param ... Unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.aridagri_split <- function(x, digits = 2, ...) {
  core <- c(x[c("anova_table", "grand_mean", "error_ms", "error_df", "cv",
                "se_cd", "factor_means")], x$core)
  core$interaction_means <- x$interaction_tables
  .aridagri_print_spd(core, x$design, x$response, x$letters_method, x$alpha,
                      digits)
  invisible(x)
}


#' Combined Split-Plot Analysis over Environments
#'
#' @description
#' Combined (pooled) analysis of a split-plot experiment repeated in
#' several environments (years, seasons or locations). Environments are
#' fixed and replications random within environments: the environment
#' effect is tested against replications within environments, main-plot
#' sources against the pooled Error(a) and sub-plot sources against the
#' pooled Error(b). SEm, SEd and CD are given for every comparison.
#'
#' @param data Data frame containing the combined data.
#' @param response Name of the response variable.
#' @param main_plot Name of the main-plot factor.
#' @param sub_plot Name of the sub-plot factor.
#' @param environment Name of the environment column.
#' @param replication Name of the replication column (nested within
#'   environment; labels may repeat across environments).
#' @param verbose Logical. If TRUE (default), prints the report.
#' @param posthoc Method for the letters (default \code{"lsd"}).
#' @param alpha Significance level. Default 0.05.
#'
#' @return An object of class \code{aridagri_spd_pooled} with
#'   \code{anova_table}, \code{cv} (by stratum), \code{se_cd},
#'   \code{factor_means} (with letters), \code{main_means},
#'   \code{sub_means}, \code{env_means} and \code{model}.
#'
#' @examples
#' data <- expand.grid(year = c("Y1", "Y2", "Y3"), rep = 1:3,
#'                     irrigation = c("I1", "I2", "I3"),
#'                     variety = c("V1", "V2", "V3", "V4"))
#' set.seed(7)
#' data$yield <- rnorm(nrow(data), 1200, 90)
#' res <- anova_spd_pooled(data, response = "yield", main_plot = "irrigation",
#'                         sub_plot = "variety", environment = "year",
#'                         replication = "rep", verbose = FALSE)
#' res$anova_table
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
anova_spd_pooled <- function(data, response, main_plot, sub_plot,
                             environment, replication, verbose = TRUE,
                             posthoc = "lsd", alpha = 0.05) {
  .check_alpha(alpha)
  methods <- setdiff(.resolve_methods(posthoc), "dunnett")
  data <- .renumber_reps(data, environment, replication)
  core <- .aridagri_spd_pooled_anova(data, response, environment, replication,
                                     split_factors = c(main_plot, sub_plot),
                                     alpha = alpha, methods = methods,
                                     verbose = verbose)
  res <- list(
    design = "Pooled Split Plot Design",
    response = response,
    anova_table = core$anova_table,
    model = core$model,
    grand_mean = core$grand_mean,
    cv = core$cv,
    alpha = alpha,
    se_cd = core$se_cd,
    factor_means = core$factor_means,
    main_means = core$factor_means[[main_plot]],
    sub_means = core$factor_means[[sub_plot]],
    env_means = core$factor_means[[environment]],
    letters_method = .letters_method(methods),
    core = core[c("r", "levels", "environment", "split_factors",
                  "error_ms", "error_df")]
  )
  class(res) <- c("aridagri_spd_pooled", "aridagri_split_pooled", "list")
  if (verbose) print(res)
  invisible(res)
}


#' Combined Split-Split-Plot Analysis over Environments
#'
#' @description
#' Combined (pooled) analysis of a split-split-plot experiment repeated in
#' several environments, with environments fixed and replications random
#' within environments. See \code{\link{anova_spd_pooled}}.
#'
#' @param data Data frame containing the combined data.
#' @param response Name of the response variable.
#' @param main_plot Name of the main-plot factor.
#' @param sub_plot Name of the sub-plot factor.
#' @param sub_sub_plot Name of the sub-sub-plot factor.
#' @param environment Name of the environment column.
#' @param replication Name of the replication column.
#' @param verbose Logical. If TRUE (default), prints the report.
#' @param posthoc Method for the letters (default \code{"lsd"}).
#' @param alpha Significance level. Default 0.05.
#'
#' @return An object of class \code{aridagri_sspd_pooled}; see
#'   \code{\link{anova_spd_pooled}}. Also contains \code{subsub_means}.
#'
#' @examples
#' data <- expand.grid(env = c("E1", "E2"), rep = 1:3, A = c("A1", "A2"),
#'                     B = c("B1", "B2"), C = c("C1", "C2"))
#' set.seed(8)
#' data$yield <- rnorm(nrow(data), 1200, 90)
#' res <- anova_sspd_pooled(data, "yield", "A", "B", "C", "env", "rep",
#'                          verbose = FALSE)
#' res$anova_table
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram
#'   Choudhary, Hemant Poonia
anova_sspd_pooled <- function(data, response, main_plot, sub_plot,
                              sub_sub_plot, environment, replication,
                              verbose = TRUE, posthoc = "lsd", alpha = 0.05) {
  .check_alpha(alpha)
  methods <- setdiff(.resolve_methods(posthoc), "dunnett")
  data <- .renumber_reps(data, environment, replication)
  fs <- c(main_plot, sub_plot, sub_sub_plot)
  core <- .aridagri_spd_pooled_anova(data, response, environment, replication,
                                     split_factors = fs, alpha = alpha,
                                     methods = methods, verbose = verbose)
  res <- list(
    design = "Pooled Split-Split Plot Design",
    response = response,
    anova_table = core$anova_table,
    model = core$model,
    grand_mean = core$grand_mean,
    cv = core$cv,
    alpha = alpha,
    se_cd = core$se_cd,
    factor_means = core$factor_means,
    main_means = core$factor_means[[main_plot]],
    sub_means = core$factor_means[[sub_plot]],
    subsub_means = core$factor_means[[sub_sub_plot]],
    env_means = core$factor_means[[environment]],
    letters_method = .letters_method(methods),
    core = core[c("r", "levels", "environment", "split_factors",
                  "error_ms", "error_df")]
  )
  class(res) <- c("aridagri_sspd_pooled", "aridagri_split_pooled", "list")
  if (verbose) print(res)
  invisible(res)
}

#' @rdname print.aridagri_split
#' @export
print.aridagri_split_pooled <- function(x, digits = 2, ...) {
  core <- c(x[c("anova_table", "grand_mean", "cv", "se_cd", "factor_means")],
            x$core)
  .aridagri_print_spd_pooled(core, toupper(x$design), x$letters_method,
                             x$alpha, digits)
  invisible(x)
}

# Number replications 1..r within each environment so that labels that are
# unique across environments (1-9 for three environments) also work.
.renumber_reps <- function(data, environment, replication) {
  .check_cols(data, c(environment, replication))
  e <- as.character(data[[environment]])
  rp <- as.character(data[[replication]])
  new <- integer(length(rp))
  for (lv in unique(e)) {
    i <- which(e == lv)
    new[i] <- as.integer(factor(rp[i], levels = unique(rp[i])))
  }
  data[[replication]] <- new
  data
}
