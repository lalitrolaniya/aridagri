#' ============================================================================
#' SPLIT PLOT DESIGN ANALYSIS FUNCTIONS (All Variations)
#' Package: aridagri
#' Author: Lalit Kumar Rolaniya
#' ICAR-Indian Institute of Pulses Research, Regional Centre, Bikaner
#' ============================================================================

# Internal: assemble a correct split-plot ANOVA table from an aov() model that
# was fitted with Error() strata. Each effect keeps the F-test computed against
# its own stratum residual (so main-plot factors are tested against the pooled
# main-plot error, not the sub-plot residual). Residual rows are relabelled, in
# stratum order, as Replication, Error(a), Error(b), Error(c), ... The named
# list of error mean squares and df is attached as attribute "error_terms".
# Not exported.
.spd_error_table <- function(fit) {
  s <- summary(fit)
  rows <- list()
  err  <- list()
  for (k in seq_along(s)) {
    tab <- s[[k]][[1]]
    rn  <- rownames(tab)
    res_label <- if (k == 1L) "Replication" else paste0("Error(", letters[k - 1L], ")")
    cn  <- colnames(tab)
    for (i in seq_len(nrow(tab))) {
      lab <- trimws(rn[i])
      if (lab == "Residuals") lab <- res_label
      rows[[length(rows) + 1L]] <- data.frame(
        Df        = tab[i, "Df"],
        `Sum Sq`  = tab[i, "Sum Sq"],
        `Mean Sq` = tab[i, "Mean Sq"],
        `F value` = if ("F value" %in% cn) tab[i, "F value"] else NA_real_,
        `Pr(>F)`  = if ("Pr(>F)"  %in% cn) tab[i, "Pr(>F)"]  else NA_real_,
        row.names = lab, check.names = FALSE, stringsAsFactors = FALSE)
    }
    ri <- which(trimws(rn) == "Residuals")
    if (length(ri)) err[[res_label]] <- list(ms = tab[ri, "Mean Sq"], df = tab[ri, "Df"])
  }
  out <- do.call(rbind, rows)
  attr(out, "error_terms") <- err
  out
}

#' Split Plot Design ANOVA (Standard)
#'
#' @description
#' Performs complete ANOVA for Split Plot Design with proper error terms for 
#' main plot and sub-plot factors. Includes all standard post-hoc comparisons.
#'
#' @param data Data frame containing the experimental data
#' @param response Name of the response variable
#' @param main_plot Name of main plot factor
#' @param sub_plot Name of sub-plot factor
#' @param replication Name of replication/block factor
#' @param posthoc Post-hoc test method
#' @param alpha Significance level
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return List containing ANOVA table, means, and post-hoc results
#'
#' @examples
#' data <- data.frame(
#'   rep = rep(1:3, each = 12),
#'   irrigation = rep(rep(c("I1", "I2", "I3"), each = 4), 3),
#'   variety = rep(c("V1", "V2", "V3", "V4"), 9),
#'   yield = rnorm(36, 1200, 150)
#' )
#' anova_spd(data, response = "yield", main_plot = "irrigation", 
#'           sub_plot = "variety", replication = "rep")
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
anova_spd <- function(data, response, main_plot, sub_plot, replication, 
                       posthoc = "lsd", alpha = 0.05,
                            verbose = TRUE) {
  
  # Convert to factors
  data[[main_plot]] <- as.factor(data[[main_plot]])
  data[[sub_plot]] <- as.factor(data[[sub_plot]])
  data[[replication]] <- as.factor(data[[replication]])
  
  # Get dimensions
  a <- nlevels(data[[main_plot]])      # Main plot levels
  b <- nlevels(data[[sub_plot]])       # Sub-plot levels
  r <- nlevels(data[[replication]])    # Replications
  N <- nrow(data)
  
  # Check for balanced design
  if (N != a * b * r) {
    warning("Unbalanced design detected: expected ", a * b * r, " observations (", 
            a, " main-plot x ", b, " sub-plot x ", r, " reps) but found ", N, ". ",
            "Results assume balanced data and may be unreliable. ",
            "Consider using lme4::lmer() for unbalanced data.", call. = FALSE)
  }
  
  # Fit full model
  formula_spd <- as.formula(paste(response, "~", replication, "+", 
                                   main_plot, "+", replication, ":", main_plot, "+",
                                   sub_plot, "+", main_plot, ":", sub_plot))
  
  model <- aov(formula_spd, data = data)
  anova_full <- anova(model)

  # Extract sums of squares BY TERM NAME. R sorts main-effect terms before
  # interaction terms in the ANOVA table regardless of formula order, so fixed
  # positions [3]/[4] would swap Error(a) (rep:main) with the Sub-plot effect.
  rn  <- rownames(anova_full)
  SSq <- anova_full$`Sum Sq`
  SS_rep         <- SSq[which(rn == replication)]
  SS_main        <- SSq[which(rn == main_plot)]
  SS_error_a     <- SSq[which(rn %in% c(paste0(replication, ":", main_plot),
                                        paste0(main_plot, ":", replication)))]
  SS_sub         <- SSq[which(rn == sub_plot)]
  SS_interaction <- SSq[which(rn %in% c(paste0(main_plot, ":", sub_plot),
                                        paste0(sub_plot, ":", main_plot)))]
  SS_error_b     <- SSq[length(SSq)]        # Residuals is always the last row
  SS_total <- sum(SSq)
  
  # Degrees of freedom
  df_rep <- r - 1
  df_main <- a - 1
  df_error_a <- (r - 1) * (a - 1)
  df_sub <- b - 1
  df_interaction <- (a - 1) * (b - 1)
  df_error_b <- a * (r - 1) * (b - 1)
  df_total <- N - 1
  
  # Mean squares
  MS_rep <- SS_rep / df_rep
  MS_main <- SS_main / df_main
  MS_error_a <- SS_error_a / df_error_a
  MS_sub <- SS_sub / df_sub
  MS_interaction <- SS_interaction / df_interaction
  MS_error_b <- SS_error_b / df_error_b
  
  # F-values and p-values
  F_main <- MS_main / MS_error_a
  F_sub <- MS_sub / MS_error_b
  F_interaction <- MS_interaction / MS_error_b
  
  p_main <- 1 - pf(F_main, df_main, df_error_a)
  p_sub <- 1 - pf(F_sub, df_sub, df_error_b)
  p_interaction <- 1 - pf(F_interaction, df_interaction, df_error_b)
  
  # Statistics
  grand_mean <- mean(data[[response]], na.rm = TRUE)
  cv_a <- sqrt(MS_error_a) / grand_mean * 100
  cv_b <- sqrt(MS_error_b) / grand_mean * 100
  
  # SE and CD calculations
  # For Main Plot
  se_main <- sqrt(MS_error_a / (b * r))
  cd_main_5 <- qt(0.975, df_error_a) * sqrt(2 * MS_error_a / (b * r))
  
  # For Sub Plot
  se_sub <- sqrt(MS_error_b / (a * r))
  cd_sub_5 <- qt(0.975, df_error_b) * sqrt(2 * MS_error_b / (a * r))
  
  # For comparing sub-plot means at same main plot level
  se_sub_same <- sqrt(MS_error_b / r)
  cd_sub_same_5 <- qt(0.975, df_error_b) * sqrt(2 * MS_error_b / r)
  
  # For comparing main plot means at same or different sub-plot levels
  # Using pooled error
  MS_pooled <- ((df_error_a * MS_error_a) + (a - 1) * (b - 1) * MS_error_b) / 
               (df_error_a + (a - 1) * (b - 1))
  se_main_sub <- sqrt((MS_error_a + (b - 1) * MS_error_b) / (b * r))
  
  # Print ANOVA table
  if (verbose) {
    cat("\n")
    cat("\n")
    cat("                    SPLIT PLOT DESIGN - ANOVA                                 \n")
    cat("\n")
    cat(" Response Variable  :", sprintf("%-55s", response), "\n")
    cat(" Main Plot Factor (A):", sprintf("%-53s", main_plot), "\n")
    cat(" Sub-Plot Factor (B) :", sprintf("%-53s", sub_plot), "\n")
    cat(" Replications        :", sprintf("%-53d", r), "\n")
    cat(" Main Plot Levels    :", sprintf("%-53d", a), "\n")
    cat(" Sub-Plot Levels     :", sprintf("%-53d", b), "\n")
    cat(" Total Observations  :", sprintf("%-53d", N), "\n")
    cat("\n")

    cat("\n\n")
    cat("                           ANALYSIS OF VARIANCE                               \n")
    cat("\n")
    cat(" Source                   df        SS           MS         F     Pr>F   \n")
    cat("\n")
    cat(sprintf(" Replication            %6d  %11.2f  %11.2f                 \n",
                df_rep, SS_rep, MS_rep))
    cat(sprintf(" Main Plot (A)          %6d  %11.2f  %11.2f  %5.2f  %7.4f \n",
                df_main, SS_main, MS_main, F_main, p_main))
    cat(sprintf(" Error (a)              %6d  %11.2f  %11.2f                 \n",
                df_error_a, SS_error_a, MS_error_a))
    cat("\n")
    cat(sprintf(" Sub-Plot (B)           %6d  %11.2f  %11.2f  %5.2f  %7.4f \n",
                df_sub, SS_sub, MS_sub, F_sub, p_sub))
    cat(sprintf(" A  B                  %6d  %11.2f  %11.2f  %5.2f  %7.4f \n",
                df_interaction, SS_interaction, MS_interaction, F_interaction, p_interaction))
    cat(sprintf(" Error (b)              %6d  %11.2f  %11.2f                 \n",
                df_error_b, SS_error_b, MS_error_b))
    cat("\n")
    cat(sprintf(" Total                  %6d  %11.2f                              \n",
                df_total, SS_total))
    cat("\n")
  }
  
  # Significance symbols
  sig_main <- get_sig_symbol(p_main)
  sig_sub <- get_sig_symbol(p_sub)
  sig_int <- get_sig_symbol(p_interaction)
  
  if (verbose) {
    cat("\n\n")
    cat("                           SUMMARY STATISTICS                                 \n")
    cat("\n")
    cat(sprintf(" Grand Mean                    : %10.2f                                   \n", grand_mean))
    cat(sprintf(" CV (a) for Main Plot          : %10.2f%%                                  \n", cv_a))
    cat(sprintf(" CV (b) for Sub-Plot           : %10.2f%%                                  \n", cv_b))
    cat("\n")
    cat(sprintf(" Main Plot (A)                 : %10s                                   \n", sig_main))
    cat(sprintf(" Sub-Plot (B)                  : %10s                                   \n", sig_sub))
    cat(sprintf(" Interaction (AB)             : %10s                                   \n", sig_int))
    cat("\n")

    cat("\n\n")
    cat("                    SE AND CD VALUES FOR COMPARISONS                          \n")
    cat("\n")
    cat(" For MAIN PLOT factor (A):                                                    \n")
    cat(sprintf("   SE(m) = %8.3f    SE(d) = %8.3f    CD(5%%) = %8.3f                 \n",
                se_main, se_main * sqrt(2), cd_main_5))
    cat("\n")
    cat(" For SUB-PLOT factor (B):                                                     \n")
    cat(sprintf("   SE(m) = %8.3f    SE(d) = %8.3f    CD(5%%) = %8.3f                 \n",
                se_sub, se_sub * sqrt(2), cd_sub_5))
    cat("\n")
    cat(" For B at same level of A:                                                    \n")
    cat(sprintf("   SE(m) = %8.3f    SE(d) = %8.3f    CD(5%%) = %8.3f                 \n",
                se_sub_same, se_sub_same * sqrt(2), cd_sub_same_5))
    cat("\n")
    cat(" For A at same or different level of B:                                       \n")
    cat(sprintf("   SE(m) = %8.3f    SE(d) = %8.3f                                     \n",
                se_main_sub, se_main_sub * sqrt(2)))
    cat("\n")
  }
  
  # Factor means
  main_means <- aggregate(data[[response]], by = list(data[[main_plot]]), 
                           FUN = mean, na.rm = TRUE)
  names(main_means) <- c("Level", "Mean")
  
  sub_means <- aggregate(data[[response]], by = list(data[[sub_plot]]), 
                          FUN = mean, na.rm = TRUE)
  names(sub_means) <- c("Level", "Mean")
  
  interaction_means <- aggregate(data[[response]], 
                                  by = list(data[[main_plot]], data[[sub_plot]]), 
                                  FUN = mean, na.rm = TRUE)
  names(interaction_means) <- c("MainPlot", "SubPlot", "Mean")
  
  if (verbose) {
    cat("\n\n")
    cat("                           MAIN PLOT MEANS                                    \n")
    cat("\n")
    for (i in 1:nrow(main_means)) {
      cat(sprintf(" %-31s  %42.2f \n", 
                  as.character(main_means$Level[i]), main_means$Mean[i]))
    }
    cat("\n")

    cat("\n\n")
    cat("                           SUB-PLOT MEANS                                     \n")
    cat("\n")
    for (i in 1:nrow(sub_means)) {
      cat(sprintf(" %-31s  %42.2f \n", 
                  as.character(sub_means$Level[i]), sub_means$Mean[i]))
    }
    cat("\n")
  }
  
  # Interaction table
  int_table <- reshape(interaction_means, idvar = "MainPlot", 
                        timevar = "SubPlot", direction = "wide")
  
  if (verbose) {
    cat("\n\n")
    cat("                        INTERACTION MEANS (A  B)                             \n")
    cat("\n")
    print(int_table, row.names = FALSE)
  }
  
  # Return results
  result <- list(
    design = "Split Plot Design",
    anova_table = anova_full,
    model = model,
    grand_mean = grand_mean,
    cv_a = cv_a,
    cv_b = cv_b,
    ms_error_a = MS_error_a,
    ms_error_b = MS_error_b,
    df_error_a = df_error_a,
    df_error_b = df_error_b,
    se_main = se_main,
    cd_main = cd_main_5,
    se_sub = se_sub,
    cd_sub = cd_sub_5,
    se_sub_same = se_sub_same,
    cd_sub_same = cd_sub_same_5,
    main_means = main_means,
    sub_means = sub_means,
    interaction_means = interaction_means,
    significance = list(main = sig_main, sub = sig_sub, interaction = sig_int)
  )
  
  class(result) <- c("aridagri_spd", "list")
  return(invisible(result))
}


#' Split Plot Design with AB in Main Plot
#'
#' @description
#' Performs ANOVA for Split Plot Design where main plot contains factorial 
#' combination of two factors (AB) and sub-plot contains factor C.
#' Common in irrigation  variety as main plot and nitrogen as sub-plot.
#'
#' @param data Data frame containing the experimental data
#' @param response Name of the response variable
#' @param main_factor1 First factor in main plot (A)
#' @param main_factor2 Second factor in main plot (B)
#' @param sub_plot Sub-plot factor (C)
#' @param replication Name of replication factor
#' @param posthoc Post-hoc test method
#' @param alpha Significance level
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return List containing ANOVA table and means
#'
#' @details
#' Design structure:
#' \itemize{
#'   \item Main Plot: A  B factorial
#'   \item Sub-Plot: C
#'   \item Error (a): For testing A, B, and AB
#'   \item Error (b): For testing C and all interactions with C
#' }
#'
#' @examples
#' data <- data.frame(
#'   rep = rep(1:3, each = 24),
#'   irrigation = rep(rep(c("I1", "I2"), each = 12), 3),
#'   variety = rep(rep(c("V1", "V2", "V3"), each = 4), 6),
#'   nitrogen = rep(c("N0", "N1", "N2", "N3"), 18),
#'   yield = rnorm(72, 1200, 150)
#' )
#' anova_spd_ab_main(data, response = "yield", 
#'                   main_factor1 = "irrigation", main_factor2 = "variety",
#'                   sub_plot = "nitrogen", replication = "rep")
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
anova_spd_ab_main <- function(data, response, main_factor1, main_factor2, 
                               sub_plot, replication, posthoc = "lsd", alpha = 0.05,
                            verbose = TRUE) {
  groups <- list(main = c(main_factor1, main_factor2), sub = sub_plot)
  core <- .aridagri_spd_anova(data, response, replication, groups, alpha)
  if (verbose)
    .aridagri_print_spd(core, "Split-Plot Design: A x B factorial in main plots")
  result <- list(
    design       = "Split-Plot (A x B main, sub)",
    anova_table  = core$anova_table,
    grand_mean   = core$grand_mean,
    cv           = core$cv,
    ms_error_a   = core$error_ms[["main"]], df_error_a = core$error_df[["main"]],
    ms_error_b   = core$error_ms[["sub"]],  df_error_b = core$error_df[["sub"]],
    factor_means = core$factor_means,
    comparisons  = core$comparisons,
    model        = core$model
  )
  class(result) <- c("aridagri_spd_ab", "list")
  return(invisible(result))
}


#' Split Plot Design with C in Main Plot, AB in Sub-Plot
#'
#' @description
#' Performs ANOVA for Split Plot Design where main plot contains single factor C 
#' and sub-plot contains factorial combination of AB.
#'
#' @param data Data frame containing the experimental data
#' @param response Name of the response variable
#' @param main_plot Main plot factor (C)
#' @param sub_factor1 First factor in sub-plot (A)
#' @param sub_factor2 Second factor in sub-plot (B)
#' @param replication Name of replication factor
#' @param posthoc Post-hoc test method
#' @param alpha Significance level
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return List containing ANOVA table and means
#'
#' @examples
#' data <- data.frame(
#'   rep = rep(1:3, each = 24),
#'   irrigation = rep(rep(c("I1", "I2", "I3"), each = 8), 3),
#'   variety = rep(rep(c("V1", "V2"), each = 4), 9),
#'   nitrogen = rep(c("N1", "N2", "N3", "N4"), 18),
#'   yield = rnorm(72, 1200, 150)
#' )
#' anova_spd_c_main_ab_sub(data, response = "yield",
#'                          main_plot = "irrigation",
#'                          sub_factor1 = "variety", sub_factor2 = "nitrogen",
#'                          replication = "rep")
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
anova_spd_c_main_ab_sub <- function(data, response, main_plot, sub_factor1, 
                                     sub_factor2, replication, posthoc = "lsd", 
                                     alpha = 0.05,
                            verbose = TRUE) {
  groups <- list(main = main_plot, sub = c(sub_factor1, sub_factor2))
  core <- .aridagri_spd_anova(data, response, replication, groups, alpha)
  if (verbose)
    .aridagri_print_spd(core, "Split-Plot Design: C in main plots, A x B factorial in sub plots")
  result <- list(
    design       = "Split-Plot (C main, A x B sub)",
    anova_table  = core$anova_table,
    grand_mean   = core$grand_mean,
    cv           = core$cv,
    ms_error_a   = core$error_ms[["main"]], df_error_a = core$error_df[["main"]],
    ms_error_b   = core$error_ms[["sub"]],  df_error_b = core$error_df[["sub"]],
    factor_means = core$factor_means,
    comparisons  = core$comparisons,
    model        = core$model
  )
  class(result) <- c("aridagri_spd_cab", "list")
  return(invisible(result))
}


#' Split Plot Design with (AB) Main and (CD) Sub
#'
#' @description
#' Performs ANOVA for Split Plot Design where main plot contains factorial 
#' combination of AB and sub-plot contains factorial combination of CD.
#' Complex design for multi-factor experiments.
#'
#' @param data Data frame containing the experimental data
#' @param response Name of the response variable
#' @param main_factor1 First factor in main plot (A)
#' @param main_factor2 Second factor in main plot (B)
#' @param sub_factor1 First factor in sub-plot (C)
#' @param sub_factor2 Second factor in sub-plot (D)
#' @param replication Name of replication factor
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return List containing ANOVA table and means
#'
#' @examples
#' # Example: Irrigation  Tillage (main), Variety  Nitrogen (sub)
#' data <- expand.grid(
#'   rep = 1:3,
#'   irrigation = c("I1", "I2"),
#'   tillage = c("CT", "ZT"),
#'   variety = c("V1", "V2"),
#'   nitrogen = c("N1", "N2", "N3")
#' )
#' data$yield <- rnorm(nrow(data), 1200, 150)
#' 
#' anova_spd_ab_cd(data, response = "yield",
#'                 main_factor1 = "irrigation", main_factor2 = "tillage",
#'                 sub_factor1 = "variety", sub_factor2 = "nitrogen",
#'                 replication = "rep")
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
anova_spd_ab_cd <- function(data, response, main_factor1, main_factor2,
                             sub_factor1, sub_factor2, replication,
                            verbose = TRUE) {
  groups <- list(main = c(main_factor1, main_factor2),
                 sub  = c(sub_factor1, sub_factor2))
  core <- .aridagri_spd_anova(data, response, replication, groups, 0.05)
  if (verbose)
    .aridagri_print_spd(core, "Split-Plot Design: A x B factorial main, C x D factorial sub")
  result <- list(
    design       = "Split-Plot (A x B main, C x D sub)",
    anova_table  = core$anova_table,
    grand_mean   = core$grand_mean,
    cv           = core$cv,
    ms_error_a   = core$error_ms[["main"]], df_error_a = core$error_df[["main"]],
    ms_error_b   = core$error_ms[["sub"]],  df_error_b = core$error_df[["sub"]],
    factor_means = core$factor_means,
    comparisons  = core$comparisons,
    model        = core$model
  )
  class(result) <- c("aridagri_spd_abcd", "list")
  return(invisible(result))
}


#' Pooled Split Plot Design Analysis
#'
#' @description
#' Performs pooled analysis of Split Plot Design experiments conducted 
#' across multiple environments/years/locations.
#'
#' @param data Data frame containing combined data
#' @param response Name of the response variable
#' @param main_plot Name of main plot factor
#' @param sub_plot Name of sub-plot factor
#' @param environment Name of environment factor
#' @param replication Name of replication factor (nested within environment)
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return List containing pooled ANOVA and component analyses
#'
#' @examples
#' data <- expand.grid(
#'   year = c("Y1", "Y2", "Y3"),
#'   rep = 1:3,
#'   irrigation = c("I1", "I2", "I3"),
#'   variety = c("V1", "V2", "V3", "V4")
#' )
#' data$yield <- rnorm(nrow(data), 1200, 180)
#' 
#' anova_spd_pooled(data, response = "yield", main_plot = "irrigation",
#'                  sub_plot = "variety", environment = "year", replication = "rep")
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
anova_spd_pooled <- function(data, response, main_plot, sub_plot, 
                              environment, replication,
                            verbose = TRUE) {

  data[[main_plot]]   <- as.factor(data[[main_plot]])
  data[[sub_plot]]    <- as.factor(data[[sub_plot]])
  data[[environment]] <- as.factor(data[[environment]])
  data[[replication]] <- as.factor(data[[replication]])

  a <- nlevels(data[[main_plot]]); b <- nlevels(data[[sub_plot]])
  e <- nlevels(data[[environment]]); r <- nlevels(data[[replication]])
  N <- nrow(data)
  if (N != a * b * e * r) {
    warning("Unbalanced design detected: expected ", a * b * e * r,
            " observations but found ", N, ". ",
            "Results assume balanced data and may be unreliable.", call. = FALSE)
  }

  core <- .aridagri_spd_pooled_anova(data, response, environment, replication,
                                     split_factors = c(main_plot, sub_plot))

  main_means <- aggregate(data[[response]], by = list(data[[main_plot]]), FUN = mean)
  names(main_means) <- c("Level", "Mean")
  sub_means <- aggregate(data[[response]], by = list(data[[sub_plot]]), FUN = mean)
  names(sub_means) <- c("Level", "Mean")
  env_means <- aggregate(data[[response]], by = list(data[[environment]]), FUN = mean)
  names(env_means) <- c("Level", "Mean")

  if (verbose) {
    .aridagri_print_spd_pooled(core, "POOLED SPLIT PLOT DESIGN ANALYSIS")
    cat("\n--- Factor Means ---\n")
    cat("\nMain Plot Means:\n"); print(main_means, row.names = FALSE)
    cat("\nSub-Plot Means:\n");  print(sub_means,  row.names = FALSE)
    cat("\nEnvironment Means:\n"); print(env_means, row.names = FALSE)
  }

  result <- list(
    design      = "Pooled Split Plot Design",
    anova_table = core$anova_table,
    model       = core$model,
    grand_mean  = core$grand_mean,
    cv          = core$cv,
    main_means  = main_means,
    sub_means   = sub_means,
    env_means   = env_means
  )
  class(result) <- c("aridagri_spd_pooled", "list")
  return(invisible(result))
}
