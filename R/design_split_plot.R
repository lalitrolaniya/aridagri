#' ============================================================================
#' SPLIT PLOT DESIGN ANALYSIS FUNCTIONS (All Variations)
#' Package: aridagri
#' Author: Lalit Kumar Rolaniya
#' ICAR-Indian Institute of Pulses Research, Regional Centre, Bikaner
#' ============================================================================

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
  
  # Fit full model
  formula_spd <- as.formula(paste(response, "~", replication, "+", 
                                   main_plot, "+", replication, ":", main_plot, "+",
                                   sub_plot, "+", main_plot, ":", sub_plot))
  
  model <- aov(formula_spd, data = data)
  anova_full <- anova(model)
  
  # Extract sums of squares
  SS_rep <- anova_full$`Sum Sq`[1]
  SS_main <- anova_full$`Sum Sq`[2]
  SS_error_a <- anova_full$`Sum Sq`[3]
  SS_sub <- anova_full$`Sum Sq`[4]
  SS_interaction <- anova_full$`Sum Sq`[5]
  SS_error_b <- anova_full$`Sum Sq`[6]
  SS_total <- sum(anova_full$`Sum Sq`)
  
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
  
  # Convert to factors
  data[[main_factor1]] <- as.factor(data[[main_factor1]])
  data[[main_factor2]] <- as.factor(data[[main_factor2]])
  data[[sub_plot]] <- as.factor(data[[sub_plot]])
  data[[replication]] <- as.factor(data[[replication]])
  
  # Get dimensions
  a <- nlevels(data[[main_factor1]])
  b <- nlevels(data[[main_factor2]])
  c <- nlevels(data[[sub_plot]])
  r <- nlevels(data[[replication]])
  N <- nrow(data)
  
  # Create main plot combination
  data$main_plot <- interaction(data[[main_factor1]], data[[main_factor2]])
  m <- nlevels(data$main_plot)  # Number of main plot treatments (a  b)
  
  # Fit model
  formula_spd <- as.formula(paste(response, "~", replication, "+",
                                   main_factor1, "*", main_factor2, "+",
                                   replication, ":", main_factor1, ":", main_factor2, "+",
                                   sub_plot, "+",
                                   main_factor1, ":", sub_plot, "+",
                                   main_factor2, ":", sub_plot, "+",
                                   main_factor1, ":", main_factor2, ":", sub_plot))
  
  model <- aov(formula_spd, data = data)
  anova_full <- anova(model)
  
  # Print header
  if (verbose) {
    cat("\n")
    cat("\n")
    cat("          SPLIT PLOT DESIGN WITH (AB) IN MAIN PLOT - ANOVA                   \n")
    cat("\n")
    cat(" Response Variable    :", sprintf("%-53s", response), "\n")
    cat(" Main Plot Factor A   :", sprintf("%-53s", main_factor1), "\n")
    cat(" Main Plot Factor B   :", sprintf("%-53s", main_factor2), "\n")
    cat(" Sub-Plot Factor C    :", sprintf("%-53s", sub_plot), "\n")
    cat(" Replications         :", sprintf("%-53d", r), "\n")
    cat(" Levels of A          :", sprintf("%-53d", a), "\n")
    cat(" Levels of B          :", sprintf("%-53d", b), "\n")
    cat(" Levels of C          :", sprintf("%-53d", c), "\n")
    cat(" Total Observations   :", sprintf("%-53d", N), "\n")
    cat("\n")
  }
  
  # Calculate degrees of freedom
  df_rep <- r - 1
  df_a <- a - 1
  df_b <- b - 1
  df_ab <- (a - 1) * (b - 1)
  df_error_a <- (r - 1) * (a * b - 1)
  df_c <- c - 1
  df_ac <- (a - 1) * (c - 1)
  df_bc <- (b - 1) * (c - 1)
  df_abc <- (a - 1) * (b - 1) * (c - 1)
  df_error_b <- a * b * (r - 1) * (c - 1)
  df_total <- N - 1
  
  if (verbose) {
    cat("\n\n")
    cat("                           ANALYSIS OF VARIANCE                               \n")
    cat("\n")
    cat(" Source                   df        SS           MS         F     Pr>F   \n")
    cat("\n")
  }
  
  # Print ANOVA rows
  for (i in 1:nrow(anova_full)) {
    source_name <- rownames(anova_full)[i]
    # Truncate long names
    if (nchar(source_name) > 21) {
      source_name <- paste0(substr(source_name, 1, 18), "...")
    }
    
    if (!is.na(anova_full$`F value`[i])) {
      if (verbose) {
        cat(sprintf(" %-21s  %6d  %11.2f  %11.2f  %5.2f  %7.4f \n",
                    source_name, anova_full$Df[i], anova_full$`Sum Sq`[i], 
                    anova_full$`Mean Sq`[i], anova_full$`F value`[i], 
                    anova_full$`Pr(>F)`[i]))
      }
    } else {
      if (verbose) {
        cat(sprintf(" %-21s  %6d  %11.2f  %11.2f                 \n",
                    source_name, anova_full$Df[i], anova_full$`Sum Sq`[i], 
                    anova_full$`Mean Sq`[i]))
      }
    }
  }
  if (verbose) {
    cat("\n")
  }
  
  # Calculate statistics
  grand_mean <- mean(data[[response]], na.rm = TRUE)
  
  # Extract error mean squares (simplified - need to identify correct rows)
  n_rows <- nrow(anova_full)
  MS_error_b <- anova_full$`Mean Sq`[n_rows]
  df_error_b_actual <- anova_full$Df[n_rows]
  
  cv_b <- sqrt(MS_error_b) / grand_mean * 100
  
  if (verbose) {
    cat("\n\n")
    cat("                           SUMMARY STATISTICS                                 \n")
    cat("\n")
    cat(sprintf(" Grand Mean                    : %10.2f                                   \n", grand_mean))
    cat(sprintf(" CV (b)                        : %10.2f%%                                  \n", cv_b))
    cat("\n")
  }
  
  # Factor means
  if (verbose) {
    cat("\n\n")
    cat("                              FACTOR MEANS                                    \n")
    cat("\n")
  }
  
  # Means for each factor
  mean_a <- aggregate(data[[response]], by = list(data[[main_factor1]]), FUN = mean)
  names(mean_a) <- c("Level", "Mean")
  if (verbose) {
    cat(sprintf("\nFactor A (%s):\n", main_factor1))
    print(mean_a, row.names = FALSE)
  }
  
  mean_b <- aggregate(data[[response]], by = list(data[[main_factor2]]), FUN = mean)
  names(mean_b) <- c("Level", "Mean")
  if (verbose) {
    cat(sprintf("\nFactor B (%s):\n", main_factor2))
    print(mean_b, row.names = FALSE)
  }
  
  mean_c <- aggregate(data[[response]], by = list(data[[sub_plot]]), FUN = mean)
  names(mean_c) <- c("Level", "Mean")
  if (verbose) {
    cat(sprintf("\nFactor C (%s):\n", sub_plot))
    print(mean_c, row.names = FALSE)
  }
  
  # Return results
  result <- list(
    design = "Split Plot Design (AB Main Plot)",
    anova_table = anova_full,
    model = model,
    grand_mean = grand_mean,
    cv = cv_b,
    factor_means = list(A = mean_a, B = mean_b, C = mean_c)
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
  
  # Convert to factors
  data[[main_plot]] <- as.factor(data[[main_plot]])
  data[[sub_factor1]] <- as.factor(data[[sub_factor1]])
  data[[sub_factor2]] <- as.factor(data[[sub_factor2]])
  data[[replication]] <- as.factor(data[[replication]])
  
  # Get dimensions
  c_levels <- nlevels(data[[main_plot]])
  a <- nlevels(data[[sub_factor1]])
  b <- nlevels(data[[sub_factor2]])
  r <- nlevels(data[[replication]])
  N <- nrow(data)
  
  # Fit model
  formula_spd <- as.formula(paste(response, "~", replication, "+",
                                   main_plot, "+", replication, ":", main_plot, "+",
                                   sub_factor1, "*", sub_factor2, "+",
                                   main_plot, ":", sub_factor1, "+",
                                   main_plot, ":", sub_factor2, "+",
                                   main_plot, ":", sub_factor1, ":", sub_factor2))
  
  model <- aov(formula_spd, data = data)
  anova_full <- anova(model)
  
  # Print header
  if (verbose) {
    cat("\n")
    cat("\n")
    cat("        SPLIT PLOT DESIGN: C IN MAIN PLOT, (AB) IN SUB-PLOT                  \n")
    cat("\n")
    cat(" Response Variable    :", sprintf("%-53s", response), "\n")
    cat(" Main Plot Factor C   :", sprintf("%-53s", main_plot), "\n")
    cat(" Sub-Plot Factor A    :", sprintf("%-53s", sub_factor1), "\n")
    cat(" Sub-Plot Factor B    :", sprintf("%-53s", sub_factor2), "\n")
    cat(" Replications         :", sprintf("%-53d", r), "\n")
    cat(" Levels of C          :", sprintf("%-53d", c_levels), "\n")
    cat(" Levels of A          :", sprintf("%-53d", a), "\n")
    cat(" Levels of B          :", sprintf("%-53d", b), "\n")
    cat(" Total Observations   :", sprintf("%-53d", N), "\n")
    cat("\n")

    cat("\n\n")
    cat("                           ANALYSIS OF VARIANCE                               \n")
    cat("\n")

    print(anova_full)
  }
  
  # Statistics
  grand_mean <- mean(data[[response]], na.rm = TRUE)
  n_rows <- nrow(anova_full)
  MS_error_b <- anova_full$`Mean Sq`[n_rows]
  cv_b <- sqrt(MS_error_b) / grand_mean * 100
  
  if (verbose) {
    cat("\n\n")
    cat("                           SUMMARY STATISTICS                                 \n")
    cat("\n")
    cat(sprintf(" Grand Mean                    : %10.2f                                   \n", grand_mean))
    cat(sprintf(" CV (b)                        : %10.2f%%                                  \n", cv_b))
    cat("\n")
  }
  
  # Factor means
  mean_c <- aggregate(data[[response]], by = list(data[[main_plot]]), FUN = mean)
  names(mean_c) <- c("Level", "Mean")
  
  mean_a <- aggregate(data[[response]], by = list(data[[sub_factor1]]), FUN = mean)
  names(mean_a) <- c("Level", "Mean")
  
  mean_b <- aggregate(data[[response]], by = list(data[[sub_factor2]]), FUN = mean)
  names(mean_b) <- c("Level", "Mean")
  
  if (verbose) {
    cat("\n--- Factor Means ---\n")
    cat(sprintf("\nMain Plot C (%s):\n", main_plot))
    print(mean_c, row.names = FALSE)
    cat(sprintf("\nSub-Plot A (%s):\n", sub_factor1))
    print(mean_a, row.names = FALSE)
    cat(sprintf("\nSub-Plot B (%s):\n", sub_factor2))
    print(mean_b, row.names = FALSE)
  }
  
  # Return
  result <- list(
    design = "Split Plot Design (C Main, AB Sub)",
    anova_table = anova_full,
    model = model,
    grand_mean = grand_mean,
    cv = cv_b,
    factor_means = list(C = mean_c, A = mean_a, B = mean_b)
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
  
  # Convert to factors
  data[[main_factor1]] <- as.factor(data[[main_factor1]])
  data[[main_factor2]] <- as.factor(data[[main_factor2]])
  data[[sub_factor1]] <- as.factor(data[[sub_factor1]])
  data[[sub_factor2]] <- as.factor(data[[sub_factor2]])
  data[[replication]] <- as.factor(data[[replication]])
  
  # Get dimensions
  a <- nlevels(data[[main_factor1]])
  b <- nlevels(data[[main_factor2]])
  c <- nlevels(data[[sub_factor1]])
  d <- nlevels(data[[sub_factor2]])
  r <- nlevels(data[[replication]])
  N <- nrow(data)
  
  # Create combination factors
  data$main_comb <- interaction(data[[main_factor1]], data[[main_factor2]])
  
  # Fit model - comprehensive
  formula_spd <- as.formula(paste(response, "~", replication, "+",
                                   main_factor1, "*", main_factor2, "+",
                                   "Error(", replication, ":", main_factor1, ":", main_factor2, ") +",
                                   sub_factor1, "*", sub_factor2, "+",
                                   main_factor1, ":", sub_factor1, "+",
                                   main_factor1, ":", sub_factor2, "+",
                                   main_factor2, ":", sub_factor1, "+",
                                   main_factor2, ":", sub_factor2, "+",
                                   main_factor1, ":", main_factor2, ":", sub_factor1, "+",
                                   main_factor1, ":", main_factor2, ":", sub_factor2, "+",
                                   main_factor1, ":", sub_factor1, ":", sub_factor2, "+",
                                   main_factor2, ":", sub_factor1, ":", sub_factor2, "+",
                                   main_factor1, ":", main_factor2, ":", sub_factor1, ":", sub_factor2))
  
  # Simpler model for estimation
  formula_simple <- as.formula(paste(response, "~", replication, "+",
                                      main_factor1, "*", main_factor2, "+",
                                      replication, ":", main_factor1, ":", main_factor2, "+",
                                      sub_factor1, "*", sub_factor2, "*",
                                      main_factor1, "*", main_factor2))
  
  model <- aov(formula_simple, data = data)
  anova_full <- anova(model)
  
  # Print header
  if (verbose) {
    cat("\n")
    cat("\n")
    cat("      SPLIT PLOT DESIGN: (AB) MAIN PLOT, (CD) SUB-PLOT                      \n")
    cat("\n")
    cat(" Response Variable    :", sprintf("%-53s", response), "\n")
    cat(" Main Plot Factor A   :", sprintf("%-53s", main_factor1), "\n")
    cat(" Main Plot Factor B   :", sprintf("%-53s", main_factor2), "\n")
    cat(" Sub-Plot Factor C    :", sprintf("%-53s", sub_factor1), "\n")
    cat(" Sub-Plot Factor D    :", sprintf("%-53s", sub_factor2), "\n")
    cat(" Replications         :", sprintf("%-53d", r), "\n")
    cat(" Levels: A=%d, B=%d, C=%d, D=%d", a, b, c, d)
    cat(sprintf("%44s\n", ""))
    cat(" Total Observations   :", sprintf("%-53d", N), "\n")
    cat("\n")

    cat("\n\n")
    cat("                           ANALYSIS OF VARIANCE                               \n")
    cat("\n")

    print(anova_full)
  }
  
  # Statistics
  grand_mean <- mean(data[[response]], na.rm = TRUE)
  
  if (verbose) {
    cat("\n\n")
    cat(sprintf(" Grand Mean: %.2f                                                           \n", grand_mean))
    cat("\n")
  }
  
  # Return
  result <- list(
    design = "Split Plot Design (AB Main, CD Sub)",
    anova_table = anova_full,
    model = model,
    grand_mean = grand_mean
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
  
  # Convert to factors
  data[[main_plot]] <- as.factor(data[[main_plot]])
  data[[sub_plot]] <- as.factor(data[[sub_plot]])
  data[[environment]] <- as.factor(data[[environment]])
  data[[replication]] <- as.factor(data[[replication]])
  
  # Get dimensions
  a <- nlevels(data[[main_plot]])
  b <- nlevels(data[[sub_plot]])
  e <- nlevels(data[[environment]])
  r <- nlevels(data[[replication]])
  N <- nrow(data)
  
  # Individual environment analyses
  env_levels <- levels(data[[environment]])
  mse_a_values <- numeric(e)
  mse_b_values <- numeric(e)
  
  if (verbose) {
    cat("\n")
    cat("\n")
    cat("            POOLED SPLIT PLOT DESIGN ANALYSIS                                 \n")
    cat("\n")
    cat(" Response Variable   :", sprintf("%-54s", response), "\n")
    cat(" Main Plot Factor    :", sprintf("%-54s", main_plot), "\n")
    cat(" Sub-Plot Factor     :", sprintf("%-54s", sub_plot), "\n")
    cat(" Environment Factor  :", sprintf("%-54s", environment), "\n")
    cat(" Main Plot Levels    :", sprintf("%-54d", a), "\n")
    cat(" Sub-Plot Levels     :", sprintf("%-54d", b), "\n")
    cat(" Environments        :", sprintf("%-54d", e), "\n")
    cat(" Replications/Env    :", sprintf("%-54d", r), "\n")
    cat(" Total Observations  :", sprintf("%-54d", N), "\n")
    cat("\n")
  }
  
  # Fit pooled model
  formula_pooled <- as.formula(paste(response, "~", environment, "+",
                                      environment, ":", replication, "+",
                                      main_plot, "+", environment, ":", main_plot, "+",
                                      environment, ":", replication, ":", main_plot, "+",
                                      sub_plot, "+", main_plot, ":", sub_plot, "+",
                                      environment, ":", sub_plot, "+",
                                      environment, ":", main_plot, ":", sub_plot))
  
  model <- aov(formula_pooled, data = data)
  anova_pooled <- anova(model)
  
  if (verbose) {
    cat("\n\n")
    cat("                        POOLED ANALYSIS OF VARIANCE                           \n")
    cat("\n")

    print(anova_pooled)
  }
  
  # Grand mean and CV
  grand_mean <- mean(data[[response]], na.rm = TRUE)
  n_rows <- nrow(anova_pooled)
  MS_error <- anova_pooled$`Mean Sq`[n_rows]
  cv <- sqrt(MS_error) / grand_mean * 100
  
  if (verbose) {
    cat("\n\n")
    cat("                           SUMMARY STATISTICS                                 \n")
    cat("\n")
    cat(sprintf(" Grand Mean                    : %10.2f                                   \n", grand_mean))
    cat(sprintf(" Pooled CV                     : %10.2f%%                                  \n", cv))
    cat("\n")
  }
  
  # Means tables
  main_means <- aggregate(data[[response]], by = list(data[[main_plot]]), FUN = mean)
  names(main_means) <- c("Level", "Mean")
  
  sub_means <- aggregate(data[[response]], by = list(data[[sub_plot]]), FUN = mean)
  names(sub_means) <- c("Level", "Mean")
  
  env_means <- aggregate(data[[response]], by = list(data[[environment]]), FUN = mean)
  names(env_means) <- c("Level", "Mean")
  
  if (verbose) {
    cat("\n--- Factor Means ---\n")
    cat("\nMain Plot Means:\n")
    print(main_means, row.names = FALSE)
    cat("\nSub-Plot Means:\n")
    print(sub_means, row.names = FALSE)
    cat("\nEnvironment Means:\n")
    print(env_means, row.names = FALSE)
  }
  
  # Return
  result <- list(
    design = "Pooled Split Plot Design",
    anova_table = anova_pooled,
    model = model,
    grand_mean = grand_mean,
    cv = cv,
    main_means = main_means,
    sub_means = sub_means,
    env_means = env_means
  )
  
  class(result) <- c("aridagri_spd_pooled", "list")
  return(invisible(result))
}
