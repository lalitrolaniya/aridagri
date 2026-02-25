#' ============================================================================
#' COMPREHENSIVE EXPERIMENTAL DESIGN ANALYSIS FUNCTIONS
#' Package: aridagri
#' Author: Lalit Kumar Rolaniya
#' ICAR-Indian Institute of Pulses Research, Regional Centre, Bikaner
#' ============================================================================

#' Completely Randomized Design (CRD) ANOVA
#'
#' @description
#' Performs complete ANOVA for Completely Randomized Design with post-hoc tests,
#' assumptions checking, and publication-ready output.
#'
#' @param data Data frame containing the experimental data
#' @param response Name of the response variable (as string)
#' @param treatment Name of treatment factor
#' @param posthoc Post-hoc test: "lsd", "duncan", "tukey", "snk", "scheffe", or "all"
#' @param alpha Significance level (default 0.05)
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return List containing ANOVA table, means, post-hoc results, and diagnostics
#'
#' @examples
#' data <- data.frame(
#'   treatment = rep(c("T1", "T2", "T3", "T4"), each = 5),
#'   yield = c(rnorm(5, 1200, 50), rnorm(5, 1350, 60), 
#'             rnorm(5, 1100, 55), rnorm(5, 1450, 65))
#' )
#' anova_crd(data, response = "yield", treatment = "treatment", posthoc = "all")
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
anova_crd <- function(data, response, treatment, posthoc = "lsd", alpha = 0.05,
                            verbose = TRUE) {
  
  # Input validation
  if (!response %in% names(data)) stop(paste("Response variable", response, "not found"))
  if (!treatment %in% names(data)) stop(paste("Treatment variable", treatment, "not found"))
  
  # Convert to factor
  data[[treatment]] <- as.factor(data[[treatment]])
  
  # Get dimensions
  t <- nlevels(data[[treatment]])
  N <- nrow(data)
  
  # Fit model
  formula_crd <- as.formula(paste(response, "~", treatment))
  model <- aov(formula_crd, data = data)
  anova_table <- anova(model)
  
  # Extract values
  SS_treatment <- anova_table$`Sum Sq`[1]
  SS_error <- anova_table$`Sum Sq`[2]
  SS_total <- SS_treatment + SS_error
  
  df_treatment <- anova_table$Df[1]
  df_error <- anova_table$Df[2]
  df_total <- df_treatment + df_error
  
  MS_treatment <- anova_table$`Mean Sq`[1]
  MS_error <- anova_table$`Mean Sq`[2]
  
  F_value <- anova_table$`F value`[1]
  p_value <- anova_table$`Pr(>F)`[1]
  
  # Calculate statistics
  grand_mean <- mean(data[[response]], na.rm = TRUE)
  cv <- sqrt(MS_error) / grand_mean * 100
  se_mean <- sqrt(MS_error / (N / t))
  cd_5 <- qt(0.975, df_error) * sqrt(2 * MS_error / (N / t))
  cd_1 <- qt(0.995, df_error) * sqrt(2 * MS_error / (N / t))
  
  # Print ANOVA table
  if (verbose) {
    cat("\n")
    cat("\n")
    cat("         COMPLETELY RANDOMIZED DESIGN (CRD) - ANOVA                   \n")
    cat("\n")
    cat(" Response Variable:", sprintf("%-50s", response), "\n")
    cat(" Treatment Factor:", sprintf("%-51s", treatment), "\n")
    cat(" Number of Treatments:", sprintf("%-47d", t), "\n")
    cat(" Total Observations:", sprintf("%-49d", N), "\n")
    cat("\n")

    cat("\n\n")
    cat("                        ANALYSIS OF VARIANCE                          \n")
    cat("\n")
    cat(" Source             df        SS           MS         F    Pr>F  \n")
    cat("\n")
    cat(sprintf(" Treatment        %6d  %11.2f  %11.2f  %5.2f  %5.4f\n",
                df_treatment, SS_treatment, MS_treatment, F_value, p_value))
    cat(sprintf(" Error            %6d  %11.2f  %11.2f               \n",
                df_error, SS_error, MS_error))
    cat("\n")
    cat(sprintf(" Total            %6d  %11.2f                            \n",
                df_total, SS_total))
    cat("\n")
  }
  
  # Significance
  sig <- ifelse(p_value < 0.001, "***", ifelse(p_value < 0.01, "**", 
         ifelse(p_value < 0.05, "*", "NS")))
  
  if (verbose) {
    cat("\n\n")
    cat("                        SUMMARY STATISTICS                            \n")
    cat("\n")
    cat(sprintf(" Grand Mean          : %10.2f                                     \n", grand_mean))
    cat(sprintf(" CV (%%)              : %10.2f                                     \n", cv))
    cat(sprintf(" SE(m)               : %10.3f                                     \n", se_mean))
    cat(sprintf(" CD (5%%)             : %10.3f                                     \n", cd_5))
    cat(sprintf(" CD (1%%)             : %10.3f                                     \n", cd_1))
    cat(sprintf(" Significance        : %10s                                     \n", sig))
    cat("\n")
  }
  
  # Treatment means
  means <- aggregate(data[[response]], by = list(data[[treatment]]), 
                     FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                         sd = sd(x, na.rm = TRUE),
                                         n = length(x)))
  means_df <- data.frame(
    Treatment = means$Group.1,
    Mean = means$x[, "mean"],
    SD = means$x[, "sd"],
    N = means$x[, "n"]
  )
  means_df <- means_df[order(-means_df$Mean), ]
  
  if (verbose) {
    cat("\n\n")
    cat("                        TREATMENT MEANS                               \n")
    cat("\n")
    cat(" Treatment           Mean          SD            N         Rank   \n")
    cat("\n")
    for (i in 1:nrow(means_df)) {
      cat(sprintf(" %-15s  %11.2f  %11.2f  %11.0f  %8d \n",
                  as.character(means_df$Treatment[i]), means_df$Mean[i], 
                  means_df$SD[i], means_df$N[i], i))
    }
    cat("\n")
  }
  
  # Post-hoc tests
  posthoc_results <- perform_posthoc(model, data, response, treatment, 
                                      MS_error, df_error, posthoc, alpha)
  
  # Model diagnostics
  diagnostics <- check_assumptions(model)
  
  # Return results
  result <- list(
    design = "CRD",
    anova_table = anova_table,
    model = model,
    grand_mean = grand_mean,
    cv = cv,
    se_mean = se_mean,
    cd_5 = cd_5,
    cd_1 = cd_1,
    treatment_means = means_df,
    posthoc = posthoc_results,
    diagnostics = diagnostics,
    significance = sig
  )
  
  class(result) <- c("aridagri_anova", "list")
  return(invisible(result))
}


#' Randomized Block Design (RBD) ANOVA
#'
#' @description
#' Performs complete ANOVA for Randomized Block Design (RCBD) with post-hoc tests,
#' assumptions checking, and publication-ready output.
#'
#' @param data Data frame containing the experimental data
#' @param response Name of the response variable
#' @param treatment Name of treatment factor
#' @param block Name of block/replication factor
#' @param posthoc Post-hoc test: "lsd", "duncan", "tukey", "snk", "scheffe", or "all"
#' @param alpha Significance level (default 0.05)
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return List containing ANOVA table, means, post-hoc results
#'
#' @examples
#' data <- data.frame(
#'   rep = rep(1:4, each = 5),
#'   treatment = rep(c("T1", "T2", "T3", "T4", "T5"), 4),
#'   yield = c(rnorm(5, 1200, 50), rnorm(5, 1250, 55),
#'             rnorm(5, 1180, 45), rnorm(5, 1270, 60))
#' )
#' anova_rbd(data, response = "yield", treatment = "treatment", 
#'           block = "rep", posthoc = "duncan")
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
anova_rbd <- function(data, response, treatment, block, posthoc = "lsd", alpha = 0.05,
                            verbose = TRUE) {
  
  # Input validation
  if (!response %in% names(data)) stop(paste("Response variable", response, "not found"))
  if (!treatment %in% names(data)) stop(paste("Treatment variable", treatment, "not found"))
  if (!block %in% names(data)) stop(paste("Block variable", block, "not found"))
  
  # Convert to factors
  data[[treatment]] <- as.factor(data[[treatment]])
  data[[block]] <- as.factor(data[[block]])
  
  # Get dimensions
  t <- nlevels(data[[treatment]])
  r <- nlevels(data[[block]])
  N <- nrow(data)
  
  # Fit model
  formula_rbd <- as.formula(paste(response, "~", block, "+", treatment))
  model <- aov(formula_rbd, data = data)
  anova_table <- anova(model)
  
  # Extract values
  SS_block <- anova_table$`Sum Sq`[1]
  SS_treatment <- anova_table$`Sum Sq`[2]
  SS_error <- anova_table$`Sum Sq`[3]
  SS_total <- SS_block + SS_treatment + SS_error
  
  df_block <- anova_table$Df[1]
  df_treatment <- anova_table$Df[2]
  df_error <- anova_table$Df[3]
  df_total <- df_block + df_treatment + df_error
  
  MS_block <- anova_table$`Mean Sq`[1]
  MS_treatment <- anova_table$`Mean Sq`[2]
  MS_error <- anova_table$`Mean Sq`[3]
  
  F_block <- anova_table$`F value`[1]
  F_treatment <- anova_table$`F value`[2]
  p_block <- anova_table$`Pr(>F)`[1]
  p_treatment <- anova_table$`Pr(>F)`[2]
  
  # Calculate statistics
  grand_mean <- mean(data[[response]], na.rm = TRUE)
  cv <- sqrt(MS_error) / grand_mean * 100
  se_mean <- sqrt(MS_error / r)
  cd_5 <- qt(0.975, df_error) * sqrt(2 * MS_error / r)
  cd_1 <- qt(0.995, df_error) * sqrt(2 * MS_error / r)
  
  # Relative efficiency of RBD over CRD
  MS_block_adj <- max(MS_block, MS_error)
  re_crd <- ((df_block * MS_block_adj + (df_treatment + df_block) * MS_error) /
              ((df_block + df_treatment + df_error) * MS_error)) * 100
  
  # Print ANOVA table
  if (verbose) {
    cat("\n")
    cat("\n")
    cat("       RANDOMIZED BLOCK DESIGN (RBD/RCBD) - ANOVA                     \n")
    cat("\n")
    cat(" Response Variable:", sprintf("%-50s", response), "\n")
    cat(" Treatment Factor:", sprintf("%-51s", treatment), "\n")
    cat(" Block Factor:", sprintf("%-55s", block), "\n")
    cat(" Number of Treatments:", sprintf("%-47d", t), "\n")
    cat(" Number of Blocks:", sprintf("%-51d", r), "\n")
    cat(" Total Observations:", sprintf("%-49d", N), "\n")
    cat("\n")

    cat("\n\n")
    cat("                        ANALYSIS OF VARIANCE                          \n")
    cat("\n")
    cat(" Source             df        SS           MS         F    Pr>F  \n")
    cat("\n")
    cat(sprintf(" Block            %6d  %11.2f  %11.2f  %5.2f  %5.4f\n",
                df_block, SS_block, MS_block, F_block, p_block))
    cat(sprintf(" Treatment        %6d  %11.2f  %11.2f  %5.2f  %5.4f\n",
                df_treatment, SS_treatment, MS_treatment, F_treatment, p_treatment))
    cat(sprintf(" Error            %6d  %11.2f  %11.2f               \n",
                df_error, SS_error, MS_error))
    cat("\n")
    cat(sprintf(" Total            %6d  %11.2f                            \n",
                df_total, SS_total))
    cat("\n")
  }
  
  # Significance
  sig <- ifelse(p_treatment < 0.001, "***", ifelse(p_treatment < 0.01, "**", 
         ifelse(p_treatment < 0.05, "*", "NS")))
  
  if (verbose) {
    cat("\n\n")
    cat("                        SUMMARY STATISTICS                            \n")
    cat("\n")
    cat(sprintf(" Grand Mean                    : %10.2f                         \n", grand_mean))
    cat(sprintf(" CV (%%)                        : %10.2f                         \n", cv))
    cat(sprintf(" SE(m)                         : %10.3f                         \n", se_mean))
    cat(sprintf(" SE(d)                         : %10.3f                         \n", se_mean * sqrt(2)))
    cat(sprintf(" CD (5%%)                       : %10.3f                         \n", cd_5))
    cat(sprintf(" CD (1%%)                       : %10.3f                         \n", cd_1))
    cat(sprintf(" Treatment Significance        : %10s                         \n", sig))
    cat(sprintf(" Relative Efficiency over CRD  : %10.1f%%                        \n", re_crd))
    cat("\n")
  }
  
  # Treatment means
  means <- aggregate(data[[response]], by = list(data[[treatment]]), 
                     FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                         sd = sd(x, na.rm = TRUE),
                                         n = length(x)))
  means_df <- data.frame(
    Treatment = means$Group.1,
    Mean = means$x[, "mean"],
    SD = means$x[, "sd"],
    N = means$x[, "n"]
  )
  means_df <- means_df[order(-means_df$Mean), ]
  
  if (verbose) {
    cat("\n\n")
    cat("                        TREATMENT MEANS                               \n")
    cat("\n")
    cat(" Treatment           Mean          SD            N         Rank   \n")
    cat("\n")
    for (i in 1:nrow(means_df)) {
      cat(sprintf(" %-15s  %11.2f  %11.2f  %11.0f  %8d \n",
                  as.character(means_df$Treatment[i]), means_df$Mean[i], 
                  means_df$SD[i], means_df$N[i], i))
    }
    cat("\n")
  }
  
  # Post-hoc tests
  posthoc_results <- perform_posthoc(model, data, response, treatment, 
                                      MS_error, df_error, posthoc, alpha)
  
  # Return results
  result <- list(
    design = "RBD",
    anova_table = anova_table,
    model = model,
    grand_mean = grand_mean,
    cv = cv,
    se_mean = se_mean,
    se_diff = se_mean * sqrt(2),
    cd_5 = cd_5,
    cd_1 = cd_1,
    relative_efficiency = re_crd,
    treatment_means = means_df,
    posthoc = posthoc_results,
    significance = sig
  )
  
  class(result) <- c("aridagri_anova", "list")
  return(invisible(result))
}


#' Pooled Analysis of RBD Experiments (Multi-Environment/Year)
#'
#' @description
#' Performs pooled ANOVA for RBD experiments conducted across multiple 
#' environments, years, or locations. Tests homogeneity of error variances
#' using Bartlett's test before pooling.
#'
#' @param data Data frame containing combined data from all environments
#' @param response Name of the response variable
#' @param treatment Name of treatment factor
#' @param environment Name of environment/year/location factor
#' @param block Name of block factor (nested within environment)
#' @param posthoc Post-hoc test method
#' @param alpha Significance level
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return List containing pooled ANOVA, individual ANOVAs, and interaction analysis
#'
#' @examples
#' # Data from 3 years
#' data <- data.frame(
#'   year = rep(c("Y1", "Y2", "Y3"), each = 20),
#'   rep = rep(rep(1:4, each = 5), 3),
#'   treatment = rep(c("T1", "T2", "T3", "T4", "T5"), 12),
#'   yield = rnorm(60, 1200, 150)
#' )
#' anova_rbd_pooled(data, response = "yield", treatment = "treatment",
#'                  environment = "year", block = "rep")
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
anova_rbd_pooled <- function(data, response, treatment, environment, block, 
                              posthoc = "lsd", alpha = 0.05,
                            verbose = TRUE) {
  
  # Convert to factors
  data[[treatment]] <- as.factor(data[[treatment]])
  data[[environment]] <- as.factor(data[[environment]])
  data[[block]] <- as.factor(data[[block]])
  
  # Get dimensions
  t <- nlevels(data[[treatment]])
  e <- nlevels(data[[environment]])
  r <- nlevels(data[[block]])
  N <- nrow(data)
  
  # Individual ANOVAs for each environment
  env_levels <- levels(data[[environment]])
  individual_anovas <- list()
  mse_values <- numeric(e)
  df_error_values <- numeric(e)
  
  if (verbose) {
    cat("\n")
    cat("\n")
    cat("              POOLED RBD ANALYSIS (Multi-Environment)                 \n")
    cat("\n")
    cat(" Response Variable:", sprintf("%-50s", response), "\n")
    cat(" Treatment Factor:", sprintf("%-51s", treatment), "\n")
    cat(" Environment Factor:", sprintf("%-49s", environment), "\n")
    cat(" Block Factor:", sprintf("%-55s", block), "\n")
    cat(" Number of Treatments:", sprintf("%-47d", t), "\n")
    cat(" Number of Environments:", sprintf("%-45d", e), "\n")
    cat(" Replications per Environment:", sprintf("%-39d", r), "\n")
    cat(" Total Observations:", sprintf("%-49d", N), "\n")
    cat("\n")

    cat("\n\n")
    cat("              INDIVIDUAL ENVIRONMENT ANALYSES                         \n")
    cat("\n")

    for (i in 1:e) {
      env_data <- data[data[[environment]] == env_levels[i], ]
      formula_rbd <- as.formula(paste(response, "~", block, "+", treatment))
      model_i <- aov(formula_rbd, data = env_data)
      anova_i <- anova(model_i)

      individual_anovas[[env_levels[i]]] <- anova_i
      mse_values[i] <- anova_i$`Mean Sq`[3]
      df_error_values[i] <- anova_i$Df[3]

      mean_i <- mean(env_data[[response]], na.rm = TRUE)
      cv_i <- sqrt(mse_values[i]) / mean_i * 100

      cat(sprintf("\n--- %s ---\n", env_levels[i]))
      cat(sprintf("Mean: %.2f | MSE: %.2f | CV: %.2f%%\n", mean_i, mse_values[i], cv_i))
    }
  }
  
  # Bartlett's test for homogeneity of error variances
  if (verbose) {
    cat("\n\n")
    cat("           HOMOGENEITY OF ERROR VARIANCES (Bartlett's Test)           \n")
    cat("\n")
  }
  
  # Calculate Bartlett's test manually
  mse_pooled <- sum(df_error_values * mse_values) / sum(df_error_values)
  
  numerator <- sum(df_error_values) * log(mse_pooled) - sum(df_error_values * log(mse_values))
  c_factor <- 1 + (1 / (3 * (e - 1))) * (sum(1 / df_error_values) - 1 / sum(df_error_values))
  chi_sq <- numerator / c_factor
  p_bartlett <- 1 - pchisq(chi_sq, e - 1)
  
  if (verbose) {
    cat(sprintf("\nBartlett's Chi-square: %.4f\n", chi_sq))
    cat(sprintf("Degrees of freedom: %d\n", e - 1))
    cat(sprintf("P-value: %.4f\n", p_bartlett))

    if (p_bartlett < 0.05) {
      cat("WARNING: Error variances are heterogeneous (p < 0.05)\n")
      cat("Consider data transformation or weighted analysis.\n")
    } else {
      cat("Error variances are homogeneous. Pooling is appropriate.\n")
    }
  }
  
  # Pooled ANOVA
  if (verbose) {
    cat("\n\n")
    cat("                        POOLED ANALYSIS OF VARIANCE                   \n")
    cat("\n")
  }
  
  # Create nested block factor
  data$block_nested <- interaction(data[[environment]], data[[block]])
  
  # Fit pooled model
  formula_pooled <- as.formula(paste(response, "~", environment, "+", 
                                      environment, ":", block, "+",
                                      treatment, "+", environment, ":", treatment))
  model_pooled <- aov(formula_pooled, data = data)
  anova_pooled <- anova(model_pooled)
  
  if (verbose) {
    print(anova_pooled)
  }
  
  # Calculate pooled statistics
  pooled_error_df <- anova_pooled$Df[length(anova_pooled$Df)]
  pooled_error_ms <- anova_pooled$`Mean Sq`[length(anova_pooled$`Mean Sq`)]
  
  grand_mean <- mean(data[[response]], na.rm = TRUE)
  cv_pooled <- sqrt(pooled_error_ms) / grand_mean * 100
  
  # SE and CD for different comparisons
  se_treatment <- sqrt(pooled_error_ms / (e * r))
  cd_treatment_5 <- qt(0.975, pooled_error_df) * sqrt(2 * pooled_error_ms / (e * r))
  
  se_env <- sqrt(pooled_error_ms / (t * r))
  cd_env_5 <- qt(0.975, pooled_error_df) * sqrt(2 * pooled_error_ms / (t * r))
  
  se_interaction <- sqrt(pooled_error_ms / r)
  cd_interaction_5 <- qt(0.975, pooled_error_df) * sqrt(2 * pooled_error_ms / r)
  
  if (verbose) {
    cat("\n\n")
    cat("                        POOLED SUMMARY STATISTICS                     \n")
    cat("\n")
    cat(sprintf(" Grand Mean                    : %10.2f                         \n", grand_mean))
    cat(sprintf(" Pooled CV (%%)                 : %10.2f                         \n", cv_pooled))
    cat("\n")
    cat(" For comparing TREATMENT means:                                       \n")
    cat(sprintf("   SE(m)  : %10.3f    CD(5%%) : %10.3f                        \n", se_treatment, cd_treatment_5))
    cat("\n")
    cat(" For comparing ENVIRONMENT means:                                     \n")
    cat(sprintf("   SE(m)  : %10.3f    CD(5%%) : %10.3f                        \n", se_env, cd_env_5))
    cat("\n")
    cat(" For comparing INTERACTION means:                                     \n")
    cat(sprintf("   SE(m)  : %10.3f    CD(5%%) : %10.3f                        \n", se_interaction, cd_interaction_5))
    cat("\n")
  }
  
  # Treatment means across environments
  treatment_means <- aggregate(data[[response]], by = list(data[[treatment]]), 
                                FUN = mean, na.rm = TRUE)
  names(treatment_means) <- c("Treatment", "Mean")
  treatment_means <- treatment_means[order(-treatment_means$Mean), ]
  
  if (verbose) {
    cat("\n\n")
    cat("                   TREATMENT MEANS (Across Environments)              \n")
    cat("\n")
    for (i in 1:nrow(treatment_means)) {
      cat(sprintf(" %-31s  %34.2f \n", 
                  as.character(treatment_means$Treatment[i]), treatment_means$Mean[i]))
    }
    cat("\n")
  }
  
  # Interaction means table
  interaction_means <- aggregate(data[[response]], 
                                  by = list(data[[treatment]], data[[environment]]), 
                                  FUN = mean, na.rm = TRUE)
  names(interaction_means) <- c("Treatment", "Environment", "Mean")
  interaction_table <- reshape(interaction_means, idvar = "Treatment", 
                                timevar = "Environment", direction = "wide")
  
  if (verbose) {
    cat("\n\n")
    cat("                   TREATMENT  ENVIRONMENT MEANS                      \n")
    cat("\n")
    print(interaction_table, row.names = FALSE)
  }
  
  # Return results
  result <- list(
    design = "Pooled RBD",
    pooled_anova = anova_pooled,
    individual_anovas = individual_anovas,
    model = model_pooled,
    grand_mean = grand_mean,
    cv = cv_pooled,
    bartlett_test = list(chi_sq = chi_sq, df = e - 1, p_value = p_bartlett),
    se_treatment = se_treatment,
    cd_treatment = cd_treatment_5,
    se_environment = se_env,
    cd_environment = cd_env_5,
    se_interaction = se_interaction,
    cd_interaction = cd_interaction_5,
    treatment_means = treatment_means,
    interaction_means = interaction_table,
    mse_values = mse_values
  )
  
  class(result) <- c("aridagri_pooled", "list")
  return(invisible(result))
}
