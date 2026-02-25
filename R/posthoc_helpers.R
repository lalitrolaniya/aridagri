#' ============================================================================
#' POST-HOC TESTS AND HELPER FUNCTIONS
#' Package: aridagri
#' Author: Lalit Kumar Rolaniya
#' ICAR-Indian Institute of Pulses Research, Regional Centre, Bikaner
#' ============================================================================

#' Get Significance Symbol
#'
#' @description
#' Returns significance symbol based on p-value.
#'
#' @param p_value P-value from statistical test
#' @return Character string with significance symbol
#' @keywords internal
get_sig_symbol <- function(p_value) {
  if (is.na(p_value)) return("")
  if (p_value < 0.001) return("***")
  if (p_value < 0.01) return("**")
  if (p_value < 0.05) return("*")
  return("NS")
}


#' Perform Post-Hoc Tests
#'
#' @description
#' Performs multiple comparison tests after ANOVA.
#'
#' @param model ANOVA model object
#' @param data Data frame
#' @param response Response variable name
#' @param treatment Treatment factor name
#' @param mse Mean square error
#' @param df_error Error degrees of freedom
#' @param posthoc Test type: "lsd", "duncan", "tukey", "snk", "scheffe", or "all"
#' @param alpha Significance level
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return List containing post-hoc test results
#' @export
perform_posthoc <- function(model, data, response, treatment, mse, df_error, 
                            posthoc = "lsd", alpha = 0.05,
                            verbose = TRUE) {
  
  # Get number of replications per treatment
  n_per_trt <- table(data[[treatment]])
  r <- mean(n_per_trt)  # Assumes balanced design
  
  results <- list()
  
  if (verbose) {
    cat("\n\n")
    cat("                        POST-HOC COMPARISONS                          \n")
    cat("\n")
  }
  
  # LSD Test
  if (posthoc %in% c("lsd", "all")) {
    if (verbose) {
      cat("\n--- Fisher's Least Significant Difference (LSD) Test ---\n")

      if (requireNamespace("agricolae", quietly = TRUE)) {
        lsd_result <- agricolae::LSD.test(model, treatment, alpha = alpha)

        cat(sprintf("\nLSD value ( = %.2f): %.3f\n", alpha, lsd_result$statistics$LSD))
        cat("\nTreatment Groups:\n")

        # Create formatted output
        groups_df <- lsd_result$groups
        groups_df$treatment <- rownames(groups_df)
        groups_df <- groups_df[order(-groups_df[, 1]), ]

        cat("\n")
        print(groups_df[, c("treatment", names(groups_df)[1], "groups")], row.names = FALSE)

        results$lsd <- lsd_result
      } else {
        # Manual LSD calculation
        t_value <- qt(1 - alpha/2, df_error)
        lsd_value <- t_value * sqrt(2 * mse / r)

        cat(sprintf("\nLSD value ( = %.2f): %.3f\n", alpha, lsd_value))
        cat("(Install 'agricolae' package for grouping letters)\n")

        results$lsd <- list(value = lsd_value, t_value = t_value)
      }
    }
  }
  
  # Duncan's Multiple Range Test
  if (posthoc %in% c("duncan", "all")) {
    if (verbose) {
      cat("\n--- Duncan's Multiple Range Test (DMRT) ---\n")

      if (requireNamespace("agricolae", quietly = TRUE)) {
        duncan_result <- agricolae::duncan.test(model, treatment, alpha = alpha)

        cat("\nTreatment Groups:\n")
        groups_df <- duncan_result$groups
        groups_df$treatment <- rownames(groups_df)
        groups_df <- groups_df[order(-groups_df[, 1]), ]

        cat("\n")
        print(groups_df[, c("treatment", names(groups_df)[1], "groups")], row.names = FALSE)

        results$duncan <- duncan_result
      } else {
        cat("(Install 'agricolae' package for Duncan's test)\n")
      }
    }
  }
  
  # Tukey's HSD Test
  if (posthoc %in% c("tukey", "all")) {
    if (verbose) {
      cat("\n--- Tukey's Honest Significant Difference (HSD) Test ---\n")

      if (requireNamespace("agricolae", quietly = TRUE)) {
        hsd_result <- agricolae::HSD.test(model, treatment, alpha = alpha)

        cat(sprintf("\nHSD value ( = %.2f): %.3f\n", alpha, hsd_result$statistics$HSD))
        cat("\nTreatment Groups:\n")

        groups_df <- hsd_result$groups
        groups_df$treatment <- rownames(groups_df)
        groups_df <- groups_df[order(-groups_df[, 1]), ]

        cat("\n")
        print(groups_df[, c("treatment", names(groups_df)[1], "groups")], row.names = FALSE)

        results$tukey <- hsd_result
      } else {
        # Manual Tukey calculation
        q_value <- qtukey(1 - alpha, nlevels(data[[treatment]]), df_error)
        hsd_value <- q_value * sqrt(mse / r)

        cat(sprintf("\nHSD value ( = %.2f): %.3f\n", alpha, hsd_value))
        cat("(Install 'agricolae' package for grouping letters)\n")

        results$tukey <- list(value = hsd_value, q_value = q_value)
      }
    }
  }
  
  # Student-Newman-Keuls (SNK) Test
  if (posthoc %in% c("snk", "all")) {
    if (verbose) {
      cat("\n--- Student-Newman-Keuls (SNK) Test ---\n")

      if (requireNamespace("agricolae", quietly = TRUE)) {
        snk_result <- agricolae::SNK.test(model, treatment, alpha = alpha)

        cat("\nTreatment Groups:\n")
        groups_df <- snk_result$groups
        groups_df$treatment <- rownames(groups_df)
        groups_df <- groups_df[order(-groups_df[, 1]), ]

        cat("\n")
        print(groups_df[, c("treatment", names(groups_df)[1], "groups")], row.names = FALSE)

        results$snk <- snk_result
      } else {
        cat("(Install 'agricolae' package for SNK test)\n")
      }
    }
  }
  
  # Scheffe's Test
  if (posthoc %in% c("scheffe", "all")) {
    if (verbose) {
      cat("\n--- Scheff's Test ---\n")

      if (requireNamespace("agricolae", quietly = TRUE)) {
        scheffe_result <- agricolae::scheffe.test(model, treatment, alpha = alpha)

        cat("\nTreatment Groups:\n")
        groups_df <- scheffe_result$groups
        groups_df$treatment <- rownames(groups_df)
        groups_df <- groups_df[order(-groups_df[, 1]), ]

        cat("\n")
        print(groups_df[, c("treatment", names(groups_df)[1], "groups")], row.names = FALSE)

        results$scheffe <- scheffe_result
      } else {
        cat("(Install 'agricolae' package for Scheff test)\n")
      }
    }
  }
  
  # Dunnett's Test (comparison with control)
  if (posthoc %in% c("dunnett", "all")) {
    if (verbose) {
      cat("\n--- Dunnett's Test (vs Control) ---\n")

      if (requireNamespace("multcomp", quietly = TRUE)) {
        # Assuming first level is control
        dunnett_result <- multcomp::glht(model, linfct = multcomp::mcp(
          treatment = "Dunnett"))
        dunnett_summary <- summary(dunnett_result)
        print(dunnett_summary)
        results$dunnett <- dunnett_result
      } else {
        cat("(Install 'multcomp' package for Dunnett's test)\n")
      }
    }
  }
  
  # Bonferroni Correction
  if (posthoc %in% c("bonferroni", "all")) {
    if (verbose) {
      cat("\n--- Bonferroni-Corrected Comparisons ---\n")
    }
    
    t <- nlevels(data[[treatment]])
    n_comparisons <- t * (t - 1) / 2
    alpha_bonf <- alpha / n_comparisons
    t_bonf <- qt(1 - alpha_bonf/2, df_error)
    bonf_value <- t_bonf * sqrt(2 * mse / r)
    
    if (verbose) {
      cat(sprintf("Number of pairwise comparisons: %d\n", n_comparisons))
      cat(sprintf("Bonferroni-corrected : %.5f\n", alpha_bonf))
      cat(sprintf("Critical t-value: %.3f\n", t_bonf))
      cat(sprintf("Bonferroni CD: %.3f\n", bonf_value))
    }
    
    results$bonferroni <- list(
      n_comparisons = n_comparisons,
      alpha_corrected = alpha_bonf,
      t_value = t_bonf,
      cd_value = bonf_value
    )
  }
  
  return(invisible(results))
}


#' Check ANOVA Assumptions
#'
#' @description
#' Tests assumptions of ANOVA: normality and homogeneity of variances.
#'
#' @param model ANOVA model object
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return List containing diagnostic test results
#' @export
check_assumptions <- function(model,
                            verbose = TRUE) {
  
  residuals <- residuals(model)
  fitted_values <- fitted(model)
  
  if (verbose) {
    cat("\n\n")
    cat("                    ANOVA ASSUMPTION CHECKS                           \n")
    cat("\n")
  }
  
  results <- list()
  
  # 1. Normality Test (Shapiro-Wilk)
  if (verbose) {
    cat("\n--- Normality of Residuals (Shapiro-Wilk Test) ---\n")

    if (length(residuals) <= 5000) {
      shapiro_test <- shapiro.test(residuals)
      cat(sprintf("W statistic: %.4f\n", shapiro_test$statistic))
      cat(sprintf("P-value: %.4f\n", shapiro_test$p.value))

      if (shapiro_test$p.value >= 0.05) {
        cat("Conclusion: Residuals are normally distributed (p >= 0.05)\n")
      } else {
        cat("WARNING: Residuals may not be normally distributed (p < 0.05)\n")
        cat("Consider data transformation (log, sqrt, arcsin)\n")
      }

      results$shapiro <- shapiro_test
    } else {
      cat("Sample size too large for Shapiro-Wilk test (n > 5000)\n")
      cat("Using Anderson-Darling test instead\n")
    }
  }
  
  # 2. Homogeneity of Variances (Levene's Test)
  if (verbose) {
    cat("\n--- Homogeneity of Variances (Bartlett's Test) ---\n")
  }
  
  model_data <- model$model
  response_name <- names(model_data)[1]
  
  # Find factor columns
  factor_cols <- sapply(model_data, is.factor)
  if (sum(factor_cols) > 0) {
    treatment_col <- names(model_data)[factor_cols][1]
    
    bartlett_test <- bartlett.test(model_data[[response_name]] ~ model_data[[treatment_col]])
    
    if (verbose) {
      cat(sprintf("Bartlett's K-squared: %.4f\n", bartlett_test$statistic))
      cat(sprintf("df: %d\n", bartlett_test$parameter))
      cat(sprintf("P-value: %.4f\n", bartlett_test$p.value))

      if (bartlett_test$p.value >= 0.05) {
        cat("Conclusion: Variances are homogeneous (p >= 0.05)\n")
      } else {
        cat("WARNING: Variances may be heterogeneous (p < 0.05)\n")
        cat("Consider Welch's ANOVA or data transformation\n")
      }
    }
    
    results$bartlett <- bartlett_test
  }
  
  # 3. Independence check (Durbin-Watson)
  if (verbose) {
    cat("\n--- Independence of Residuals (Durbin-Watson) ---\n")
  }
  
  # Calculate Durbin-Watson statistic manually
  n <- length(residuals)
  dw_num <- sum((residuals[-1] - residuals[-n])^2)
  dw_den <- sum(residuals^2)
  dw_stat <- dw_num / dw_den
  
  if (verbose) {
    cat(sprintf("Durbin-Watson statistic: %.4f\n", dw_stat))

    if (dw_stat > 1.5 && dw_stat < 2.5) {
      cat("Conclusion: No significant autocorrelation detected\n")
    } else {
      cat("NOTE: Possible autocorrelation in residuals\n")
    }
  }
  
  results$durbin_watson <- dw_stat
  
  # 4. Outlier detection
  if (verbose) {
    cat("\n--- Outlier Detection ---\n")
  }
  
  # Standardized residuals
  std_residuals <- scale(residuals)
  outliers <- which(abs(std_residuals) > 3)
  
  if (length(outliers) > 0) {
    if (verbose) {
      cat(sprintf("Potential outliers (|standardized residual| > 3): %d observations\n", 
                  length(outliers)))
      cat("Observation indices:", paste(outliers, collapse = ", "), "\n")
    }
  } else {
    if (verbose) {
      cat("No extreme outliers detected (all |standardized residuals| <= 3)\n")
    }
  }
  
  results$outliers <- outliers
  
  return(invisible(results))
}


#' Factorial ANOVA (Two-Factor)
#'
#' @description
#' Performs two-factor factorial ANOVA with interaction analysis.
#'
#' @param data Data frame containing the data
#' @param response Name of the response variable
#' @param factor1 Name of first factor (A)
#' @param factor2 Name of second factor (B)
#' @param replication Name of replication factor (optional)
#' @param posthoc Post-hoc test method
#' @param alpha Significance level
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return List containing ANOVA results
#'
#' @examples
#' data <- expand.grid(
#'   rep = 1:4,
#'   nitrogen = c("N0", "N40", "N80", "N120"),
#'   phosphorus = c("P0", "P30", "P60")
#' )
#' data$yield <- rnorm(nrow(data), 1200, 150)
#' 
#' anova_factorial(data, response = "yield", 
#'                 factor1 = "nitrogen", factor2 = "phosphorus",
#'                 replication = "rep")
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
anova_factorial <- function(data, response, factor1, factor2, 
                            replication = NULL, posthoc = "lsd", alpha = 0.05,
                            verbose = TRUE) {
  
  # Convert to factors
  data[[factor1]] <- as.factor(data[[factor1]])
  data[[factor2]] <- as.factor(data[[factor2]])
  
  # Get dimensions
  a <- nlevels(data[[factor1]])
  b <- nlevels(data[[factor2]])
  N <- nrow(data)
  
  # Build formula
  if (!is.null(replication)) {
    data[[replication]] <- as.factor(data[[replication]])
    r <- nlevels(data[[replication]])
    formula_fact <- as.formula(paste(response, "~", replication, "+", 
                                      factor1, "*", factor2))
  } else {
    r <- N / (a * b)
    formula_fact <- as.formula(paste(response, "~", factor1, "*", factor2))
  }
  
  # Fit model
  model <- aov(formula_fact, data = data)
  anova_table <- anova(model)
  
  # Print header
  if (verbose) {
    cat("\n")
    cat("\n")
    cat("                TWO-FACTOR FACTORIAL ANOVA                            \n")
    cat("\n")
    cat(" Response Variable    :", sprintf("%-46s", response), "\n")
    cat(" Factor A             :", sprintf("%-46s", factor1), "\n")
    cat(" Factor B             :", sprintf("%-46s", factor2), "\n")
    cat(" Levels of A          :", sprintf("%-46d", a), "\n")
    cat(" Levels of B          :", sprintf("%-46d", b), "\n")
    cat(" Replications         :", sprintf("%-46.0f", r), "\n")
    cat(" Total Observations   :", sprintf("%-46d", N), "\n")
    cat("\n")

    cat("\n\n")
    cat("                        ANALYSIS OF VARIANCE                          \n")
    cat("\n")

    print(anova_table)
  }
  
  # Statistics
  grand_mean <- mean(data[[response]], na.rm = TRUE)
  ms_error <- anova_table$`Mean Sq`[nrow(anova_table)]
  df_error <- anova_table$Df[nrow(anova_table)]
  cv <- sqrt(ms_error) / grand_mean * 100
  
  # SE and CD
  se_a <- sqrt(ms_error / (b * r))
  cd_a <- qt(0.975, df_error) * sqrt(2 * ms_error / (b * r))
  
  se_b <- sqrt(ms_error / (a * r))
  cd_b <- qt(0.975, df_error) * sqrt(2 * ms_error / (a * r))
  
  se_ab <- sqrt(ms_error / r)
  cd_ab <- qt(0.975, df_error) * sqrt(2 * ms_error / r)
  
  if (verbose) {
    cat("\n\n")
    cat("                        SUMMARY STATISTICS                            \n")
    cat("\n")
    cat(sprintf(" Grand Mean          : %10.2f                                     \n", grand_mean))
    cat(sprintf(" CV (%%)              : %10.2f                                     \n", cv))
    cat("\n")
    cat(sprintf(" Factor A: SE(m)=%7.3f  SE(d)=%7.3f  CD(5%%)=%7.3f           \n",
                se_a, se_a*sqrt(2), cd_a))
    cat(sprintf(" Factor B: SE(m)=%7.3f  SE(d)=%7.3f  CD(5%%)=%7.3f           \n",
                se_b, se_b*sqrt(2), cd_b))
    cat(sprintf(" A  B   : SE(m)=%7.3f  SE(d)=%7.3f  CD(5%%)=%7.3f           \n",
                se_ab, se_ab*sqrt(2), cd_ab))
    cat("\n")
  }
  
  # Means
  mean_a <- aggregate(data[[response]], by = list(data[[factor1]]), FUN = mean)
  names(mean_a) <- c("Level", "Mean")
  
  mean_b <- aggregate(data[[response]], by = list(data[[factor2]]), FUN = mean)
  names(mean_b) <- c("Level", "Mean")
  
  mean_ab <- aggregate(data[[response]], 
                        by = list(data[[factor1]], data[[factor2]]), 
                        FUN = mean)
  names(mean_ab) <- c("A", "B", "Mean")
  
  if (verbose) {
    cat("\n--- Factor A Means ---\n")
    print(mean_a, row.names = FALSE)

    cat("\n--- Factor B Means ---\n")
    print(mean_b, row.names = FALSE)

    cat("\n--- A  B Interaction Means ---\n")
  }
  ab_table <- reshape(mean_ab, idvar = "A", timevar = "B", direction = "wide")
  if (verbose) {
    print(ab_table, row.names = FALSE)
  }
  
  # Return
  result <- list(
    design = "Two-Factor Factorial",
    anova_table = anova_table,
    model = model,
    grand_mean = grand_mean,
    cv = cv,
    se_a = se_a, cd_a = cd_a,
    se_b = se_b, cd_b = cd_b,
    se_ab = se_ab, cd_ab = cd_ab,
    mean_a = mean_a,
    mean_b = mean_b,
    mean_ab = ab_table
  )
  
  class(result) <- c("aridagri_factorial", "list")
  return(invisible(result))
}

# NOTE: anova_factorial_3way is defined in design_additional.R


#' Latin Square Design ANOVA
#'
#' @description
#' Performs ANOVA for Latin Square Design with row and column blocking.
#'
#' @param data Data frame containing the data
#' @param response Name of the response variable
#' @param treatment Name of treatment factor
#' @param row Name of row factor
#' @param column Name of column factor
#' @param posthoc Post-hoc test method
#' @param alpha Significance level
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return List containing ANOVA results
#'
#' @examples
#' # 5x5 Latin Square
#' data <- data.frame(
#'   row = rep(1:5, each = 5),
#'   col = rep(1:5, 5),
#'   treatment = c("A","B","C","D","E", "B","C","D","E","A",
#'                 "C","D","E","A","B", "D","E","A","B","C",
#'                 "E","A","B","C","D"),
#'   yield = rnorm(25, 1200, 100)
#' )
#' anova_latin(data, response = "yield", treatment = "treatment",
#'             row = "row", column = "col")
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
anova_latin <- function(data, response, treatment, row, column,
                         posthoc = "lsd", alpha = 0.05,
                            verbose = TRUE) {
  
  # Convert to factors
  data[[treatment]] <- as.factor(data[[treatment]])
  data[[row]] <- as.factor(data[[row]])
  data[[column]] <- as.factor(data[[column]])
  
  # Get dimensions
  t <- nlevels(data[[treatment]])
  N <- nrow(data)
  
  # Fit model
  formula_lsd <- as.formula(paste(response, "~", row, "+", column, "+", treatment))
  model <- aov(formula_lsd, data = data)
  anova_table <- anova(model)
  
  # Print header
  if (verbose) {
    cat("\n")
    cat("\n")
    cat("                   LATIN SQUARE DESIGN - ANOVA                        \n")
    cat("\n")
    cat(" Response Variable    :", sprintf("%-46s", response), "\n")
    cat(" Treatment Factor     :", sprintf("%-46s", treatment), "\n")
    cat(" Row Factor           :", sprintf("%-46s", row), "\n")
    cat(" Column Factor        :", sprintf("%-46s", column), "\n")
    cat(" Size (t  t)         :", sprintf("%-46s", paste0(t, "  ", t)), "\n")
    cat(" Total Observations   :", sprintf("%-46d", N), "\n")
    cat("\n")

    cat("\n\n")
    cat("                        ANALYSIS OF VARIANCE                          \n")
    cat("\n")

    print(anova_table)
  }
  
  # Statistics
  grand_mean <- mean(data[[response]], na.rm = TRUE)
  ms_error <- anova_table$`Mean Sq`[4]
  df_error <- anova_table$Df[4]
  cv <- sqrt(ms_error) / grand_mean * 100
  
  se_mean <- sqrt(ms_error / t)
  cd_5 <- qt(0.975, df_error) * sqrt(2 * ms_error / t)
  
  if (verbose) {
    cat("\n\n")
    cat("                        SUMMARY STATISTICS                            \n")
    cat("\n")
    cat(sprintf(" Grand Mean          : %10.2f                                     \n", grand_mean))
    cat(sprintf(" CV (%%)              : %10.2f                                     \n", cv))
    cat(sprintf(" SE(m)               : %10.3f                                     \n", se_mean))
    cat(sprintf(" SE(d)               : %10.3f                                     \n", se_mean * sqrt(2)))
    cat(sprintf(" CD (5%%)             : %10.3f                                     \n", cd_5))
    cat("\n")
  }
  
  # Treatment means
  trt_means <- aggregate(data[[response]], by = list(data[[treatment]]), FUN = mean)
  names(trt_means) <- c("Treatment", "Mean")
  trt_means <- trt_means[order(-trt_means$Mean), ]
  
  if (verbose) {
    cat("\n--- Treatment Means ---\n")
    print(trt_means, row.names = FALSE)
  }
  
  # Post-hoc
  posthoc_results <- perform_posthoc(model, data, response, treatment,
                                      ms_error, df_error, posthoc, alpha)
  
  # Return
  result <- list(
    design = "Latin Square Design",
    anova_table = anova_table,
    model = model,
    grand_mean = grand_mean,
    cv = cv,
    se_mean = se_mean,
    cd = cd_5,
    treatment_means = trt_means,
    posthoc = posthoc_results
  )
  
  class(result) <- c("aridagri_latin", "list")
  return(invisible(result))
}
