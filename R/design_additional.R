#' ============================================================================
#' ADDITIONAL EXPERIMENTAL DESIGN FUNCTIONS
#' Package: aridagri
#' Authors: Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram Choudhary
#' ============================================================================

#' Strip Plot Design ANOVA
#'
#' @description
#' Performs ANOVA for Strip Plot (Strip-Split) Design where two factors
#' are applied in horizontal and vertical strips.
#'
#' @param data Data frame containing the experimental data
#' @param response Name of the response variable
#' @param horizontal_factor Factor applied in horizontal strips (A)
#' @param vertical_factor Factor applied in vertical strips (B)
#' @param replication Name of replication factor
#' @param alpha Significance level
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return List containing ANOVA results
#'
#' @examples
#' data <- expand.grid(
#'   rep = 1:4,
#'   irrigation = c("I1", "I2", "I3"),
#'   tillage = c("CT", "MT", "ZT")
#' )
#' data$yield <- rnorm(nrow(data), 1200, 150)
#' anova_strip(data, response = "yield", horizontal_factor = "irrigation",
#'             vertical_factor = "tillage", replication = "rep")
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram Choudhary
anova_strip <- function(data, response, horizontal_factor, vertical_factor,
                         replication, alpha = 0.05,
                            verbose = TRUE) {
  
  # Convert to factors
  data[[horizontal_factor]] <- as.factor(data[[horizontal_factor]])
  data[[vertical_factor]] <- as.factor(data[[vertical_factor]])
  data[[replication]] <- as.factor(data[[replication]])
  
  # Get dimensions
  a <- nlevels(data[[horizontal_factor]])
  b <- nlevels(data[[vertical_factor]])
  r <- nlevels(data[[replication]])
  N <- nrow(data)
  
  # Grand mean
  grand_mean <- mean(data[[response]], na.rm = TRUE)
  
  # Fit model for strip plot
  formula_strip <- as.formula(paste(response, "~", replication, "+",
                                     horizontal_factor, "+", replication, ":", horizontal_factor, "+",
                                     vertical_factor, "+", replication, ":", vertical_factor, "+",
                                     horizontal_factor, ":", vertical_factor))
  
  model <- aov(formula_strip, data = data)
  anova_table <- anova(model)
  
  # Print header
  if (verbose) {
    cat("\n")
    cat("\n")
    cat("                      STRIP PLOT DESIGN - ANOVA                               \n")
    cat("\n")
    cat(" Response Variable     :", sprintf("%-52s", response), "\n")
    cat(" Horizontal Factor (A) :", sprintf("%-52s", horizontal_factor), "\n")
    cat(" Vertical Factor (B)   :", sprintf("%-52s", vertical_factor), "\n")
    cat(" Replications          :", sprintf("%-52d", r), "\n")
    cat(" Levels of A           :", sprintf("%-52d", a), "\n")
    cat(" Levels of B           :", sprintf("%-52d", b), "\n")
    cat("\n")

    cat("\n\n")
    cat("                        ANALYSIS OF VARIANCE                                  \n")
    cat("\n")

    print(anova_table)
  }
  
  # Statistics
  n_rows <- nrow(anova_table)
  MS_error <- anova_table$`Mean Sq`[n_rows]
  cv <- sqrt(MS_error) / grand_mean * 100
  
  if (verbose) {
    cat(sprintf("\nGrand Mean: %.2f\n", grand_mean))
    cat(sprintf("CV: %.2f%%\n", cv))
  }
  
  # Factor means
  mean_a <- aggregate(data[[response]], by = list(data[[horizontal_factor]]), FUN = mean)
  names(mean_a) <- c("Level", "Mean")
  
  mean_b <- aggregate(data[[response]], by = list(data[[vertical_factor]]), FUN = mean)
  names(mean_b) <- c("Level", "Mean")
  
  if (verbose) {
    cat("\n--- Horizontal Factor Means ---\n")
    print(mean_a, row.names = FALSE)
    cat("\n--- Vertical Factor Means ---\n")
    print(mean_b, row.names = FALSE)
  }
  
  # Return
  result <- list(
    design = "Strip Plot Design",
    anova_table = anova_table,
    model = model,
    grand_mean = grand_mean,
    cv = cv,
    horizontal_means = mean_a,
    vertical_means = mean_b
  )
  
  class(result) <- c("aridagri_strip", "list")
  return(invisible(result))
}


#' Three-Factor Factorial ANOVA
#'
#' @description
#' Performs three-factor factorial ANOVA with all interactions.
#'
#' @param data Data frame containing the data
#' @param response Name of the response variable
#' @param factor1 Name of first factor (A)
#' @param factor2 Name of second factor (B)
#' @param factor3 Name of third factor (C)
#' @param replication Name of replication factor (optional)
#' @param alpha Significance level
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return List containing ANOVA results
#'
#' @examples
#' data <- expand.grid(rep = 1:3, A = c("A1", "A2"), B = c("B1", "B2"), C = c("C1", "C2"))
#' data$yield <- rnorm(nrow(data), 1200, 150)
#' anova_factorial_3way(data, "yield", "A", "B", "C", "rep")
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram Choudhary
anova_factorial_3way <- function(data, response, factor1, factor2, factor3,
                                  replication = NULL, alpha = 0.05,
                            verbose = TRUE) {
  
  # Convert to factors
  data[[factor1]] <- as.factor(data[[factor1]])
  data[[factor2]] <- as.factor(data[[factor2]])
  data[[factor3]] <- as.factor(data[[factor3]])
  
  # Get dimensions
  a <- nlevels(data[[factor1]])
  b <- nlevels(data[[factor2]])
  c <- nlevels(data[[factor3]])
  N <- nrow(data)
  
  # Build formula
  if (!is.null(replication)) {
    data[[replication]] <- as.factor(data[[replication]])
    r <- nlevels(data[[replication]])
    formula_fact <- as.formula(paste(response, "~", replication, "+", 
                                      factor1, "*", factor2, "*", factor3))
  } else {
    r <- N / (a * b * c)
    formula_fact <- as.formula(paste(response, "~", factor1, "*", factor2, "*", factor3))
  }
  
  # Fit model
  model <- aov(formula_fact, data = data)
  anova_table <- anova(model)
  
  # Print header
  if (verbose) {
    cat("\n")
    cat("\n")
    cat("              THREE-FACTOR FACTORIAL ANOVA                            \n")
    cat("\n")
    cat(" Response Variable    :", sprintf("%-46s", response), "\n")
    cat(" Factor A             :", sprintf("%-46s", factor1), "\n")
    cat(" Factor B             :", sprintf("%-46s", factor2), "\n")
    cat(" Factor C             :", sprintf("%-46s", factor3), "\n")
    cat(sprintf(" Levels: A=%d, B=%d, C=%d", a, b, c))
    cat(sprintf("%36s\n", ""))
    cat(" Total Observations   :", sprintf("%-46d", N), "\n")
    cat("\n")

    cat("\n\n")
    cat("                        ANALYSIS OF VARIANCE                          \n")
    cat("\n")

    print(anova_table)
  }
  
  # Grand mean and CV
  grand_mean <- mean(data[[response]], na.rm = TRUE)
  ms_error <- anova_table$`Mean Sq`[nrow(anova_table)]
  cv <- sqrt(ms_error) / grand_mean * 100
  
  if (verbose) {
    cat(sprintf("\nGrand Mean: %.2f\n", grand_mean))
    cat(sprintf("CV: %.2f%%\n", cv))
  }
  
  # Return
  result <- list(
    design = "Three-Factor Factorial",
    anova_table = anova_table,
    model = model,
    grand_mean = grand_mean,
    cv = cv
  )
  
  class(result) <- c("aridagri_factorial3", "list")
  return(invisible(result))
}


#' Augmented Block Design ANOVA
#'
#' @description
#' Performs ANOVA for Augmented Randomized Block Design where checks are
#' replicated and test entries appear once.
#'
#' @param data Data frame containing the data
#' @param response Name of response variable
#' @param genotype Name of genotype/entry column
#' @param block Name of block column
#' @param check_names Vector of check variety names
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return List containing ANOVA results and adjusted means
#'
#' @examples
#' data <- data.frame(
#'   block = rep(1:5, each = 7),
#'   genotype = c(rep(c("C1","C2","C3"), 5), paste0("T", 1:20)),
#'   yield = rnorm(35, 1200, 150)
#' )
#' # Note: This is simplified example
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram Choudhary
anova_augmented <- function(data, response, genotype, block, check_names,
                            verbose = TRUE) {
  
  # Convert to factors
  data[[genotype]] <- as.factor(data[[genotype]])
  data[[block]] <- as.factor(data[[block]])
  
  # Identify checks and test entries
  data$Entry_Type <- ifelse(data[[genotype]] %in% check_names, "Check", "Test")
  
  # Get dimensions
  n_checks <- length(check_names)
  n_blocks <- nlevels(data[[block]])
  n_test <- nlevels(data[[genotype]]) - n_checks
  N <- nrow(data)
  
  if (verbose) {
    cat("\n")
    cat("\n")
    cat("              AUGMENTED RANDOMIZED BLOCK DESIGN - ANOVA               \n")
    cat("\n")
    cat(" Response Variable   :", sprintf("%-47s", response), "\n")
    cat(" Number of Blocks    :", sprintf("%-47d", n_blocks), "\n")
    cat(" Number of Checks    :", sprintf("%-47d", n_checks), "\n")
    cat(" Number of Test Entries:", sprintf("%-44d", n_test), "\n")
    cat(" Total Observations  :", sprintf("%-47d", N), "\n")
    cat("\n")
  }
  
  # Separate check data
  check_data <- data[data$Entry_Type == "Check", ]
  
  # Analyze checks (RBD)
  formula_checks <- as.formula(paste(response, "~", block, "+", genotype))
  model_checks <- aov(formula_checks, data = check_data)
  anova_checks <- anova(model_checks)
  
  if (verbose) {
    cat("\n\n")
    cat("                  ANOVA FOR CHECK VARIETIES (RBD)                     \n")
    cat("\n")
    print(anova_checks)
  }
  
  # Error mean square
  mse_checks <- anova_checks$`Mean Sq`[nrow(anova_checks)]
  df_error <- anova_checks$Df[nrow(anova_checks)]
  
  # Block effects
  block_means <- tapply(check_data[[response]], check_data[[block]], mean)
  grand_mean_checks <- mean(check_data[[response]])
  block_effects <- block_means - grand_mean_checks
  
  # Calculate CD
  cd_check <- qt(0.975, df_error) * sqrt(2 * mse_checks / n_blocks)
  
  if (verbose) {
    cat("\n\n")
    cat("                        SUMMARY                                       \n")
    cat("\n")
    cat(sprintf(" Error MSE                          : %10.2f                      \n", mse_checks))
    cat(sprintf(" CD for comparing two checks        : %10.3f                      \n", cd_check))
    cat("\n")
  }
  
  # Return
  result <- list(
    design = "Augmented Block Design",
    anova_checks = anova_checks,
    mse = mse_checks,
    df_error = df_error,
    block_effects = block_effects,
    cd_check = cd_check
  )
  
  class(result) <- c("aridagri_augmented", "list")
  return(invisible(result))
}


#' Alpha Lattice Design ANOVA
#'
#' @description
#' Performs ANOVA for Alpha Lattice (Resolvable Incomplete Block) Design.
#'
#' @param data Data frame containing the data
#' @param response Name of response variable
#' @param genotype Name of genotype column
#' @param replication Name of replication column
#' @param block Name of incomplete block column (nested within replication)
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return List containing ANOVA results
#'
#' @examples
#' data <- expand.grid(rep = 1:2, block = 1:5, entry = 1:4)
#' data$genotype <- paste0("G", 1:nrow(data))
#' data$yield <- rnorm(nrow(data), 1200, 150)
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram Choudhary
anova_alpha_lattice <- function(data, response, genotype, replication, block,
                            verbose = TRUE) {
  
  # Convert to factors
  data[[genotype]] <- as.factor(data[[genotype]])
  data[[replication]] <- as.factor(data[[replication]])
  data[[block]] <- as.factor(data[[block]])
  
  # Get dimensions
  g <- nlevels(data[[genotype]])
  r <- nlevels(data[[replication]])
  N <- nrow(data)
  
  if (verbose) {
    cat("\n")
    cat("\n")
    cat("                  ALPHA LATTICE DESIGN - ANOVA                        \n")
    cat("\n")
    cat(" Response Variable   :", sprintf("%-47s", response), "\n")
    cat(" Number of Genotypes :", sprintf("%-47d", g), "\n")
    cat(" Number of Replications:", sprintf("%-44d", r), "\n")
    cat(" Total Observations  :", sprintf("%-47d", N), "\n")
    cat("\n")
  }
  
  # Fit model: Block nested within replication
  formula_alpha <- as.formula(paste(response, "~", replication, "+",
                                     replication, ":", block, "+", genotype))
  
  model <- aov(formula_alpha, data = data)
  anova_table <- anova(model)
  
  if (verbose) {
    cat("\n\n")
    cat("                        ANALYSIS OF VARIANCE                          \n")
    cat("\n")
    print(anova_table)
  }
  
  # Statistics
  grand_mean <- mean(data[[response]], na.rm = TRUE)
  mse <- anova_table$`Mean Sq`[nrow(anova_table)]
  cv <- sqrt(mse) / grand_mean * 100
  
  if (verbose) {
    cat(sprintf("\nGrand Mean: %.2f\n", grand_mean))
    cat(sprintf("CV: %.2f%%\n", cv))
  }
  
  # Genotype means
  g_means <- aggregate(data[[response]], by = list(data[[genotype]]), FUN = mean)
  names(g_means) <- c("Genotype", "Mean")
  g_means <- g_means[order(-g_means$Mean), ]
  
  if (verbose) {
    cat("\n--- Genotype Means (Top 10) ---\n")
    print(head(g_means, 10), row.names = FALSE)
  }
  
  # Return
  result <- list(
    design = "Alpha Lattice Design",
    anova_table = anova_table,
    model = model,
    grand_mean = grand_mean,
    cv = cv,
    genotype_means = g_means
  )
  
  class(result) <- c("aridagri_alpha", "list")
  return(invisible(result))
}


#' Pooled Split-Split Plot Design ANOVA
#'
#' @description
#' Performs pooled ANOVA for SSPD experiments across multiple environments.
#'
#' @param data Data frame containing combined data
#' @param response Name of the response variable
#' @param main_plot Name of main plot factor
#' @param sub_plot Name of sub-plot factor
#' @param sub_sub_plot Name of sub-sub-plot factor
#' @param environment Name of environment factor
#' @param replication Name of replication factor
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return List containing pooled ANOVA results
#'
#' @examples
#' data <- expand.grid(env = c("E1","E2"), rep = 1:3, A = c("A1","A2"), 
#'                     B = c("B1","B2"), C = c("C1","C2"))
#' data$yield <- rnorm(nrow(data), 1200, 150)
#' anova_sspd_pooled(data, "yield", "A", "B", "C", "env", "rep")
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram Choudhary
anova_sspd_pooled <- function(data, response, main_plot, sub_plot, sub_sub_plot,
                               environment, replication,
                            verbose = TRUE) {
  
  # Convert to factors
  data[[main_plot]] <- as.factor(data[[main_plot]])
  data[[sub_plot]] <- as.factor(data[[sub_plot]])
  data[[sub_sub_plot]] <- as.factor(data[[sub_sub_plot]])
  data[[environment]] <- as.factor(data[[environment]])
  data[[replication]] <- as.factor(data[[replication]])
  
  # Get dimensions
  a <- nlevels(data[[main_plot]])
  b <- nlevels(data[[sub_plot]])
  c <- nlevels(data[[sub_sub_plot]])
  e <- nlevels(data[[environment]])
  r <- nlevels(data[[replication]])
  N <- nrow(data)
  
  if (verbose) {
    cat("\n")
    cat("\n")
    cat("               POOLED SPLIT-SPLIT PLOT DESIGN ANALYSIS                              \n")
    cat("\n")
    cat(" Main Plot (A)         :", sprintf("%-59s", main_plot), "\n")
    cat(" Sub-Plot (B)          :", sprintf("%-59s", sub_plot), "\n")
    cat(" Sub-Sub-Plot (C)      :", sprintf("%-59s", sub_sub_plot), "\n")
    cat(" Environment (E)       :", sprintf("%-59s", environment), "\n")
    cat(sprintf(" Levels: A=%d, B=%d, C=%d, E=%d, R=%d", a, b, c, e, r))
    cat(sprintf("%39s\n", ""))
    cat(" Total Observations    :", sprintf("%-59d", N), "\n")
    cat("\n")
  }
  
  # Fit pooled model
  formula_pooled <- as.formula(paste(response, "~", 
                                      environment, "*", main_plot, "*", sub_plot, "*", sub_sub_plot, "+",
                                      environment, ":", replication))
  
  model <- aov(formula_pooled, data = data)
  anova_table <- anova(model)
  
  if (verbose) {
    cat("\n\n")
    cat("                        POOLED ANALYSIS OF VARIANCE                                 \n")
    cat("\n")

    print(anova_table)
  }
  
  # Grand mean
  grand_mean <- mean(data[[response]], na.rm = TRUE)
  if (verbose) {
    cat(sprintf("\nGrand Mean: %.2f\n", grand_mean))
  }
  
  # Return
  result <- list(
    design = "Pooled Split-Split Plot Design",
    anova_table = anova_table,
    model = model,
    grand_mean = grand_mean
  )
  
  class(result) <- c("aridagri_sspd_pooled", "list")
  return(invisible(result))
}
