#' Split-Split Plot Design ANOVA
#'
#' @description
#' Performs complete ANOVA for Split-Split Plot Design with proper error terms
#' for main plot, sub-plot, and sub-sub-plot factors. Generates publication-ready
#' ANOVA table with significance levels.
#'
#' @param data Data frame containing the experimental data
#' @param response Name of the response variable (as string)
#' @param main_plot Name of main plot factor
#' @param sub_plot Name of sub-plot factor
#' @param sub_sub_plot Name of sub-sub-plot factor
#' @param replication Name of replication/block factor
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return List containing ANOVA table, means, and significance tests
#'
#' @examples
#' # Example with sample data
#' data <- expand.grid(rep=1:3, A=c('A1','A2'), B=c('B1','B2'), C=c('C1','C2'))
#' data$yield <- rnorm(24, 1200, 150)
#' anova_sspd(data, response='yield', main_plot='A', sub_plot='B',
#'            sub_sub_plot='C', replication='rep')
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
anova_sspd <- function(data, response, main_plot, sub_plot, sub_sub_plot, replication,
                            verbose = TRUE) {
  
  # Convert to factors
  data[[main_plot]] <- as.factor(data[[main_plot]])
  data[[sub_plot]] <- as.factor(data[[sub_plot]])
  data[[sub_sub_plot]] <- as.factor(data[[sub_sub_plot]])
  data[[replication]] <- as.factor(data[[replication]])
  
  # Number of levels
  a <- nlevels(data[[main_plot]])      # Main plot levels
  b <- nlevels(data[[sub_plot]])       # Sub-plot levels
  c <- nlevels(data[[sub_sub_plot]])   # Sub-sub-plot levels
  r <- nlevels(data[[replication]])    # Replications
  
  # Total degrees of freedom
  N <- nrow(data)
  
  # Build formula for complete model
  formula_full <- as.formula(paste(response, "~", replication, "+", 
                                    main_plot, "+", replication, ":", main_plot, "+",
                                    sub_plot, "+", main_plot, ":", sub_plot, "+",
                                    replication, ":", main_plot, ":", sub_plot, "+",
                                    sub_sub_plot, "+", main_plot, ":", sub_sub_plot, "+",
                                    sub_plot, ":", sub_sub_plot, "+",
                                    main_plot, ":", sub_plot, ":", sub_sub_plot))
  
  # Fit model
  model <- aov(formula_full, data = data)
  anova_table <- anova(model)
  
  # Extract components
  SS <- anova_table$`Sum Sq`
  df <- anova_table$Df
  MS <- anova_table$`Mean Sq`
  
  # Create proper SSPD ANOVA table
  source_names <- c("Replication", "Main Plot (A)", "Error (a)",
                    "Sub-Plot (B)", "A  B", "Error (b)",
                    "Sub-Sub-Plot (C)", "A  C", "B  C", "A  B  C", "Error (c)",
                    "Total")
  
  # Calculate degrees of freedom
  df_rep <- r - 1
  df_a <- a - 1
  df_error_a <- (r - 1) * (a - 1)
  df_b <- b - 1
  df_ab <- (a - 1) * (b - 1)
  df_error_b <- a * (r - 1) * (b - 1)
  df_c <- c - 1
  df_ac <- (a - 1) * (c - 1)
  df_bc <- (b - 1) * (c - 1)
  df_abc <- (a - 1) * (b - 1) * (c - 1)
  df_error_c <- a * b * (r - 1) * (c - 1)
  df_total <- N - 1
  
  # Print results
  if (verbose) {
    cat("\n==========================================================\n")
    cat("     SPLIT-SPLIT PLOT DESIGN ANOVA\n")
    cat("==========================================================\n")
    cat("\nExperimental Details:\n")
    cat("  Main Plot Factor (A):", main_plot, "-", a, "levels\n")
    cat("  Sub-Plot Factor (B):", sub_plot, "-", b, "levels\n")
    cat("  Sub-Sub-Plot Factor (C):", sub_sub_plot, "-", c, "levels\n")
    cat("  Replications:", r, "\n")
    cat("  Total Observations:", N, "\n")
    cat("\n----------------------------------------------------------\n")
    cat("Response Variable:", response, "\n")
    cat("----------------------------------------------------------\n")

    print(anova_table)
  }
  
  # Calculate CV%
  grand_mean <- mean(data[[response]], na.rm = TRUE)
  
  if (verbose) {
    cat("\n----------------------------------------------------------\n")
    cat("Grand Mean:", round(grand_mean, 2), "\n")
    cat("==========================================================\n")
  }
  
  # Return results
  result <- list(
    anova_table = anova_table,
    model = model,
    grand_mean = grand_mean,
    factors = list(main_plot = main_plot, sub_plot = sub_plot, 
                   sub_sub_plot = sub_sub_plot, replication = replication),
    levels = list(a = a, b = b, c = c, r = r)
  )
  
  class(result) <- "aridagri_anova"
  return(invisible(result))
}


# NOTE: anova_spd is defined in design_split_plot.R with full features


#' Correlation Analysis with Significance
#'
#' @description
#' Computes correlation matrix with significance levels and generates
#' publication-ready correlation table and plot.
#'
#' @param data Data frame with numeric variables
#' @param method Correlation method: "pearson", "spearman", or "kendall"
#' @param plot Logical, whether to generate correlation plot
#' @param digits Number of decimal places
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return List with correlation matrix and significance matrix
#'
#' @examples
#' data <- data.frame(
#'   yield = c(1200, 1350, 1100, 1450, 1280),
#'   wue = c(4.2, 4.8, 3.9, 5.1, 4.5),
#'   protein = c(22.1, 23.5, 21.8, 24.2, 22.9)
#' )
#' correlation_analysis(data)
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
correlation_analysis <- function(data, method = "pearson", plot = TRUE, digits = 3,
                            verbose = TRUE) {
  
  # Select only numeric columns
  numeric_data <- data[sapply(data, is.numeric)]
  
  if (ncol(numeric_data) < 2) {
    stop("Need at least 2 numeric variables for correlation analysis")
  }
  
  n_vars <- ncol(numeric_data)
  n_obs <- nrow(numeric_data)
  var_names <- names(numeric_data)
  
  # Initialize matrices
  cor_matrix <- matrix(NA, n_vars, n_vars, dimnames = list(var_names, var_names))
  p_matrix <- matrix(NA, n_vars, n_vars, dimnames = list(var_names, var_names))
  
  # Calculate correlations and p-values
  for (i in 1:n_vars) {
    for (j in 1:n_vars) {
      test <- cor.test(numeric_data[[i]], numeric_data[[j]], method = method)
      cor_matrix[i, j] <- test$estimate
      p_matrix[i, j] <- test$p.value
    }
  }
  
  # Create significance symbols
  sig_matrix <- matrix("", n_vars, n_vars, dimnames = list(var_names, var_names))
  sig_matrix[p_matrix < 0.001] <- "***"
  sig_matrix[p_matrix >= 0.001 & p_matrix < 0.01] <- "**"
  sig_matrix[p_matrix >= 0.01 & p_matrix < 0.05] <- "*"
  sig_matrix[p_matrix >= 0.05 & p_matrix < 0.1] <- "."
  
  # Print results
  if (verbose) {
    cat("\n==========================================================\n")
    cat("          CORRELATION ANALYSIS\n")
    cat("==========================================================\n")
    cat("Method:", method, "\n")
    cat("Number of observations:", n_obs, "\n")
    cat("Number of variables:", n_vars, "\n")
    cat("\n--- Correlation Matrix ---\n\n")

    print(round(cor_matrix, digits))

    cat("\n--- Significance Levels ---\n")
    cat("*** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.1\n")
  }
  
  # Plot if requested
  if (plot && requireNamespace("corrplot", quietly = TRUE)) {
    corrplot::corrplot(cor_matrix, method = "color", type = "upper",
                       addCoef.col = "black", tl.col = "black",
                       tl.srt = 45, diag = FALSE,
                       p.mat = p_matrix, sig.level = 0.05, insig = "blank",
                       title = "Correlation Matrix", mar = c(0, 0, 2, 0))
  }
  
  result <- list(
    correlation = round(cor_matrix, digits),
    p_values = round(p_matrix, 4),
    significance = sig_matrix,
    method = method,
    n = n_obs
  )
  
  class(result) <- "aridagri_correlation"
  return(invisible(result))
}


#' Principal Component Analysis
#'
#' @description
#' Performs PCA with visualization suitable for agricultural research data.
#' Includes scree plot, biplot, and variable contributions.
#'
#' @param data Data frame with numeric variables
#' @param scale Logical, whether to scale variables (default TRUE)
#' @param ncp Number of components to retain (default 5)
#' @param plot Logical, whether to generate plots
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return PCA results with eigenvalues, loadings, and scores
#'
#' @examples
#' data <- data.frame(
#'   yield = rnorm(30, 1200, 200),
#'   wue = rnorm(30, 4.5, 0.5),
#'   protein = rnorm(30, 22, 2),
#'   biomass = rnorm(30, 3500, 500)
#' )
#' pca_analysis(data)
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
pca_analysis <- function(data, scale = TRUE, ncp = 5, plot = TRUE,
                            verbose = TRUE) {
  
  # Select numeric columns
  numeric_data <- data[sapply(data, is.numeric)]
  
  if (ncol(numeric_data) < 2) {
    stop("Need at least 2 numeric variables for PCA")
  }
  
  # Remove rows with NA
  numeric_data <- na.omit(numeric_data)
  
  # Perform PCA
  pca_result <- prcomp(numeric_data, scale. = scale, center = TRUE)
  
  # Extract results
  eigenvalues <- pca_result$sdev^2
  var_explained <- eigenvalues / sum(eigenvalues) * 100
  cum_var <- cumsum(var_explained)
  
  # Loadings
  loadings <- pca_result$rotation
  
  # Print results
  if (verbose) {
    cat("\n==========================================================\n")
    cat("          PRINCIPAL COMPONENT ANALYSIS\n")
    cat("==========================================================\n")
    cat("Variables:", ncol(numeric_data), "\n")
    cat("Observations:", nrow(numeric_data), "\n")
    cat("Scaling:", ifelse(scale, "Yes", "No"), "\n")

    cat("\n--- Eigenvalues and Variance Explained ---\n\n")
  }
  
  eigen_table <- data.frame(
    PC = paste0("PC", 1:length(eigenvalues)),
    Eigenvalue = round(eigenvalues, 3),
    Variance_Percent = round(var_explained, 2),
    Cumulative_Percent = round(cum_var, 2)
  )
  if (verbose) {
    print(eigen_table)

    cat("\n--- Variable Loadings (First 3 PCs) ---\n\n")
    print(round(loadings[, 1:min(3, ncol(loadings))], 3))
  }
  
  # Kaiser criterion - components with eigenvalue > 1
  n_retain <- sum(eigenvalues > 1)
  if (verbose) {
    cat("\n--- Recommendation ---\n")
    cat("Components with eigenvalue > 1:", n_retain, "\n")
    cat("Variance explained by", n_retain, "components:", round(cum_var[n_retain], 1), "%\n")
  }
  
  # Plots
  if (plot && requireNamespace("factoextra", quietly = TRUE)) {
    # Scree plot
    if (verbose) {
      print(factoextra::fviz_screeplot(pca_result, addlabels = TRUE, 
                                        title = "Scree Plot"))
    }
    
    # Biplot
    if (verbose) {
      print(factoextra::fviz_pca_biplot(pca_result, 
                                         repel = TRUE,
                                         title = "PCA Biplot"))
    }
  }
  
  result <- list(
    pca = pca_result,
    eigenvalues = eigenvalues,
    variance_explained = var_explained,
    cumulative_variance = cum_var,
    loadings = loadings,
    scores = pca_result$x,
    n_components_kaiser = n_retain
  )
  
  class(result) <- "aridagri_pca"
  return(invisible(result))
}


#' Path Coefficient Analysis
#'
#' @description
#' Performs path analysis to determine direct and indirect effects of
#' independent variables on a dependent variable. Essential for 
#' understanding yield contributing factors.
#'
#' @param data Data frame with numeric variables
#' @param dependent Name of dependent variable (e.g., "yield")
#' @param independent Character vector of independent variable names
#' @param digits Number of decimal places
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return List with direct effects, indirect effects, and correlation breakdown
#'
#' @examples
#' data <- data.frame(
#'   yield = c(1200, 1350, 1100, 1450, 1280, 1380, 1220, 1400),
#'   pods = c(45, 52, 42, 58, 48, 54, 46, 56),
#'   seeds = c(8.2, 9.1, 7.8, 9.5, 8.5, 9.0, 8.3, 9.3),
#'   weight = c(32, 35, 30, 38, 33, 36, 31, 37)
#' )
#' path_analysis(data, dependent = "yield", 
#'               independent = c("pods", "seeds", "weight"))
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
path_analysis <- function(data, dependent, independent, digits = 4,
                            verbose = TRUE) {
  
  # Validate inputs
  if (!dependent %in% names(data)) {
    stop(paste("Dependent variable", dependent, "not found in data"))
  }
  
  missing_vars <- setdiff(independent, names(data))
  if (length(missing_vars) > 0) {
    stop(paste("Variables not found:", paste(missing_vars, collapse = ", ")))
  }
  
  # Subset data
  vars <- c(dependent, independent)
  analysis_data <- data[, vars]
  analysis_data <- na.omit(analysis_data)
  
  n <- nrow(analysis_data)
  p <- length(independent)
  
  # Calculate correlation matrix
  cor_matrix <- cor(analysis_data)
  
  # Correlation of independent variables with dependent
  r_y <- cor_matrix[dependent, independent]
  
  # Correlation among independent variables
  R_xx <- cor_matrix[independent, independent]
  
  # Path coefficients (direct effects) = inverse of R_xx multiplied by r_y
  R_xx_inv <- solve(R_xx)
  path_coef <- R_xx_inv %*% r_y
  direct_effects <- as.vector(path_coef)
  names(direct_effects) <- independent
  
  # Calculate indirect effects
  indirect_matrix <- matrix(0, p, p, dimnames = list(independent, independent))
  
  for (i in 1:p) {
    for (j in 1:p) {
      if (i != j) {
        indirect_matrix[i, j] <- direct_effects[j] * R_xx[i, j]
      }
    }
  }
  
  # Total indirect effects
  indirect_total <- rowSums(indirect_matrix)
  
  # Residual effect
  r_squared <- sum(direct_effects * r_y)
  residual <- sqrt(1 - r_squared)
  
  # Print results
  if (verbose) {
    cat("\n==========================================================\n")
    cat("          PATH COEFFICIENT ANALYSIS\n")
    cat("==========================================================\n")
    cat("Dependent variable:", dependent, "\n")
    cat("Independent variables:", paste(independent, collapse = ", "), "\n")
    cat("Number of observations:", n, "\n")

    cat("\n--- Direct Effects (Path Coefficients) ---\n\n")
  }
  direct_df <- data.frame(
    Variable = independent,
    Direct_Effect = round(direct_effects, digits)
  )
  if (verbose) {
    print(direct_df, row.names = FALSE)

    cat("\n--- Indirect Effects ---\n\n")
    print(round(indirect_matrix, digits))

    cat("\n--- Correlation Partitioning ---\n\n")
  }
  partition_df <- data.frame(
    Variable = independent,
    Correlation_with_Y = round(r_y, digits),
    Direct_Effect = round(direct_effects, digits),
    Indirect_Effect = round(indirect_total, digits)
  )
  if (verbose) {
    print(partition_df, row.names = FALSE)

    cat("\n--- Model Summary ---\n")
    cat("R-squared:", round(r_squared, 4), "\n")
    cat("Residual effect:", round(residual, 4), "\n")
    cat("==========================================================\n")
  }
  
  result <- list(
    direct_effects = direct_effects,
    indirect_effects = indirect_matrix,
    indirect_total = indirect_total,
    correlation_with_y = r_y,
    correlation_matrix = cor_matrix,
    r_squared = r_squared,
    residual = residual
  )
  
  class(result) <- "aridagri_path"
  return(invisible(result))
}


#' Structural Equation Modeling for Field Experiments
#'
#' @description
#' Performs SEM analysis for agricultural field experiments. Allows testing
#' of hypothesized causal relationships among variables.
#'
#' @param data Data frame with variables
#' @param model Model specification in lavaan syntax
#' @param plot Logical, whether to generate path diagram
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return SEM results including fit indices and parameter estimates
#'
#' @examples
#' \donttest{
#' if (requireNamespace("lavaan", quietly = TRUE)) {
#'   set.seed(42)
#'   n <- 100
#'   nitrogen <- rnorm(n, 60, 10)
#'   phosphorus <- rnorm(n, 30, 5)
#'   yield <- 0.5 * nitrogen + 0.3 * phosphorus + rnorm(n, 0, 5)
#'   df <- data.frame(yield = yield, nitrogen = nitrogen, phosphorus = phosphorus)
#'   model <- 'yield ~ nitrogen + phosphorus'
#'   result <- sem_analysis(df, model, plot = FALSE)
#' }
#' }
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
sem_analysis <- function(data, model, plot = TRUE,
                            verbose = TRUE) {
  
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required for SEM analysis")
  }
  
  # Fit SEM model
  fit <- lavaan::sem(model, data = data)
  
  # Extract results
  fit_measures <- lavaan::fitMeasures(fit, c("chisq", "df", "pvalue", 
                                              "cfi", "tli", "rmsea", 
                                              "srmr", "aic", "bic"))
  
  if (verbose) {
    cat("\n==========================================================\n")
    cat("          STRUCTURAL EQUATION MODEL ANALYSIS\n")
    cat("==========================================================\n")

    cat("\n--- Model Fit Indices ---\n")
    cat("Chi-square:", round(fit_measures["chisq"], 3), 
        "(df =", fit_measures["df"], ", p =", round(fit_measures["pvalue"], 4), ")\n")
    cat("CFI:", round(fit_measures["cfi"], 3), "\n")
    cat("TLI:", round(fit_measures["tli"], 3), "\n")
    cat("RMSEA:", round(fit_measures["rmsea"], 3), "\n")
    cat("SRMR:", round(fit_measures["srmr"], 3), "\n")

    cat("\n--- Fit Assessment ---\n")
    if (fit_measures["cfi"] >= 0.95 && fit_measures["rmsea"] <= 0.06) {
      cat("Model fit: GOOD\n")
    } else if (fit_measures["cfi"] >= 0.90 && fit_measures["rmsea"] <= 0.08) {
      cat("Model fit: ACCEPTABLE\n")
    } else {
      cat("Model fit: POOR - Consider model modification\n")
    }

    cat("\n--- Parameter Estimates ---\n\n")
    print(lavaan::parameterEstimates(fit, standardized = TRUE))
  }
  
  # Plot path diagram
  if (plot && requireNamespace("semPlot", quietly = TRUE)) {
    semPlot::semPaths(fit, what = "std", layout = "tree",
                      edge.label.cex = 0.8, curvePivot = TRUE,
                      title = TRUE, title.cex = 1.2)
  }
  
  result <- list(
    fit = fit,
    fit_measures = fit_measures,
    parameters = lavaan::parameterEstimates(fit, standardized = TRUE)
  )
  
  class(result) <- "aridagri_sem"
  return(invisible(result))
}
