#' aridagri: Comprehensive Statistical Tools for Agricultural Research
#'
#' @description
#' The aridagri package provides comprehensive statistical and analytical tools 
#' for agricultural research. It includes complete ANOVA functions for all 
#' experimental designs, multiple post-hoc tests, stability analysis methods,
#' thermal indices, crop growth analysis, and advanced statistical methods.
#'
#' @section Experimental Design ANOVA Functions:
#' \itemize{
#'   \item \code{\link{anova_crd}}: Completely Randomized Design
#'   \item \code{\link{anova_rbd}}: Randomized Block Design
#'   \item \code{\link{anova_rbd_pooled}}: Pooled RBD (Multi-Environment)
#'   \item \code{\link{anova_latin}}: Latin Square Design
#'   \item \code{\link{anova_factorial}}: Two-Factor Factorial
#'   \item \code{\link{anova_factorial_3way}}: Three-Factor Factorial
#'   \item \code{\link{anova_spd}}: Split Plot Design
#'   \item \code{\link{anova_sspd}}: Split-Split Plot Design
#'   \item \code{\link{anova_strip}}: Strip Plot Design
#'   \item \code{\link{anova_augmented}}: Augmented Block Design
#'   \item \code{\link{anova_alpha_lattice}}: Alpha Lattice Design
#' }
#'
#' @section Post-Hoc Tests:
#' \itemize{
#'   \item \code{\link{perform_posthoc}}: Multiple comparison tests (LSD, Duncan, Tukey, SNK, Scheffe, Bonferroni, Dunnett)
#'   \item \code{\link{check_assumptions}}: ANOVA assumption checking
#' }
#'
#' @section Agronomic Analysis Functions:
#' \itemize{
#'   \item \code{\link{stability_analysis}}: Multi-method stability analysis (Eberhart-Russell, AMMI, Finlay-Wilkinson, Shukla, Wricke, CV, Superiority)
#'   \item \code{\link{thermal_indices}}: GDD, HTU, PTU, Heat Use Efficiency
#'   \item \code{\link{crop_growth_analysis}}: CGR, RGR, NAR, LAI
#'   \item \code{\link{harvest_index}}: Harvest index and partitioning
#'   \item \code{\link{yield_gap_analysis}}: Yield gap calculations
#'   \item \code{\link{economic_indices}}: B:C ratio, net returns
#' }
#'
#' @section Statistical Functions:
#' \itemize{
#'   \item \code{\link{correlation_analysis}}: Correlation matrix with significance
#'   \item \code{\link{pca_analysis}}: Principal component analysis
#'   \item \code{\link{path_analysis}}: Path coefficient analysis
#'   \item \code{\link{sem_analysis}}: Structural equation modeling
#' }
#'
#' @section Nutrient Analysis Functions:
#' \itemize{
#'   \item \code{\link{nue_calculate}}: Nutrient use efficiency calculations
#'   \item \code{\link{nutrient_response}}: Response curve analysis
#'   \item \code{\link{economic_analysis}}: Economic viability assessment
#' }
#'
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram Choudhary
#' @docType package
#' @name aridagri-package
#' @aliases aridagri
#'
#' @keywords package
"_PACKAGE"


#' Visualization Functions for aridagri
#'
#' @description
#' Generate publication-quality plots for aridagri analyses.
#'
#' @param x An object from aridagri analysis functions
#' @param type Plot type: "bar", "line", "interaction", "boxplot"
#' @param ... Additional arguments passed to plotting functions
#'
#' @return A ggplot2 object
#'
#' @examples
#' df <- data.frame(
#'   yield = c(1200, 1350, 1100, 1450, 1280),
#'   wue = c(4.2, 4.8, 3.9, 5.1, 4.5),
#'   protein = c(22.1, 23.5, 21.8, 24.2, 22.9)
#' )
#' result <- correlation_analysis(df, plot = FALSE)
#' arid_plot(result)
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram Choudhary
arid_plot <- function(x, type = "bar", ...) {
  UseMethod("arid_plot")
}


#' @export
arid_plot.default <- function(x, type = "bar", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for plotting")
  }
  message("Plotting method not available for this object type.")
  invisible(NULL)
}


#' Export Results to Publication Format
#'
#' @description
#' Exports aridagri analysis results to Excel format for publication.
#'
#' @param x An object from aridagri analysis functions
#' @param file Output file path
#' @param format Output format: "xlsx"
#' @param digits Number of decimal places
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return Invisibly returns the file path
#'
#' @examples
#' \donttest{
#' df <- data.frame(
#'   yield = c(1200, 1350, 1100, 1450, 1280),
#'   wue = c(4.2, 4.8, 3.9, 5.1, 4.5)
#' )
#' result <- correlation_analysis(df, plot = FALSE)
#' export_results(result, tempfile(fileext = ".xlsx"))
#' }
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram Choudhary
export_results <- function(x, file, format = "xlsx", digits = 3,
                            verbose = TRUE) {
  
  # Extract data based on object class
  if (inherits(x, "aridagri_correlation")) {
    export_data <- round(x$correlation, digits)
  } else if (inherits(x, "data.frame")) {
    export_data <- x
  } else if (is.list(x) && !is.null(x$anova_table)) {
    export_data <- x$anova_table
  } else {
    stop("Object type not supported for export")
  }
  
  if (format == "xlsx") {
    if (!requireNamespace("writexl", quietly = TRUE)) {
      stop("Package 'writexl' required for Excel export. Install with: install.packages('writexl')")
    }
    writexl::write_xlsx(as.data.frame(export_data), file)
    message("Results exported to: ", file)
  }
  
  invisible(file)
}


#' Print Method for PCA Results
#'
#' @description
#' Prints a formatted summary of Principal Component Analysis (PCA) results.
#'
#' @param x An object of class 'aridagri_pca' from \code{\link{pca_analysis}}
#' @param ... Additional arguments (currently unused)
#'
#' @return No return value, called for side effects. Prints the number of
#'   components retained by Kaiser criterion and cumulative variance explained
#'   to the console. The input object is returned invisibly.
#'
#' @export
print.aridagri_pca <- function(x, ...) {
  cat("\n=== PCA Summary ===\n")
  cat("Components retained (Kaiser):", x$n_components_kaiser, "\n")
  cat("Variance explained:", round(x$cumulative_variance[x$n_components_kaiser], 1), "%\n")
  invisible(x)
}

#' @export
print.aridagri_path <- function(x, ...) {
  cat("\n=== Path Analysis ===\n")
  cat("R-squared:", round(x$r_squared, 4), "\n")
  cat("Residual:", round(x$residual, 4), "\n\n")
  cat("Direct Effects:\n")
  print(round(x$direct_effects, 4))
  invisible(x)
}

#' @export
print.aridagri_stability <- function(x, ...) {
  cat("\n=== Stability Analysis ===\n")
  cat("Grand Mean:", round(x$grand_mean, 2), "\n")
  if (!is.null(x$integrated)) {
    cat("\nIntegrated Ranking:\n")
    print(head(x$integrated, 10))
  }
  invisible(x)
}
