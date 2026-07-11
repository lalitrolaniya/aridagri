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
#' @return No return value; called for its side effect of drawing a plot on the
#'   current graphics device. The plot produced depends on the class of \code{x}:
#'   a correlation heatmap for \code{aridagri_correlation} objects, a bar chart of
#'   treatment means for \code{aridagri_anova} objects, and a bar chart of the
#'   integrated stability ranking for \code{aridagri_stability} objects. Uses only
#'   base graphics, so no additional packages are required.
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

  # --- Correlation objects: draw a heatmap of the correlation matrix ---------
  if (inherits(x, "aridagri_correlation")) {
    cormat <- as.matrix(x$correlation)
    p <- ncol(cormat)
    vars <- colnames(cormat)
    if (is.null(vars)) vars <- paste0("V", seq_len(p))

    # Blue (negative) - white (0) - red (positive) palette
    pal <- grDevices::colorRampPalette(c("#2166AC", "#FFFFFF", "#B2182B"))(101)
    brks <- seq(-1, 1, length.out = 102)

    op <- graphics::par(mar = c(6, 6, 3, 2))
    on.exit(graphics::par(op), add = TRUE)

    # image() plots columns bottom-to-top; reverse rows so it reads top-to-bottom
    graphics::image(seq_len(p), seq_len(p), t(cormat[p:1, , drop = FALSE]),
                    col = pal, breaks = brks, axes = FALSE, xlab = "", ylab = "",
                    main = paste0("Correlation matrix (",
                                  if (!is.null(x$method)) x$method else "pearson", ")"))
    graphics::axis(1, at = seq_len(p), labels = vars, las = 2, cex.axis = 0.8)
    graphics::axis(2, at = seq_len(p), labels = rev(vars), las = 2, cex.axis = 0.8)
    graphics::box()

    # Overlay the coefficients
    for (i in seq_len(p)) {
      for (j in seq_len(p)) {
        val <- cormat[p - i + 1, j]
        graphics::text(j, i, formatC(val, format = "f", digits = 2),
                       cex = 0.7,
                       col = if (abs(val) > 0.6) "white" else "black")
      }
    }
    return(invisible(NULL))
  }

  # --- ANOVA-family objects: bar chart(s) of factor / treatment means --------
  anova_classes <- c("aridagri_anova", "aridagri_factorial", "aridagri_factorial3",
                     "aridagri_latin", "aridagri_spd", "aridagri_spd_ab",
                     "aridagri_spd_cab", "aridagri_spd_abcd", "aridagri_strip",
                     "aridagri_alpha", "aridagri_augmented", "aridagri_pooled",
                     "aridagri_spd_pooled", "aridagri_sspd_pooled")
  means_fields <- c("treatment_means", "genotype_means", "main_means", "sub_means",
                    "env_means", "horizontal_means", "vertical_means",
                    "mean_a", "mean_b", "mean_c", "factor_means")
  if (inherits(x, anova_classes) || any(means_fields %in% names(x))) {

    # normalise a means data.frame to list(labs, mean, sd)
    norm_means <- function(df) {
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(NULL)
      num <- vapply(df, is.numeric, logical(1))
      mcol <- if ("Mean" %in% names(df)) "Mean" else if (any(num)) names(df)[num][1] else NA
      if (is.na(mcol)) return(NULL)
      lcol <- setdiff(names(df), c(mcol, "SD"))
      labs <- if (length(lcol)) as.character(df[[lcol[1]]]) else as.character(seq_len(nrow(df)))
      list(labs = labs, mean = as.numeric(df[[mcol]]),
           sd = if ("SD" %in% names(df)) as.numeric(df$SD) else NULL)
    }

    panels <- list()
    push <- function(title, df) {
      m <- norm_means(df)
      if (!is.null(m)) panels[[length(panels) + 1L]] <<- c(list(title = title), m)
    }
    if (!is.null(x$treatment_means))  push("Treatment means", x$treatment_means)
    if (!is.null(x$genotype_means))   push("Genotype means",  x$genotype_means)
    if (!is.null(x$main_means))       push("Main-plot means", x$main_means)
    if (!is.null(x$sub_means))        push("Sub-plot means",  x$sub_means)
    if (!is.null(x$env_means))        push("Environment means", x$env_means)
    if (!is.null(x$horizontal_means)) push("Horizontal factor", x$horizontal_means)
    if (!is.null(x$vertical_means))   push("Vertical factor",   x$vertical_means)
    if (!is.null(x$mean_a))           push("Factor A", x$mean_a)
    if (!is.null(x$mean_b))           push("Factor B", x$mean_b)
    if (!is.null(x$mean_c))           push("Factor C", x$mean_c)
    if (is.list(x$factor_means) && !is.data.frame(x$factor_means)) {
      fm <- x$factor_means
      nms <- if (!is.null(names(fm))) names(fm) else paste0("Factor ", seq_along(fm))
      for (k in seq_along(fm)) push(paste0("Factor ", nms[k]), fm[[k]])
    }

    # fall back: derive means from the fitted model's first treatment factor
    if (length(panels) == 0L && !is.null(x$model)) {
      mf <- tryCatch(stats::model.frame(x$model), error = function(e) NULL)
      if (!is.null(mf) && ncol(mf) >= 2L) {
        is_fac <- vapply(mf, function(z) is.factor(z) || is.character(z), logical(1))
        fac <- names(mf)[is_fac]
        fac <- fac[!tolower(fac) %in% c("rep", "reps", "replication", "block",
                                        "blocks", "rep_block")]
        if (!length(fac)) fac <- names(mf)[is_fac]
        if (length(fac)) {
          ag <- aggregate(mf[[1]], by = list(mf[[fac[1]]]), FUN = mean, na.rm = TRUE)
          names(ag) <- c(fac[1], "Mean")
          push(paste0(fac[1], " means"), ag)
        }
      }
    }

    if (length(panels) == 0L) {
      message("arid_plot(): no factor means available to plot for class ",
              paste(class(x), collapse = "/"), ".")
      return(invisible(NULL))
    }

    np <- length(panels)
    nc <- if (np <= 1L) 1L else if (np <= 4L) 2L else 3L
    nr <- ceiling(np / nc)
    op <- graphics::par(mfrow = c(nr, nc), mar = c(7, 5, 3, 2))
    on.exit(graphics::par(op), add = TRUE)

    dlab <- if (!is.null(x$design)) x$design else "ANOVA"
    cols <- c("#4292C6", "#41AB5D", "#EF6548", "#8C6BB1", "#FDAE6B", "#66C2A4")
    for (k in seq_len(np)) {
      pnl  <- panels[[k]]
      ord  <- order(-pnl$mean)
      mval <- pnl$mean[ord]; labs <- pnl$labs[ord]
      sdv  <- if (!is.null(pnl$sd)) pnl$sd[ord] else NULL
      ytop <- max(mval, if (!is.null(sdv)) mval + sdv else mval, na.rm = TRUE) * 1.15
      bp <- graphics::barplot(mval, names.arg = labs, las = 2,
                              col = cols[(k - 1L) %% length(cols) + 1L], border = "grey20",
                              ylim = c(0, ytop), ylab = "Mean",
                              main = paste0(dlab, ": ", pnl$title))
      if (!is.null(sdv)) {
        ok <- !is.na(sdv) & sdv > 0
        if (any(ok))
          graphics::arrows(bp[ok], mval[ok] - sdv[ok], bp[ok], mval[ok] + sdv[ok],
                           angle = 90, code = 3, length = 0.04, col = "grey30")
      }
      graphics::box()
    }
    return(invisible(NULL))
  }

  # --- Stability objects: bar chart of the integrated stability ranking ------
  if (inherits(x, "aridagri_stability")) {
    tab <- x$integrated
    if (!is.null(tab) && all(c("Genotype", "Rank") %in% names(tab))) {
      tab <- tab[order(tab$Rank), , drop = FALSE]
      op <- graphics::par(mar = c(7, 5, 3, 2))
      on.exit(graphics::par(op), add = TRUE)
      graphics::barplot(tab$Rank, names.arg = as.character(tab$Genotype),
                        las = 2, col = "#41AB5D", border = "#006D2C",
                        ylab = "Integrated stability rank (1 = most stable)",
                        main = "Integrated stability ranking")
      graphics::box()
      return(invisible(NULL))
    }
    # Fall back to genotype means from the GE matrix
    gem <- x$ge_matrix
    if (!is.null(gem)) {
      num <- gem[, vapply(gem, is.numeric, logical(1)), drop = FALSE]
      gmean <- sort(rowMeans(num, na.rm = TRUE), decreasing = TRUE)
      op <- graphics::par(mar = c(7, 5, 3, 2))
      on.exit(graphics::par(op), add = TRUE)
      graphics::barplot(gmean, las = 2, col = "#41AB5D", border = "#006D2C",
                        ylab = "Mean", main = "Genotype means")
      graphics::box()
      return(invisible(NULL))
    }
    message("arid_plot(): no stability summary available to plot.")
    return(invisible(NULL))
  }

  message("arid_plot(): no plotting method available for an object of class ",
          paste(class(x), collapse = "/"), ".")
  invisible(NULL)
}


#' Export Results to Publication Format
#'
#' @description
#' Exports aridagri analysis results to Excel format for publication.
#'
#' @param x An object from aridagri analysis functions
#' @param file Output file path
#' @param format Output format: \code{"xlsx"} (default) or \code{"csv"}
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
  
  format <- match.arg(tolower(format), c("xlsx", "csv"))

  if (format == "xlsx") {
    if (!requireNamespace("writexl", quietly = TRUE)) {
      stop("Package 'writexl' required for Excel export. Install with: install.packages('writexl')")
    }
    writexl::write_xlsx(as.data.frame(export_data), file)
  } else if (format == "csv") {
    utils::write.csv(as.data.frame(export_data), file, row.names = TRUE)
  }

  if (verbose) message("Results exported to: ", file)

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
