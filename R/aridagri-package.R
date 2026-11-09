#' aridagri: Statistical Tools for Agricultural Research
#'
#' @description
#' Analysis of variance for the common field designs with correct error
#' terms, standard errors, critical differences and letter groupings, plus
#' stability analysis, repeated measures, thermal indices, crop growth
#' analysis, nutrient use efficiency, response curves and economics.
#'
#' @section Experimental designs:
#' \itemize{
#'   \item \code{\link{anova_crd}}, \code{\link{anova_rbd}},
#'     \code{\link{anova_latin}}: single-factor designs (missing plots and
#'     unequal replication handled exactly)
#'   \item \code{\link{anova_factorial}}, \code{\link{anova_factorial_3way}}:
#'     factorial experiments
#'   \item \code{\link{anova_spd}}, \code{\link{anova_spd_ab_main}},
#'     \code{\link{anova_spd_c_main_ab_sub}}, \code{\link{anova_spd_ab_cd}},
#'     \code{\link{anova_sspd}}, \code{\link{anova_strip}}: split-plot
#'     family with SE and CD for every comparison
#'   \item \code{\link{anova_rbd_pooled}}, \code{\link{anova_spd_pooled}},
#'     \code{\link{anova_sspd_pooled}}: combined analysis over environments
#'   \item \code{\link{anova_augmented}}, \code{\link{anova_alpha_lattice}}:
#'     designs for germplasm and breeding trials
#'   \item \code{\link{anova_repeated}}: repeated measures with sphericity
#'     corrections and covariance-structure models
#'   \item \code{\link{anova_traits}}: publication table for many traits
#' }
#'
#' @section Mean comparisons and diagnostics:
#' \code{\link{perform_posthoc}}, \code{\link{compare_means}},
#' \code{\link{check_assumptions}}.
#'
#' @section Genotype x environment analysis:
#' \code{\link{stability_analysis}} (Eberhart-Russell, AMMI,
#' Finlay-Wilkinson, Shukla, Wricke, CV, Lin-Binns) and
#' \code{\link{gge_biplot}}.
#'
#' @section Agronomic indices:
#' \code{\link{thermal_indices}}, \code{\link{day_length}},
#' \code{\link{crop_growth_analysis}}, \code{\link{harvest_index}},
#' \code{\link{yield_gap_analysis}}, \code{\link{economic_indices}},
#' \code{\link{economic_analysis}}, \code{\link{nue_calculate}},
#' \code{\link{nutrient_response}}.
#'
#' @section Multivariate methods:
#' \code{\link{correlation_analysis}}, \code{\link{pca_analysis}},
#' \code{\link{path_analysis}}, \code{\link{sem_analysis}}.
#'
#' @section Figures and export:
#' \code{\link{plot_means}}, \code{\link{plot_interaction}},
#' \code{\link{plot_stability}}, \code{\link{plot_response}},
#' \code{plot()} methods, \code{\link{theme_aridagri}} (\pkg{ggplot2});
#' \code{\link{arid_plot}} (base graphics); \code{\link{export_results}}.
#'
#' @author Lalit Kumar Rolaniya (maintainer), Ram Lal Jat, Monika Punia
#'   (deceased), Raja Ram Choudhary and Hemant Poonia (deceased)
#' @docType package
#' @name aridagri-package
#' @aliases aridagri
#' @keywords package
"_PACKAGE"


#' Base-Graphics Plots for aridagri Results
#'
#' @description
#' Quick plots with base graphics (no extra packages): bar charts of
#' treatment or factor means with their letters for ANOVA results, a
#' heatmap for correlation results and the integrated ranking for stability
#' results. For publication figures with \pkg{ggplot2} see
#' \code{\link{plot_means}} and the \code{plot()} methods.
#'
#' @param x An object from an aridagri analysis function.
#' @param type Kept for compatibility; bar charts are drawn.
#' @param ... Unused.
#'
#' @return No return value; called for its side effect of drawing a plot.
#'
#' @examples
#' data(mothbean_rbd)
#' res <- anova_rbd(mothbean_rbd, "seed_yield", "treatment", "block",
#'                  verbose = FALSE)
#' arid_plot(res)
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram
#'   Choudhary, Hemant Poonia
arid_plot <- function(x, type = "bar", ...) {
  UseMethod("arid_plot")
}

#' @export
arid_plot.default <- function(x, type = "bar", ...) {
  ink <- .arid_ink
  if (inherits(x, "aridagri_correlation")) {
    cormat <- as.matrix(x$correlation)
    p <- ncol(cormat)
    vars <- colnames(cormat)
    if (is.null(vars)) vars <- paste0("V", seq_len(p))
    pal <- grDevices::colorRampPalette(c("#e34948", "#f0efec", "#2a78d6"))(101)
    brks <- seq(-1, 1, length.out = 102)
    op <- graphics::par(mar = c(6, 6, 3, 2))
    on.exit(graphics::par(op), add = TRUE)
    graphics::image(seq_len(p), seq_len(p), t(cormat[p:1, , drop = FALSE]),
                    col = pal, breaks = brks, axes = FALSE, xlab = "",
                    ylab = "",
                    main = paste0("Correlation matrix (", x$method %||%
                                    "pearson", ")"))
    graphics::axis(1, at = seq_len(p), labels = vars, las = 2, cex.axis = 0.8)
    graphics::axis(2, at = seq_len(p), labels = rev(vars), las = 2,
                   cex.axis = 0.8)
    graphics::box()
    for (i in seq_len(p)) {
      for (j in seq_len(p)) {
        val <- cormat[p - i + 1, j]
        graphics::text(j, i, formatC(val, format = "f", digits = 2),
                       cex = 0.7, col = if (abs(val) > 0.6) "white" else
                         ink[["primary"]])
      }
    }
    return(invisible(NULL))
  }
  if (inherits(x, "aridagri_stability")) {
    tab <- x$integrated
    if (!is.null(tab)) {
      tab <- tab[order(tab$Rank), , drop = FALSE]
      op <- graphics::par(mar = c(7, 5, 3, 2))
      on.exit(graphics::par(op), add = TRUE)
      graphics::barplot(tab$Rank, names.arg = as.character(tab$Genotype),
                        las = 2, col = .arid_cols[1], border = NA,
                        ylab = "Integrated rank (1 = most stable)",
                        main = "Integrated stability ranking")
      return(invisible(NULL))
    }
    gmean <- sort(rowMeans(x$ge_matrix), decreasing = TRUE)
    op <- graphics::par(mar = c(7, 5, 3, 2))
    on.exit(graphics::par(op), add = TRUE)
    graphics::barplot(gmean, las = 2, col = .arid_cols[1], border = NA,
                      ylab = "Mean", main = "Genotype means")
    return(invisible(NULL))
  }
  panels <- list()
  push <- function(title, labs, m, g = NULL, se = NULL) {
    panels[[length(panels) + 1L]] <<- list(title = title, labs = labs,
                                           mean = m, group = g, se = se)
  }
  if (!is.null(x$treatment_means)) {
    tm <- x$treatment_means
    push("Treatment means", tm$Treatment, tm$Mean, tm$Group, tm$SE)
  }
  if (!is.null(x$genotype_means) && is.data.frame(x$genotype_means)) {
    gm <- utils::head(x$genotype_means, 25)
    mcol <- if ("Adjusted_mean" %in% names(gm)) "Adjusted_mean" else "Mean"
    push("Genotype means", gm[[1]], gm[[mcol]], gm$Group)
  }
  if (is.list(x$factor_means) && !is.data.frame(x$factor_means)) {
    for (nm in names(x$factor_means)) {
      fm <- x$factor_means[[nm]]
      push(paste("Means of", nm), fm$Level, fm$Mean, fm$Group)
    }
  } else {
    for (fld in c("mean_a", "mean_b", "horizontal_means", "vertical_means")) {
      fm <- x[[fld]]
      if (is.data.frame(fm)) push(fld, fm$Level, fm$Mean, fm$Group)
    }
  }
  if (inherits(x, "aridagri_augmented")) {
    am <- x$adjusted_means
    am <- utils::head(am[order(-am$Adjusted_mean), ], 25)
    push("Adjusted means (top 25)", am$Entry, am$Adjusted_mean)
  }
  if (length(panels) == 0L) {
    message("arid_plot(): nothing to plot for class ",
            paste(class(x), collapse = "/"), ".")
    return(invisible(NULL))
  }
  np <- length(panels)
  nc <- if (np <= 1L) 1L else if (np <= 4L) 2L else 3L
  nr <- ceiling(np / nc)
  op <- graphics::par(mfrow = c(nr, nc), mar = c(7, 5, 3, 1))
  on.exit(graphics::par(op), add = TRUE)
  ttl <- x$response %||% x$design %||% "ANOVA"
  for (k in seq_len(np)) {
    pn <- panels[[k]]
    ord <- order(-pn$mean)
    m <- pn$mean[ord]
    se <- if (!is.null(pn$se)) pn$se[ord] else NULL
    top <- max(m + if (!is.null(se)) se else 0, na.rm = TRUE)
    bp <- graphics::barplot(m, names.arg = pn$labs[ord], las = 2,
                            col = .arid_cols[1], border = NA,
                            ylim = c(0, top * 1.15), ylab = "Mean",
                            main = paste0(ttl, ": ", pn$title))
    if (!is.null(se) && any(se > 0, na.rm = TRUE)) {
      graphics::arrows(bp, m - se, bp, m + se, angle = 90, code = 3,
                       length = 0.03, col = ink[["secondary"]])
    }
    if (!is.null(pn$group)) {
      graphics::text(bp, (m + if (!is.null(se)) se else 0) + 0.04 * top,
                     pn$group[ord], cex = 0.85, col = ink[["primary"]])
    }
  }
  invisible(NULL)
}


#' Export Results to Excel or CSV
#'
#' @description
#' Writes the main tables of an aridagri result to an Excel workbook (one
#' sheet per table) or a CSV file. Multi-trait tables from
#' \code{\link{anova_traits}} are written as they appear in the report.
#'
#' @param x An aridagri result object or a data frame.
#' @param file Output file path.
#' @param format \code{"xlsx"} (default; needs \pkg{writexl}) or
#'   \code{"csv"} (the first table only).
#' @param digits Number of decimal places for numeric columns.
#' @param verbose Logical. If TRUE (default), prints the file name.
#'
#' @return Invisibly returns the file path.
#'
#' @examples
#' data(mothbean_rbd)
#' res <- anova_rbd(mothbean_rbd, "seed_yield", "treatment", "block",
#'                  verbose = FALSE)
#' export_results(res, tempfile(fileext = ".csv"), format = "csv")
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram
#'   Choudhary, Hemant Poonia
export_results <- function(x, file, format = "xlsx", digits = 3,
                           verbose = TRUE) {
  format <- match.arg(tolower(format), c("xlsx", "csv"))
  rnd <- function(df) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    for (cl in names(df)) if (is.numeric(df[[cl]])) df[[cl]] <- round(
      df[[cl]], digits)
    df
  }
  with_source <- function(tab) {
    tab <- rnd(tab)
    cbind(Source = rownames(tab), tab, stringsAsFactors = FALSE)
  }
  sheets <- list()
  if (inherits(x, "aridagri_traits")) {
    sheets$Table <- x$table
    sheets$Statistics <- rnd(x$stats)
  } else if (inherits(x, "aridagri_correlation")) {
    sheets$Correlation <- cbind(Variable = rownames(x$correlation),
                                rnd(x$correlation))
    sheets$P_values <- cbind(Variable = rownames(x$p_values),
                             rnd(x$p_values))
  } else if (is.data.frame(x)) {
    sheets$Data <- rnd(x)
  } else if (is.list(x)) {
    tab <- x$anova_table %||% x$pooled_anova %||% x$anova_adjusted %||%
      x$anova
    if (!is.null(tab)) sheets$ANOVA <- with_source(tab)
    if (!is.null(x$anova_partitioned)) {
      sheets$ANOVA_partitioned <- with_source(x$anova_partitioned)
    }
    for (fld in c("treatment_means", "genotype_means", "adjusted_means",
                  "test_vs_checks", "environment_summary", "se_cd",
                  "comparisons", "integrated", "eberhart", "by_time")) {
      if (is.data.frame(x[[fld]])) sheets[[fld]] <- rnd(x[[fld]])
    }
    if (is.list(x$factor_means) && !is.data.frame(x$factor_means)) {
      for (nm in names(x$factor_means)) {
        sheets[[substr(paste0("means_", make.names(nm)), 1, 31)]] <-
          rnd(x$factor_means[[nm]])
      }
    }
  }
  if (length(sheets) == 0) stop("Object type not supported for export",
                                call. = FALSE)
  if (format == "xlsx") {
    if (!requireNamespace("writexl", quietly = TRUE)) {
      stop("Package 'writexl' required for Excel export. Install with: ",
           "install.packages('writexl')", call. = FALSE)
    }
    writexl::write_xlsx(sheets, file)
  } else {
    utils::write.csv(sheets[[1]], file, row.names = FALSE)
  }
  if (verbose) message("Results exported to: ", file)
  invisible(file)
}


#' Print Method for PCA Results
#'
#' @description
#' Prints the number of components retained and the variance they explain.
#'
#' @param x An object of class \code{aridagri_pca} from
#'   \code{\link{pca_analysis}}.
#' @param ... Additional arguments (unused).
#'
#' @return The input object, invisibly.
#'
#' @export
print.aridagri_pca <- function(x, ...) {
  cat("\n=== PCA Summary ===\n")
  k <- x$n_components_kaiser
  cat("Components retained (Kaiser):", k, "\n")
  if (k > 0) {
    cat("Variance explained:", round(x$cumulative_variance[k], 1), "%\n")
  }
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
