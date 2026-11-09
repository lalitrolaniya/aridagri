# =============================================================================
# PUBLICATION-READY FIGURES (ggplot2)
# Package: aridagri
# =============================================================================

# Column names used inside aes() (ggplot2 is in Suggests).
utils::globalVariables(c("level", "mean", "se", "group", "lab", "xf", "yv",
                         "series", "x", "y", "xend", "yend", "type", "label",
                         "fitted", "model", "PC1", "PC2", "lo", "hi", "hj"))

# Colour-blind-checked categorical palette (fixed order, never cycled) and
# text inks.
.arid_cols <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
                "#008300", "#6250d6", "#e34948")
.arid_ink <- c(primary = "#0b0b0b", secondary = "#52514e", muted = "#898781",
               grid = "#e1e0d9", axis = "#c3c2b7")

.need_ggplot <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this plot. Install it with ",
         "install.packages(\"ggplot2\"), or use arid_plot() for base ",
         "graphics.", call. = FALSE)
  }
}

# Up to 8 groups get the fixed palette; more groups share the first colour
# and are identified by labels.
.group_cols <- function(n) {
  if (n <= length(.arid_cols)) .arid_cols[seq_len(n)] else
    rep(.arid_cols[1], n)
}

#' A Clean ggplot2 Theme for Agricultural Figures
#'
#' @description
#' The theme used by the plotting functions of \pkg{aridagri}: white
#' background, light grid, dark text and a bold left-aligned title. Add it
#' to any ggplot2 figure.
#'
#' @param base_size Base font size in points. Default 11.
#'
#' @return A ggplot2 theme object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
#'     ggplot2::geom_point() + theme_aridagri()
#' }
#'
#' @export
theme_aridagri <- function(base_size = 11) {
  .need_ggplot()
  ink <- .arid_ink
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = ink[["primary"]]),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0,
                                         size = ggplot2::rel(1.12)),
      plot.title.position = "plot",
      plot.subtitle = ggplot2::element_text(colour = ink[["secondary"]],
                                            hjust = 0),
      plot.caption = ggplot2::element_text(colour = ink[["muted"]], hjust = 0,
                                           size = ggplot2::rel(0.8)),
      plot.caption.position = "plot",
      axis.text = ggplot2::element_text(colour = ink[["secondary"]]),
      axis.title = ggplot2::element_text(colour = ink[["secondary"]]),
      panel.grid.major = ggplot2::element_line(colour = ink[["grid"]],
                                               linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0,
                                         colour = ink[["primary"]]),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(colour = ink[["secondary"]]),
      legend.text = ggplot2::element_text(colour = ink[["primary"]]),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(10, 14, 10, 10))
}

# Bar chart of means with SE bars and letters. Long labels or many bars
# switch to horizontal bars so that labels stay readable.
.bar_means <- function(df, ylab, title, sort = TRUE, caption = NULL,
                       fill = .arid_cols[1], horizontal = NULL) {
  if (is.null(horizontal)) {
    horizontal <- max(nchar(as.character(df$level))) > 10 || nrow(df) > 12
  }
  lev <- if (sort) df$level[order(-df$mean)] else unique(df$level)
  if (horizontal) lev <- rev(lev)
  df$level <- factor(df$level, levels = lev)
  top <- max(df$mean + ifelse(is.na(df$se), 0, df$se), na.rm = TRUE)
  df$hi <- df$mean + ifelse(is.na(df$se), 0, df$se) + 0.02 * top
  p <- ggplot2::ggplot(df, ggplot2::aes(x = level, y = mean)) +
    ggplot2::geom_col(width = if (horizontal) 0.68 else 0.62, fill = fill) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(
      mult = c(0, if (horizontal) 0.12 else 0.1)))
  if (any(!is.na(df$se))) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean - se, ymax = mean + se),
      width = if (horizontal) 0.25 else 0.18, linewidth = 0.45,
      colour = .arid_ink[["secondary"]], na.rm = TRUE)
  }
  if (!is.null(df$group) && any(nzchar(df$group))) {
    p <- p + ggplot2::geom_text(ggplot2::aes(y = hi, label = group),
                                vjust = if (horizontal) 0.5 else 0,
                                hjust = if (horizontal) 0 else 0.5,
                                size = 3.4, colour = .arid_ink[["primary"]])
  }
  p <- p + ggplot2::labs(title = title, x = NULL, y = ylab,
                         caption = caption) + theme_aridagri()
  if (horizontal) {
    p + ggplot2::coord_flip() +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
  } else {
    p + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
  }
}

#' Bar Chart of Treatment Means with Letters
#'
#' @description
#' Publication-ready bar chart of treatment means with standard-error bars
#' and the letter groupings from an \pkg{aridagri} analysis. Works with the
#' results of \code{\link{anova_crd}}, \code{\link{anova_rbd}},
#' \code{\link{anova_latin}}, \code{\link{anova_rbd_pooled}},
#' \code{\link{anova_factorial}} (one factor or the combinations),
#' \code{\link{anova_spd}} and the other split-plot functions (one factor),
#' \code{\link{anova_alpha_lattice}} and \code{\link{anova_augmented}}
#' (adjusted means).
#'
#' @param x A result object from one of the functions above.
#' @param which For factorial and split-plot results, the factor to plot
#'   (name), or \code{"interaction"} for the treatment combinations of a
#'   two-factor factorial. Default: the first factor.
#' @param sort Logical. Order bars from highest to lowest mean. Default
#'   TRUE.
#' @param top For alpha lattice and augmented results, the number of
#'   entries to show. Default 20.
#' @param title Plot title. Default: the response name.
#' @param ylab Y-axis label. Default: the response name.
#' @param horizontal Logical. Draw horizontal bars. Default \code{NULL}
#'   chooses horizontal bars when labels are long or there are more than
#'   12 bars.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   data(mothbean_rbd)
#'   res <- anova_rbd(mothbean_rbd, "seed_yield", "treatment", "block",
#'                    posthoc = "duncan", verbose = FALSE)
#'   plot_means(res, ylab = "Seed yield (kg/ha)")
#' }
#'
#' @export
plot_means <- function(x, which = NULL, sort = TRUE, top = 20, title = NULL,
                       ylab = NULL, horizontal = NULL) {
  .need_ggplot()
  resp <- x$response %||% "Mean"
  if (is.null(ylab)) ylab <- resp
  cap <- if (!is.null(x$letters_method) && !is.na(x$letters_method)) {
    paste0("Error bars: SE. Bars sharing a letter do not differ (",
           x$letters_method, ", p = ", x$alpha %||% 0.05, ").")
  } else "Error bars: SE."
  if (inherits(x, c("aridagri_anova", "aridagri_pooled"))) {
    tm <- x$treatment_means
    se <- if (!is.null(tm$SE)) tm$SE else rep(x$se_treatment %||% NA,
                                              nrow(tm))
    df <- data.frame(level = tm$Treatment, mean = tm$Mean, se = se,
                     group = tm$Group %||% rep("", nrow(tm)),
                     stringsAsFactors = FALSE)
    return(.bar_means(df, ylab, title %||% resp, sort, cap,
                      horizontal = horizontal))
  }
  if (inherits(x, "aridagri_factorial")) {
    f <- x$factors
    which <- which %||% f[["A"]]
    if (identical(which, "interaction")) {
      im <- x$interaction_means
      df <- data.frame(level = paste(im[[f[["A"]]]], im[[f[["B"]]]],
                                     sep = " : "),
                       mean = im$Mean, se = im$SE,
                       group = im$Group %||% rep("", nrow(im)),
                       stringsAsFactors = FALSE)
    } else {
      tb <- if (which == f[["A"]]) x$mean_a else if (which == f[["B"]])
        x$mean_b else stop("'which' must be a factor name or \"interaction\".",
                           call. = FALSE)
      df <- data.frame(level = tb$Level, mean = tb$Mean, se = tb$SE,
                       group = tb$Group %||% rep("", nrow(tb)),
                       stringsAsFactors = FALSE)
    }
    return(.bar_means(df, ylab, title %||% resp, sort, cap,
                      horizontal = horizontal))
  }
  if (inherits(x, c("aridagri_split", "aridagri_split_pooled"))) {
    fm <- x$factor_means
    which <- which %||% names(fm)[1]
    if (!which %in% names(fm)) {
      stop("'which' must be one of: ", paste(names(fm), collapse = ", "),
           call. = FALSE)
    }
    tb <- fm[[which]]
    sem <- x$se_cd$SEm[x$se_cd$Effect == which][1]
    df <- data.frame(level = tb$Level, mean = tb$Mean, se = sem,
                     group = tb$Group %||% rep("", nrow(tb)),
                     stringsAsFactors = FALSE)
    return(.bar_means(df, ylab, title %||% paste(resp, "by", which), sort,
                      "Error bars: SEm from the ANOVA.",
                      horizontal = horizontal))
  }
  if (inherits(x, "aridagri_alpha")) {
    gm <- utils::head(x$genotype_means, top)
    df <- data.frame(level = gm$Genotype, mean = gm$Adjusted_mean,
                     se = gm$SE, group = gm$Group %||% rep("", nrow(gm)),
                     stringsAsFactors = FALSE)
    return(.bar_means(df, ylab, title %||% paste("Adjusted means:", resp),
                      sort, cap, horizontal = horizontal %||% TRUE))
  }
  if (inherits(x, "aridagri_augmented")) {
    am <- x$adjusted_means
    tv <- utils::head(am[am$Type == "Test", ][order(-am$Adjusted_mean[
      am$Type == "Test"]), ], top)
    ck <- am[am$Type == "Check", ]
    df <- rbind(data.frame(level = tv$Entry, mean = tv$Adjusted_mean,
                           type = "Test", stringsAsFactors = FALSE),
                data.frame(level = ck$Entry, mean = ck$Adjusted_mean,
                           type = "Check", stringsAsFactors = FALSE))
    df$level <- factor(df$level, levels = df$level[order(df$mean)])
    cols <- c(Test = .arid_cols[1], Check = .arid_cols[2])
    return(ggplot2::ggplot(df, ggplot2::aes(x = mean, y = level,
                                            colour = type)) +
             ggplot2::geom_vline(xintercept = mean(ck$Adjusted_mean),
                                 colour = .arid_ink[["muted"]],
                                 linetype = "dashed", linewidth = 0.4) +
             ggplot2::geom_point(size = 2.8) +
             ggplot2::scale_colour_manual(values = cols, name = NULL) +
             ggplot2::labs(title = title %||% paste("Adjusted means:", resp),
                           x = ylab, y = NULL,
                           caption = "Dashed line: mean of checks.") +
             theme_aridagri())
  }
  stop("No plot_means() method for an object of class ",
       paste(class(x), collapse = "/"), ".", call. = FALSE)
}


#' Interaction Plot of Two Factors
#'
#' @description
#' Line chart of the means of every combination of two factors, with the
#' first factor on the x-axis and one line per level of the second factor.
#' Works with \code{\link{anova_factorial}}, \code{\link{anova_spd}} and
#' the other split-plot functions (first two factors, or the two named in
#' \code{factors}), \code{\link{anova_strip}} and
#' \code{\link{anova_repeated}} (time course).
#'
#' @param x A result object.
#' @param factors Optional character vector of two factor names (split-plot
#'   results with more than two factors).
#' @param title Plot title.
#' @param ylab Y-axis label.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   data(clusterbean_spd)
#'   res <- anova_spd(clusterbean_spd, "seed_yield", "irrigation",
#'                    "nitrogen", "rep", verbose = FALSE)
#'   plot_interaction(res, ylab = "Seed yield (kg/ha)")
#' }
#'
#' @export
plot_interaction <- function(x, factors = NULL, title = NULL, ylab = NULL) {
  .need_ggplot()
  resp <- x$response %||% "Mean"
  ylab <- ylab %||% resp
  sem <- NA_real_
  if (inherits(x, "aridagri_factorial")) {
    im <- x$interaction_means
    fa <- x$factors[["A"]]
    fb <- x$factors[["B"]]
    df <- data.frame(xf = im[[fa]], series = im[[fb]], mean = im$Mean,
                     se = im$SE, stringsAsFactors = FALSE)
    df$xf <- factor(df$xf, levels = unique(im[[fa]]))
    xl <- fa
    sl <- fb
  } else if (inherits(x, c("aridagri_split", "aridagri_strip"))) {
    tabs <- if (inherits(x, "aridagri_strip")) {
      stats::setNames(list(x$interaction_means),
                      paste(x$factors[1], "x", x$factors[2]))
    } else x$interaction_tables
    nm <- if (is.null(factors)) names(tabs)[1] else {
      k <- c(paste(factors, collapse = " x "),
             paste(rev(factors), collapse = " x "))
      k[k %in% names(tabs)][1]
    }
    if (is.na(nm) || is.null(tabs[[nm]])) {
      stop("Interaction not found; give two factor names in 'factors'.",
           call. = FALSE)
    }
    m <- tabs[[nm]]
    parts <- strsplit(nm, " x ", fixed = TRUE)[[1]]
    df <- data.frame(xf = rep(rownames(m), ncol(m)),
                     series = rep(colnames(m), each = nrow(m)),
                     mean = as.vector(m), stringsAsFactors = FALSE)
    df$xf <- factor(df$xf, levels = rownames(m))
    df$se <- NA_real_
    xl <- parts[1]
    sl <- parts[2]
  } else if (inherits(x, "aridagri_repeated")) {
    m <- x$means
    se <- x$se_table
    df <- data.frame(xf = rep(colnames(m), each = nrow(m)),
                     series = rep(rownames(m), ncol(m)),
                     mean = as.vector(m), se = as.vector(se),
                     stringsAsFactors = FALSE)
    df$xf <- factor(df$xf, levels = colnames(m))
    xl <- x$factors[["time"]]
    sl <- x$factors[["treatment"]]
  } else {
    stop("No plot_interaction() method for an object of class ",
         paste(class(x), collapse = "/"), ".", call. = FALSE)
  }
  df$series <- factor(df$series, levels = unique(df$series))
  cols <- stats::setNames(.group_cols(nlevels(df$series)), levels(df$series))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = xf, y = mean, colour = series,
                                        group = series))
  if (any(!is.na(df$se))) {
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se,
                                                 ymax = mean + se),
                                    width = 0.12, linewidth = 0.4,
                                    na.rm = TRUE)
  }
  p + ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::scale_colour_manual(values = cols, name = sl) +
    ggplot2::labs(title = title %||% paste(resp, ":", xl, "x", sl), x = xl,
                  y = ylab,
                  caption = if (any(!is.na(df$se))) "Error bars: SE." else
                    NULL) +
    theme_aridagri()
}


#' Stability Plots for Multi-Environment Trials
#'
#' @description
#' Figures for \code{\link{stability_analysis}} results:
#' \itemize{
#'   \item \code{"eberhart"}: genotype mean against the regression
#'     coefficient \eqn{b_i}, with reference lines at the grand mean and
#'     \eqn{b_i = 1}.
#'   \item \code{"ammi1"}: AMMI1 biplot of mean against IPCA1 scores for
#'     genotypes and environments.
#'   \item \code{"ammi2"}: AMMI2 biplot of IPCA1 against IPCA2, with
#'     environments as vectors.
#' }
#'
#' @param x An object of class \code{aridagri_stability}.
#' @param type \code{"eberhart"} (default), \code{"ammi1"} or
#'   \code{"ammi2"}.
#' @param title Plot title.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   data(pearlmillet_met)
#'   st <- stability_analysis(pearlmillet_met, "genotype", "environment",
#'                            "rep", "grain_yield", verbose = FALSE)
#'   plot_stability(st, "eberhart")
#'   plot_stability(st, "ammi2")
#' }
#'
#' @export
plot_stability <- function(x, type = c("eberhart", "ammi1", "ammi2"),
                           title = NULL) {
  .need_ggplot()
  type <- match.arg(type)
  if (!inherits(x, "aridagri_stability")) {
    stop("'x' must come from stability_analysis().", call. = FALSE)
  }
  ink <- .arid_ink
  if (type == "eberhart") {
    if (is.null(x$eberhart)) stop("Run stability_analysis() with the ",
                                  "Eberhart-Russell method.", call. = FALSE)
    er <- x$eberhart
    df <- data.frame(x = er$Mean, y = er$bi, label = er$Genotype,
                     stringsAsFactors = FALSE)
    return(ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
             ggplot2::geom_hline(yintercept = 1, colour = ink[["muted"]],
                                 linetype = "dashed", linewidth = 0.4) +
             ggplot2::geom_vline(xintercept = x$grand_mean,
                                 colour = ink[["muted"]], linetype = "dashed",
                                 linewidth = 0.4) +
             ggplot2::geom_point(size = 2.8, colour = .arid_cols[1]) +
             ggplot2::geom_text(ggplot2::aes(label = label), vjust = -0.8,
                                size = 3.1, colour = ink[["primary"]]) +
             ggplot2::labs(title = title %||% paste("Eberhart-Russell:",
                                                    x$trait),
                           x = "Mean", y = expression(b[i]),
                           caption = paste0("Dashed lines: grand mean and ",
                                            "bi = 1.\nHigh mean with bi near 1 ",
                                            "indicates wide adaptation.")) +
             theme_aridagri())
  }
  if (is.null(x$ammi)) stop("Run stability_analysis() with the AMMI method.",
                            call. = FALSE)
  gs <- x$ammi$genotype_scores
  es <- x$ammi$environment_scores
  pct <- x$ammi$components$Variance_Percent
  if (type == "ammi1") {
    df <- rbind(data.frame(x = rowMeans(x$ge_matrix), y = gs[, 1],
                           label = rownames(gs), type = "Genotype",
                           stringsAsFactors = FALSE),
                data.frame(x = colMeans(x$ge_matrix), y = es[, 1],
                           label = rownames(es), type = "Environment",
                           stringsAsFactors = FALSE))
    cols <- c(Genotype = .arid_cols[1], Environment = .arid_cols[2])
    return(ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = type)) +
             ggplot2::geom_hline(yintercept = 0, colour = ink[["muted"]],
                                 linewidth = 0.4) +
             ggplot2::geom_vline(xintercept = x$grand_mean,
                                 colour = ink[["muted"]], linetype = "dashed",
                                 linewidth = 0.4) +
             ggplot2::geom_point(size = 2.6) +
             ggplot2::geom_text(ggplot2::aes(label = label), vjust = -0.8,
                                size = 3, show.legend = FALSE) +
             ggplot2::scale_colour_manual(values = cols, name = NULL) +
             ggplot2::labs(title = title %||% paste("AMMI1 biplot:", x$trait),
                           x = "Mean",
                           y = sprintf("IPCA1 (%.1f%%)", pct[1])) +
             theme_aridagri())
  }
  if (ncol(gs) < 2) stop("AMMI2 needs at least two IPCA axes.", call. = FALSE)
  dg <- data.frame(x = gs[, 1], y = gs[, 2], label = rownames(gs),
                   stringsAsFactors = FALSE)
  de <- data.frame(x = 0, y = 0, xend = es[, 1], yend = es[, 2],
                   label = rownames(es), stringsAsFactors = FALSE)
  ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, colour = ink[["grid"]],
                        linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, colour = ink[["grid"]],
                        linewidth = 0.4) +
    ggplot2::geom_segment(data = de, ggplot2::aes(x = x, y = y, xend = xend,
                                                  yend = yend),
                          colour = .arid_cols[2], linewidth = 0.6,
                          arrow = ggplot2::arrow(length = ggplot2::unit(
                            0.15, "cm"))) +
    ggplot2::geom_text(data = .end_labels(de), ggplot2::aes(
      x = x, y = y, label = label, hjust = hj), colour = .arid_cols[2],
      size = 3) +
    ggplot2::geom_point(data = dg, ggplot2::aes(x = x, y = y),
                        colour = .arid_cols[1], size = 2.6) +
    ggplot2::geom_text(data = dg, ggplot2::aes(x = x, y = y, label = label),
                       colour = ink[["primary"]], size = 3, vjust = -0.8) +
    ggplot2::coord_equal(xlim = .biplot_lim(c(dg$x, de$xend)),
                         ylim = .biplot_lim(c(dg$y, de$yend))) +
    ggplot2::labs(title = title %||% paste("AMMI2 biplot:", x$trait),
                  x = sprintf("IPCA1 (%.1f%%)", pct[1]),
                  y = sprintf("IPCA2 (%.1f%%)", pct[2]),
                  caption = paste0("Points: genotypes. Arrows: environments.",
                                   "\nGenotypes near the origin are stable.")) +
    theme_aridagri()
}



# Label positions just beyond arrow heads, left-aligned on the right side
# of the origin and right-aligned on the left side.
.end_labels <- function(de, push = 0.06) {
  r <- sqrt(de$xend^2 + de$yend^2)
  rmax <- max(r)
  k <- ifelse(r > 0, (r + push * rmax) / r, 1)
  data.frame(x = de$xend * k, y = de$yend * k, label = de$label,
             hj = ifelse(de$xend < 0, 1, 0), stringsAsFactors = FALSE)
}

# Axis limits with room for labels.
.biplot_lim <- function(v, pad = 0.22) {
  r <- range(c(v, 0))
  w <- diff(r)
  if (w == 0) w <- 1
  c(r[1] - pad * w, r[2] + pad * w)
}

#' Plot a Nutrient Response Curve
#'
#' @description
#' Observed yields and the fitted response curve(s) from
#' \code{\link{nutrient_response}}, with the economic optimum dose marked.
#'
#' @param x An object of class \code{aridagri_response} or
#'   \code{aridagri_response_set}.
#' @param dose,yield The data used to fit the model.
#' @param title Plot title.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   dose <- c(0, 20, 40, 60, 80, 100)
#'   yield <- c(1100, 1350, 1520, 1610, 1650, 1660)
#'   nr <- nutrient_response(dose, yield, model = "all", verbose = FALSE)
#'   plot_response(nr, dose, yield)
#' }
#'
#' @export
plot_response <- function(x, dose, yield, title = NULL) {
  .need_ggplot()
  fits <- if (inherits(x, "aridagri_response_set")) x$models else
    stats::setNames(list(x), x$model)
  grid_x <- seq(min(dose), max(dose), length.out = 200)
  curves <- do.call(rbind, lapply(names(fits), function(m) {
    data.frame(x = grid_x, fitted = fits[[m]]$predict(grid_x),
               model = gsub("_", " ", m), stringsAsFactors = FALSE)
  }))
  curves$model <- factor(curves$model, levels = unique(curves$model))
  opt <- do.call(rbind, lapply(names(fits), function(m) {
    data.frame(x = fits[[m]]$dose_economic, y = fits[[m]]$yield_economic,
               model = gsub("_", " ", m), stringsAsFactors = FALSE)
  }))
  opt <- opt[is.finite(opt$x), , drop = FALSE]
  opt$model <- factor(opt$model, levels = levels(curves$model))
  cols <- stats::setNames(.group_cols(nlevels(curves$model)),
                          levels(curves$model))
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = data.frame(x = dose, y = yield),
                        ggplot2::aes(x = x, y = y), size = 2.4,
                        colour = .arid_ink[["secondary"]]) +
    ggplot2::geom_line(data = curves, ggplot2::aes(x = x, y = fitted,
                                                   colour = model),
                       linewidth = 0.8)
  if (nrow(opt) > 0) {
    p <- p + ggplot2::geom_point(data = opt, ggplot2::aes(x = x, y = y,
                                                          colour = model),
                                 shape = 21, fill = "white", size = 3.2,
                                 stroke = 1.2)
  }
  p + ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::labs(title = title %||% "Nutrient response",
                  x = "Nutrient dose (kg/ha)", y = "Yield (kg/ha)",
                  caption = "Open circles: economic optimum dose.") +
    theme_aridagri()
}


#' Plot Methods for aridagri Results
#'
#' @description
#' \code{plot()} methods that draw the standard publication figure for each
#' kind of result with \pkg{ggplot2}: a bar chart of means with letters
#' (CRD, RCBD, Latin square, pooled RCBD, alpha lattice, augmented), an
#' interaction plot or bar chart (factorial and split-plot designs), a
#' time course (repeated measures) or a stability plot. For base graphics
#' without \pkg{ggplot2}, use \code{\link{arid_plot}}.
#'
#' @param x A result object.
#' @param type For factorial and split-plot results: \code{"interaction"}
#'   (default) or \code{"means"}. For stability results:
#'   \code{"eberhart"} (default), \code{"ammi1"} or \code{"ammi2"}.
#' @param ... Passed to \code{\link{plot_means}},
#'   \code{\link{plot_interaction}} or \code{\link{plot_stability}}.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   data(mothbean_rbd)
#'   res <- anova_rbd(mothbean_rbd, "seed_yield", "treatment", "block",
#'                    verbose = FALSE)
#'   plot(res)
#' }
#'
#' @name plot.aridagri
#' @export
plot.aridagri_anova <- function(x, ...) plot_means(x, ...)

#' @rdname plot.aridagri
#' @export
plot.aridagri_pooled <- function(x, ...) plot_means(x, ...)

#' @rdname plot.aridagri
#' @export
plot.aridagri_alpha <- function(x, ...) plot_means(x, ...)

#' @rdname plot.aridagri
#' @export
plot.aridagri_augmented <- function(x, ...) plot_means(x, ...)

#' @rdname plot.aridagri
#' @export
plot.aridagri_factorial <- function(x, type = c("interaction", "means"),
                                    ...) {
  type <- match.arg(type)
  if (type == "interaction") plot_interaction(x, ...) else plot_means(x, ...)
}

#' @rdname plot.aridagri
#' @export
plot.aridagri_split <- function(x, type = c("interaction", "means"), ...) {
  type <- match.arg(type)
  if (type == "interaction") plot_interaction(x, ...) else plot_means(x, ...)
}

#' @rdname plot.aridagri
#' @export
plot.aridagri_strip <- function(x, ...) plot_interaction(x, ...)

#' @rdname plot.aridagri
#' @export
plot.aridagri_repeated <- function(x, ...) plot_interaction(x, ...)

#' @rdname plot.aridagri
#' @export
plot.aridagri_stability <- function(x, type = c("eberhart", "ammi1", "ammi2"),
                                    ...) {
  plot_stability(x, type = match.arg(type), ...)
}
