# =============================================================================
# GGE BIPLOT ANALYSIS
# Package: aridagri
# =============================================================================

#' GGE Biplot Analysis (Which-Won-Where and Mean vs Stability)
#'
#' @description
#' Genotype plus genotype-by-environment (GGE) biplot analysis of a
#' multi-environment trial (Yan et al., 2000). Gives the genotype and
#' environment scores on the first two principal components, the winning
#' genotype in each environment and the resulting mega-environments
#' (which-won-where), and the average-environment coordination (AEC)
#' ranking of genotypes for mean performance and stability.
#'
#' @param data Data frame with genotype, environment and trait columns
#'   (plot data or genotype x environment means).
#' @param genotype Name of the genotype column.
#' @param environment Name of the environment column.
#' @param trait Name of the trait column.
#' @param scaling \code{"none"} (default) or \code{"sd"} (each environment
#'   divided by its standard deviation, giving environments equal weight).
#' @param svp Singular value partitioning: \code{"symmetric"} (default,
#'   suited to which-won-where), \code{"genotype"} (genotype-focused,
#'   suited to comparing genotypes and the mean vs stability view) or
#'   \code{"environment"} (environment-focused).
#' @param verbose Logical. If TRUE (default), prints a summary.
#'
#' @details
#' The genotype x environment table of means is centred on the environment
#' means (and optionally scaled) and decomposed by singular value
#' decomposition. With singular values \eqn{d_k}, genotype scores are
#' \eqn{u_{ik} d_k^f} and environment scores \eqn{v_{jk} d_k^{1-f}}, with
#' \eqn{f = 0.5, 1, 0} for symmetric, genotype- and environment-focused
#' partitioning. The winner in an environment is the genotype with the
#' highest value predicted by the two-component biplot; environments with
#' the same winner form a mega-environment. On the AEC view the abscissa
#' (projection on the average-environment axis) ranks mean performance
#' and the ordinate (distance from that axis) measures instability.
#'
#' @return An object of class \code{aridagri_gge} with
#'   \code{genotype_scores}, \code{environment_scores},
#'   \code{variance_percent}, \code{which_won_where} (environment and
#'   winner), \code{mega_environments}, \code{aec} (genotype AEC abscissa,
#'   ordinate and ranks) and the settings used. Plot it with
#'   \code{plot()}.
#'
#' @references
#' Yan, W., Hunt, L.A., Sheng, Q. and Szlavnics, Z. (2000). Cultivar
#' evaluation and mega-environment investigation based on the GGE biplot.
#' \emph{Crop Science}, 40(3), 597-605. \doi{10.2135/cropsci2000.403597x}
#'
#' Yan, W. and Kang, M.S. (2003). \emph{GGE Biplot Analysis: A Graphical
#' Tool for Breeders, Geneticists, and Agronomists}. CRC Press, Boca Raton.
#'
#' @examples
#' data(pearlmillet_met)
#' gg <- gge_biplot(pearlmillet_met, "genotype", "environment", "grain_yield")
#' gg$which_won_where
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot(gg, type = "which_won_where")
#' }
#'
#' @export
gge_biplot <- function(data, genotype, environment, trait,
                       scaling = c("none", "sd"),
                       svp = c("symmetric", "genotype", "environment"),
                       verbose = TRUE) {
  scaling <- match.arg(scaling)
  svp <- match.arg(svp)
  d <- .prep_data(data, trait, c(genotype, environment), verbose)
  G <- d[[genotype]]
  E <- d[[environment]]
  if (nlevels(G) < 3 || nlevels(E) < 2) {
    stop("Need at least 3 genotypes and 2 environments.", call. = FALSE)
  }
  if (any(table(G, E) == 0)) {
    stop("Every genotype must be tested in every environment.",
         call. = FALSE)
  }
  M <- tapply(d[[trait]], list(G, E), mean)
  X <- sweep(M, 2, colMeans(M))
  if (scaling == "sd") {
    s <- apply(M, 2, stats::sd)
    if (any(s <= 0)) stop("An environment has no variation.", call. = FALSE)
    X <- sweep(X, 2, s, "/")
  }
  sv <- svd(X)
  f <- switch(svp, symmetric = 0.5, genotype = 1, environment = 0)
  k <- min(2, length(sv$d))
  gs <- sv$u[, 1:k, drop = FALSE] %*% diag(sv$d[1:k]^f, k)
  es <- sv$v[, 1:k, drop = FALSE] %*% diag(sv$d[1:k]^(1 - f), k)
  dimnames(gs) <- list(rownames(M), paste0("PC", 1:k))
  dimnames(es) <- list(colnames(M), paste0("PC", 1:k))
  pct <- 100 * sv$d^2 / sum(sv$d^2)
  pred <- gs %*% t(es)
  win <- rownames(M)[apply(pred, 2, which.max)]
  wwww <- data.frame(Environment = colnames(M), Winner = win,
                     stringsAsFactors = FALSE)
  mega <- split(wwww$Environment, wwww$Winner)
  a <- colMeans(es)
  an <- sqrt(sum(a^2))
  abscissa <- as.numeric(gs %*% a) / an
  ordinate <- if (k == 2) abs(gs[, 1] * a[2] - gs[, 2] * a[1]) / an else
    rep(0, nrow(gs))
  aec <- data.frame(Genotype = rownames(M), Mean = rowMeans(M),
                    AEC_abscissa = abscissa, AEC_ordinate = ordinate,
                    Mean_rank = rank(-abscissa), Stability_rank = rank(ordinate),
                    stringsAsFactors = FALSE)
  aec <- aec[order(aec$Mean_rank), ]
  rownames(aec) <- NULL
  res <- list(genotype_scores = gs, environment_scores = es,
              variance_percent = pct, which_won_where = wwww,
              mega_environments = mega, aec = aec, trait = trait,
              scaling = scaling, svp = svp, ge_matrix = M)
  class(res) <- "aridagri_gge"
  if (verbose) print(res)
  invisible(res)
}

#' Print Method for GGE Biplot Results
#'
#' @param x An object of class \code{aridagri_gge}.
#' @param ... Unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.aridagri_gge <- function(x, ...) {
  .print_info("GGE BIPLOT ANALYSIS",
              c(Trait = x$trait, Genotypes = nrow(x$genotype_scores),
                Environments = nrow(x$environment_scores),
                Scaling = x$scaling, `Singular value partitioning` = x$svp,
                `PC1 + PC2` = sprintf("%.1f%% of G + GE", sum(
                  x$variance_percent[1:min(2, length(x$variance_percent))]))))
  cat("\n Which won where\n")
  print(x$which_won_where, row.names = FALSE)
  cat("\n Mega-environments (environments sharing a winner)\n")
  for (w in names(x$mega_environments)) {
    cat("  ", w, ": ", paste(x$mega_environments[[w]], collapse = ", "),
        "\n", sep = "")
  }
  cat("\n Mean vs stability (average-environment coordination)\n")
  a <- x$aec
  a[, 2:4] <- round(a[, 2:4], 3)
  print(a, row.names = FALSE)
  invisible(x)
}

#' Plot Method for GGE Biplot Results
#'
#' @param x An object of class \code{aridagri_gge}.
#' @param type \code{"which_won_where"} (default), \code{"mean_stability"}
#'   or \code{"basic"}.
#' @param title Plot title.
#' @param ... Unused.
#'
#' @return A ggplot2 object.
#'
#' @export
plot.aridagri_gge <- function(x, type = c("which_won_where", "mean_stability",
                                          "basic"), title = NULL, ...) {
  .need_ggplot()
  type <- match.arg(type)
  gs <- x$genotype_scores
  es <- x$environment_scores
  if (ncol(gs) < 2) stop("The biplot needs two principal components.",
                         call. = FALSE)
  ink <- .arid_ink
  pct <- x$variance_percent
  dg <- data.frame(x = gs[, 1], y = gs[, 2], label = rownames(gs),
                   stringsAsFactors = FALSE)
  de <- data.frame(x = 0, y = 0, xend = es[, 1], yend = es[, 2],
                   label = rownames(es), stringsAsFactors = FALSE)
  R <- 1.35 * max(abs(c(gs, es)))
  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, colour = ink[["grid"]],
                        linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, colour = ink[["grid"]],
                        linewidth = 0.4)
  cap <- "Points: genotypes. Arrows: environments."
  if (type == "which_won_where") {
    hull <- grDevices::chull(gs[, 1], gs[, 2])
    hp <- gs[hull, , drop = FALSE]
    poly <- data.frame(x = hp[, 1], y = hp[, 2])
    nxt <- c(seq_len(nrow(hp))[-1], 1)
    segs <- do.call(rbind, lapply(seq_len(nrow(hp)), function(i) {
      v <- hp[nxt[i], ] - hp[i, ]
      nrm <- c(v[2], -v[1])
      if (sum(nrm * hp[i, ]) < 0) nrm <- -nrm
      nrm <- nrm / sqrt(sum(nrm^2))
      data.frame(x = 0, y = 0, xend = R * nrm[1], yend = R * nrm[2])
    }))
    p <- p + ggplot2::geom_polygon(data = poly, ggplot2::aes(x = x, y = y),
                                   fill = NA, colour = ink[["secondary"]],
                                   linewidth = 0.5) +
      ggplot2::geom_segment(data = segs, ggplot2::aes(x = x, y = y,
                                                      xend = xend,
                                                      yend = yend),
                            colour = ink[["muted"]], linetype = "dashed",
                            linewidth = 0.4)
    cap <- paste0(cap, "\nPolygon vertices are the best genotypes in the ",
                  "sectors between the dashed lines.")
  }
  if (type == "mean_stability") {
    a <- colMeans(es)
    a <- a / sqrt(sum(a^2))
    axis_df <- data.frame(x = -R * a[1], y = -R * a[2], xend = R * a[1],
                          yend = R * a[2])
    proj <- as.numeric(gs %*% a)
    foot <- data.frame(x = gs[, 1], y = gs[, 2], xend = proj * a[1],
                       yend = proj * a[2])
    p <- p + ggplot2::geom_segment(data = axis_df,
                                   ggplot2::aes(x = x, y = y, xend = xend,
                                                yend = yend),
                                   colour = ink[["secondary"]],
                                   linewidth = 0.6,
                                   arrow = ggplot2::arrow(
                                     length = ggplot2::unit(0.2, "cm"))) +
      ggplot2::geom_segment(data = foot, ggplot2::aes(x = x, y = y,
                                                      xend = xend,
                                                      yend = yend),
                            colour = ink[["muted"]], linetype = "dashed",
                            linewidth = 0.35)
    cap <- paste0("Arrow: average-environment axis (higher mean towards ",
                  "the head).\nDashed lines: distance from the axis ",
                  "(instability).")
  }
  p + ggplot2::geom_segment(data = de, ggplot2::aes(x = x, y = y, xend = xend,
                                                    yend = yend),
                            colour = .arid_cols[2], linewidth = 0.55,
                            arrow = ggplot2::arrow(length = ggplot2::unit(
                              0.15, "cm"))) +
    ggplot2::geom_text(data = .end_labels(de), ggplot2::aes(
      x = x, y = y, label = label, hjust = hj), colour = .arid_cols[2],
      size = 3) +
    ggplot2::geom_point(data = dg, ggplot2::aes(x = x, y = y),
                        colour = .arid_cols[1], size = 2.5) +
    ggplot2::geom_text(data = dg, ggplot2::aes(x = x, y = y, label = label),
                       colour = ink[["primary"]], size = 3, vjust = -0.8) +
    ggplot2::coord_equal(xlim = c(-R, R), ylim = c(-R, R)) +
    ggplot2::labs(title = title %||% paste0("GGE biplot (",
                                            gsub("_", "-", type), "): ",
                                            x$trait),
                  x = sprintf("PC1 (%.1f%%)", pct[1]),
                  y = sprintf("PC2 (%.1f%%)", pct[2]), caption = cap) +
    theme_aridagri()
}
