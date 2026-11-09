# =============================================================================
# PUBLICATION TABLE FOR MANY TRAITS
# Package: aridagri
# =============================================================================

#' Publication Table of Treatment Means for Many Traits
#'
#' @description
#' Runs the analysis of variance for several traits at once and assembles
#' the table that agronomy journals and theses expect: treatment means
#' followed by letters, then SEm, CD and CV for each trait (with "NS" when
#' the F-test is not significant). Works for completely randomized,
#' randomized block, two-factor factorial and split-plot designs.
#'
#' @param data Data frame with one row per plot.
#' @param traits Character vector of response columns.
#' @param design One of \code{"rbd"} (default), \code{"crd"},
#'   \code{"factorial"} or \code{"spd"}.
#' @param treatment Treatment column (CRD and RBD).
#' @param block Block (replication) column (RBD, factorial in blocks and
#'   split-plot).
#' @param factor1,factor2 The two factors of a factorial design.
#' @param main_plot,sub_plot Main-plot and sub-plot factors of a
#'   split-plot design.
#' @param posthoc Method for the letters: \code{"lsd"} (default),
#'   \code{"duncan"}, \code{"tukey"}, \code{"snk"}, \code{"scheffe"} or
#'   \code{"bonferroni"}.
#' @param alpha Significance level. Default 0.05.
#' @param digits Decimal places for means, SEm and CD. A single number, or
#'   one per trait. Default 2.
#' @param letters_if_ns Logical. If FALSE (default), letters are omitted
#'   when the F-test of that factor is not significant.
#' @param verbose Logical. If TRUE (default), prints the table.
#'
#' @return An object of class \code{aridagri_traits} with \code{table}
#'   (character data frame ready to export), \code{stats} (F, p, SEm, SEd,
#'   CD and CV for each trait and factor) and \code{results} (the full
#'   analysis for each trait). Export with \code{\link{export_results}}.
#'
#' @examples
#' data(mothbean_rbd)
#' tab <- anova_traits(mothbean_rbd,
#'                     traits = c("plant_height", "pods_per_plant",
#'                                "seed_yield", "harvest_index"),
#'                     design = "rbd", treatment = "treatment",
#'                     block = "block", posthoc = "duncan")
#' tab$table
#'
#' data(clusterbean_spd)
#' anova_traits(clusterbean_spd, traits = c("seed_yield", "wue"),
#'              design = "spd", main_plot = "irrigation",
#'              sub_plot = "nitrogen", block = "rep")
#'
#' @export
anova_traits <- function(data, traits, design = c("rbd", "crd", "factorial",
                                                  "spd"),
                         treatment = NULL, block = NULL, factor1 = NULL,
                         factor2 = NULL, main_plot = NULL, sub_plot = NULL,
                         posthoc = "lsd", alpha = 0.05, digits = 2,
                         letters_if_ns = FALSE, verbose = TRUE) {
  design <- match.arg(design)
  .check_alpha(alpha)
  if (!is.character(traits) || length(traits) < 1) {
    stop("'traits' must name at least one response column.", call. = FALSE)
  }
  .check_cols(data, traits, "Trait")
  method <- setdiff(.resolve_methods(posthoc), "dunnett")
  if (length(method) != 1) {
    stop("Choose one post-hoc method for the letters.", call. = FALSE)
  }
  need <- switch(design, crd = list(treatment = treatment),
                 rbd = list(treatment = treatment, block = block),
                 factorial = list(factor1 = factor1, factor2 = factor2),
                 spd = list(main_plot = main_plot, sub_plot = sub_plot,
                            block = block))
  miss <- names(need)[vapply(need, is.null, logical(1))]
  if (length(miss) > 0) {
    stop("For design = \"", design, "\" supply: ",
         paste(miss, collapse = ", "), ".", call. = FALSE)
  }
  digs <- rep_len(digits, length(traits))
  names(digs) <- traits
  pm <- .pm()
  cd_lab <- sprintf("CD (p = %s)", alpha)
  results <- list()
  stats_rows <- list()
  sections <- list()

  fmt <- function(v, dg) formatC(v, format = "f", digits = dg)
  for (tr in traits) {
    dg <- digs[[tr]]
    if (design %in% c("crd", "rbd")) {
      res <- if (design == "crd") {
        anova_crd(data, tr, treatment, posthoc = method, alpha = alpha,
                  verbose = FALSE)
      } else {
        anova_rbd(data, tr, treatment, block, posthoc = method,
                  alpha = alpha, verbose = FALSE)
      }
      i <- if (design == "crd") 1 else 2
      p <- res$anova_table[["Pr(>F)"]][i]
      tm <- res$treatment_means
      lv <- levels(stats::model.frame(res$model)[[treatment]])
      tm <- tm[match(lv, tm$Treatment), , drop = FALSE]
      sig <- !is.na(p) && p < alpha
      let <- if (sig || letters_if_ns) tm$Group else rep("", nrow(tm))
      sec <- list(list(title = NULL, factor = treatment,
                       levels = tm$Treatment,
                       cells = trimws(paste(fmt(tm$Mean, dg), let)),
                       sem = res$se_mean, cd = res$cd, sig = sig))
      cv <- c(`CV (%)` = res$cv)
      stats_rows[[tr]] <- data.frame(Trait = tr, Factor = treatment,
                                     F_value = res$anova_table[["F value"]][i],
                                     p_value = p, SEm = res$se_mean,
                                     SEd = res$se_diff, CD = res$cd,
                                     CV = res$cv, stringsAsFactors = FALSE)
    } else if (design == "factorial") {
      res <- anova_factorial(data, tr, factor1, factor2, replication = block,
                             posthoc = method, alpha = alpha, verbose = FALSE)
      at <- res$anova_table
      rn <- trimws(gsub("`", "", rownames(at)))
      pa <- at[["Pr(>F)"]][rn == factor1]
      pb <- at[["Pr(>F)"]][rn == factor2]
      pab <- at[["Pr(>F)"]][rn %in% c(paste0(factor1, ":", factor2),
                                       paste0(factor2, ":", factor1))]
      mk <- function(tb, p, cmp_row, fac, title) {
        sig <- !is.na(p) && p < alpha
        let <- if ((sig || letters_if_ns) && !is.null(tb$Group)) tb$Group else
          rep("", nrow(tb))
        list(title = title, factor = fac, levels = tb$Level,
             cells = trimws(paste(fmt(tb$Mean, dg), let)),
             sem = res$comparisons$SEm[cmp_row],
             cd = res$comparisons$CD[cmp_row], sig = sig)
      }
      sig_ab <- !is.na(pab) && pab < alpha
      sec <- list(mk(res$mean_a, pa, 1, factor1, paste("Factor A:", factor1)),
                  mk(res$mean_b, pb, 2, factor2, paste("Factor B:", factor2)),
                  list(title = paste(factor1, "x", factor2), factor = "AxB",
                       levels = character(0), cells = character(0),
                       sem = res$comparisons$SEm[3],
                       cd = res$comparisons$CD[3], sig = sig_ab))
      cv <- c(`CV (%)` = res$cv)
      stats_rows[[tr]] <- data.frame(
        Trait = tr, Factor = c(factor1, factor2, paste(factor1, "x", factor2)),
        F_value = at[["F value"]][match(c(factor1, factor2),
                                        rn)][c(1, 2, NA)],
        p_value = c(pa, pb, pab), SEm = res$comparisons$SEm,
        SEd = res$comparisons$SEd, CD = res$comparisons$CD, CV = res$cv,
        stringsAsFactors = FALSE)
      stats_rows[[tr]]$F_value[3] <- at[["F value"]][rn %in% c(
        paste0(factor1, ":", factor2), paste0(factor2, ":", factor1))]
    } else {
      res <- anova_spd(data, tr, main_plot, sub_plot, block, posthoc = method,
                       alpha = alpha, verbose = FALSE)
      at <- res$anova_table
      rn <- trimws(gsub("`", "", rownames(at)))
      pa <- at[["Pr(>F)"]][rn == main_plot]
      pb <- at[["Pr(>F)"]][rn == sub_plot]
      pab <- at[["Pr(>F)"]][rn %in% c(paste0(main_plot, ":", sub_plot),
                                       paste0(sub_plot, ":", main_plot))]
      mk2 <- function(tb, p, sem, cd, fac, title) {
        sig <- !is.na(p) && p < alpha
        let <- if ((sig || letters_if_ns) && !is.null(tb$Group)) tb$Group else
          rep("", nrow(tb))
        list(title = title, factor = fac, levels = tb$Level,
             cells = trimws(paste(fmt(tb$Mean, dg), let)), sem = sem,
             cd = cd, sig = sig)
      }
      sig_ab <- !is.na(pab) && pab < alpha
      sec <- list(
        mk2(res$main_means, pa, res$se_main, res$cd_main, main_plot,
            paste("Main plot:", main_plot)),
        mk2(res$sub_means, pb, res$se_sub, res$cd_sub, sub_plot,
            paste("Sub-plot:", sub_plot)),
        list(title = paste(main_plot, "x", sub_plot), factor = "AxB",
             levels = character(0), cells = character(0),
             sem = res$se_sub_same, cd = res$cd_sub_same, sig = sig_ab,
             cd2 = res$cd_main_same_sub,
             cd_labels = c(paste(sub_plot, "at same", main_plot),
                           paste(main_plot, "at same", sub_plot))))
      cv <- c(`CV(a) (%)` = res$cv_a, `CV(b) (%)` = res$cv_b)
      stats_rows[[tr]] <- data.frame(
        Trait = tr, Factor = c(main_plot, sub_plot,
                               paste(main_plot, "x", sub_plot)),
        F_value = c(at[rn == main_plot, "F value"],
                    at[rn == sub_plot, "F value"],
                    at[rn %in% c(paste0(main_plot, ":", sub_plot),
                                 paste0(sub_plot, ":", main_plot)),
                       "F value"]),
        p_value = c(pa, pb, pab),
        SEm = c(res$se_main, res$se_sub, res$se_sub_same),
        SEd = sqrt(2) * c(res$se_main, res$se_sub, res$se_sub_same),
        CD = c(res$cd_main, res$cd_sub, res$cd_sub_same),
        CV = c(res$cv_a, res$cv_b, res$cv_b), stringsAsFactors = FALSE)
    }
    results[[tr]] <- res
    sections[[tr]] <- list(sec = sec, cv = cv, dg = dg)
  }

  # Assemble the character table
  first <- sections[[1]]
  rows_lab <- character(0)
  build <- function(s, dg) {
    out <- character(0)
    for (blk in s$sec) {
      if (!is.null(blk$title)) out <- c(out, "")
      out <- c(out, blk$cells)
      if (identical(blk$factor, "AxB")) {
        if (!is.null(blk$cd2)) {
          out <- c(out, if (blk$sig) fmt(blk$cd, dg) else "NS",
                   if (blk$sig) fmt(blk$cd2, dg) else "NS")
        } else {
          out <- c(out, fmt(blk$sem, dg), if (blk$sig) fmt(blk$cd, dg) else
            "NS")
        }
      } else {
        out <- c(out, fmt(blk$sem, dg), if (blk$sig) fmt(blk$cd, dg) else
          "NS")
      }
    }
    c(out, formatC(s$cv, format = "f", digits = 2))
  }
  for (blk in first$sec) {
    if (!is.null(blk$title)) rows_lab <- c(rows_lab, blk$title)
    rows_lab <- c(rows_lab, blk$levels)
    if (identical(blk$factor, "AxB")) {
      if (!is.null(blk$cd2)) {
        rows_lab <- c(rows_lab, paste0(cd_lab, ": ", blk$cd_labels))
      } else {
        rows_lab <- c(rows_lab, paste("SEm", pm), cd_lab)
      }
    } else {
      rows_lab <- c(rows_lab, paste("SEm", pm), cd_lab)
    }
  }
  rows_lab <- c(rows_lab, names(first$cv))
  first_col <- switch(design, crd = treatment, rbd = treatment,
                      factorial = "Treatment", spd = "Treatment")
  tbl <- data.frame(rows_lab, stringsAsFactors = FALSE)
  names(tbl) <- first_col
  for (tr in traits) {
    tbl[[tr]] <- build(sections[[tr]], sections[[tr]]$dg)
  }
  st <- do.call(rbind, stats_rows)
  rownames(st) <- NULL
  out <- list(table = tbl, stats = st, results = results, design = design,
              posthoc = method, alpha = alpha,
              letters_method = .letters_method(method))
  class(out) <- "aridagri_traits"
  if (verbose) print(out)
  invisible(out)
}

#' Print Method for Multi-Trait Tables
#'
#' @param x An object of class \code{aridagri_traits}.
#' @param ... Unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.aridagri_traits <- function(x, ...) {
  .heading(paste0("TREATMENT MEANS FOR ", ncol(x$table) - 1, " TRAIT(S) (",
                  toupper(x$design), ")"))
  print(x$table, row.names = FALSE, right = FALSE)
  cat("\n Means followed by the same letter do not differ significantly (",
      x$letters_method, ", p = ", x$alpha, ").\n", sep = "")
  cat(" NS: not significant. Full results for each trait are in $results.\n")
  invisible(x)
}
