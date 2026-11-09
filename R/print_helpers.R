# =============================================================================
# PRINTING HELPERS
# Plain-text report layout shared by the analysis functions and their print
# methods. Not exported.
# =============================================================================

# Print a header line with design information.
.print_info <- function(title, info) {
  .heading(title)
  lab <- names(info)
  w <- max(nchar(lab))
  for (i in seq_along(info)) {
    cat(" ", formatC(lab[i], width = -w), " : ", info[[i]], "\n", sep = "")
  }
}

# Print an ANOVA data frame (columns Df, Sum Sq, Mean Sq, F value, Pr(>F))
# with readable source labels and significance codes.
.print_anova_df <- function(tab, labels = rownames(tab), total = TRUE,
                            title = "Analysis of variance") {
  cat("\n ", title, "\n", sep = "")
  df <- tab[["Df"]]
  ss <- tab[["Sum Sq"]]
  ms <- tab[["Mean Sq"]]
  fv <- if ("F value" %in% names(tab)) tab[["F value"]] else rep(NA, nrow(tab))
  pv <- if ("Pr(>F)" %in% names(tab)) tab[["Pr(>F)"]] else rep(NA, nrow(tab))
  src <- c(labels, if (total) "Total")
  src <- formatC(src, width = -max(nchar(src)))
  dfs <- c(formatC(df, format = "fg"), if (total) formatC(sum(df),
                                                           format = "fg"))
  dig <- if (max(abs(ss), na.rm = TRUE) >= 1000) 2 else 4
  sss <- c(.fmt(ss, dig), if (total) .fmt(sum(ss), dig))
  mss <- c(.fmt(ms, dig), if (total) "")
  fvs <- c(ifelse(is.na(fv), "", formatC(fv, format = "f", digits = 3)),
           if (total) "")
  pvs <- c(ifelse(is.na(pv), "", formatC(pv, format = "f", digits = 4)),
           if (total) "")
  sig <- c(ifelse(is.na(pv), "", .sig_code(pv)), if (total) "")
  out <- data.frame(Source = src, Df = dfs, `Sum Sq` = sss, `Mean Sq` = mss,
                    `F value` = fvs, `Pr(>F)` = pvs, ` ` = sig,
                    check.names = FALSE, stringsAsFactors = FALSE)
  print(out, row.names = FALSE, right = TRUE)
  cat(" Signif. codes: *** p < 0.001, ** p < 0.01, * p < 0.05, NS not",
      "significant\n")
  invisible(NULL)
}

# Print a means table with letters.
.print_means <- function(tbl, title, digits = 2) {
  cat("\n ", title, "\n", sep = "")
  out <- tbl
  for (cl in names(out)) {
    if (is.numeric(out[[cl]])) {
      d <- if (cl %in% c("SE", "SD", "SEm")) max(digits, 2) + 1 else digits
      if (cl %in% c("N", "n", "Obs", "Rank")) d <- 0
      out[[cl]] <- formatC(out[[cl]], format = "f", digits = d)
    }
  }
  print(out, row.names = FALSE, right = TRUE)
  invisible(NULL)
}

# Print a block of "label = value" statistics on one or more lines.
.print_stats_line <- function(x, digits = 3) {
  parts <- vapply(names(x), function(n) {
    v <- x[[n]]
    if (is.numeric(v)) v <- formatC(v, format = "f", digits = digits)
    paste0(n, " = ", v)
  }, character(1))
  cat("\n ", paste(parts, collapse = "   "), "\n", sep = "")
}
