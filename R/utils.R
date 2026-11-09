# =============================================================================
# INTERNAL UTILITIES
# Input checks, factor handling and small formatting helpers shared by the
# analysis functions. None of these are exported.
# =============================================================================

# Stop unless every name in 'cols' is a column of 'data'.
.check_cols <- function(data, cols, what = "Column") {
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  cols <- unlist(cols, use.names = FALSE)
  cols <- cols[!is.null(cols)]
  if (!is.character(cols) || anyNA(cols)) {
    stop("Column names must be given as character strings.", call. = FALSE)
  }
  miss <- setdiff(cols, names(data))
  if (length(miss) > 0) {
    stop(what, " not found in data: ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  invisible(TRUE)
}

# Stop unless 'x' is numeric.
.check_numeric <- function(x, name) {
  if (!is.numeric(x)) {
    stop("'", name, "' must be numeric.", call. = FALSE)
  }
  invisible(TRUE)
}

# Stop unless 0 < alpha < 1.
.check_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single number between 0 and 1.", call. = FALSE)
  }
  invisible(TRUE)
}

# Prepare an analysis data frame: keep the named columns, convert factors,
# drop rows with a missing response (with a message) and drop unused levels.
.prep_data <- function(data, response, factors, verbose = TRUE) {
  .check_cols(data, c(response, factors))
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  .check_numeric(data[[response]], response)
  d <- data[, unique(c(response, factors)), drop = FALSE]
  for (f in factors) {
    if (anyNA(d[[f]])) {
      stop("Factor '", f, "' contains missing values.", call. = FALSE)
    }
    d[[f]] <- if (is.factor(d[[f]])) d[[f]] else
      factor(d[[f]], levels = unique(d[[f]]))
  }
  miss <- is.na(d[[response]])
  if (any(miss)) {
    if (verbose) {
      message(sum(miss), " row(s) with missing '", response,
              "' removed before analysis.")
    }
    d <- d[!miss, , drop = FALSE]
  }
  for (f in factors) d[[f]] <- droplevels(d[[f]])
  rownames(d) <- NULL
  d
}

# Build a formula from character pieces, back-quoting non-syntactic names.
.bq <- function(x) {
  ifelse(make.names(x) == x, x, paste0("`", x, "`"))
}
.formula <- function(response, rhs) {
  stats::as.formula(paste(.bq(response), "~", rhs))
}

# Significance code for a p-value.
.sig_code <- function(p) {
  vapply(p, function(v) {
    if (is.na(v)) "" else if (v < 0.001) "***" else if (v < 0.01) "**" else
      if (v < 0.05) "*" else "NS"
  }, character(1))
}

# Plus-minus sign that degrades to "+/-" when the session is not UTF-8.
.pm <- function() {
  if (isTRUE(l10n_info()[["UTF-8"]])) "\u00b1" else "+/-"
}

# Multiplication sign for printed labels (ASCII "x" when not UTF-8).
.times <- function() {
  if (isTRUE(l10n_info()[["UTF-8"]])) "\u00d7" else "x"
}

# Number of complete replications in a balanced layout, or NA.
.reps_per_cell <- function(d, factors) {
  tab <- table(d[, factors, drop = FALSE])
  if (all(tab == tab[1])) as.numeric(tab[1]) else NA_real_
}

# Every combination of 'factors' present exactly once per level of 'rep'?
.is_balanced <- function(d, rep, factors) {
  tab <- table(d[, c(rep, factors), drop = FALSE])
  all(tab == 1)
}

# Format numbers for printed tables.
.fmt <- function(x, digits = 2) {
  ifelse(is.na(x), "-", formatC(x, format = "f", digits = digits))
}

# Print a rule of a given width.
.rule <- function(width = 70, char = "-") {
  cat(strrep(char, width), "\n", sep = "")
}

# Heading for printed reports.
.heading <- function(title, width = 70) {
  cat("\n")
  .rule(width, "=")
  cat(" ", title, "\n", sep = "")
  .rule(width, "=")
}

# Round the numeric columns of an ANOVA-type data frame for printing.
.round_anova <- function(tab, digits = 3) {
  out <- as.data.frame(tab, stringsAsFactors = FALSE)
  for (cl in intersect(c("Sum Sq", "Mean Sq", "F value"), names(out))) {
    out[[cl]] <- round(out[[cl]], digits)
  }
  if ("Pr(>F)" %in% names(out)) out[["Pr(>F)"]] <- signif(out[["Pr(>F)"]], 4)
  out
}

# Add a significance-code column to an ANOVA data frame.
.add_sig <- function(tab) {
  tab <- as.data.frame(tab, stringsAsFactors = FALSE)
  if ("Pr(>F)" %in% names(tab)) {
    s <- .sig_code(tab[["Pr(>F)"]])
    tab$Sig <- s
  }
  tab
}

# Cochran-Cox weighted t for a variance that is a positive combination of
# several error mean squares: t' = sum(w_k t_k) / sum(w_k), w_k = c_k E_k.
.weighted_t <- function(coef, ms, df, alpha) {
  w <- coef * ms
  tk <- stats::qt(1 - alpha / 2, df)
  sum(w * tk) / sum(w)
}

# Satterthwaite degrees of freedom for sum(c_k E_k).
.satterthwaite <- function(coef, ms, df) {
  w <- coef * ms
  sum(w)^2 / sum(w^2 / df)
}

# Remove back-quotes that anova() adds to non-syntactic term names.
.clean_rownames <- function(tab) {
  rownames(tab) <- trimws(gsub("`", "", rownames(tab)))
  tab
}

# Null-coalescing helper.
`%||%` <- function(a, b) if (is.null(a)) b else a
