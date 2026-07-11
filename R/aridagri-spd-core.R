# Internal helpers for split-plot / split-split-plot ANOVA.
# Not exported. These compute the correct error strata by pooling the
# replication-by-treatment interactions to the proper nesting level, so that
# every treatment source is tested against the error term at the level of its
# most-nested factor (main-plot factors against the main-plot error, sub-plot
# factors against the sub-plot error, and so on).

.aridagri_spd_anova <- function(data, response, replication, groups, alpha = 0.05) {
  facs <- unlist(groups, use.names = FALSE)
  nlev <- length(groups)
  lab_names <- names(groups)
  fac_level <- integer(0)
  for (g in seq_along(groups)) for (f in groups[[g]]) fac_level[f] <- g

  data[[replication]] <- as.factor(data[[replication]])
  for (f in facs) data[[f]] <- as.factor(data[[f]])
  y <- data[[response]]
  N <- length(y)

  form <- stats::as.formula(paste0("`", response, "` ~ `", replication, "` * ",
                                   paste0("`", facs, "`", collapse = " * ")))
  full <- suppressWarnings(stats::aov(form, data = data))
  ## The full rep*factor crossing is saturated (one observation per cell), so
  ## anova() would warn that the full-model residual is an essentially perfect
  ## fit. That residual is not used here: each effect is tested against a pooled
  ## error stratum built below. Suppress the spurious warning.
  af <- suppressWarnings(stats::anova(full))
  rn <- trimws(gsub("`", "", rownames(af)))
  SS <- af[["Sum Sq"]]; DF <- af[["Df"]]
  terms <- strsplit(rn, ":")
  elevel <- function(tf) {
    nz <- setdiff(tf, c(replication, "Residuals"))
    if (!length(nz)) NA_integer_ else max(fac_level[nz])
  }
  has_rep  <- vapply(terms, function(tf) replication %in% tf, logical(1))
  is_resid <- rn == "Residuals"

  # pool replication-interactions and the residual into per-level error terms
  errSS <- numeric(nlev); errDF <- numeric(nlev)
  for (i in seq_along(rn)) {
    if (is_resid[i]) {
      errSS[nlev] <- errSS[nlev] + SS[i]; errDF[nlev] <- errDF[nlev] + DF[i]
    } else if (has_rep[i]) {
      L <- elevel(terms[[i]]); errSS[L] <- errSS[L] + SS[i]; errDF[L] <- errDF[L] + DF[i]
    }
  }
  errMS <- errSS / errDF
  errMS[!is.finite(errMS)] <- NA

  trt_idx <- which(!has_rep & !is_resid & rn != replication)
  trt_lvl <- vapply(trt_idx, function(i) elevel(terms[[i]]), integer(1))
  grand <- mean(y)

  rows <- list()
  add <- function(s, df, ss, ms, Fv, p)
    rows[[length(rows) + 1L]] <<- data.frame(
      .s = s, Df = df, `Sum Sq` = ss, `Mean Sq` = ms,
      `F value` = Fv, `Pr(>F)` = p, check.names = FALSE, stringsAsFactors = FALSE)
  ri <- which(rn == replication)
  add(replication, DF[ri], SS[ri], SS[ri] / DF[ri], NA, NA)
  for (L in seq_len(nlev)) {
    for (i in trt_idx[trt_lvl == L]) {
      ms <- SS[i] / DF[i]; Fv <- ms / errMS[L]
      add(rn[i], DF[i], SS[i], ms, Fv, stats::pf(Fv, DF[i], errDF[L], lower.tail = FALSE))
    }
    add(paste0("Error(", letters[L], ")"), errDF[L], errSS[L], errMS[L], NA, NA)
  }
  tab <- do.call(rbind, rows)
  rownames(tab) <- make.unique(tab$.s); tab$.s <- NULL

  factor_means <- lapply(facs, function(f) {
    m <- stats::aggregate(y, by = list(data[[f]]), FUN = mean)
    names(m) <- c("Level", "Mean"); m
  })
  names(factor_means) <- facs
  cv <- stats::setNames(sqrt(errMS) / grand * 100, lab_names)
  comp <- do.call(rbind, lapply(facs, function(f) {
    L <- fac_level[f]; k <- nlevels(data[[f]]); n_f <- N / k
    se <- sqrt(2 * errMS[L] / n_f); cd <- stats::qt(1 - alpha / 2, errDF[L]) * se
    data.frame(Factor = f, Levels = k, SE_diff = se, CD = cd,
               Error = paste0("Error(", letters[L], ")"), stringsAsFactors = FALSE)
  }))

  list(anova_table = tab, grand_mean = grand,
       error_ms = stats::setNames(errMS, lab_names),
       error_df = stats::setNames(errDF, lab_names),
       cv = cv, factor_means = factor_means, comparisons = comp, model = full)
}

.aridagri_print_spd <- function(core, design) {
  cat("\n", design, "\n", sep = "")
  cat(strrep("=", nchar(design)), "\n\n", sep = "")
  cat("Analysis of Variance\n")
  print(round(core$anova_table, 4))
  cat("\nGrand mean:", round(core$grand_mean, 3), "\n\n")
  cat("Error terms\n")
  nms <- names(core$error_ms)
  for (k in seq_along(nms))
    cat(sprintf("  Error(%s) [%s]: df = %d, MS = %.4f, CV = %.2f%%\n",
                letters[k], nms[k], core$error_df[[k]], core$error_ms[[k]], core$cv[[k]]))
  cat("\nStandard error of difference and critical difference (alpha = 0.05)\n")
  cmp <- core$comparisons[, c("Factor", "SE_diff", "CD")]
  cmp$SE_diff <- round(cmp$SE_diff, 3)
  cmp$CD <- round(cmp$CD, 3)
  print(cmp, row.names = FALSE)
  invisible(NULL)
}


## ---------------------------------------------------------------------------
## Internal helper for the combined ("pooled") split-plot family over
## environments. Environments are treated as fixed and replications as random
## within environment. Each treatment source is tested against the error term
## at the level of its deepest split factor; the environment main effect is
## tested against the environment/replication error (E:Rep). Balanced data are
## assumed (the calling functions warn otherwise). Sums of squares are taken
## from a single flat aov() model (which, for balanced data, are identical to
## the corresponding aov() Error()-strata sums of squares) and each effect is
## then assigned its correct denominator. Not exported.
.aridagri_spd_pooled_anova <- function(data, response, environment, replication,
                                       split_factors, alpha = 0.05) {
  E <- environment; R <- replication; FS <- split_factors; K <- length(FS)
  data[[E]] <- as.factor(data[[E]]); data[[R]] <- as.factor(data[[R]])
  for (f in FS) data[[f]] <- as.factor(data[[f]])

  trt  <- paste(c(E, FS), collapse = " * ")
  errs <- paste(E, R, sep = ":")
  if (K >= 2) for (k in 1:(K - 1))
    errs <- c(errs, paste(c(E, R, FS[1:k]), collapse = ":"))
  form  <- stats::as.formula(paste(response, "~", trt, "+", paste(errs, collapse = " + ")))
  model <- stats::aov(form, data = data)
  at    <- suppressWarnings(stats::anova(model))

  terms <- rownames(at)
  MS <- stats::setNames(at$`Mean Sq`, terms)
  DF <- stats::setNames(at$Df,        terms)
  SS <- stats::setNames(at$`Sum Sq`,  terms)

  keyof    <- function(v) paste(sort(v), collapse = "|")
  term_key <- vapply(terms, function(tm) keyof(strsplit(tm, ":", fixed = TRUE)[[1]]), character(1))
  find_term <- function(fset) { hit <- terms[term_key == keyof(fset)]; if (length(hit) == 1) hit else NA_character_ }

  level_error_name <- function(maxlev) {
    if (maxlev == 0) return(find_term(c(E, R)))
    if (maxlev == K) return("Residuals")
    find_term(c(E, R, FS[1:maxlev]))
  }
  split_level <- function(tm) {
    parts <- strsplit(tm, ":", fixed = TRUE)[[1]]
    lv <- 0L; for (i in seq_along(FS)) if (FS[i] %in% parts) lv <- i
    lv
  }
  err_rows <- c(vapply(0:(K - 1), function(k)
                  if (k == 0) find_term(c(E, R)) else find_term(c(E, R, FS[1:k])),
                  character(1)), "Residuals")

  Fv <- stats::setNames(rep(NA_real_, length(terms)), terms); Pv <- Fv
  for (tm in terms) {
    if (tm %in% err_rows) next
    den <- level_error_name(split_level(tm))
    Fv[tm] <- MS[[tm]] / MS[[den]]
    Pv[tm] <- stats::pf(Fv[tm], DF[[tm]], DF[[den]], lower.tail = FALSE)
  }

  out <- data.frame(Df = DF, `Sum Sq` = SS, `Mean Sq` = MS,
                    `F value` = Fv, `Pr(>F)` = Pv,
                    check.names = FALSE, row.names = terms)

  gm  <- mean(data[[response]], na.rm = TRUE)
  err_labels <- c("a (main-plot)", "b (sub-plot)", "c (sub-sub-plot)", "d", "e")
  cvs <- c(); cvs["Environment (E:Rep)"] <- sqrt(MS[[find_term(c(E, R))]]) / gm * 100
  for (k in 1:K) {
    den <- if (k == K) "Residuals" else find_term(c(E, R, FS[1:k]))
    cvs[paste0("Error ", err_labels[k])] <- sqrt(MS[[den]]) / gm * 100
  }

  list(anova_table = out, model = model, grand_mean = gm, cv = cvs,
       error_terms = err_rows, split_factors = FS, environment = E)
}

## Internal printer for the pooled split-plot family. Not exported.
.aridagri_print_spd_pooled <- function(core, design) {
  cat("\n", design, "\n", sep = "")
  cat(rep("=", 70), "\n", sep = "")
  cat("Environments fixed; replications random within environment.\n")
  cat("Each source tested against the error at the level of its deepest\n")
  cat("split factor; environment tested against E:Rep.\n\n")
  at <- core$anova_table
  print(round_df_signif(at))
  cat(sprintf("\nGrand Mean : %.4f\n", core$grand_mean))
  cat("\nCoefficient of variation by stratum (%):\n")
  cvp <- data.frame(Stratum = names(core$cv), CV = round(as.numeric(core$cv), 2))
  print(cvp, row.names = FALSE)
  invisible(NULL)
}

## small formatter shared by the pooled printer
round_df_signif <- function(df) {
  d <- df
  for (col in c("Sum Sq", "Mean Sq", "F value")) if (col %in% names(d)) d[[col]] <- round(d[[col]], 3)
  if ("Pr(>F)" %in% names(d)) d[["Pr(>F)"]] <- signif(d[["Pr(>F)"]], 4)
  d
}
