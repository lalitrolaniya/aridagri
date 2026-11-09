# =============================================================================
# INTERNAL ENGINE FOR SPLIT-PLOT TYPE DESIGNS
# Correct error strata, SE/CD for every type of comparison, and letters.
# Not exported.
# =============================================================================

# Stop unless every combination of 'facs' occurs once in every level of
# 'rep' (complete, balanced split-plot data).
.check_split_balance <- function(d, rep, facs, extra = NULL) {
  tab <- table(d[, c(extra, rep, facs), drop = FALSE])
  if (any(tab != 1)) {
    n_miss <- sum(tab == 0)
    n_dup <- sum(tab > 1)
    stop("Split-plot analysis needs complete data: every treatment ",
         "combination exactly once in every replication",
         if (!is.null(extra)) " of every environment" else "", ". Found ",
         n_miss, " missing and ", n_dup, " duplicated cell(s). Estimate ",
         "missing plots first, or fit a mixed model (for example ",
         "lme4::lmer or nlme::lme).", call. = FALSE)
  }
  invisible(TRUE)
}

# Standard errors of difference for every type of comparison in a
# hierarchical (split-plot type) design.
#   facs  : factor names
#   lev_of: stratum (1 = top) of each factor, named by factor
#   nl    : number of levels of each factor, named by factor
#   r     : replications (units at the top stratum per factor level)
#   E, DFE: error mean squares and df by stratum (length K)
# For two cell means of effect F that first differ at stratum s,
#   Var = (2/r) sum_{l >= s} sigma2_l / u_l,
# with sigma2_l = (E_l - E_{l+1}) / n_l (n_l = observations per stratum-l
# unit, E_{K+1} = 0) and u_l = product of the levels of factors not in F
# at strata <= l. This reproduces the standard formulas of Gomez and Gomez
# (1984, tables 3.10 and 4.?) and extends them to any number of strata.
.split_sed <- function(F, s, facs, lev_of, nl, r, E) {
  K <- length(E)
  n_unit <- vapply(seq_len(K), function(l) prod(nl[lev_of > l]), numeric(1))
  u <- vapply(seq_len(K), function(l) {
    prod(nl[setdiff(facs[lev_of <= l], F)])
  }, numeric(1))
  co <- numeric(K)
  for (m in s:K) {
    co[m] <- 1 / (n_unit[m] * u[m])
    if (m - 1 >= s) co[m] <- co[m] - 1 / (n_unit[m - 1] * u[m - 1])
  }
  co <- (2 / r) * co
  co[abs(co) < 1e-14] <- 0
  co
}

# Table of SEd and CD for main effects and interactions up to 'max_order'.
.split_comparisons <- function(facs, lev_of, nl, r, E, DFE, alpha,
                               max_order = 3, labels = NULL) {
  rows <- list()
  for (ord in seq_len(min(max_order, length(facs)))) {
    sets <- utils::combn(facs, ord, simplify = FALSE)
    for (F in sets) {
      S <- sort(unique(lev_of[F]))
      for (s in S) {
        co <- .split_sed(F, s, facs, lev_of, nl, r, E)
        v <- sum(co * E)
        sed <- sqrt(v)
        use <- co > 0
        tw <- .weighted_t(co[use], E[use], DFE[use], alpha)
        dfs <- .satterthwaite(co[use], E[use], DFE[use])
        at_s <- F[lev_of[F] == s]
        above <- F[lev_of[F] < s]
        below <- F[lev_of[F] > s]
        if (ord == 1) {
          lab <- paste(F, "means")
        } else {
          lab <- paste(paste(at_s, collapse = " x "), "levels")
          if (length(above)) {
            lab <- paste0(lab, " at the same ", paste(above, collapse = " x "))
          }
          if (length(below)) {
            lab <- paste0(lab, if (length(above)) "," else "",
                          " at the same or different ",
                          paste(below, collapse = " x "))
          }
        }
        err <- paste(paste0("E", letters[seq_along(E)][use]), collapse = "+")
        rows[[length(rows) + 1L]] <- data.frame(
          Effect = paste(F, collapse = " x "), Comparison = lab,
          SEm = sed / sqrt(2), SEd = sed, t_value = tw, df = dfs,
          CD = tw * sed, Error = err, stratum = s,
          stringsAsFactors = FALSE)
      }
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (!is.null(labels)) out$Error <- out$Error
  out
}

# Letters for the levels of factor X within each level of factor Y (Y may
# be NULL for a main effect). Means are cell means of the balanced data.
.split_letters <- function(d, response, X, Y, sed, df, tcrit, methods,
                           alpha, mse) {
  if (length(methods) == 0) return(NULL)
  m1 <- methods[1]
  y <- d[[response]]
  if (is.null(Y)) {
    mx <- tapply(y, d[[X]], mean)
    k <- length(mx)
    V <- diag(sed^2 / 2, k)
    dimnames(V) <- list(names(mx), names(mx))
    tc <- if (m1 == "lsd") tcrit else NULL
    cm <- .compare(mx, V, df, mse, m1, alpha, t_crit = tc,
                   response = response)
    return(data.frame(Level = names(mx), Mean = as.numeric(mx),
                      Group = cm$letters[names(mx)],
                      stringsAsFactors = FALSE))
  }
  out <- list()
  for (yl in levels(d[[Y]])) {
    sub <- d[d[[Y]] == yl, , drop = FALSE]
    mx <- tapply(sub[[response]], droplevels(sub[[X]]), mean)
    k <- length(mx)
    V <- diag(sed^2 / 2, k)
    dimnames(V) <- list(names(mx), names(mx))
    tc <- if (m1 == "lsd") tcrit else NULL
    cm <- .compare(mx, V, df, mse, m1, alpha, t_crit = tc,
                   response = response)
    out[[yl]] <- data.frame(Y = yl, X = names(mx), Mean = as.numeric(mx),
                            Group = cm$letters[names(mx)],
                            stringsAsFactors = FALSE)
  }
  res <- do.call(rbind, out)
  names(res)[1:2] <- c(Y, X)
  rownames(res) <- NULL
  res
}

# Split-plot family ANOVA. 'groups' is a list of character vectors, one
# per stratum from the main plot down, e.g. list(main = "A", sub = "B").
.aridagri_spd_anova <- function(data, response, replication, groups,
                                alpha = 0.05, methods = "lsd",
                                verbose = TRUE) {
  facs <- unlist(groups, use.names = FALSE)
  nlev <- length(groups)
  lab_names <- names(groups)
  lev_of <- integer(0)
  for (g in seq_along(groups)) for (f in groups[[g]]) lev_of[f] <- g
  d <- .prep_data(data, response, c(replication, facs), verbose)
  .check_split_balance(d, replication, facs)
  nl <- vapply(facs, function(f) nlevels(d[[f]]), numeric(1))
  names(nl) <- facs
  if (any(nl < 2)) stop("Each factor needs at least 2 levels.", call. = FALSE)
  r <- nlevels(d[[replication]])
  if (r < 2) stop("Need at least 2 replications.", call. = FALSE)
  y <- d[[response]]
  N <- length(y)

  form <- .formula(response, paste(.bq(replication), "*",
                                   paste(.bq(facs), collapse = " * ")))
  full <- suppressWarnings(stats::aov(form, data = d))
  # The full rep x factor crossing is saturated; its residual is not used.
  af <- suppressWarnings(stats::anova(full))
  rn <- trimws(gsub("`", "", rownames(af)))
  SS <- af[["Sum Sq"]]
  DF <- af[["Df"]]
  terms <- strsplit(rn, ":", fixed = TRUE)
  elevel <- function(tf) {
    nz <- setdiff(tf, c(replication, "Residuals"))
    if (!length(nz)) NA_integer_ else max(lev_of[nz])
  }
  has_rep <- vapply(terms, function(tf) replication %in% tf, logical(1))
  is_resid <- rn == "Residuals"
  errSS <- numeric(nlev)
  errDF <- numeric(nlev)
  for (i in seq_along(rn)) {
    if (is_resid[i]) {
      errSS[nlev] <- errSS[nlev] + SS[i]
      errDF[nlev] <- errDF[nlev] + DF[i]
    } else if (has_rep[i] && rn[i] != replication) {
      L <- elevel(terms[[i]])
      errSS[L] <- errSS[L] + SS[i]
      errDF[L] <- errDF[L] + DF[i]
    }
  }
  errMS <- errSS / errDF
  errMS[!is.finite(errMS)] <- NA
  if (any(errDF < 1)) stop("An error stratum has no degrees of freedom.",
                           call. = FALSE)

  trt_idx <- which(!has_rep & !is_resid & rn != replication)
  trt_lvl <- vapply(trt_idx, function(i) elevel(terms[[i]]), integer(1))
  grand <- mean(y)
  rows <- list()
  add <- function(s, df, ss, ms, Fv, p) {
    rows[[length(rows) + 1L]] <<- data.frame(
      .s = s, Df = df, `Sum Sq` = ss, `Mean Sq` = ms, `F value` = Fv,
      `Pr(>F)` = p, check.names = FALSE, stringsAsFactors = FALSE)
  }
  ri <- which(rn == replication)
  add(replication, DF[ri], SS[ri], SS[ri] / DF[ri], SS[ri] / DF[ri] / errMS[1],
      stats::pf(SS[ri] / DF[ri] / errMS[1], DF[ri], errDF[1],
                lower.tail = FALSE))
  for (L in seq_len(nlev)) {
    for (i in trt_idx[trt_lvl == L]) {
      ms <- SS[i] / DF[i]
      Fv <- ms / errMS[L]
      add(rn[i], DF[i], SS[i], ms, Fv,
          stats::pf(Fv, DF[i], errDF[L], lower.tail = FALSE))
    }
    add(paste0("Error(", letters[L], ")"), errDF[L], errSS[L], errMS[L], NA,
        NA)
  }
  tab <- do.call(rbind, rows)
  rownames(tab) <- make.unique(tab$.s)
  tab$.s <- NULL

  factor_means <- lapply(facs, function(f) {
    m <- tapply(y, d[[f]], mean)
    data.frame(Level = names(m), Mean = as.numeric(m),
               stringsAsFactors = FALSE)
  })
  names(factor_means) <- facs
  cv <- stats::setNames(sqrt(errMS) / grand * 100, lab_names)

  comps <- .split_comparisons(facs, lev_of, nl, r, errMS, errDF, alpha)
  # Main-effect summary kept for compatibility with earlier versions.
  main_rows <- comps[comps$Effect %in% facs, ]
  comp <- data.frame(Factor = main_rows$Effect,
                     Levels = as.numeric(nl[main_rows$Effect]),
                     SE_diff = main_rows$SEd, CD = main_rows$CD,
                     Error = paste0("Error(", letters[main_rows$stratum], ")"),
                     stringsAsFactors = FALSE)

  # Letters: main effects, then each two-factor interaction both ways.
  letters_out <- list()
  if (length(methods) > 0) {
    for (f in facs) {
      cr <- main_rows[main_rows$Effect == f, ]
      lt <- .split_letters(d, response, f, NULL, cr$SEd, cr$df, cr$t_value,
                           methods, alpha, errMS[lev_of[f]])
      factor_means[[f]]$Group <- lt$Group[match(factor_means[[f]]$Level,
                                                lt$Level)]
    }
    if (length(facs) >= 2) {
      for (pr in utils::combn(facs, 2, simplify = FALSE)) {
        for (dirn in 1:2) {
          X <- pr[dirn]
          Y <- pr[3 - dirn]
          cr <- comps[comps$Effect == paste(pr, collapse = " x ") &
                        comps$stratum == lev_of[X], ]
          lt <- .split_letters(d, response, X, Y, cr$SEd, cr$df, cr$t_value,
                               methods, alpha, errMS[lev_of[X]])
          letters_out[[paste(X, "within", Y)]] <- lt
        }
      }
    }
  }
  inter_means <- list()
  if (length(facs) >= 2) {
    for (pr in utils::combn(facs, 2, simplify = FALSE)) {
      m <- tapply(y, list(d[[pr[1]]], d[[pr[2]]]), mean)
      inter_means[[paste(pr, collapse = " x ")]] <- m
    }
  }
  list(anova_table = tab, grand_mean = grand,
       error_ms = stats::setNames(errMS, lab_names),
       error_df = stats::setNames(errDF, lab_names),
       cv = cv, factor_means = factor_means, comparisons = comp,
       se_cd = comps, letters = letters_out, interaction_means = inter_means,
       model = full, data = d, r = r, levels = nl, lev_of = lev_of)
}

# Printer for the split-plot family.
.aridagri_print_spd <- function(core, design, response = NULL,
                                letters_method = NA, alpha = 0.05,
                                digits = 2) {
  info <- c(Design = design)
  if (!is.null(response)) info <- c(Response = response, info)
  info <- c(info, Replications = core$r,
            Factors = paste(paste0(names(core$levels), " (", core$levels, ")"),
                            collapse = ", "))
  .print_info("SPLIT-PLOT TYPE DESIGN", info)
  tab <- core$anova_table
  .print_anova_df(tab, rownames(tab))
  cat(sprintf("\n Grand mean = %.3f\n", core$grand_mean))
  nms <- names(core$error_ms)
  for (k in seq_along(nms)) {
    cat(sprintf("  Error(%s) [%s]: df = %d, MS = %.4f, CV = %.2f%%\n",
                letters[k], nms[k], as.integer(core$error_df[[k]]),
                core$error_ms[[k]], core$cv[[k]]))
  }
  cat("\n SEm, SEd and CD (alpha = ", alpha, ")\n", sep = "")
  sc <- core$se_cd[, c("Comparison", "SEm", "SEd", "CD", "Error")]
  sc$SEm <- round(sc$SEm, 3)
  sc$SEd <- round(sc$SEd, 3)
  sc$CD <- round(sc$CD, 3)
  print(sc, row.names = FALSE, right = FALSE)
  cat(" CD for comparisons mixing error terms uses the weighted t of",
      "Cochran and Cox.\n")
  lm_lab <- if (is.na(letters_method)) "" else paste0(" (", letters_method, ")")
  for (f in names(core$factor_means)) {
    .print_means(core$factor_means[[f]], paste0("Means of ", f, lm_lab),
                 digits)
  }
  if (length(core$interaction_means) > 0) {
    for (nm in names(core$interaction_means)) {
      cat("\n Interaction means: ", nm, "\n", sep = "")
      print(round(core$interaction_means[[nm]], digits))
    }
  }
  invisible(NULL)
}

## ---------------------------------------------------------------------------
## Combined ("pooled") split-plot family over environments. Environments
## fixed, replications random within environment. Each treatment source is
## tested against the error at the level of its deepest split factor; the
## environment main effect against the environment/replication error
## (E:Rep). Balanced data are required.
.aridagri_spd_pooled_anova <- function(data, response, environment,
                                       replication, split_factors,
                                       alpha = 0.05, methods = "lsd",
                                       verbose = TRUE) {
  E <- environment
  R <- replication
  FS <- split_factors
  K <- length(FS)
  d <- .prep_data(data, response, c(E, R, FS), verbose)
  .check_split_balance(d, R, FS, extra = E)
  nl <- vapply(c(E, FS), function(f) nlevels(d[[f]]), numeric(1))
  names(nl) <- c(E, FS)
  if (nl[[E]] < 2) stop("Need at least 2 environments.", call. = FALSE)
  r <- nlevels(droplevels(d[[R]]))
  if (r < 2) stop("Need at least 2 replications per environment.",
                  call. = FALSE)

  trt <- paste(.bq(c(E, FS)), collapse = " * ")
  errs <- paste(.bq(E), .bq(R), sep = ":")
  if (K >= 2) {
    for (k in 1:(K - 1)) {
      errs <- c(errs, paste(.bq(c(E, R, FS[1:k])), collapse = ":"))
    }
  }
  form <- .formula(response, paste(trt, "+", paste(errs, collapse = " + ")))
  model <- stats::aov(form, data = d)
  at <- suppressWarnings(stats::anova(model))
  terms <- trimws(gsub("`", "", rownames(at)))
  MS <- stats::setNames(at[["Mean Sq"]], terms)
  DF <- stats::setNames(at[["Df"]], terms)
  SS <- stats::setNames(at[["Sum Sq"]], terms)
  keyof <- function(v) paste(sort(v), collapse = "|")
  term_key <- vapply(terms, function(tm) keyof(strsplit(tm, ":",
                                                        fixed = TRUE)[[1]]),
                     character(1))
  find_term <- function(fset) {
    hit <- terms[term_key == keyof(fset)]
    if (length(hit) == 1) hit else NA_character_
  }
  level_error_name <- function(maxlev) {
    if (maxlev == 0) return(find_term(c(E, R)))
    if (maxlev == K) return("Residuals")
    find_term(c(E, R, FS[1:maxlev]))
  }
  split_level <- function(tm) {
    parts <- strsplit(tm, ":", fixed = TRUE)[[1]]
    lv <- 0L
    for (i in seq_along(FS)) if (FS[i] %in% parts) lv <- i
    lv
  }
  err_rows <- c(vapply(0:(K - 1), function(k) {
    if (k == 0) find_term(c(E, R)) else find_term(c(E, R, FS[1:k]))
  }, character(1)), "Residuals")
  Fv <- stats::setNames(rep(NA_real_, length(terms)), terms)
  Pv <- Fv
  for (tm in terms) {
    if (tm %in% err_rows) next
    den <- level_error_name(split_level(tm))
    Fv[tm] <- MS[[tm]] / MS[[den]]
    Pv[tm] <- stats::pf(Fv[tm], DF[[tm]], DF[[den]], lower.tail = FALSE)
  }
  # order: treatment terms by stratum, each followed by its error
  ord <- character(0)
  for (lv in 0:K) {
    tt <- terms[!terms %in% err_rows & vapply(terms, split_level,
                                              integer(1)) == lv]
    ord <- c(ord, tt, err_rows[lv + 1])
  }
  out <- data.frame(Df = DF[ord], `Sum Sq` = SS[ord], `Mean Sq` = MS[ord],
                    `F value` = Fv[ord], `Pr(>F)` = Pv[ord],
                    check.names = FALSE, row.names = ord)
  lab <- rownames(out)
  lab[lab == err_rows[1]] <- paste0(R, "(", E, ")")
  for (k in seq_len(K - 1)) {
    lab[lab == err_rows[k + 1]] <- paste0("Error(", letters[k], ")")
  }
  lab[lab == "Residuals"] <- paste0("Error(", letters[K], ")")
  rownames(out) <- make.unique(lab)

  gm <- mean(d[[response]])
  Ems <- vapply(err_rows, function(t) MS[[t]], numeric(1))
  Edf <- vapply(err_rows, function(t) DF[[t]], numeric(1))
  err_labels <- c("a (main-plot)", "b (sub-plot)", "c (sub-sub-plot)", "d",
                  "e")
  cvs <- c()
  cvs["Environment (E:Rep)"] <- sqrt(Ems[1]) / gm * 100
  for (k in 1:K) {
    cvs[paste0("Error ", err_labels[k])] <- sqrt(Ems[k + 1]) / gm * 100
  }
  # SE / CD: environment is a stratum-1 factor; split factors follow.
  facs <- c(E, FS)
  lev_of <- stats::setNames(c(1L, seq_len(K) + 1L), facs)
  comps <- .split_comparisons(facs, lev_of, nl, r, Ems, Edf, alpha,
                              max_order = 3)
  comps$Error <- vapply(strsplit(comps$Error, "+", fixed = TRUE), function(z) {
    z <- sub("^Ea$", "E(R:E)", z)
    for (k in seq_len(K)) z <- sub(paste0("^E", letters[k + 1], "$"),
                                   paste0("E", letters[k]), z)
    paste(z, collapse = "+")
  }, character(1))
  # Letters for the split factors (environment fixed).
  factor_means <- lapply(c(FS, E), function(f) {
    m <- tapply(d[[response]], d[[f]], mean)
    tb <- data.frame(Level = names(m), Mean = as.numeric(m),
                     stringsAsFactors = FALSE)
    if (length(methods) > 0) {
      cr <- comps[comps$Effect == f, ]
      lt <- .split_letters(d, response, f, NULL, cr$SEd, cr$df, cr$t_value,
                           methods, alpha, Ems[lev_of[f]])
      tb$Group <- lt$Group[match(tb$Level, lt$Level)]
    }
    tb
  })
  names(factor_means) <- c(FS, E)
  list(anova_table = out, model = model, grand_mean = gm, cv = cvs,
       error_terms = err_rows, split_factors = FS, environment = E,
       se_cd = comps, factor_means = factor_means, r = r, levels = nl,
       error_ms = Ems, error_df = Edf, data = d)
}

# Printer for the pooled split-plot family.
.aridagri_print_spd_pooled <- function(core, design, letters_method = NA,
                                       alpha = 0.05, digits = 2) {
  .print_info(design, c(Environments = core$levels[[core$environment]],
                        `Reps per environment` = core$r,
                        Factors = paste(core$split_factors, collapse = ", ")))
  cat(" Environments fixed; replications random within environment.\n")
  cat(" Each source is tested against the error at the level of its\n")
  cat(" deepest split factor; environment is tested against Rep(E).\n")
  at <- core$anova_table
  .print_anova_df(at, rownames(at))
  cat(sprintf("\n Grand mean = %.4f\n", core$grand_mean))
  cat("\n Coefficient of variation by stratum (%):\n")
  cvp <- data.frame(Stratum = names(core$cv), CV = round(as.numeric(core$cv),
                                                          2))
  print(cvp, row.names = FALSE)
  cat("\n SEm, SEd and CD (alpha = ", alpha, ")\n", sep = "")
  sc <- core$se_cd[, c("Comparison", "SEm", "SEd", "CD", "Error")]
  sc$SEm <- round(sc$SEm, 3)
  sc$SEd <- round(sc$SEd, 3)
  sc$CD <- round(sc$CD, 3)
  print(sc, row.names = FALSE, right = FALSE)
  lm_lab <- if (is.na(letters_method)) "" else paste0(" (", letters_method, ")")
  for (f in names(core$factor_means)) {
    .print_means(core$factor_means[[f]], paste0("Means of ", f, lm_lab),
                 digits)
  }
  invisible(NULL)
}
