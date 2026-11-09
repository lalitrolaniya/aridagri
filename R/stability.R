# =============================================================================
# STABILITY ANALYSIS OF MULTI-ENVIRONMENT TRIALS
# Package: aridagri
# =============================================================================

#' Stability Analysis for Multi-Environment Trials (Seven Methods)
#'
#' @description
#' Stability analysis of genotypes (or treatments) tested in several
#' environments, with the combined ANOVA and seven established methods:
#' \itemize{
#'   \item Eberhart and Russell (1966) joint regression: \eqn{b_i},
#'     \eqn{S^2_{di}}, tests of \eqn{b_i = 1} and \eqn{S^2_{di} = 0}, and
#'     the Eberhart-Russell ANOVA
#'   \item AMMI: ANOVA with Gollob's degrees of freedom for each IPCA,
#'     IPCA scores, AMMI stability value (ASV) and yield stability index
#'   \item Finlay and Wilkinson (1963) regression on the environment mean
#'   \item Shukla (1972) stability variance with its F-test
#'   \item Wricke (1962) ecovalence
#'   \item Coefficient of variation (Francis and Kannenberg, 1978)
#'   \item Cultivar superiority index (Lin and Binns, 1988)
#' }
#' An integrated ranking averages the stability ranks of the seven methods.
#'
#' @param data Data frame with genotype, environment, replication and trait
#'   columns (one row per plot). Genotype x environment means (one row per
#'   combination) can be used with \code{replication = NULL}; tests that
#'   need the pooled error are then omitted.
#' @param genotype Name of the genotype (or treatment) column.
#' @param environment Name of the environment (location, year or season)
#'   column.
#' @param replication Name of the replication column, or \code{NULL}.
#' @param trait Name of the trait (response) column.
#' @param method \code{"all"} (default) or one or more of
#'   \code{"eberhart"}, \code{"ammi"}, \code{"finlay"}, \code{"shukla"},
#'   \code{"wricke"}, \code{"cv"}, \code{"superiority"}.
#' @param verbose Logical. If TRUE (default), prints the report.
#' @param alpha Significance level for the tests. Default 0.05.
#' @param fw_log Logical. If TRUE, the Finlay-Wilkinson regression uses
#'   log10-transformed means as in the original paper. Default FALSE (on
#'   the original scale, where the slope equals the Eberhart-Russell
#'   \eqn{b_i}).
#'
#' @details
#' Calculations use the genotype x environment table of means
#' (\eqn{g} genotypes, \eqn{e} environments, \eqn{r} replications). The
#' environmental index is \eqn{I_j = \bar{y}_{.j} - \bar{y}_{..}}.
#'
#' \strong{Eberhart-Russell.} \eqn{b_i = \sum_j y_{ij} I_j / \sum_j I_j^2};
#' \eqn{S^2_{di} = \sum_j \delta^2_{ij}/(e - 2) - MSE/r}, where
#' \eqn{\delta_{ij}} are deviations from the regression and MSE is the
#' pooled error. \eqn{b_i = 1} is tested with
#' \eqn{t = (b_i - 1)/SE(b_i)} on \eqn{e - 2} df, and \eqn{S^2_{di} = 0}
#' with \eqn{F = [\sum_j \delta^2_{ij}/(e-2)] / (MSE/r)}.
#'
#' \strong{Shukla.} \eqn{\sigma^2_i = [g(g-1)W_i - \sum W] /
#' [(e-1)(g-1)(g-2)]}, where \eqn{W_i} is Wricke's ecovalence
#' \eqn{\sum_j (\bar{y}_{ij} - \bar{y}_{i.} - \bar{y}_{.j} + \bar{y}_{..})^2};
#' tested with \eqn{F = \sigma^2_i / (MSE/r)} on \eqn{e - 1} and the error
#' df. Values are on the scale of the means.
#'
#' \strong{AMMI.} The interaction matrix is decomposed by singular value
#' decomposition. IPCA \eqn{k} has sum of squares \eqn{r \lambda_k^2} and
#' \eqn{g + e - 1 - 2k} df (Gollob, 1968). Scores are
#' \eqn{u_{ik}\sqrt{\lambda_k}}. ASV = \eqn{\sqrt{[(SS_1/SS_2)
#' IPCA1]^2 + IPCA2^2}} (Purchase et al., 2000) and the yield stability
#' index is the sum of the ranks of mean (high = 1) and ASV.
#'
#' \strong{Lin-Binns.} \eqn{P_i = \sum_j (y_{ij} - M_j)^2 / (2e)}, with
#' \eqn{M_j} the highest mean in environment \eqn{j}.
#'
#' @return An object of class \code{aridagri_stability} containing
#'   \code{anova} (combined ANOVA, when replications are given),
#'   \code{genotype_means}, \code{environment_means}, \code{ge_matrix}
#'   and, depending on \code{method}, \code{eberhart} (with
#'   \code{eberhart_anova}), \code{ammi} (\code{anova}, \code{components},
#'   \code{stability}, \code{genotype_scores}, \code{environment_scores}),
#'   \code{finlay}, \code{shukla}, \code{wricke}, \code{cv_stability},
#'   \code{superiority}, \code{integrated} and \code{concordance}
#'   (Kendall's W among the seven stability rankings).
#'
#' @references
#' Eberhart, S.A. and Russell, W.A. (1966). Stability parameters for
#' comparing varieties. \emph{Crop Science}, 6(1), 36-40.
#' \doi{10.2135/cropsci1966.0011183X000600010011x}
#'
#' Finlay, K.W. and Wilkinson, G.N. (1963). The analysis of adaptation in
#' a plant-breeding programme. \emph{Australian Journal of Agricultural
#' Research}, 14(6), 742-754. \doi{10.1071/AR9630742}
#'
#' Shukla, G.K. (1972). Some statistical aspects of partitioning
#' genotype-environmental components of variability. \emph{Heredity},
#' 29(2), 237-245. \doi{10.1038/hdy.1972.87}
#'
#' Wricke, G. (1962). Uber eine Methode zur Erfassung der okologischen
#' Streubreite in Feldversuchen. \emph{Zeitschrift fur Pflanzenzuchtung},
#' 47, 92-96.
#'
#' Gollob, H.F. (1968). A statistical model which combines features of
#' factor analytic and analysis of variance techniques.
#' \emph{Psychometrika}, 33(1), 73-115. \doi{10.1007/BF02289676}
#'
#' Purchase, J.L., Hatting, H. and van Deventer, C.S. (2000).
#' Genotype x environment interaction of winter wheat (Triticum aestivum
#' L.) in South Africa: II. Stability analysis of yield performance.
#' \emph{South African Journal of Plant and Soil}, 17(3), 101-107.
#' \doi{10.1080/02571862.2000.10634878}
#'
#' Francis, T.R. and Kannenberg, L.W. (1978). Yield stability studies in
#' short-season maize. I. A descriptive method for grouping genotypes.
#' \emph{Canadian Journal of Plant Science}, 58(4), 1029-1034.
#' \doi{10.4141/cjps78-157}
#'
#' Lin, C.S. and Binns, M.R. (1988). A superiority measure of cultivar
#' performance for cultivar x location data. \emph{Canadian Journal of
#' Plant Science}, 68(1), 193-198. \doi{10.4141/cjps88-018}
#'
#' @examples
#' data(pearlmillet_met)
#' st <- stability_analysis(pearlmillet_met, genotype = "genotype",
#'                          environment = "environment", replication = "rep",
#'                          trait = "grain_yield")
#' st$eberhart
#' st$ammi$anova
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram
#'   Choudhary, Hemant Poonia
stability_analysis <- function(data, genotype, environment, replication,
                               trait, method = "all", verbose = TRUE,
                               alpha = 0.05, fw_log = FALSE) {
  .check_alpha(alpha)
  all_m <- c("eberhart", "ammi", "finlay", "shukla", "wricke", "cv",
             "superiority")
  method <- tolower(method)
  if ("all" %in% method) method <- all_m
  bad <- setdiff(method, all_m)
  if (length(bad) > 0) {
    stop("Unknown method: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  facs <- c(genotype, environment, replication)
  d <- .prep_data(data, trait, facs, verbose)
  G <- d[[genotype]]
  Ev <- d[[environment]]
  g <- nlevels(G)
  e <- nlevels(Ev)
  if (g < 3) stop("Need at least 3 genotypes.", call. = FALSE)
  if (e < 3) stop("Need at least 3 environments.", call. = FALSE)
  cnt <- table(G, Ev)
  if (any(cnt == 0)) {
    stop("Every genotype must be tested in every environment (",
         sum(cnt == 0), " genotype x environment combination(s) missing).",
         call. = FALSE)
  }
  y <- d[[trait]]
  M <- tapply(y, list(G, Ev), mean)
  r <- if (is.null(replication)) 1 else stats::median(cnt)
  if (!is.null(replication) && any(cnt != r)) {
    warning("Unequal replication across genotype x environment cells; ",
            "the harmonic mean replication is used.", call. = FALSE)
    r <- 1 / mean(1 / cnt)
  }

  # Combined ANOVA and pooled error
  anova_tab <- NULL
  mse <- NA_real_
  df_e <- NA_real_
  if (!is.null(replication) && all(cnt >= 2)) {
    d$.rep <- interaction(Ev, d[[replication]], drop = TRUE)
    fit <- stats::aov(.formula(trait, paste(.bq(environment), "+ .rep +",
                                            .bq(genotype), "+",
                                            paste0(.bq(environment), ":",
                                                   .bq(genotype)))),
                      data = d)
    at <- stats::anova(fit)
    ms <- at[["Mean Sq"]]
    dfv <- at[["Df"]]
    mse <- ms[5]
    df_e <- dfv[5]
    fv <- c(ms[1] / ms[2], ms[2] / mse, ms[3] / mse, ms[4] / mse)
    pv <- stats::pf(fv, dfv[1:4], c(dfv[2], df_e, df_e, df_e),
                    lower.tail = FALSE)
    anova_tab <- data.frame(Df = dfv, `Sum Sq` = at[["Sum Sq"]],
                            `Mean Sq` = ms, `F value` = c(fv, NA),
                            `Pr(>F)` = c(pv, NA), check.names = FALSE)
    rownames(anova_tab) <- c("Environment", "Rep(Environment)", "Genotype",
                             "Genotype x Environment", "Pooled error")
  }
  err_mean <- if (is.na(mse)) NA_real_ else mse / r

  gm_i <- rowMeans(M)
  em_j <- colMeans(M)
  grand <- mean(M)
  I <- em_j - grand
  sI <- sum(I^2)
  Z <- M - outer(gm_i, rep(1, e)) - outer(rep(1, g), em_j) + grand
  W <- rowSums(Z^2)

  results <- list(
    trait = trait, anova = anova_tab, mse = mse, df_error = df_e, r = r,
    grand_mean = grand,
    genotype_means = data.frame(Genotype = rownames(M), Mean = gm_i,
                                row.names = NULL, stringsAsFactors = FALSE),
    environment_means = data.frame(Environment = colnames(M), Mean = em_j,
                                   Index = I, row.names = NULL,
                                   stringsAsFactors = FALSE),
    ge_matrix = M, alpha = alpha, methods = method)

  # Eberhart and Russell
  b <- as.numeric(M %*% I) / sI
  dev <- M - outer(gm_i, rep(1, e)) - outer(b, I)
  ssd <- rowSums(dev^2)
  msd <- ssd / (e - 2)
  if ("eberhart" %in% method) {
    se_b <- sqrt(msd / sI)
    t_b <- (b - 1) / se_b
    p_b <- 2 * stats::pt(-abs(t_b), e - 2)
    s2d <- if (is.na(err_mean)) msd else msd - err_mean
    p_s2d <- if (is.na(err_mean)) NA_real_ else
      stats::pf(msd / err_mean, e - 2, df_e, lower.tail = FALSE)
    sst <- rowSums((M - gm_i)^2)
    er <- data.frame(Genotype = rownames(M), Mean = gm_i, bi = b, SE_bi = se_b,
                     t_bi_eq_1 = t_b, p_bi = p_b, S2di = s2d, p_S2di = p_s2d,
                     R2 = 1 - ssd / sst, stringsAsFactors = FALSE)
    er$Response <- ifelse(er$p_bi < alpha, ifelse(er$bi > 1,
                                                  "above average (bi > 1)",
                                                  "below average (bi < 1)"),
                          "average (bi = 1)")
    er$Stable <- if (is.na(err_mean)) NA else er$p_S2di >= alpha
    er <- er[order(-er$Mean), ]
    rownames(er) <- NULL
    results$eberhart <- er
    # Eberhart-Russell ANOVA on the scale of means
    ss_g <- e * sum((gm_i - grand)^2)
    ss_eg <- sum((M - gm_i)^2)
    ss_el <- g * sI
    ss_gel <- sum(as.numeric(M %*% I)^2) / sI - ss_el
    ss_pd <- sum(ssd)
    df_v <- c(g * e - 1, g - 1, g * (e - 1), 1, g - 1, g * (e - 2))
    ss_v <- c(ss_g + ss_eg, ss_g, ss_eg, ss_el, ss_gel, ss_pd)
    ms_v <- ss_v / df_v
    ms_pd <- ms_v[6]
    f_v <- c(NA, ms_v[2] / ms_pd, NA, NA, ms_v[5] / ms_pd,
             if (is.na(err_mean)) NA else ms_pd / err_mean)
    p_v <- c(NA, stats::pf(f_v[2], df_v[2], df_v[6], lower.tail = FALSE), NA,
             NA, stats::pf(f_v[5], df_v[5], df_v[6], lower.tail = FALSE),
             if (is.na(err_mean)) NA else
               stats::pf(f_v[6], df_v[6], df_e, lower.tail = FALSE))
    era <- data.frame(Df = df_v, `Sum Sq` = ss_v, `Mean Sq` = ms_v,
                      `F value` = f_v, `Pr(>F)` = p_v, check.names = FALSE)
    rownames(era) <- c("Total", "Genotypes", "Env + (G x E)",
                       "Environment (linear)", "G x E (linear)",
                       "Pooled deviation")
    if (!is.na(err_mean)) {
      era <- rbind(era, data.frame(Df = df_e, `Sum Sq` = err_mean * df_e,
                                   `Mean Sq` = err_mean, `F value` = NA,
                                   `Pr(>F)` = NA, check.names = FALSE,
                                   row.names = "Pooled error (mean basis)"))
    }
    results$eberhart_anova <- era
  }

  # Finlay and Wilkinson
  if ("finlay" %in% method) {
    Mf <- M
    if (fw_log) {
      if (any(M <= 0)) stop("fw_log = TRUE needs positive means.",
                            call. = FALSE)
      Mf <- log10(M)
    }
    xm <- colMeans(Mf)
    xc <- xm - mean(xm)
    bf <- as.numeric(Mf %*% xc) / sum(xc^2)
    resf <- Mf - rowMeans(Mf) - outer(bf, xc)
    sef <- sqrt(rowSums(resf^2) / (e - 2) / sum(xc^2))
    fw <- data.frame(Genotype = rownames(M), Mean = gm_i, bi = bf, SE_bi = sef,
                     p_bi_eq_1 = 2 * stats::pt(-abs((bf - 1) / sef), e - 2),
                     R2 = 1 - rowSums(resf^2) / rowSums((Mf - rowMeans(Mf))^2),
                     stringsAsFactors = FALSE)
    fw <- fw[order(-fw$Mean), ]
    rownames(fw) <- NULL
    results$finlay <- fw
    results$finlay_scale <- if (fw_log) "log10" else "original"
  }

  # Shukla's stability variance
  sig2 <- (g * (g - 1) * W - sum(W)) / ((e - 1) * (g - 1) * (g - 2))
  if ("shukla" %in% method) {
    ps <- if (is.na(err_mean)) NA_real_ else
      stats::pf(sig2 / err_mean, e - 1, df_e, lower.tail = FALSE)
    sh <- data.frame(Genotype = rownames(M), Mean = gm_i, Sigma2_i = sig2,
                     F_value = if (is.na(err_mean)) NA_real_ else
                       sig2 / err_mean,
                     p_value = ps, Rank = rank(sig2),
                     stringsAsFactors = FALSE)
    sh <- sh[order(sh$Sigma2_i), ]
    rownames(sh) <- NULL
    results$shukla <- sh
  }

  # Wricke's ecovalence
  if ("wricke" %in% method) {
    wr <- data.frame(Genotype = rownames(M), Mean = gm_i, Wi = W,
                     Wi_percent = 100 * W / sum(W), Rank = rank(W),
                     stringsAsFactors = FALSE)
    wr <- wr[order(wr$Wi), ]
    rownames(wr) <- NULL
    results$wricke <- wr
  }

  # AMMI
  sv <- svd(Z)
  n_pc <- min(g - 1, e - 1)
  lam <- sv$d[seq_len(n_pc)]
  ss_pc <- r * lam^2
  gs <- sv$u[, seq_len(n_pc), drop = FALSE] %*% diag(sqrt(lam), n_pc)
  es <- sv$v[, seq_len(n_pc), drop = FALSE] %*% diag(sqrt(lam), n_pc)
  dimnames(gs) <- list(rownames(M), paste0("IPCA", seq_len(n_pc)))
  dimnames(es) <- list(colnames(M), paste0("IPCA", seq_len(n_pc)))
  asv <- if (n_pc >= 2 && ss_pc[2] > 0) {
    sqrt((ss_pc[1] / ss_pc[2] * gs[, 1])^2 + gs[, 2]^2)
  } else abs(gs[, 1])
  if ("ammi" %in% method) {
    df_pc <- g + e - 1 - 2 * seq_len(n_pc)
    keep <- df_pc > 0
    comp <- data.frame(IPCA = paste0("IPCA", seq_len(n_pc))[keep],
                       Eigenvalue = lam[keep]^2,
                       Variance_Percent = 100 * lam[keep]^2 / sum(lam^2),
                       Cumulative_Percent = cumsum(100 * lam^2 /
                                                     sum(lam^2))[keep],
                       Df = df_pc[keep], Sum_Sq = ss_pc[keep],
                       Mean_Sq = ss_pc[keep] / df_pc[keep],
                       stringsAsFactors = FALSE)
    if (!is.na(mse)) {
      comp$F_value <- comp$Mean_Sq / mse
      comp$p_value <- stats::pf(comp$F_value, comp$Df, df_e,
                                lower.tail = FALSE)
    }
    am_tab <- NULL
    if (!is.null(anova_tab)) {
      ss_ge <- anova_tab["Genotype x Environment", "Sum Sq"]
      df_ge <- anova_tab["Genotype x Environment", "Df"]
      rows <- anova_tab[1:4, ]
      ip <- data.frame(Df = comp$Df, `Sum Sq` = comp$Sum_Sq,
                       `Mean Sq` = comp$Mean_Sq, `F value` = comp$F_value,
                       `Pr(>F)` = comp$p_value, check.names = FALSE,
                       row.names = paste0("  ", comp$IPCA))
      res_df <- df_ge - sum(comp$Df)
      res_row <- NULL
      if (res_df > 0) {
        res_ss <- ss_ge - sum(comp$Sum_Sq)
        res_row <- data.frame(Df = res_df, `Sum Sq` = res_ss,
                              `Mean Sq` = res_ss / res_df,
                              `F value` = res_ss / res_df / mse,
                              `Pr(>F)` = stats::pf(res_ss / res_df / mse,
                                                   res_df, df_e,
                                                   lower.tail = FALSE),
                              check.names = FALSE, row.names = "  Residual")
      }
      am_tab <- rbind(rows, ip, res_row, anova_tab[5, ])
    }
    mean_rank <- rank(-gm_i)
    asv_rank <- rank(asv)
    stab <- data.frame(Genotype = rownames(M), Mean = gm_i,
                       Mean_Rank = mean_rank,
                       IPCA1 = gs[, 1],
                       IPCA2 = if (n_pc >= 2) gs[, 2] else NA_real_,
                       ASV = asv, ASV_Rank = asv_rank,
                       YSI = mean_rank + asv_rank, stringsAsFactors = FALSE)
    stab <- stab[order(stab$YSI, -stab$Mean), ]
    rownames(stab) <- NULL
    results$ammi <- list(anova = am_tab, components = comp, stability = stab,
                         genotype_scores = gs, environment_scores = es)
  }

  # Coefficient of variation
  cvv <- 100 * apply(M, 1, stats::sd) / gm_i
  if ("cv" %in% method) {
    cvs <- data.frame(Genotype = rownames(M), Mean = gm_i,
                      SD = apply(M, 1, stats::sd), CV = cvv,
                      CV_Rank = rank(cvv), stringsAsFactors = FALSE)
    cvs <- cvs[order(cvs$CV), ]
    rownames(cvs) <- NULL
    results$cv_stability <- cvs
  }

  # Lin and Binns superiority
  Mj <- apply(M, 2, max)
  Pi <- rowSums((M - outer(rep(1, g), Mj))^2) / (2 * e)
  if ("superiority" %in% method) {
    sp <- data.frame(Genotype = rownames(M), Mean = gm_i, Pi = Pi,
                     Pi_Rank = rank(Pi), stringsAsFactors = FALSE)
    sp <- sp[order(sp$Pi), ]
    rownames(sp) <- NULL
    results$superiority <- sp
  }

  # Integrated ranking over the seven methods
  if (setequal(method, all_m)) {
    integ <- data.frame(Genotype = rownames(M), Mean = gm_i,
                        ER = rank(msd), AMMI = rank(asv),
                        FW = rank(abs(b - 1)), Shukla = rank(sig2),
                        Wricke = rank(W), CV = rank(cvv), LB = rank(Pi),
                        stringsAsFactors = FALSE)
    rk <- as.matrix(integ[, 3:9])
    integ$Rank <- rank(rowMeans(rk), ties.method = "min")
    integ <- integ[order(integ$Rank, -integ$Mean), ]
    rownames(integ) <- NULL
    results$integrated <- integ
    m_r <- ncol(rk)
    rs <- rowSums(rk)
    S <- sum((rs - mean(rs))^2)
    Tj <- sum(apply(rk, 2, function(v) {
      tt <- table(v)
      sum(tt^3 - tt)
    }))
    den <- m_r^2 * (g^3 - g) - m_r * Tj
    Wk <- if (den > 0) 12 * S / den else NA_real_
    results$concordance <- list(kendall_w = Wk, rankings = m_r)
  }
  class(results) <- c("aridagri_stability", "list")
  if (verbose) print(results)
  invisible(results)
}


#' Print Method for Stability Analysis
#'
#' @param x An object of class \code{aridagri_stability}.
#' @param digits Number of decimal places. Default 3.
#' @param ... Unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.aridagri_stability <- function(x, digits = 3, ...) {
  .print_info("STABILITY ANALYSIS",
              c(Trait = x$trait %||% "",
                Genotypes = nrow(x$ge_matrix),
                Environments = ncol(x$ge_matrix),
                Replications = format(x$r, digits = 3),
                `Grand mean` = format(round(x$grand_mean, 2), nsmall = 2)))
  if (!is.null(x$anova)) .print_anova_df(x$anova, rownames(x$anova),
                                         title = "Combined analysis of variance")
  rnd <- function(df) {
    for (cl in names(df)) {
      if (is.numeric(df[[cl]])) df[[cl]] <- signif(df[[cl]], 4)
    }
    df
  }
  if (!is.null(x$eberhart)) {
    cat("\n Eberhart and Russell (1966)\n")
    print(rnd(x$eberhart[, c("Genotype", "Mean", "bi", "SE_bi", "p_bi",
                             "S2di", "p_S2di", "R2")]), row.names = FALSE)
    cat(" bi = 1 and non-significant S2di indicate average, stable response.\n")
  }
  if (!is.null(x$ammi)) {
    if (!is.null(x$ammi$anova)) {
      .print_anova_df(x$ammi$anova, rownames(x$ammi$anova), total = FALSE,
                      title = "AMMI analysis of variance")
    }
    cat("\n AMMI stability value (ASV) and yield stability index (YSI)\n")
    print(rnd(x$ammi$stability), row.names = FALSE)
  }
  if (!is.null(x$shukla)) {
    cat("\n Shukla (1972) stability variance\n")
    print(rnd(x$shukla), row.names = FALSE)
  }
  if (!is.null(x$integrated)) {
    cat("\n Integrated stability ranking (1 = most stable)\n")
    print(rnd(x$integrated), row.names = FALSE)
    if (!is.null(x$concordance$kendall_w)) {
      cat(sprintf(" Agreement among the seven rankings: Kendall's W = %.3f\n",
                  x$concordance$kendall_w))
    }
  }
  cat("\n Full tables: $eberhart, $eberhart_anova, $finlay, $shukla, $wricke,",
      "$ammi, $cv_stability, $superiority\n")
  invisible(x)
}
