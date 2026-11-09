# =============================================================================
# NUTRIENT USE EFFICIENCY, RESPONSE CURVES AND ECONOMICS
# Package: aridagri
# =============================================================================

#' Nutrient Use Efficiency Indices
#'
#' @description
#' Agronomic efficiency (AE), partial factor productivity (PFP), apparent
#' recovery efficiency (ARE), physiological efficiency (PE), internal
#' utilization efficiency (IE) and, when biomass is given, the
#' biomass-based physiological efficiency and utilization efficiency of
#' Fageria and Baligar (2005). All arguments may be vectors (one value per
#' treatment).
#'
#' @param yield_fertilized Grain yield with the nutrient (kg/ha).
#' @param yield_control Grain yield of the control without the nutrient
#'   (kg/ha).
#' @param nutrient_applied Nutrient applied (kg/ha), greater than 0.
#' @param nutrient_uptake_fert Total nutrient uptake with the nutrient
#'   (kg/ha), optional.
#' @param nutrient_uptake_ctrl Total nutrient uptake of the control
#'   (kg/ha), optional.
#' @param biomass_fert Total biomass with the nutrient (kg/ha), optional.
#' @param biomass_ctrl Total biomass of the control (kg/ha), optional.
#' @param verbose Logical. If TRUE (default), prints the table.
#'
#' @details
#' \itemize{
#'   \item AE = (Yf - Yc) / N (kg grain per kg nutrient)
#'   \item PFP = Yf / N
#'   \item ARE = 100 (Uf - Uc) / N (percent)
#'   \item PE = (Yf - Yc) / (Uf - Uc) (kg grain per kg nutrient taken up;
#'     called agro-physiological efficiency by Fageria and Baligar)
#'   \item IE = Yf / Uf
#'   \item PE_biomass = (Bf - Bc) / (Uf - Uc) and UE = PE_biomass x ARE /
#'     100 (Fageria and Baligar, 2005)
#' }
#' AE = PE x ARE / 100 holds by definition. PE is not computed when the
#' uptake difference is zero or negative.
#'
#' @return A data frame of class \code{aridagri_nue}.
#'
#' @references
#' Dobermann, A. (2005). Nitrogen use efficiency - state of the art. IFA
#' International Workshop on Enhanced-Efficiency Fertilizers, Frankfurt.
#'
#' Fageria, N.K. and Baligar, V.C. (2005). Enhancing nitrogen use
#' efficiency in crop plants. \emph{Advances in Agronomy}, 88, 97-185.
#' \doi{10.1016/S0065-2113(05)88004-6}
#'
#' @examples
#' nue_calculate(yield_fertilized = 1850, yield_control = 1200,
#'               nutrient_applied = 40)
#' # Several N doses with uptake and biomass
#' nue_calculate(yield_fertilized = c(1450, 1720, 1850),
#'               yield_control = 1200, nutrient_applied = c(20, 40, 60),
#'               nutrient_uptake_fert = c(48, 61, 68),
#'               nutrient_uptake_ctrl = 35,
#'               biomass_fert = c(4300, 4950, 5300), biomass_ctrl = 3600)
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
nue_calculate <- function(yield_fertilized, yield_control, nutrient_applied,
                          nutrient_uptake_fert = NULL,
                          nutrient_uptake_ctrl = NULL, biomass_fert = NULL,
                          biomass_ctrl = NULL, verbose = TRUE) {
  .check_vectors(yield_fertilized = yield_fertilized,
                 yield_control = yield_control,
                 nutrient_applied = nutrient_applied)
  if (any(nutrient_applied <= 0)) {
    stop("nutrient_applied must be greater than 0", call. = FALSE)
  }
  n <- max(length(yield_fertilized), length(yield_control),
           length(nutrient_applied))
  yf <- rep_len(yield_fertilized, n)
  yc <- rep_len(yield_control, n)
  N <- rep_len(nutrient_applied, n)
  res <- data.frame(Yield_Fertilized = yf, Yield_Control = yc,
                    Nutrient_Applied = N, Yield_Response = yf - yc,
                    AE_kg_kg = (yf - yc) / N, PFP_kg_kg = yf / N)
  has_u <- !is.null(nutrient_uptake_fert) && !is.null(nutrient_uptake_ctrl)
  if (has_u) {
    .check_vectors(nutrient_uptake_fert = nutrient_uptake_fert,
                   nutrient_uptake_ctrl = nutrient_uptake_ctrl)
    uf <- rep_len(nutrient_uptake_fert, n)
    uc <- rep_len(nutrient_uptake_ctrl, n)
    du <- uf - uc
    res$ARE_percent <- 100 * du / N
    res$PE_kg_kg <- ifelse(du > 0, (yf - yc) / du, NA_real_)
    res$IE_kg_kg <- ifelse(uf > 0, yf / uf, NA_real_)
    if (!is.null(biomass_fert) && !is.null(biomass_ctrl)) {
      .check_vectors(biomass_fert = biomass_fert, biomass_ctrl = biomass_ctrl)
      bf <- rep_len(biomass_fert, n)
      bc <- rep_len(biomass_ctrl, n)
      res$PE_biomass_kg_kg <- ifelse(du > 0, (bf - bc) / du, NA_real_)
      res$UE_kg_kg <- res$PE_biomass_kg_kg * res$ARE_percent / 100
    }
  }
  if (!is.null(biomass_fert)) {
    res$HI_Fertilized <- 100 * yf / rep_len(biomass_fert, n)
  }
  if (!is.null(biomass_ctrl)) {
    res$HI_Control <- 100 * yc / rep_len(biomass_ctrl, n)
  }
  if (verbose) {
    .heading("NUTRIENT USE EFFICIENCY")
    o <- res
    for (cl in names(o)) o[[cl]] <- round(o[[cl]], 2)
    print(o, row.names = FALSE)
    cat("\n AE, PFP, PE, IE: kg grain per kg nutrient; ARE: percent of",
        "applied nutrient recovered\n")
  }
  class(res) <- c("aridagri_nue", "data.frame")
  invisible(res)
}


# --- Response models ----------------------------------------------------------

.fit_quadratic <- function(x, y, pr) {
  fit <- stats::lm(y ~ x + I(x^2))
  cf <- stats::coef(fit)
  a <- cf[[1]]
  b <- cf[[2]]
  c2 <- cf[[3]]
  concave <- is.finite(c2) && c2 < 0
  xmax <- if (concave) -b / (2 * c2) else NA_real_
  xeco <- if (concave) min(max((pr - b) / (2 * c2), 0), xmax) else NA_real_
  f <- function(z) a + b * z + c2 * z^2
  list(model = "quadratic", fit = fit, coefficients = c(a = a, b = b, c = c2),
       predict = f, dose_max_yield = xmax,
       max_yield = if (concave) f(xmax) else NA_real_,
       dose_economic = xeco,
       yield_economic = if (concave) f(xeco) else NA_real_,
       rss = sum(stats::residuals(fit)^2), k = 3,
       note = if (!concave) "Quadratic term is not negative: no maximum." else
         "")
}

# Profile least squares over the join point x0.
.profile_x0 <- function(x, y, make_X) {
  u <- sort(unique(x))
  rss_at <- function(x0) {
    X <- make_X(x0)
    f <- stats::lm.fit(cbind(1, X), y)
    sum(f$residuals^2)
  }
  grid <- seq(u[2], u[length(u)], length.out = 200)
  rs <- vapply(grid, rss_at, numeric(1))
  i <- which.min(rs)
  lo <- grid[max(1, i - 1)]
  hi <- grid[min(length(grid), i + 1)]
  opt <- stats::optimize(rss_at, c(lo, hi))
  if (opt$objective <= rs[i]) opt$minimum else grid[i]
}

.fit_linear_plateau <- function(x, y, pr) {
  if (length(unique(x)) < 3) stop("Need at least 3 distinct doses.",
                                  call. = FALSE)
  x0 <- .profile_x0(x, y, function(x0) pmin(x, x0))
  xx <- pmin(x, x0)
  fit <- stats::lm(y ~ xx)
  a <- stats::coef(fit)[[1]]
  b <- stats::coef(fit)[[2]]
  f <- function(z) a + b * pmin(z, x0)
  xeco <- if (b > pr) x0 else 0
  list(model = "linear_plateau", fit = fit,
       coefficients = c(a = a, b = b, join_point = x0, plateau = f(x0)),
       predict = f, dose_max_yield = x0, max_yield = f(x0),
       dose_economic = xeco, yield_economic = f(xeco),
       rss = sum(stats::residuals(fit)^2), k = 3,
       note = if (x0 >= max(x) - 1e-8) "Join point at the highest dose: no plateau reached." else "")
}

.fit_quadratic_plateau <- function(x, y, pr) {
  if (length(unique(x)) < 3) stop("Need at least 3 distinct doses.",
                                  call. = FALSE)
  x0 <- .profile_x0(x, y, function(x0) (pmin(x, x0) - x0)^2)
  q <- (pmin(x, x0) - x0)^2
  fit <- stats::lm(y ~ q)
  P <- stats::coef(fit)[[1]]
  c2 <- stats::coef(fit)[[2]]
  a <- P + c2 * x0^2
  b <- -2 * c2 * x0
  f <- function(z) P + c2 * (pmin(z, x0) - x0)^2
  concave <- c2 < 0
  xeco <- if (concave) min(max((pr - b) / (2 * c2), 0), x0) else NA_real_
  list(model = "quadratic_plateau", fit = fit,
       coefficients = c(a = a, b = b, c = c2, join_point = x0, plateau = P),
       predict = f, dose_max_yield = x0, max_yield = P,
       dose_economic = xeco,
       yield_economic = if (concave) f(xeco) else NA_real_,
       rss = sum(stats::residuals(fit)^2), k = 3,
       note = if (!concave) "Curvature is not negative." else "")
}

.fit_mitscherlich <- function(x, y, pr) {
  if (length(unique(x)) < 4) stop("Need at least 4 distinct doses.",
                                  call. = FALSE)
  fit <- tryCatch(stats::nls(y ~ SSasymp(x, Asym, R0, lrc)),
                  error = function(e) NULL)
  if (is.null(fit)) stop("The Mitscherlich model did not converge.",
                         call. = FALSE)
  cf <- stats::coef(fit)
  A <- cf[["Asym"]]
  B <- A - cf[["R0"]]
  cc <- exp(cf[["lrc"]])
  f <- function(z) A - B * exp(-cc * z)
  xeco <- if (B * cc > pr) log(B * cc / pr) / cc else 0
  x95 <- if (B > 0.05 * A) log(B / (0.05 * A)) / cc else 0
  list(model = "mitscherlich", fit = fit,
       coefficients = c(A = A, B = B, c = cc), predict = f,
       dose_max_yield = NA_real_, max_yield = A, dose_95_percent = x95,
       dose_economic = xeco, yield_economic = f(xeco),
       rss = sum(stats::residuals(fit)^2), k = 3,
       note = "max_yield is the asymptote A.")
}

.fit_stats <- function(m, y, dose) {
  n <- length(y)
  ext <- c(m$dose_economic, m$dose_max_yield)
  if (any(is.finite(ext) & ext > max(dose) + 1e-8)) {
    m$note <- trimws(paste(m$note, "Optimum lies beyond the highest dose",
                           "tested; treat it as an extrapolation."))
  }
  tss <- sum((y - mean(y))^2)
  r2 <- 1 - m$rss / tss
  p <- m$k
  m$r_squared <- r2
  m$adj_r_squared <- 1 - (1 - r2) * (n - 1) / max(n - p, 1)
  m$rmse <- sqrt(m$rss / n)
  m$aic <- n * log(m$rss / n) + 2 * (p + 1)
  m
}


#' Nutrient Response Curves and Economic Optimum Dose
#'
#' @description
#' Fits nutrient response models and finds the dose for maximum yield and
#' the economic optimum dose, where the marginal return equals the price
#' ratio. Quadratic, linear-plateau, quadratic-plateau and Mitscherlich
#' models are available; \code{model = "all"} fits all four and compares
#' them.
#'
#' @param dose Numeric vector of nutrient doses (kg/ha).
#' @param yield Numeric vector of yields (kg/ha), one per dose (replicated
#'   doses are allowed).
#' @param model \code{"quadratic"} (default), \code{"linear_plateau"},
#'   \code{"quadratic_plateau"}, \code{"mitscherlich"} or \code{"all"}.
#' @param price_output Price of the produce (Rs/kg).
#' @param price_nutrient Price of the nutrient (Rs/kg of nutrient).
#' @param verbose Logical. If TRUE (default), prints the results.
#'
#' @details
#' \itemize{
#'   \item Quadratic: \eqn{Y = a + bX + cX^2}; maximum at \eqn{-b/(2c)};
#'     economic optimum \eqn{(P_n/P_y - b)/(2c)}.
#'   \item Linear-plateau: \eqn{Y = a + b \min(X, X_0)}; economic optimum
#'     at the join point \eqn{X_0} when \eqn{b > P_n/P_y}.
#'   \item Quadratic-plateau: quadratic up to its vertex \eqn{X_0}, then
#'     flat.
#'   \item Mitscherlich: \eqn{Y = A - B e^{-cX}}; economic optimum
#'     \eqn{\ln(B c P_y / P_n)/c}; also the dose giving 95\% of the
#'     asymptote.
#' }
#' Plateau models are fitted by profile least squares over the join
#' point; the Mitscherlich model by \code{nls} with a self-starting model.
#' All four models have three parameters (the join point counts as one), so
#' they are compared by R-squared, adjusted R-squared, RMSE and AIC on an
#' equal footing (AIC differences equal those of \code{stats::AIC()}). An
#' optimum beyond the highest dose tested is flagged as an extrapolation.
#'
#' @return For a single model, an object of class \code{aridagri_response}
#'   with \code{model}, \code{coefficients}, \code{r_squared},
#'   \code{adj_r_squared}, \code{rmse}, \code{aic}, \code{dose_max_yield},
#'   \code{max_yield}, \code{dose_economic}, \code{yield_economic} and
#'   \code{fit}. For \code{model = "all"}, a list with \code{comparison}
#'   (one row per model), \code{models} and \code{best} (lowest AIC).
#'
#' @references
#' Cerrato, M.E. and Blackmer, A.M. (1990). Comparison of models for
#' describing corn yield response to nitrogen fertilizer. \emph{Agronomy
#' Journal}, 82(1), 138-143.
#' \doi{10.2134/agronj1990.00021962008200010030x}
#'
#' @examples
#' dose <- c(0, 20, 40, 60, 80, 100)
#' yield <- c(1100, 1350, 1520, 1610, 1650, 1660)
#' nutrient_response(dose, yield, model = "quadratic",
#'                   price_output = 60, price_nutrient = 15)
#' nutrient_response(dose, yield, model = "all", price_output = 60,
#'                   price_nutrient = 15)
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
nutrient_response <- function(dose, yield, model = "quadratic",
                              price_output = 50, price_nutrient = 15,
                              verbose = TRUE) {
  .check_vectors(dose = dose, yield = yield)
  if (length(dose) != length(yield)) {
    stop("dose and yield must have same length", call. = FALSE)
  }
  if (length(unique(dose)) < 3) stop("Need at least 3 distinct doses.",
                                     call. = FALSE)
  if (price_output <= 0 || price_nutrient <= 0) {
    stop("Prices must be positive.", call. = FALSE)
  }
  model <- match.arg(model, c("quadratic", "linear_plateau",
                              "quadratic_plateau", "mitscherlich", "all"))
  pr <- price_nutrient / price_output
  fitters <- list(quadratic = .fit_quadratic,
                  linear_plateau = .fit_linear_plateau,
                  quadratic_plateau = .fit_quadratic_plateau,
                  mitscherlich = .fit_mitscherlich)
  run <- if (model == "all") names(fitters) else model
  fits <- list()
  for (m in run) {
    f <- tryCatch(fitters[[m]](dose, yield, pr), error = function(e) {
      if (model != "all") stop(conditionMessage(e), call. = FALSE)
      NULL
    })
    if (!is.null(f)) fits[[m]] <- .fit_stats(f, yield, dose)
  }
  if (model != "all") {
    res <- fits[[1]]
    res$price_ratio <- pr
    class(res) <- "aridagri_response"
    if (verbose) print(res)
    return(invisible(res))
  }
  cmp <- do.call(rbind, lapply(fits, function(m) {
    data.frame(Model = m$model, R2 = m$r_squared, Adj_R2 = m$adj_r_squared,
               RMSE = m$rmse, AIC = m$aic, Dose_max = m$dose_max_yield,
               Max_yield = m$max_yield, Dose_economic = m$dose_economic,
               Yield_economic = m$yield_economic, stringsAsFactors = FALSE)
  }))
  rownames(cmp) <- NULL
  best <- cmp$Model[which.min(cmp$AIC)]
  res <- list(comparison = cmp, models = fits, best = best, price_ratio = pr)
  class(res) <- "aridagri_response_set"
  if (verbose) print(res)
  invisible(res)
}

#' Print Methods for Nutrient Response Results
#'
#' @param x An object of class \code{aridagri_response} or
#'   \code{aridagri_response_set}.
#' @param ... Unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @name print.aridagri_response
#' @export
print.aridagri_response <- function(x, ...) {
  .heading(paste("NUTRIENT RESPONSE:", toupper(gsub("_", " ", x$model))))
  cf <- x$coefficients
  cat(" Coefficients:", paste(names(cf), signif(cf, 5), sep = " = ",
                              collapse = ", "), "\n")
  cat(sprintf(" R-squared = %.4f   Adjusted R-squared = %.4f   RMSE = %.2f\n",
              x$r_squared, x$adj_r_squared, x$rmse))
  if (!is.na(x$dose_max_yield)) {
    cat(sprintf(" Dose for maximum yield = %.1f kg/ha (yield %.0f)\n",
                x$dose_max_yield, x$max_yield))
  }
  if (!is.null(x$dose_95_percent)) {
    cat(sprintf(" Asymptote = %.0f; dose for 95%% of asymptote = %.1f kg/ha\n",
                x$max_yield, x$dose_95_percent))
  }
  if (!is.na(x$dose_economic)) {
    cat(sprintf(" Economic optimum dose = %.1f kg/ha (yield %.0f; price ratio %.3f)\n",
                x$dose_economic, x$yield_economic, x$price_ratio))
  }
  if (nzchar(x$note)) cat(" Note:", x$note, "\n")
  invisible(x)
}

#' @rdname print.aridagri_response
#' @export
print.aridagri_response_set <- function(x, ...) {
  .heading("NUTRIENT RESPONSE: MODEL COMPARISON")
  o <- x$comparison
  for (cl in names(o)[-1]) o[[cl]] <- round(o[[cl]], 3)
  print(o, row.names = FALSE)
  cat("\n Best model by AIC:", x$best, "\n")
  invisible(x)
}


#' Economic Analysis of a Crop
#'
#' @description
#' Cost of cultivation, gross and net returns, benefit-cost ratio and cost
#' of production. Arguments may be vectors (one value per treatment).
#'
#' @param yield Crop yield (kg/ha).
#' @param price Market price (Rs/kg).
#' @param cost_fixed Fixed costs (Rs/ha).
#' @param cost_variable Variable costs (Rs/ha).
#' @param byproduct_yield By-product yield (kg/ha), default 0.
#' @param byproduct_price By-product price (Rs/kg), default 0.
#' @param verbose Logical. If TRUE (default), prints the table.
#'
#' @return A data frame of class \code{aridagri_economics}.
#'
#' @examples
#' economic_analysis(yield = 1200, price = 65, cost_fixed = 15000,
#'                   cost_variable = 12000)
#' economic_analysis(yield = c(1050, 1200, 1340), price = 65,
#'                   cost_fixed = 15000, cost_variable = c(11000, 12000, 13500),
#'                   byproduct_yield = c(1600, 1800, 1950), byproduct_price = 5)
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
economic_analysis <- function(yield, price, cost_fixed, cost_variable,
                              byproduct_yield = 0, byproduct_price = 0,
                              verbose = TRUE) {
  .check_vectors(yield = yield, price = price, cost_fixed = cost_fixed,
                 cost_variable = cost_variable,
                 byproduct_yield = byproduct_yield,
                 byproduct_price = byproduct_price)
  total_cost <- cost_fixed + cost_variable
  if (any(total_cost <= 0)) stop("Total cost must be positive.", call. = FALSE)
  gross_main <- yield * price
  gross_by <- byproduct_yield * byproduct_price
  gross <- gross_main + gross_by
  net <- gross - total_cost
  res <- data.frame(Yield_kg_ha = yield, Price_Rs_kg = price,
                    Cost_Fixed = cost_fixed, Cost_Variable = cost_variable,
                    Total_Cost = total_cost, Gross_Return = gross,
                    Net_Return = net, BC_Ratio = gross / total_cost,
                    Return_per_Rupee = net / total_cost,
                    Cost_of_Production = ifelse(yield > 0, total_cost / yield,
                                                NA_real_))
  if (verbose) {
    .heading("ECONOMIC ANALYSIS")
    o <- res
    o$BC_Ratio <- round(o$BC_Ratio, 2)
    o$Return_per_Rupee <- round(o$Return_per_Rupee, 2)
    o$Cost_of_Production <- round(o$Cost_of_Production, 2)
    print(o, row.names = FALSE)
    cat("\n B:C ratio = gross return / total cost; cost of production in Rs/kg\n")
  }
  class(res) <- c("aridagri_economics", "data.frame")
  invisible(res)
}
