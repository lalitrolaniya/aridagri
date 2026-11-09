# =============================================================================
# AGRONOMIC INDICES: THERMAL TIME, CROP GROWTH, HARVEST INDEX, YIELD GAPS
# AND ECONOMICS
# Package: aridagri
# =============================================================================

# Stop unless all arguments are numeric vectors of the same length without
# missing values.
.check_vectors <- function(..., allow_na = FALSE) {
  args <- list(...)
  nm <- names(args)
  for (i in seq_along(args)) {
    if (!is.numeric(args[[i]])) {
      stop("'", nm[i], "' must be numeric.", call. = FALSE)
    }
    if (!allow_na && anyNA(args[[i]])) {
      stop("'", nm[i], "' contains missing values.", call. = FALSE)
    }
  }
  len <- vapply(args, length, integer(1))
  if (length(unique(len[len > 1])) > 1) {
    stop("Arguments must have the same length (or length 1): ",
         paste(nm, collapse = ", "), ".", call. = FALSE)
  }
  invisible(TRUE)
}

#' Day Length from Latitude and Date
#'
#' @description
#' Maximum possible daylight hours (astronomical day length) from latitude
#' and day of year, following FAO Irrigation and Drainage Paper 56
#' (equations 24, 25 and 34). Useful for photo-thermal units when measured
#' day length is not available.
#'
#' @param latitude Latitude in decimal degrees (north positive).
#' @param date A \code{Date} vector, or day of year (1-366) as numbers.
#'
#' @return Numeric vector of day length in hours.
#'
#' @references
#' Allen, R.G., Pereira, L.S., Raes, D. and Smith, M. (1998). \emph{Crop
#' Evapotranspiration: Guidelines for Computing Crop Water Requirements}.
#' FAO Irrigation and Drainage Paper 56. FAO, Rome.
#'
#' @examples
#' # Bikaner (28.0 N) on 15 March and 15 June
#' day_length(28.0, as.Date(c("2026-03-15", "2026-06-15")))
#'
#' @export
day_length <- function(latitude, date) {
  if (!is.numeric(latitude) || any(abs(latitude) > 90, na.rm = TRUE)) {
    stop("'latitude' must be in decimal degrees between -90 and 90.",
         call. = FALSE)
  }
  doy <- if (inherits(date, "Date")) as.numeric(format(date, "%j")) else date
  if (!is.numeric(doy) || any(doy < 1 | doy > 366, na.rm = TRUE)) {
    stop("'date' must be a Date or a day of year between 1 and 366.",
         call. = FALSE)
  }
  phi <- latitude * pi / 180
  delta <- 0.409 * sin(2 * pi * doy / 365 - 1.39)
  x <- pmin(pmax(-tan(phi) * tan(delta), -1), 1)
  24 / pi * acos(x)
}


#' Growing Degree Days and Thermal Indices
#'
#' @description
#' Growing degree days (GDD), helio-thermal units (HTU), photo-thermal
#' units (PTU), heat use efficiency (HUE) and helio-thermal use efficiency
#' from daily temperatures.
#'
#' @param tmax Daily maximum temperature (degrees C).
#' @param tmin Daily minimum temperature (degrees C).
#' @param base_temp Base temperature (degrees C). Default 10.
#' @param sunshine_hours Daily bright sunshine hours (optional, for HTU).
#' @param day_length Daily day length in hours (optional, for PTU); see
#'   \code{\link{day_length}}.
#' @param crop_yield Economic yield (kg/ha) for HUE (optional).
#' @param biomass Total biomass (kg/ha) for HUE (optional).
#' @param verbose Logical. If TRUE (default), prints the results.
#' @param upper_temp Optional upper threshold temperature (degrees C)
#'   above which development does not increase.
#' @param method \code{"average"} (default): GDD = max(0, (Tmax + Tmin)/2 -
#'   Tbase), with the mean capped at \code{upper_temp}. \code{"modified"}:
#'   Tmax and Tmin are each limited to the range from \code{base_temp} to
#'   \code{upper_temp} before averaging (McMaster and Wilhelm, 1997,
#'   method 2).
#'
#' @details
#' HTU = sum of GDD x bright sunshine hours; PTU = sum of GDD x day length
#' (Rajput, 1980). HUE = yield / GDD (kg/ha per degree-day); HTUE = yield /
#' HTU.
#'
#' @return An object of class \code{aridagri_thermal}: a list with
#'   \code{summary} (index, value, unit), \code{gdd_total} and
#'   \code{daily} (daily and cumulative GDD, HTU and PTU).
#'
#' @references
#' McMaster, G.S. and Wilhelm, W.W. (1997). Growing degree-days: one
#' equation, two interpretations. \emph{Agricultural and Forest
#' Meteorology}, 87(4), 291-300. \doi{10.1016/S0168-1923(97)00027-0}
#'
#' Rajput, R.P. (1980). Response of soybean crop to climatic and soil
#' environments. PhD thesis, IARI, New Delhi.
#'
#' @examples
#' tmax <- c(31, 33, 35, 36, 38, 40, 41, 39, 37, 36)
#' tmin <- c(16, 17, 19, 20, 22, 24, 25, 24, 22, 21)
#' thermal_indices(tmax, tmin, base_temp = 10, upper_temp = 40,
#'                 sunshine_hours = rep(9.5, 10),
#'                 day_length = day_length(28, 70:79), crop_yield = 1250)
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram
#'   Choudhary, Hemant Poonia
thermal_indices <- function(tmax, tmin, base_temp = 10, sunshine_hours = NULL,
                            day_length = NULL, crop_yield = NULL,
                            biomass = NULL, verbose = TRUE, upper_temp = NULL,
                            method = c("average", "modified")) {
  method <- match.arg(method)
  .check_vectors(tmax = tmax, tmin = tmin)
  if (length(tmax) != length(tmin) || length(tmax) < 1) {
    stop("'tmax' and 'tmin' must have the same, non-zero length.",
         call. = FALSE)
  }
  if (any(tmax < tmin)) {
    stop("'tmax' is below 'tmin' on day(s): ",
         paste(utils::head(which(tmax < tmin), 10), collapse = ", "), ".",
         call. = FALSE)
  }
  if (!is.null(upper_temp) && upper_temp <= base_temp) {
    stop("'upper_temp' must be above 'base_temp'.", call. = FALSE)
  }
  n_days <- length(tmax)
  if (method == "average") {
    tmean <- (tmax + tmin) / 2
    if (!is.null(upper_temp)) tmean <- pmin(tmean, upper_temp)
    gdd <- pmax(0, tmean - base_temp)
  } else {
    up <- if (is.null(upper_temp)) Inf else upper_temp
    tx <- pmin(pmax(tmax, base_temp), up)
    tn <- pmin(pmax(tmin, base_temp), up)
    tmean <- (tx + tn) / 2
    gdd <- tmean - base_temp
  }
  daily <- data.frame(Day = seq_len(n_days), Tmax = tmax, Tmin = tmin,
                      Tmean = (tmax + tmin) / 2, GDD = gdd,
                      Cum_GDD = cumsum(gdd))
  out <- data.frame(Index = c("Growing degree days (GDD)", "Mean daily GDD"),
                    Value = c(sum(gdd), mean(gdd)),
                    Unit = c("degree C day", "degree C day per day"),
                    stringsAsFactors = FALSE)
  htu_total <- NULL
  if (!is.null(sunshine_hours)) {
    .check_vectors(sunshine_hours = sunshine_hours)
    if (length(sunshine_hours) != n_days) {
      stop("'sunshine_hours' must have one value per day.", call. = FALSE)
    }
    daily$HTU <- gdd * sunshine_hours
    daily$Cum_HTU <- cumsum(daily$HTU)
    htu_total <- sum(daily$HTU)
    out <- rbind(out, data.frame(Index = "Helio-thermal units (HTU)",
                                 Value = htu_total,
                                 Unit = "degree C day hour"))
  }
  if (!is.null(day_length)) {
    .check_vectors(day_length = day_length)
    if (length(day_length) != n_days) {
      stop("'day_length' must have one value per day.", call. = FALSE)
    }
    daily$PTU <- gdd * day_length
    daily$Cum_PTU <- cumsum(daily$PTU)
    out <- rbind(out, data.frame(Index = "Photo-thermal units (PTU)",
                                 Value = sum(daily$PTU),
                                 Unit = "degree C day hour"))
  }
  gt <- sum(gdd)
  if (!is.null(crop_yield)) {
    out <- rbind(out, data.frame(Index = "Heat use efficiency (grain)",
                                 Value = if (gt > 0) crop_yield / gt else NA,
                                 Unit = "kg/ha per degree C day"))
    if (!is.null(htu_total)) {
      out <- rbind(out, data.frame(
        Index = "Helio-thermal use efficiency (grain)",
        Value = if (htu_total > 0) crop_yield / htu_total else NA,
        Unit = "kg/ha per degree C day hour"))
    }
  }
  if (!is.null(biomass)) {
    out <- rbind(out, data.frame(Index = "Heat use efficiency (biomass)",
                                 Value = if (gt > 0) biomass / gt else NA,
                                 Unit = "kg/ha per degree C day"))
  }
  if (verbose) {
    .heading("THERMAL INDICES")
    cat(sprintf(" Days = %d   Base temperature = %.1f C%s   Mean Tmax = %.1f C   Mean Tmin = %.1f C\n\n",
                n_days, base_temp,
                if (!is.null(upper_temp)) sprintf("   Upper = %.1f C", upper_temp) else "",
                mean(tmax), mean(tmin)))
    o <- out
    o$Value <- round(o$Value, 3)
    print(o, row.names = FALSE, right = FALSE)
  }
  output <- list(summary = out, gdd_total = gt, daily = daily,
                 base_temp = base_temp, upper_temp = upper_temp,
                 method = method)
  class(output) <- "aridagri_thermal"
  invisible(output)
}


#' Crop Growth Analysis (CGR, RGR, NAR, LAI, LAD)
#'
#' @description
#' Classical growth analysis from sequential harvests: crop growth rate,
#' relative growth rate, net assimilation rate, leaf area index and leaf
#' area duration for each interval between samplings.
#'
#' @param dry_weight Dry matter at each sampling (g per sample).
#' @param leaf_area Leaf area at each sampling (cm2 per sample), or
#'   \code{NULL} to compute only CGR and RGR.
#' @param days Days after sowing at each sampling (strictly increasing).
#' @param verbose Logical. If TRUE (default), prints the table.
#' @param ground_area Land area occupied by each sample (m2). When given,
#'   CGR is expressed in g m-2 day-1 and LAI = leaf area / ground area.
#'   When \code{NULL} (default), CGR is in g per sample per day and LAI is
#'   leaf area / 10000 (that is, leaf area in cm2 per m2 of ground).
#'
#' @details
#' For two samplings at \eqn{t_1} and \eqn{t_2} (Watson, 1952;
#' Radford, 1967):
#' \itemize{
#'   \item CGR = \eqn{(W_2 - W_1) / [P (t_2 - t_1)]}, P = ground area
#'   \item RGR = \eqn{(\ln W_2 - \ln W_1) / (t_2 - t_1)} (g g-1 day-1)
#'   \item NAR = \eqn{\frac{W_2 - W_1}{t_2 - t_1} \times
#'     \frac{\ln L_2 - \ln L_1}{L_2 - L_1}}, reported in g m-2 (leaf)
#'     day-1
#'   \item LAD = \eqn{(LAI_1 + LAI_2)(t_2 - t_1) / 2} (days)
#' }
#'
#' @return A data frame of class \code{aridagri_growth} with one row per
#'   interval: \code{Period}, \code{CGR}, \code{RGR} and, with leaf area,
#'   \code{NAR}, \code{LAI_start}, \code{LAI_end} and \code{LAD}.
#'
#' @references
#' Radford, P.J. (1967). Growth analysis formulae: their use and abuse.
#' \emph{Crop Science}, 7(3), 171-175.
#' \doi{10.2135/cropsci1967.0011183X000700030001x}
#'
#' Watson, D.J. (1952). The physiological basis of variation in yield.
#' \emph{Advances in Agronomy}, 4, 101-145.
#' \doi{10.1016/S0065-2113(08)60307-7}
#'
#' @examples
#' dry_weight <- c(0.5, 2.1, 8.5, 25, 45, 62, 75)
#' leaf_area <- c(15, 85, 350, 800, 950, 850, 600)
#' days <- c(15, 30, 45, 60, 75, 90, 105)
#' crop_growth_analysis(dry_weight, leaf_area, days, ground_area = 0.09)
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram
#'   Choudhary, Hemant Poonia
crop_growth_analysis <- function(dry_weight, leaf_area = NULL, days,
                                 verbose = TRUE, ground_area = NULL) {
  .check_vectors(dry_weight = dry_weight, days = days)
  n <- length(dry_weight)
  if (n < 2 || length(days) != n) {
    stop("'dry_weight' and 'days' need the same length (at least 2).",
         call. = FALSE)
  }
  if (any(diff(days) <= 0)) stop("'days' must be strictly increasing.",
                                 call. = FALSE)
  if (any(dry_weight <= 0)) stop("'dry_weight' must be positive.",
                                 call. = FALSE)
  if (!is.null(ground_area) && (!is.numeric(ground_area) ||
                                length(ground_area) != 1 ||
                                ground_area <= 0)) {
    stop("'ground_area' must be a single positive number (m2).",
         call. = FALSE)
  }
  dt <- diff(days)
  dw <- diff(dry_weight)
  P <- if (is.null(ground_area)) 1 else ground_area
  res <- data.frame(Period = paste0(days[-n], "-", days[-1], " DAS"),
                    CGR = dw / (P * dt),
                    RGR = diff(log(dry_weight)) / dt,
                    stringsAsFactors = FALSE)
  if (!is.null(leaf_area)) {
    .check_vectors(leaf_area = leaf_area)
    if (length(leaf_area) != n) {
      stop("'leaf_area' must have one value per sampling.", call. = FALSE)
    }
    if (any(leaf_area <= 0)) stop("'leaf_area' must be positive.",
                                  call. = FALSE)
    la1 <- leaf_area[-n]
    la2 <- leaf_area[-1]
    same <- abs(la2 - la1) < 1e-12
    nar <- ifelse(same, dw / (dt * la1),
                  dw / dt * (log(la2) - log(la1)) / (la2 - la1))
    lai <- if (is.null(ground_area)) leaf_area / 10000 else
      leaf_area / (ground_area * 10000)
    res$NAR <- nar * 10000
    res$LAI_start <- lai[-n]
    res$LAI_end <- lai[-1]
    res$LAD <- (lai[-n] + lai[-1]) * dt / 2
  }
  if (verbose) {
    .heading("CROP GROWTH ANALYSIS")
    o <- res
    for (cl in names(o)[-1]) o[[cl]] <- signif(o[[cl]], 4)
    print(o, row.names = FALSE)
    cat("\n CGR:", if (is.null(ground_area)) "g per sample per day" else
      "g m-2 day-1", "| RGR: g g-1 day-1 | NAR: g m-2 (leaf) day-1 |",
      "LAD: days\n")
  }
  class(res) <- c("aridagri_growth", "data.frame")
  invisible(res)
}


#' Harvest Index
#'
#' @description
#' Harvest index and grain-to-straw ratio.
#'
#' @param grain_yield Grain (economic) yield.
#' @param straw_yield Straw (stover) yield, in the same unit.
#' @param verbose Logical. If TRUE (default), prints the table.
#'
#' @return A data frame of class \code{aridagri_hi} with grain, straw and
#'   total biomass, \code{HI_percent} and \code{Grain_Straw_Ratio}.
#'
#' @examples
#' harvest_index(grain_yield = c(1200, 1350, 1100),
#'               straw_yield = c(2400, 2500, 2300))
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram
#'   Choudhary, Hemant Poonia
harvest_index <- function(grain_yield, straw_yield, verbose = TRUE) {
  .check_vectors(grain_yield = grain_yield, straw_yield = straw_yield)
  if (any(grain_yield < 0) || any(straw_yield < 0)) {
    stop("Yields must not be negative.", call. = FALSE)
  }
  total <- grain_yield + straw_yield
  res <- data.frame(Grain_Yield = grain_yield, Straw_Yield = straw_yield,
                    Total_Biomass = total,
                    HI_percent = ifelse(total > 0, 100 * grain_yield / total,
                                        NA_real_),
                    Grain_Straw_Ratio = ifelse(straw_yield > 0,
                                               grain_yield / straw_yield,
                                               NA_real_))
  if (verbose) {
    .heading("HARVEST INDEX")
    if (nrow(res) > 1) {
      cat(sprintf(" Mean HI = %.2f%%   Mean grain:straw ratio = %.3f\n\n",
                  mean(res$HI_percent, na.rm = TRUE),
                  mean(res$Grain_Straw_Ratio, na.rm = TRUE)))
    }
    o <- res
    o$HI_percent <- round(o$HI_percent, 2)
    o$Grain_Straw_Ratio <- round(o$Grain_Straw_Ratio, 3)
    print(o, row.names = FALSE)
  }
  class(res) <- c("aridagri_hi", "data.frame")
  invisible(res)
}


#' Yield Gap and Front-Line Demonstration Indices
#'
#' @description
#' Yield gap between actual and potential yields. When an attainable yield
#' (for example the yield of an improved-practice or front-line
#' demonstration plot) is supplied, the gap is split into a technology gap
#' and an extension gap and the technology index is given.
#'
#' @param actual_yield Actual (farmers' practice) yield.
#' @param potential_yield Potential yield of the crop or variety.
#' @param attainable_yield Attainable or demonstration yield (optional).
#' @param verbose Logical. If TRUE (default), prints the table.
#'
#' @details
#' Total gap = potential - actual; yield efficiency = 100 x actual /
#' potential. With attainable yield (Samui et al., 2000): technology gap =
#' potential - attainable; extension gap = attainable - actual; technology
#' index = 100 x technology gap / potential; yield increase over farmers'
#' practice = 100 x (attainable - actual) / actual.
#'
#' @return A data frame of class \code{aridagri_yield_gap}.
#'
#' @references
#' Samui, S.K., Maitra, S., Roy, D.K., Mondal, A.K. and Saha, D. (2000).
#' Evaluation of front line demonstration on groundnut (Arachis hypogaea
#' L.). \emph{Journal of the Indian Society of Coastal Agricultural
#' Research}, 18(2), 180-183.
#'
#' @examples
#' yield_gap_analysis(actual_yield = c(800, 950, 720), potential_yield = 1500)
#' # Front-line demonstrations: farmers' practice vs demonstration plots
#' yield_gap_analysis(actual_yield = c(820, 760, 905),
#'                    potential_yield = 1500,
#'                    attainable_yield = c(1080, 1010, 1150))
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram
#'   Choudhary, Hemant Poonia
yield_gap_analysis <- function(actual_yield, potential_yield,
                               attainable_yield = NULL, verbose = TRUE) {
  .check_vectors(actual_yield = actual_yield,
                 potential_yield = potential_yield)
  if (any(potential_yield <= 0)) stop("'potential_yield' must be positive.",
                                      call. = FALSE)
  gap <- potential_yield - actual_yield
  res <- data.frame(Actual_Yield = actual_yield,
                    Potential_Yield = potential_yield,
                    Yield_Gap = gap,
                    Gap_Percent = 100 * gap / potential_yield,
                    Yield_Efficiency = 100 * actual_yield / potential_yield)
  if (!is.null(attainable_yield)) {
    .check_vectors(actual_yield = actual_yield,
                   attainable_yield = attainable_yield)
    res$Attainable_Yield <- attainable_yield
    res$Technology_Gap <- potential_yield - attainable_yield
    res$Extension_Gap <- attainable_yield - actual_yield
    res$Technology_Index <- 100 * (potential_yield - attainable_yield) /
      potential_yield
    res$Increase_over_Actual_pct <- ifelse(actual_yield > 0,
                                           100 * (attainable_yield -
                                                    actual_yield) /
                                             actual_yield, NA_real_)
  }
  if (verbose) {
    .heading("YIELD GAP ANALYSIS")
    cat(sprintf(" Mean yield gap = %.1f (%.1f%% of potential)   Mean yield efficiency = %.1f%%\n",
                mean(gap), mean(res$Gap_Percent), mean(res$Yield_Efficiency)))
    if (!is.null(attainable_yield)) {
      cat(sprintf(" Mean technology gap = %.1f   Mean extension gap = %.1f   Mean technology index = %.1f%%\n",
                  mean(res$Technology_Gap), mean(res$Extension_Gap),
                  mean(res$Technology_Index)))
    }
    cat("\n")
    o <- res
    for (cl in names(o)) o[[cl]] <- round(o[[cl]], 1)
    print(o, row.names = FALSE)
  }
  class(res) <- c("aridagri_yield_gap", "data.frame")
  invisible(res)
}


#' Economic Efficiency Indices (B:C Ratio)
#'
#' @description
#' Net return, benefit-cost ratio and return per rupee invested from gross
#' return and total cost.
#'
#' @param gross_return Gross return (Rs/ha).
#' @param total_cost Total cost of cultivation (Rs/ha).
#' @param verbose Logical. If TRUE (default), prints the table.
#'
#' @return A data frame of class \code{aridagri_economics} with
#'   \code{Net_Return}, \code{BC_Ratio} (gross return / cost) and
#'   \code{Return_per_Rupee} (net return / cost).
#'
#' @examples
#' economic_indices(gross_return = c(75000, 82000), total_cost = c(35000, 38000))
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram
#'   Choudhary, Hemant Poonia
economic_indices <- function(gross_return, total_cost, verbose = TRUE) {
  .check_vectors(gross_return = gross_return, total_cost = total_cost)
  if (any(total_cost <= 0)) stop("'total_cost' must be positive.",
                                 call. = FALSE)
  net <- gross_return - total_cost
  res <- data.frame(Gross_Return = gross_return, Total_Cost = total_cost,
                    Net_Return = net, BC_Ratio = gross_return / total_cost,
                    Return_per_Rupee = net / total_cost)
  if (verbose) {
    .heading("ECONOMIC EFFICIENCY")
    cat(sprintf(" Mean net return = Rs %.0f/ha   Mean B:C ratio = %.2f\n\n",
                mean(net), mean(res$BC_Ratio)))
    o <- res
    o$BC_Ratio <- round(o$BC_Ratio, 2)
    o$Return_per_Rupee <- round(o$Return_per_Rupee, 2)
    print(o, row.names = FALSE)
  }
  class(res) <- c("aridagri_economics", "data.frame")
  invisible(res)
}
