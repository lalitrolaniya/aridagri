#' ============================================================================
#' ADVANCED AGRONOMIC STATISTICAL METHODS
#' Package: aridagri
#' Authors: Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram Choudhary
#' ICAR-Indian Institute of Pulses Research, Regional Centre, Bikaner
#' ICAR-Indian Institute of Groundnut Research, Regional Research Station, Bikaner
#' ============================================================================

#' Stability Analysis for Agronomic Traits (Multiple Methods)
#'
#' @description
#' Performs comprehensive stability analysis using multiple established methods for 
#' evaluating genotype/treatment performance across environments.
#'
#' Methods included:
#' \itemize{
#'   \item Eberhart & Russell (1966): Regression approach
#'   \item AMMI: Additive Main effects and Multiplicative Interaction
#'   \item Finlay & Wilkinson (1963): Linear regression on environmental mean
#'   \item Shukla (1972): Stability variance
#'   \item Wricke (1962): Ecovalence
#'   \item Coefficient of Variation: CV-based ranking
#'   \item Superiority Index: Lin & Binns (1988)
#' }
#'
#' @param data Data frame with genotype/treatment, environment, replication, and trait
#' @param genotype Name of genotype/treatment column
#' @param environment Name of environment/location/year column
#' @param replication Name of replication column
#' @param trait Name of trait/response variable
#' @param method Method: "eberhart", "ammi", "finlay", "shukla", "wricke", "cv", "superiority", or "all"
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return List containing stability parameters and rankings
#'
#' @examples
#' data <- expand.grid(
#'   variety = paste0("V", 1:10),
#'   location = paste0("L", 1:5),
#'   rep = 1:3
#' )
#' data$yield <- rnorm(nrow(data), 1200, 200)
#' 
#' stability_analysis(data, genotype = "variety", environment = "location",
#'                    replication = "rep", trait = "yield", method = "all")
#'
#' @references
#' Eberhart, S.A. and Russell, W.A. (1966). Crop Science, 6: 36-40.
#' 
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram Choudhary
stability_analysis <- function(data, genotype, environment, replication, 
                                trait, method = "all",
                            verbose = TRUE) {
  
  # Convert to factors
  data[[genotype]] <- as.factor(data[[genotype]])
  data[[environment]] <- as.factor(data[[environment]])
  data[[replication]] <- as.factor(data[[replication]])
  
  g <- nlevels(data[[genotype]])
  e <- nlevels(data[[environment]])
  r <- nlevels(data[[replication]])
  
  # Calculate means
  ge_means <- aggregate(data[[trait]], 
                         by = list(data[[genotype]], data[[environment]]), 
                         FUN = mean, na.rm = TRUE)
  names(ge_means) <- c("Genotype", "Environment", "Mean")
  
  g_means <- aggregate(ge_means$Mean, by = list(ge_means$Genotype), FUN = mean)
  names(g_means) <- c("Genotype", "Mean")
  
  e_means <- aggregate(ge_means$Mean, by = list(ge_means$Environment), FUN = mean)
  names(e_means) <- c("Environment", "Mean")
  
  grand_mean <- mean(ge_means$Mean)
  e_means$Index <- e_means$Mean - grand_mean
  
  if (verbose) {
    cat("\n")
    cat("\n")
    cat("                    STABILITY ANALYSIS FOR AGRONOMIC TRAITS                         \n")
    cat("\n")
    cat(" Number of Genotypes/Treatments:", sprintf("%-51d", g), "\n")
    cat(" Number of Environments        :", sprintf("%-51d", e), "\n")
    cat(" Replications/Environment      :", sprintf("%-51d", r), "\n")
    cat(" Grand Mean                    :", sprintf("%-51.2f", grand_mean), "\n")
    cat("\n")
  }
  
  results <- list()
  results$grand_mean <- grand_mean
  results$genotype_means <- g_means
  results$environment_means <- e_means
  
  # Create GE matrix
  ge_matrix <- reshape(ge_means, idvar = "Genotype", timevar = "Environment", 
                        direction = "wide")
  rownames(ge_matrix) <- ge_matrix$Genotype
  ge_matrix <- ge_matrix[, -1]
  names(ge_matrix) <- gsub("Mean.", "", names(ge_matrix))
  ge_matrix <- as.matrix(ge_matrix)
  
  # EBERHART AND RUSSELL
  if (method %in% c("eberhart", "all")) {
    if (verbose) {
      cat("\n\n")
      cat("              EBERHART & RUSSELL (1966) STABILITY ANALYSIS                         \n")
      cat("\n")
    }
    
    stability_er <- data.frame(Genotype = rownames(ge_matrix))
    stability_er$Mean <- rowMeans(ge_matrix)
    
    env_index <- colMeans(ge_matrix) - grand_mean
    bi <- numeric(g)
    s2di <- numeric(g)
    r_squared <- numeric(g)
    
    for (i in 1:g) {
      y <- ge_matrix[i, ]
      reg <- lm(y ~ env_index)
      bi[i] <- coef(reg)[2]
      predicted <- predict(reg)
      s2di[i] <- sum((y - predicted)^2) / (e - 2)
      r_squared[i] <- summary(reg)$r.squared
    }
    
    stability_er$bi <- round(bi, 3)
    stability_er$S2di <- round(s2di, 3)
    stability_er$R2 <- round(r_squared, 3)
    stability_er <- stability_er[order(-stability_er$Mean), ]
    
    if (verbose) {
      print(stability_er, row.names = FALSE)
      cat("\n bi  1.0, low Sdi: Wide adaptation | bi > 1.0: Favorable env | bi < 1.0: Poor env\n")
    }
    
    results$eberhart <- stability_er
  }
  
  # FINLAY & WILKINSON
  if (method %in% c("finlay", "all")) {
    if (verbose) {
      cat("\n\n")
      cat("                  FINLAY & WILKINSON (1963) STABILITY ANALYSIS                     \n")
      cat("\n")
    }
    
    stability_fw <- data.frame(Genotype = rownames(ge_matrix))
    stability_fw$Mean <- rowMeans(ge_matrix)
    
    env_means <- colMeans(ge_matrix)
    bi_fw <- numeric(g)
    for (i in 1:g) {
      reg <- lm(ge_matrix[i, ] ~ env_means)
      bi_fw[i] <- coef(reg)[2]
    }
    
    stability_fw$bi <- round(bi_fw, 3)
    stability_fw <- stability_fw[order(-stability_fw$Mean), ]
    if (verbose) {
      print(stability_fw, row.names = FALSE)
    }
    
    results$finlay <- stability_fw
  }
  
  # SHUKLA'S STABILITY VARIANCE
  if (method %in% c("shukla", "all")) {
    if (verbose) {
      cat("\n\n")
      cat("                    SHUKLA'S (1972) STABILITY VARIANCE (i)                       \n")
      cat("\n")
    }
    
    stability_shukla <- data.frame(Genotype = rownames(ge_matrix))
    stability_shukla$Mean <- rowMeans(ge_matrix)
    
    sigma2_i <- numeric(g)
    for (i in 1:g) {
      interaction <- ge_matrix[i, ] - rowMeans(ge_matrix)[i] - colMeans(ge_matrix) + grand_mean
      sigma2_i[i] <- sum(interaction^2) / (e - 1)
    }
    
    stability_shukla$Sigma2_i <- round(sigma2_i, 2)
    stability_shukla$Rank <- rank(sigma2_i)
    stability_shukla <- stability_shukla[order(stability_shukla$Sigma2_i), ]
    if (verbose) {
      print(stability_shukla, row.names = FALSE)
      cat("\n Lower i indicates higher stability\n")
    }
    
    results$shukla <- stability_shukla
  }
  
  # WRICKE'S ECOVALENCE
  if (method %in% c("wricke", "all")) {
    if (verbose) {
      cat("\n\n")
      cat("                       WRICKE'S (1962) ECOVALENCE (Wi)                             \n")
      cat("\n")
    }
    
    stability_wricke <- data.frame(Genotype = rownames(ge_matrix))
    stability_wricke$Mean <- rowMeans(ge_matrix)
    
    Wi <- numeric(g)
    for (i in 1:g) {
      interaction <- ge_matrix[i, ] - rowMeans(ge_matrix)[i] - colMeans(ge_matrix) + grand_mean
      Wi[i] <- sum(interaction^2)
    }
    
    stability_wricke$Wi <- round(Wi, 2)
    stability_wricke$Wi_percent <- round(Wi / sum(Wi) * 100, 2)
    stability_wricke$Rank <- rank(Wi)
    stability_wricke <- stability_wricke[order(stability_wricke$Wi), ]
    if (verbose) {
      print(stability_wricke, row.names = FALSE)
      cat("\n Lower Wi indicates higher stability\n")
    }
    
    results$wricke <- stability_wricke
  }
  
  # AMMI ANALYSIS
  if (method %in% c("ammi", "all")) {
    if (verbose) {
      cat("\n\n")
      cat("            AMMI (Additive Main Effects and Multiplicative Interaction)            \n")
      cat("\n")
    }
    
    interaction_matrix <- ge_matrix
    for (i in 1:g) {
      for (j in 1:e) {
        interaction_matrix[i, j] <- ge_matrix[i, j] - 
          rowMeans(ge_matrix)[i] - colMeans(ge_matrix)[j] + grand_mean
      }
    }
    
    svd_result <- svd(interaction_matrix)
    eigenvalues <- svd_result$d^2
    var_explained <- eigenvalues / sum(eigenvalues) * 100
    cum_var <- cumsum(var_explained)
    n_pc <- min(g-1, e-1)
    
    if (verbose) {
      cat("\n--- AMMI Components ---\n\n")
    }
    ammi_components <- data.frame(
      IPCA = paste0("IPCA", 1:n_pc),
      Eigenvalue = round(eigenvalues[1:n_pc], 3),
      Variance_Percent = round(var_explained[1:n_pc], 2),
      Cumulative_Percent = round(cum_var[1:n_pc], 2)
    )
    if (verbose) {
      print(ammi_components, row.names = FALSE)
    }
    
    # AMMI Stability Value (ASV)
    pc1_scores <- svd_result$u[, 1] * svd_result$d[1]
    pc2_scores <- svd_result$u[, 2] * svd_result$d[2]
    weight <- eigenvalues[1] / eigenvalues[2]
    asv <- sqrt(weight * pc1_scores^2 + pc2_scores^2)
    
    mean_rank <- rank(-rowMeans(ge_matrix))
    asv_rank <- rank(asv)
    ysi <- mean_rank + asv_rank
    
    ammi_stability <- data.frame(
      Genotype = rownames(ge_matrix),
      Mean = round(rowMeans(ge_matrix), 2),
      Mean_Rank = mean_rank,
      IPCA1 = round(pc1_scores, 3),
      IPCA2 = round(pc2_scores, 3),
      ASV = round(asv, 3),
      ASV_Rank = asv_rank,
      YSI = ysi
    )
    ammi_stability <- ammi_stability[order(ammi_stability$YSI), ]
    
    if (verbose) {
      cat("\n--- AMMI Stability Values ---\n\n")
      print(ammi_stability, row.names = FALSE)
      cat("\n ASV: Lower = More Stable | YSI: Lower = Better\n")
    }
    
    results$ammi <- list(components = ammi_components, stability = ammi_stability)
  }
  
  # CV STABILITY
  if (method %in% c("cv", "all")) {
    if (verbose) {
      cat("\n\n")
      cat("                    COEFFICIENT OF VARIATION (CV) STABILITY                        \n")
      cat("\n")
    }
    
    stability_cv <- data.frame(Genotype = rownames(ge_matrix))
    stability_cv$Mean <- rowMeans(ge_matrix)
    stability_cv$SD <- apply(ge_matrix, 1, sd)
    stability_cv$CV <- round((stability_cv$SD / stability_cv$Mean) * 100, 2)
    stability_cv$CV_Rank <- rank(stability_cv$CV)
    stability_cv <- stability_cv[order(stability_cv$CV), ]
    if (verbose) {
      print(stability_cv, row.names = FALSE)
      cat("\n Lower CV indicates higher stability\n")
    }
    
    results$cv_stability <- stability_cv
  }
  
  # SUPERIORITY INDEX
  if (method %in% c("superiority", "all")) {
    if (verbose) {
      cat("\n\n")
      cat("                  CULTIVAR SUPERIORITY INDEX (Lin & Binns, 1988)                   \n")
      cat("\n")
    }
    
    stability_sup <- data.frame(Genotype = rownames(ge_matrix))
    stability_sup$Mean <- rowMeans(ge_matrix)
    
    max_env <- apply(ge_matrix, 2, max)
    Pi <- numeric(g)
    for (i in 1:g) {
      Pi[i] <- sum((ge_matrix[i, ] - max_env)^2) / (2 * e)
    }
    
    stability_sup$Pi <- round(Pi, 2)
    stability_sup$Pi_Rank <- rank(Pi)
    stability_sup <- stability_sup[order(stability_sup$Pi), ]
    if (verbose) {
      print(stability_sup, row.names = FALSE)
      cat("\n Lower Pi indicates superior performance\n")
    }
    
    results$superiority <- stability_sup
  }
  
  # INTEGRATED RANKING
  if (verbose) {
    cat("\n\n")
    cat("                        INTEGRATED STABILITY RANKING                                \n")
    cat("\n")
  }
  
  integrated <- data.frame(Genotype = rownames(ge_matrix))
  integrated$Mean <- round(rowMeans(ge_matrix), 2)
  integrated$Mean_Rank <- rank(-integrated$Mean)
  
  env_index <- colMeans(ge_matrix) - grand_mean
  s2di_temp <- numeric(g)
  for (i in 1:g) {
    y <- ge_matrix[i, ]
    reg <- lm(y ~ env_index)
    s2di_temp[i] <- sum((y - predict(reg))^2) / (e - 2)
  }
  integrated$ER_Rank <- rank(s2di_temp)
  
  sigma2_temp <- numeric(g)
  for (i in 1:g) {
    interaction <- ge_matrix[i, ] - rowMeans(ge_matrix)[i] - colMeans(ge_matrix) + grand_mean
    sigma2_temp[i] <- sum(interaction^2) / (e - 1)
  }
  integrated$Shukla_Rank <- rank(sigma2_temp)
  
  cv_temp <- apply(ge_matrix, 1, sd) / rowMeans(ge_matrix) * 100
  integrated$CV_Rank <- rank(cv_temp)
  
  integrated$Overall_Rank <- rowMeans(integrated[, c("Mean_Rank", "ER_Rank", "Shukla_Rank", "CV_Rank")])
  integrated$Final_Rank <- rank(integrated$Overall_Rank)
  integrated <- integrated[order(integrated$Final_Rank), ]
  if (verbose) {
    print(integrated, row.names = FALSE)
  }
  
  results$integrated <- integrated
  results$ge_matrix <- ge_matrix
  
  class(results) <- c("aridagri_stability", "list")
  return(invisible(results))
}


#' Growing Degree Days (GDD) and Thermal Indices
#'
#' @description
#' Calculates Growing Degree Days, Helio-thermal Units, Photo-thermal Units,
#' and Heat Use Efficiency from temperature data.
#'
#' @param tmax Vector of daily maximum temperatures (C)
#' @param tmin Vector of daily minimum temperatures (C)
#' @param base_temp Base temperature (C)
#' @param sunshine_hours Vector of daily sunshine hours (optional)
#' @param day_length Vector of day length in hours (optional)
#' @param crop_yield Crop yield for HUE calculation (optional)
#' @param biomass Biomass for HUE calculation (optional)
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return List with thermal indices
#'
#' @examples
#' tmax <- runif(90, 30, 42)
#' tmin <- runif(90, 18, 28)
#' thermal_indices(tmax, tmin, base_temp = 10)
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram Choudhary
thermal_indices <- function(tmax, tmin, base_temp = 10, sunshine_hours = NULL,
                            day_length = NULL, crop_yield = NULL, biomass = NULL,
                            verbose = TRUE) {
  
  n_days <- length(tmax)
  tmean <- (tmax + tmin) / 2
  gdd_daily <- pmax(0, tmean - base_temp)
  gdd_cumulative <- cumsum(gdd_daily)
  gdd_total <- sum(gdd_daily)
  
  if (verbose) {
    cat("\n")
    cat("\n")
    cat("                    THERMAL INDICES ANALYSIS                                        \n")
    cat("\n")
    cat(sprintf(" Number of Days         : %10d                                               \n", n_days))
    cat(sprintf(" Base Temperature       : %10.1f C                                          \n", base_temp))
    cat(sprintf(" Mean Tmax              : %10.1f C                                          \n", mean(tmax)))
    cat(sprintf(" Mean Tmin              : %10.1f C                                          \n", mean(tmin)))
    cat("\n")
  }
  
  results <- data.frame(Index = character(), Value = numeric(), Unit = character())
  
  results <- rbind(results, data.frame(Index = "Growing Degree Days (GDD)", 
                                        Value = round(gdd_total, 1), Unit = "C-days"))
  results <- rbind(results, data.frame(Index = "Mean Daily GDD", 
                                        Value = round(mean(gdd_daily), 2), Unit = "C-days/day"))
  
  if (!is.null(sunshine_hours)) {
    htu_total <- sum(gdd_daily * sunshine_hours)
    results <- rbind(results, data.frame(Index = "Helio-thermal Units (HTU)", 
                                          Value = round(htu_total, 1), Unit = "C-days-hours"))
  }
  
  if (!is.null(day_length)) {
    ptu_total <- sum(gdd_daily * day_length)
    results <- rbind(results, data.frame(Index = "Photo-thermal Units (PTU)", 
                                          Value = round(ptu_total, 1), Unit = "C-days-hours"))
  }
  
  if (!is.null(crop_yield)) {
    hue <- crop_yield / gdd_total
    results <- rbind(results, data.frame(Index = "Heat Use Efficiency (Grain)", 
                                          Value = round(hue, 3), Unit = "kg/ha/C-day"))
  }
  
  if (!is.null(biomass)) {
    hue_b <- biomass / gdd_total
    results <- rbind(results, data.frame(Index = "Heat Use Efficiency (Biomass)", 
                                          Value = round(hue_b, 3), Unit = "kg/ha/C-day"))
  }
  
  if (verbose) {
    cat("\n")
    print(results, row.names = FALSE)
  }
  
  output <- list(summary = results, gdd_total = gdd_total, 
                  daily = data.frame(Day = 1:n_days, Tmax = tmax, Tmin = tmin, 
                                      Tmean = tmean, GDD = gdd_daily, Cum_GDD = gdd_cumulative))
  class(output) <- "aridagri_thermal"
  return(invisible(output))
}


#' Crop Growth Analysis (CGR, RGR, NAR, LAI)
#'
#' @description
#' Calculates crop growth parameters from sequential harvest data.
#'
#' @param dry_weight Vector of dry matter at different stages
#' @param leaf_area Vector of leaf area at different stages
#' @param days Vector of days after sowing
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return Data frame with growth parameters
#'
#' @examples
#' dry_weight <- c(0.5, 2.1, 8.5, 25, 45, 62, 75)
#' leaf_area <- c(15, 85, 350, 800, 950, 850, 600)
#' days <- c(15, 30, 45, 60, 75, 90, 105)
#' crop_growth_analysis(dry_weight, leaf_area, days)
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram Choudhary
crop_growth_analysis <- function(dry_weight, leaf_area, days,
                            verbose = TRUE) {
  
  n <- length(dry_weight)
  
  if (verbose) {
    cat("\n")
    cat("\n")
    cat("                         CROP GROWTH ANALYSIS                                       \n")
    cat("\n")
  }
  
  cgr <- rgr <- nar <- numeric(n - 1)
  lai <- leaf_area / 10000
  
  for (i in 1:(n-1)) {
    dt <- days[i+1] - days[i]
    dw <- dry_weight[i+1] - dry_weight[i]
    cgr[i] <- dw / dt
    rgr[i] <- (log(dry_weight[i+1]) - log(dry_weight[i])) / dt
    mean_la <- (leaf_area[i] + leaf_area[i+1]) / 2
    if (leaf_area[i+1] != leaf_area[i]) {
      nar[i] <- dw * (log(leaf_area[i+1]) - log(leaf_area[i])) / (dt * (leaf_area[i+1] - leaf_area[i]))
    } else {
      nar[i] <- dw / (dt * mean_la)
    }
  }
  
  results <- data.frame(
    Period = paste0(days[-n], "-", days[-1], " DAS"),
    CGR = round(cgr, 4),
    RGR = round(rgr, 5),
    NAR = round(nar * 10000, 4)
  )
  
  if (verbose) {
    print(results, row.names = FALSE)
    cat("\nCGR = Crop Growth Rate | RGR = Relative Growth Rate | NAR = Net Assimilation Rate\n")
  }
  
  class(results) <- c("aridagri_growth", "data.frame")
  return(invisible(results))
}


#' Harvest Index Calculation
#'
#' @description
#' Calculates Harvest Index and related partitioning indices.
#'
#' @param grain_yield Grain/economic yield
#' @param straw_yield Straw/stover yield
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return Data frame with harvest indices
#'
#' @examples
#' harvest_index(grain_yield = c(1200, 1350, 1100), straw_yield = c(2400, 2500, 2300))
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram Choudhary
harvest_index <- function(grain_yield, straw_yield,
                            verbose = TRUE) {
  
  total_biomass <- grain_yield + straw_yield
  hi <- grain_yield / total_biomass * 100
  gsr <- grain_yield / straw_yield
  
  results <- data.frame(
    Grain_Yield = grain_yield,
    Straw_Yield = straw_yield,
    Total_Biomass = total_biomass,
    HI_percent = round(hi, 2),
    Grain_Straw_Ratio = round(gsr, 3)
  )
  
  if (verbose) {
    cat("\n")
    cat("\n")
    cat("                         HARVEST INDEX ANALYSIS                                     \n")
    cat("\n")

    if (length(grain_yield) > 1) {
      cat(sprintf("\nMean HI: %.2f%% | Mean Grain:Straw Ratio: %.3f\n\n", mean(hi), mean(gsr)))
    }
    print(results, row.names = FALSE)
  }
  
  class(results) <- c("aridagri_hi", "data.frame")
  return(invisible(results))
}


#' Yield Gap Analysis
#'
#' @description
#' Calculates yield gaps comparing actual with potential yields.
#'
#' @param actual_yield Actual yield (kg/ha)
#' @param potential_yield Potential yield (kg/ha)
#' @param attainable_yield Attainable yield (kg/ha, optional)
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return Data frame with yield gap analysis
#'
#' @examples
#' yield_gap_analysis(actual_yield = c(800, 950, 720), potential_yield = 1500)
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram Choudhary
yield_gap_analysis <- function(actual_yield, potential_yield, attainable_yield = NULL,
                            verbose = TRUE) {
  
  gap <- potential_yield - actual_yield
  gap_pct <- (gap / potential_yield) * 100
  efficiency <- (actual_yield / potential_yield) * 100
  
  results <- data.frame(
    Actual_Yield = actual_yield,
    Potential_Yield = potential_yield,
    Yield_Gap = round(gap, 1),
    Gap_Percent = round(gap_pct, 1),
    Yield_Efficiency = round(efficiency, 1)
  )
  
  if (verbose) {
    cat("\n")
    cat("\n")
    cat("                           YIELD GAP ANALYSIS                                       \n")
    cat("\n")
    cat(sprintf(" Mean Yield Gap          : %10.1f kg/ha (%.1f%%)                                \n", mean(gap), mean(gap_pct)))
    cat(sprintf(" Mean Yield Efficiency   : %10.1f %%                                           \n", mean(efficiency)))
    cat("\n\n")

    print(results, row.names = FALSE)
  }
  
  class(results) <- c("aridagri_yield_gap", "data.frame")
  return(invisible(results))
}


#' Economic Efficiency Indices (B:C Ratio)
#'
#' @description
#' Calculates B:C ratio and economic efficiency metrics.
#'
#' @param gross_return Gross returns (Rs/ha)
#' @param total_cost Total cost (Rs/ha)
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return Data frame with economic indices
#'
#' @examples
#' economic_indices(gross_return = c(75000, 82000), total_cost = c(35000, 38000))
#'
#' @export
#' @author Lalit Kumar Rolaniya, Ram Lal Jat, Monika Punia, Raja Ram Choudhary
economic_indices <- function(gross_return, total_cost,
                            verbose = TRUE) {
  
  net_return <- gross_return - total_cost
  bc_ratio <- gross_return / total_cost
  return_per_rupee <- net_return / total_cost
  
  results <- data.frame(
    Gross_Return = gross_return,
    Total_Cost = total_cost,
    Net_Return = round(net_return, 0),
    BC_Ratio = round(bc_ratio, 2),
    Return_per_Rupee = round(return_per_rupee, 2)
  )
  
  if (verbose) {
    cat("\n")
    cat("\n")
    cat("                         ECONOMIC EFFICIENCY ANALYSIS                               \n")
    cat("\n")
    cat(sprintf(" Mean Net Return         : Rs. %10.0f /ha                                      \n", mean(net_return)))
    cat(sprintf(" Mean B:C Ratio          : %14.2f                                              \n", mean(bc_ratio)))
    cat("\n\n")

    print(results, row.names = FALSE)
  }
  
  class(results) <- c("aridagri_economics", "data.frame")
  return(invisible(results))
}
