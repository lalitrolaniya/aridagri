#' Nutrient Use Efficiency Calculations
#'
#' @description
#' Comprehensive nutrient use efficiency calculations including Agronomic Efficiency (AE),
#' Physiological Efficiency (PE), Apparent Recovery Efficiency (ARE), and Partial Factor
#' Productivity (PFP). Essential for INM research in arid regions.
#'
#' @param yield_fertilized Yield with fertilizer application (kg/ha)
#' @param yield_control Yield in control/unfertilized plot (kg/ha)
#' @param nutrient_applied Amount of nutrient applied (kg/ha)
#' @param nutrient_uptake_fert Nutrient uptake in fertilized plot (kg/ha), optional
#' @param nutrient_uptake_ctrl Nutrient uptake in control plot (kg/ha), optional
#' @param biomass_fert Total biomass in fertilized plot (kg/ha), optional
#' @param biomass_ctrl Total biomass in control plot (kg/ha), optional
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return Data frame with efficiency indices
#'
#' @details
#' Efficiency calculations:
#' \itemize{
#'   \item AE (Agronomic Efficiency) = (Yield_fert - Yield_ctrl) / Nutrient_applied
#'   \item PFP (Partial Factor Productivity) = Yield_fert / Nutrient_applied
#'   \item ARE (Apparent Recovery Efficiency) = (Uptake_fert - Uptake_ctrl) / Nutrient_applied  100
#'   \item PE (Physiological Efficiency) = (Yield_fert - Yield_ctrl) / (Uptake_fert - Uptake_ctrl)
#' }
#'
#' @examples
#' # Basic NUE calculation
#' nue_calculate(yield_fertilized = 1850, yield_control = 1200, nutrient_applied = 40)
#'
#' # Complete NUE with uptake data
#' nue_calculate(
#'   yield_fertilized = 1850, 
#'   yield_control = 1200,
#'   nutrient_applied = 40,
#'   nutrient_uptake_fert = 65,
#'   nutrient_uptake_ctrl = 35
#' )
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
nue_calculate <- function(yield_fertilized, yield_control, nutrient_applied,
                          nutrient_uptake_fert = NULL, nutrient_uptake_ctrl = NULL,
                          biomass_fert = NULL, biomass_ctrl = NULL,
                            verbose = TRUE) {
  
  # Input validation
  if (nutrient_applied <= 0) {
    stop("nutrient_applied must be greater than 0")
  }
  
  # Initialize results
  results <- data.frame(
    Yield_Fertilized = yield_fertilized,
    Yield_Control = yield_control,
    Nutrient_Applied = nutrient_applied
  )
  
  # Agronomic Efficiency (AE) - kg grain per kg nutrient
  ae <- (yield_fertilized - yield_control) / nutrient_applied
  results$AE_kg_kg <- round(ae, 2)
  
  # Partial Factor Productivity (PFP) - kg grain per kg nutrient
  pfp <- yield_fertilized / nutrient_applied
  results$PFP_kg_kg <- round(pfp, 2)
  
  # Apparent Recovery Efficiency (ARE) - if uptake data available
  if (!is.null(nutrient_uptake_fert) && !is.null(nutrient_uptake_ctrl)) {
    are <- ((nutrient_uptake_fert - nutrient_uptake_ctrl) / nutrient_applied) * 100
    results$ARE_percent <- round(are, 1)
    
    # Physiological Efficiency (PE)
    uptake_diff <- nutrient_uptake_fert - nutrient_uptake_ctrl
    if (uptake_diff > 0) {
      pe <- (yield_fertilized - yield_control) / uptake_diff
      results$PE_kg_kg <- round(pe, 2)
    }
    
    # Utilization Efficiency (UE) = AE  ARE
    ue <- ae * (are / 100)
    results$UE_kg_kg <- round(ue, 2)
  }
  
  # Harvest Index if biomass available
  if (!is.null(biomass_fert)) {
    hi_fert <- (yield_fertilized / biomass_fert) * 100
    results$HI_Fertilized <- round(hi_fert, 1)
  }
  if (!is.null(biomass_ctrl)) {
    hi_ctrl <- (yield_control / biomass_ctrl) * 100
    results$HI_Control <- round(hi_ctrl, 1)
  }
  
  # Print summary
  if (verbose) {
    cat("\n=== Nutrient Use Efficiency Analysis ===\n")
    cat("Yield (Fertilized):", yield_fertilized, "kg/ha\n")
    cat("Yield (Control):", yield_control, "kg/ha\n")
    cat("Yield Response:", yield_fertilized - yield_control, "kg/ha\n")
    cat("Nutrient Applied:", nutrient_applied, "kg/ha\n")
    cat("\n--- Efficiency Indices ---\n")
    cat("Agronomic Efficiency (AE):", round(ae, 2), "kg grain/kg nutrient\n")
    cat("Partial Factor Productivity (PFP):", round(pfp, 2), "kg grain/kg nutrient\n")

  }
  if (!is.null(nutrient_uptake_fert) && !is.null(nutrient_uptake_ctrl)) {
    if (verbose) {
      cat("Apparent Recovery Efficiency (ARE):", round(are, 1), "%\n")
      if (exists("pe")) {
        cat("Physiological Efficiency (PE):", round(pe, 2), "kg grain/kg uptake\n")
      }
      cat("Utilization Efficiency (UE):", round(ue, 2), "kg grain/kg nutrient\n")
    }
  }
  
  class(results) <- c("aridagri_nue", "data.frame")
  return(invisible(results))
}


#' Nutrient Response Curve Analysis
#'
#' @description
#' Fits nutrient response curves using quadratic, linear-plateau, or 
#' Mitscherlich models to determine economic optimum dose.
#'
#' @param dose Numeric vector of nutrient doses (kg/ha)
#' @param yield Numeric vector of corresponding yields (kg/ha)
#' @param model Model type: "quadratic", "linear_plateau", or "mitscherlich"
#' @param price_output Price of output (Rs/kg)
#' @param price_nutrient Price of nutrient (Rs/kg)
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return List with model parameters and economic optimum
#'
#' @examples
#' dose <- c(0, 20, 40, 60, 80, 100)
#' yield <- c(1100, 1350, 1520, 1610, 1650, 1660)
#' nutrient_response(dose, yield, model = "quadratic", 
#'                   price_output = 60, price_nutrient = 15)
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
nutrient_response <- function(dose, yield, model = "quadratic",
                               price_output = 50, price_nutrient = 15,
                            verbose = TRUE) {
  
  # Input validation
  if (length(dose) != length(yield)) {
    stop("dose and yield must have same length")
  }
  
  if (model == "quadratic") {
    # Fit quadratic model: Y = a + bX + cX
    fit <- lm(yield ~ dose + I(dose^2))
    coef_vals <- coef(fit)
    a <- coef_vals[1]
    b <- coef_vals[2]
    c <- coef_vals[3]
    
    # Maximum yield dose (dy/dx = 0)
    dose_max <- -b / (2 * c)
    yield_max <- predict(fit, newdata = data.frame(dose = dose_max))
    
    # Economic optimum dose (dy/dx = price ratio)
    price_ratio <- price_nutrient / price_output
    dose_econ <- (price_ratio - b) / (2 * c)
    yield_econ <- predict(fit, newdata = data.frame(dose = dose_econ))
    
    # Model statistics
    r_squared <- summary(fit)$r.squared
    
    if (verbose) {
      cat("\n=== Nutrient Response Curve Analysis ===\n")
      cat("Model: Quadratic (Y = a + bX + cX)\n")
      cat("\n--- Model Parameters ---\n")
      cat("Intercept (a):", round(a, 2), "\n")
      cat("Linear coefficient (b):", round(b, 4), "\n")
      cat("Quadratic coefficient (c):", round(c, 6), "\n")
      cat("R-squared:", round(r_squared, 4), "\n")
      cat("\n--- Optimum Doses ---\n")
      cat("Dose for maximum yield:", round(dose_max, 1), "kg/ha\n")
      cat("Maximum yield:", round(yield_max, 0), "kg/ha\n")
      cat("Economic optimum dose:", round(dose_econ, 1), "kg/ha\n")
      cat("Yield at economic optimum:", round(yield_econ, 0), "kg/ha\n")
      cat("\n--- Price Parameters ---\n")
      cat("Output price: Rs", price_output, "/kg\n")
      cat("Nutrient price: Rs", price_nutrient, "/kg\n")
    }
    
    result <- list(
      model = "quadratic",
      coefficients = c(a = a, b = b, c = c),
      r_squared = r_squared,
      dose_max_yield = round(dose_max, 1),
      max_yield = round(yield_max, 0),
      dose_economic = round(dose_econ, 1),
      yield_economic = round(yield_econ, 0),
      fit = fit
    )
    
  } else {
    stop("Currently only 'quadratic' model is implemented")
  }
  
  class(result) <- "aridagri_response"
  return(invisible(result))
}


#' Economic Analysis for Arid Agriculture
#'
#' @description
#' Calculates economic indicators including Cost of Cultivation, Gross Returns,
#' Net Returns, B:C Ratio, and profitability indices for arid farming systems.
#'
#' @param yield Crop yield (kg/ha)
#' @param price Market price (Rs/kg)
#' @param cost_fixed Fixed costs (Rs/ha)
#' @param cost_variable Variable costs (Rs/ha)
#' @param byproduct_yield Byproduct yield (kg/ha), optional
#' @param byproduct_price Byproduct price (Rs/kg), optional
#'
#' @param verbose Logical. If TRUE (default), prints formatted output to console.
#'
#' @return Data frame with economic analysis
#'
#' @examples
#' economic_analysis(yield = 1200, price = 65, cost_fixed = 15000, cost_variable = 12000)
#'
#' # With byproduct
#' economic_analysis(yield = 1200, price = 65, cost_fixed = 15000, 
#'                   cost_variable = 12000, byproduct_yield = 1800, byproduct_price = 5)
#'
#' @export
#' @author Lalit Kumar Rolaniya, ICAR-IIPR, Bikaner
economic_analysis <- function(yield, price, cost_fixed, cost_variable,
                               byproduct_yield = 0, byproduct_price = 0,
                            verbose = TRUE) {
  
  # Calculate costs
  total_cost <- cost_fixed + cost_variable
  
  # Calculate returns
  gross_return_main <- yield * price
  gross_return_byproduct <- byproduct_yield * byproduct_price
  gross_return <- gross_return_main + gross_return_byproduct
  
  # Net returns
  net_return <- gross_return - total_cost
  
  # Profitability ratios
  bc_ratio <- gross_return / total_cost
  return_per_rupee <- net_return / total_cost
  
  # Per kg cost of production
  cost_of_production <- total_cost / yield
  
  if (verbose) {
    cat("\n=== Economic Analysis ===\n")
    cat("\n--- Costs (Rs/ha) ---\n")
    cat("Fixed Costs:", format(cost_fixed, big.mark = ","), "\n")
    cat("Variable Costs:", format(cost_variable, big.mark = ","), "\n")
    cat("Total Cost:", format(total_cost, big.mark = ","), "\n")
    cat("\n--- Returns (Rs/ha) ---\n")
    cat("Gross Return (main):", format(gross_return_main, big.mark = ","), "\n")
    if (byproduct_yield > 0) {
      cat("Gross Return (byproduct):", format(gross_return_byproduct, big.mark = ","), "\n")
    }
    cat("Total Gross Return:", format(gross_return, big.mark = ","), "\n")
    cat("Net Return:", format(net_return, big.mark = ","), "\n")
    cat("\n--- Profitability Indices ---\n")
    cat("B:C Ratio:", round(bc_ratio, 2), "\n")
    cat("Return per Rupee Invested:", round(return_per_rupee, 2), "\n")
    cat("Cost of Production:", round(cost_of_production, 2), "Rs/kg\n")
  }
  
  result <- data.frame(
    Yield_kg_ha = yield,
    Price_Rs_kg = price,
    Cost_Fixed = cost_fixed,
    Cost_Variable = cost_variable,
    Total_Cost = total_cost,
    Gross_Return = gross_return,
    Net_Return = net_return,
    BC_Ratio = round(bc_ratio, 2),
    Cost_of_Production = round(cost_of_production, 2)
  )
  
  class(result) <- c("aridagri_economics", "data.frame")
  return(invisible(result))
}
