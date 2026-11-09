test_that("day length follows FAO-56", {
  # equator: 12 h every day; higher latitude: longer days in June
  expect_equal(day_length(0, c(1, 100, 200)), rep(12, 3), tolerance = 1e-10)
  expect_gt(day_length(28, 172), day_length(28, 355))
  expect_equal(day_length(28, as.Date("2026-06-21")), day_length(28, 172))
  expect_error(day_length(100, 10))
})

test_that("thermal indices: upper threshold and modified method", {
  tx <- c(30, 45)
  tn <- c(5, 25)
  a <- thermal_indices(tx, tn, base_temp = 10, verbose = FALSE)
  expect_equal(a$gdd_total, (17.5 - 10) + (35 - 10))
  b <- thermal_indices(tx, tn, base_temp = 10, upper_temp = 30,
                       verbose = FALSE)
  expect_equal(b$gdd_total, (17.5 - 10) + (30 - 10))
  m <- thermal_indices(tx, tn, base_temp = 10, upper_temp = 40,
                       method = "modified", verbose = FALSE)
  expect_equal(m$gdd_total, ((30 + 10) / 2 - 10) + ((40 + 25) / 2 - 10))
  h <- thermal_indices(tx, tn, 10, sunshine_hours = c(10, 9),
                       crop_yield = 1000, verbose = FALSE)
  expect_true("Helio-thermal units (HTU)" %in% h$summary$Index)
  expect_error(thermal_indices(c(20, 30), c(25, 10), verbose = FALSE),
               "below")
})

test_that("crop growth analysis returns LAI and LAD", {
  dw <- c(10, 30, 70)
  la <- c(500, 1500, 3000)
  d <- c(20, 35, 50)
  g <- crop_growth_analysis(dw, la, d, ground_area = 0.1, verbose = FALSE)
  expect_equal(g$CGR, c(20, 40) / (0.1 * 15))
  expect_equal(g$RGR, diff(log(dw)) / 15)
  lai <- la / (0.1 * 10000)
  expect_equal(g$LAI_start, lai[1:2])
  expect_equal(g$LAD, (lai[1:2] + lai[2:3]) * 15 / 2)
  nar1 <- (20 / 15) * (log(1500) - log(500)) / (1500 - 500) * 10000
  expect_equal(g$NAR[1], nar1)
  expect_error(crop_growth_analysis(dw, la, c(20, 20, 50), verbose = FALSE))
  g2 <- crop_growth_analysis(dw, days = d, verbose = FALSE)
  expect_false("NAR" %in% names(g2))
})

test_that("yield gap gives front-line demonstration indices", {
  y <- yield_gap_analysis(actual_yield = c(800, 900), potential_yield = 1500,
                          attainable_yield = c(1100, 1200), verbose = FALSE)
  expect_equal(y$Technology_Gap, c(400, 300))
  expect_equal(y$Extension_Gap, c(300, 300))
  expect_equal(y$Technology_Index, 100 * c(400, 300) / 1500)
  expect_equal(y$Increase_over_Actual_pct, 100 * c(300, 300) / c(800, 900))
})

test_that("nutrient use efficiency is vectorized and consistent", {
  n <- nue_calculate(yield_fertilized = c(1450, 1720, 1850),
                     yield_control = 1200, nutrient_applied = c(20, 40, 60),
                     nutrient_uptake_fert = c(48, 61, 68),
                     nutrient_uptake_ctrl = 35,
                     biomass_fert = c(4300, 4950, 5300), biomass_ctrl = 3600,
                     verbose = FALSE)
  expect_equal(nrow(n), 3)
  expect_equal(n$AE_kg_kg, n$PE_kg_kg * n$ARE_percent / 100)
  expect_equal(n$UE_kg_kg, (c(4300, 4950, 5300) - 3600) / c(20, 40, 60))
  expect_equal(n$IE_kg_kg, c(1450, 1720, 1850) / c(48, 61, 68))
  expect_error(nue_calculate(1500, 1200, 0, verbose = FALSE))
})

test_that("response models find the known optimum", {
  dose <- c(0, 20, 40, 60, 80, 100, 120)
  yq <- 1000 + 12 * dose - 0.06 * dose^2
  q <- nutrient_response(dose, yq, model = "quadratic", price_output = 50,
                         price_nutrient = 15, verbose = FALSE)
  expect_equal(q$dose_max_yield, 100, tolerance = 1e-6)
  expect_equal(q$dose_economic, (15 / 50 - 12) / (2 * -0.06),
               tolerance = 1e-6)
  ylp <- 900 + 10 * pmin(dose, 70)
  lp <- nutrient_response(dose, ylp, model = "linear_plateau",
                          verbose = FALSE)
  expect_equal(unname(lp$coefficients["join_point"]), 70, tolerance = 1e-3)
  ym <- 2000 - 1200 * exp(-0.03 * dose)
  mi <- nutrient_response(dose, ym + c(1, -1, 1, -1, 1, -1, 1),
                          model = "mitscherlich", verbose = FALSE)
  expect_equal(unname(mi$coefficients["A"]), 2000, tolerance = 1e-2)
  all <- nutrient_response(dose, yq, model = "all", verbose = FALSE)
  expect_s3_class(all, "aridagri_response_set")
  expect_true(all$best %in% all$comparison$Model)
  conv <- nutrient_response(dose, 1000 + 0.05 * dose^2, verbose = FALSE)
  expect_true(is.na(conv$dose_max_yield))
})

test_that("economic analysis is vectorized", {
  e <- economic_analysis(yield = c(1050, 1200), price = 65, cost_fixed = 15000,
                         cost_variable = c(11000, 12000),
                         byproduct_yield = c(1600, 1800),
                         byproduct_price = 5, verbose = FALSE)
  expect_equal(e$Gross_Return, c(1050, 1200) * 65 + c(1600, 1800) * 5)
  expect_equal(e$BC_Ratio, e$Gross_Return / e$Total_Cost)
})

test_that("path analysis stops on collinear predictors", {
  set.seed(3)
  d <- data.frame(y = rnorm(20), a = rnorm(20))
  d$b <- 2 * d$a
  expect_error(path_analysis(d, "y", c("a", "b"), verbose = FALSE),
               "correlated")
})
