
test_that("thermal_indices returns correct structure", {
  result <- thermal_indices(tmax = c(32, 34, 35), tmin = c(18, 20, 22), base_temp = 10, verbose = FALSE)
  expect_s3_class(result, "aridagri_thermal")
  expect_true("gdd_total" %in% names(result))
})

test_that("GDD calculation is correct", {
  result <- thermal_indices(tmax = c(30), tmin = c(20), base_temp = 10, verbose = FALSE)
  expect_equal(result$gdd_total, 15)
})

test_that("harvest_index returns correct structure", {
  result <- harvest_index(grain_yield = 2500, straw_yield = 4000, verbose = FALSE)
  expect_s3_class(result, "data.frame")
  expect_true("HI_percent" %in% names(result))
})

test_that("yield_gap_analysis returns correct structure", {
  result <- yield_gap_analysis(actual_yield = 2000, potential_yield = 4000, verbose = FALSE)
  expect_s3_class(result, "data.frame")
})

test_that("economic_indices returns correct BC ratio", {
  result <- economic_indices(gross_return = 80000, total_cost = 40000, verbose = FALSE)
  expect_s3_class(result, "data.frame")
  expect_true("BC_Ratio" %in% names(result))
  expect_equal(result$BC_Ratio, 2.0)
})

