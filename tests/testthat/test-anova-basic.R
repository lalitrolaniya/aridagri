
test_that("anova_crd produces correct output structure", {
  set.seed(123)
  data <- data.frame(
    treatment = rep(c("T1", "T2", "T3", "T4"), each = 5),
    yield = c(rnorm(5, 1200, 50), rnorm(5, 1350, 60),
              rnorm(5, 1100, 55), rnorm(5, 1450, 65))
  )
  result <- anova_crd(data, response = "yield", treatment = "treatment", verbose = FALSE)
  expect_type(result, "list")
  expect_true("anova_table" %in% names(result))
  expect_true("grand_mean" %in% names(result))
  expect_true("cv" %in% names(result))
  expect_true("se_mean" %in% names(result))
})

test_that("anova_crd returns correct degrees of freedom", {
  set.seed(123)
  data <- data.frame(
    treatment = rep(c("T1", "T2", "T3"), each = 4),
    yield = rnorm(12, 100, 10)
  )
  result <- anova_crd(data, response = "yield", treatment = "treatment", verbose = FALSE)
  expect_equal(result$anova_table[["Df"]][1], 2)
  expect_equal(result$anova_table[["Df"]][2], 9)
})

test_that("anova_rbd produces correct output", {
  set.seed(456)
  data <- data.frame(
    variety = rep(paste0("V", 1:5), times = 4),
    block = rep(1:4, each = 5),
    yield = rnorm(20, 2000, 200)
  )
  result <- anova_rbd(data, response = "yield", treatment = "variety", block = "block", verbose = FALSE)
  expect_type(result, "list")
  expect_true("anova_table" %in% names(result))
})

test_that("anova_rbd returns correct degrees of freedom", {
  set.seed(456)
  data <- data.frame(
    trt = rep(c("A", "B", "C", "D"), times = 3),
    block = rep(1:3, each = 4),
    yield = rnorm(12, 50, 5)
  )
  result <- anova_rbd(data, response = "yield", treatment = "trt", block = "block", verbose = FALSE)
  expect_equal(result$anova_table[["Df"]][1], 2)
  expect_equal(result$anova_table[["Df"]][2], 3)
  expect_equal(result$anova_table[["Df"]][3], 6)
})

test_that("anova_crd stops with invalid columns", {
  data <- data.frame(x = 1:10, y = rnorm(10))
  expect_error(anova_crd(data, response = "yield", treatment = "trt"))
})

test_that("anova_rbd stops with invalid columns", {
  data <- data.frame(x = 1:10, y = rnorm(10))
  expect_error(anova_rbd(data, response = "yield", treatment = "trt", block = "blk"))
})

