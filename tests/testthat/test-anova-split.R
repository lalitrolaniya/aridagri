
test_that("anova_spd produces correct output structure", {
  set.seed(789)
  data <- expand.grid(irrigation = c("I1", "I2", "I3"), variety = c("V1", "V2", "V3", "V4"), rep = 1:3)
  data$yield <- rnorm(nrow(data), 1500, 200)
  result <- anova_spd(data, response = "yield", main_plot = "irrigation", sub_plot = "variety", replication = "rep", verbose = FALSE)
  expect_type(result, "list")
  expect_true("anova_table" %in% names(result))
  expect_true("cv_a" %in% names(result))
  expect_true("cv_b" %in% names(result))
})

test_that("anova_spd returns correct error df", {
  set.seed(789)
  data <- expand.grid(main = c("A", "B"), sub = c("X", "Y", "Z"), rep = 1:4)
  data$yield <- rnorm(nrow(data), 100, 15)
  result <- anova_spd(data, response = "yield", main_plot = "main", sub_plot = "sub", replication = "rep", verbose = FALSE)
  expect_equal(result$df_error_a, 3)
})

test_that("anova_sspd produces correct output", {
  set.seed(101)
  data <- expand.grid(main = c("M1", "M2"), sub = c("S1", "S2", "S3"), subsub = c("SS1", "SS2"), rep = 1:3)
  data$yield <- rnorm(nrow(data), 200, 30)
  result <- anova_sspd(data, response = "yield", main_plot = "main", sub_plot = "sub", sub_sub_plot = "subsub", replication = "rep", verbose = FALSE)
  expect_type(result, "list")
  expect_true("anova_table" %in% names(result))
})

test_that("anova_strip produces correct output", {
  set.seed(202)
  data <- expand.grid(tillage = c("ZT", "CT", "MT"), nitrogen = c("N0", "N60", "N120"), rep = 1:3)
  data$yield <- rnorm(nrow(data), 2000, 250)
  result <- anova_strip(data, response = "yield", horizontal_factor = "tillage", vertical_factor = "nitrogen", replication = "rep", verbose = FALSE)
  expect_type(result, "list")
  expect_true("anova_table" %in% names(result))
})

