test_that("anova_crd handles unequal replication exactly", {
  set.seed(1)
  data <- data.frame(trt = c(rep("A", 5), rep("B", 3), rep("C", 4)),
                     yield = rnorm(12, 50, 5))
  res <- anova_crd(data, response = "yield", treatment = "trt",
                   verbose = FALSE)
  expect_false(res$balanced)
  expect_equal(unname(res$replications), c(5, 3, 4))
  # SE of each mean uses its own replication
  se <- res$treatment_means$SE[match(c("A", "B", "C"),
                                     res$treatment_means$Treatment)]
  expect_equal(se, sqrt(res$ms_error / c(5, 3, 4)))
})

test_that("anova_crd no warning on balanced data", {
  data <- data.frame(trt = rep(c("A", "B", "C"), each = 4),
                     yield = rnorm(12, 50, 5))
  expect_no_warning(anova_crd(data, response = "yield", treatment = "trt",
                              verbose = FALSE))
})

test_that("anova_rbd analyses missing plots by least squares", {
  set.seed(1)
  data <- data.frame(trt = rep(c("A", "B", "C"), times = 4),
                     block = rep(1:4, each = 3), yield = rnorm(12, 50, 5))
  data_unbal <- data[-5, ]
  res <- anova_rbd(data_unbal, response = "yield", treatment = "trt",
                   block = "block", verbose = FALSE)
  expect_equal(res$missing_plots, 1)
  expect_false(res$balanced)
  expect_equal(res$df_error, 5)
  expect_true(is.na(res$relative_efficiency))
  utils::capture.output(
    expect_message(anova_rbd(data_unbal, response = "yield", treatment = "trt",
                             block = "block", verbose = TRUE), "missing plot"))
})

test_that("anova_rbd no warning on balanced data", {
  set.seed(1)
  data <- data.frame(trt = rep(c("A", "B", "C"), times = 4),
                     block = rep(1:4, each = 3), yield = rnorm(12, 50, 5))
  expect_no_warning(anova_rbd(data, response = "yield", treatment = "trt",
                              block = "block", verbose = FALSE))
})

test_that("anova_spd stops on incomplete data", {
  set.seed(1)
  data <- expand.grid(main = c("M1", "M2"), sub = c("S1", "S2", "S3"),
                      rep = 1:3)
  data$yield <- rnorm(nrow(data), 100, 10)
  data_unbal <- data[-3, ]
  expect_error(anova_spd(data_unbal, response = "yield", main_plot = "main",
                         sub_plot = "sub", replication = "rep",
                         verbose = FALSE), "complete data")
})

test_that("anova_spd no warning on balanced data", {
  set.seed(1)
  data <- expand.grid(main = c("M1", "M2"), sub = c("S1", "S2", "S3"),
                      rep = 1:3)
  data$yield <- rnorm(nrow(data), 100, 10)
  expect_no_warning(anova_spd(data, response = "yield", main_plot = "main",
                              sub_plot = "sub", replication = "rep",
                              verbose = FALSE))
})

test_that("missing responses are dropped with a message", {
  data <- data.frame(trt = rep(c("A", "B", "C"), each = 4),
                     yield = c(NA, rnorm(11, 50, 5)))
  utils::capture.output(
    expect_message(anova_crd(data, "yield", "trt", verbose = TRUE),
                   "missing"))
  res <- anova_crd(data, "yield", "trt", verbose = FALSE)
  expect_equal(sum(res$replications), 11)
})
