
test_that("stability_analysis method=all produces all results", {
  set.seed(42)
  data <- expand.grid(g = paste0("G",1:5), e = paste0("E",1:4), r = 1:3)
  data$y <- rnorm(nrow(data), 1500, 250)
  result <- stability_analysis(data, "g", "e", "r", "y", method="all", verbose=FALSE)
  expect_type(result, "list")
  expect_true("eberhart" %in% names(result))
  expect_true("finlay" %in% names(result))
  expect_true("shukla" %in% names(result))
  expect_true("wricke" %in% names(result))
  expect_true("ammi" %in% names(result))
  expect_true("cv_stability" %in% names(result))
  expect_true("superiority" %in% names(result))
  expect_true("integrated" %in% names(result))
})

test_that("integrated ranking has all 7 method columns", {
  set.seed(42)
  data <- expand.grid(g = paste0("G",1:5), e = paste0("E",1:4), r = 1:3)
  data$y <- rnorm(nrow(data), 1500, 250)
  result <- stability_analysis(data, "g", "e", "r", "y", method="all", verbose=FALSE)
  expected_cols <- c("Genotype","Mean","ER","AMMI","FW","Shukla","Wricke","CV","LB","Rank")
  expect_true(all(expected_cols %in% names(result$integrated)))
  expect_equal(ncol(result$integrated), 10)
})

test_that("Shukla and Wricke produce identical rankings", {
  set.seed(42)
  data <- expand.grid(g = paste0("G",1:6), e = paste0("E",1:5), r = 1:3)
  data$y <- rnorm(nrow(data), 1500, 250)
  result <- stability_analysis(data, "g", "e", "r", "y", method="all", verbose=FALSE)
  expect_equal(result$integrated$Shukla, result$integrated$Wricke)
})

test_that("single method skips integrated ranking", {
  set.seed(42)
  data <- expand.grid(g = paste0("G",1:5), e = paste0("E",1:4), r = 1:3)
  data$y <- rnorm(nrow(data), 1500, 250)
  result <- stability_analysis(data, "g", "e", "r", "y", method="eberhart", verbose=FALSE)
  expect_true("eberhart" %in% names(result))
  expect_false("integrated" %in% names(result))
})

test_that("stability_analysis works with many genotypes", {
  set.seed(99)
  data <- expand.grid(g = paste0("G",1:15), e = paste0("E",1:6), r = 1:3)
  data$y <- rnorm(nrow(data), 2000, 300)
  result <- stability_analysis(data, "g", "e", "r", "y", method="all", verbose=FALSE)
  expect_equal(nrow(result$integrated), 15)
  expect_equal(max(result$integrated$Rank), 15)
})

