
test_that("anova_factorial produces correct output", {
  set.seed(303)
  data <- expand.grid(fert = c("F1", "F2", "F3"), variety = c("V1", "V2"), rep = 1:4)
  data$yield <- rnorm(nrow(data), 1000, 100)
  result <- anova_factorial(data, response = "yield", factor1 = "fert", factor2 = "variety", replication = "rep", verbose = FALSE)
  expect_type(result, "list")
  expect_true("anova_table" %in% names(result))
})

test_that("anova_latin produces correct output", {
  set.seed(404)
  data <- data.frame(
    trt = c("A","B","C","D","E","B","C","D","E","A","C","D","E","A","B","D","E","A","B","C","E","A","B","C","D"),
    row = rep(1:5, each = 5), col = rep(1:5, times = 5), yield = rnorm(25, 500, 50))
  result <- anova_latin(data, response = "yield", treatment = "trt", row = "row", column = "col", verbose = FALSE)
  expect_type(result, "list")
  expect_true("anova_table" %in% names(result))
})

