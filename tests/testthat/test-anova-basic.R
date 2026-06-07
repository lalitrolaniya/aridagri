
test_that("anova_crd produces correct output", {
  set.seed(123)
  data <- data.frame(treatment = rep(c("T1","T2","T3","T4"), each=5),
                     yield = c(rnorm(5,1200,50), rnorm(5,1350,60), rnorm(5,1100,55), rnorm(5,1450,65)))
  result <- anova_crd(data, response="yield", treatment="treatment", verbose=FALSE)
  expect_type(result, "list")
  expect_true("anova_table" %in% names(result))
  expect_true("grand_mean" %in% names(result))
  expect_true("cv" %in% names(result))
  expect_true("se_mean" %in% names(result))
})
test_that("anova_crd returns correct df", {
  set.seed(1)
  data <- data.frame(treatment = rep(c("T1","T2","T3"), each=4), yield = rnorm(12,100,10))
  result <- anova_crd(data, response="yield", treatment="treatment", verbose=FALSE)
  expect_equal(result$anova_table[["Df"]][1], 2)
  expect_equal(result$anova_table[["Df"]][2], 9)
})
test_that("anova_crd stops with invalid columns", {
  data <- data.frame(x=1:10, y=rnorm(10))
  expect_error(anova_crd(data, response="yield", treatment="trt"))
})
test_that("anova_rbd produces correct output", {
  set.seed(456)
  data <- data.frame(variety = rep(paste0("V",1:5), times=4), block = rep(1:4, each=5), yield = rnorm(20,2000,200))
  result <- anova_rbd(data, response="yield", treatment="variety", block="block", verbose=FALSE)
  expect_type(result, "list")
  expect_true("anova_table" %in% names(result))
  expect_true("cv" %in% names(result))
})
test_that("anova_rbd returns correct df", {
  set.seed(1)
  data <- data.frame(trt = rep(c("A","B","C","D"), times=3), block = rep(1:3, each=4), yield = rnorm(12,50,5))
  result <- anova_rbd(data, response="yield", treatment="trt", block="block", verbose=FALSE)
  expect_equal(result$anova_table[["Df"]][1], 2)
  expect_equal(result$anova_table[["Df"]][2], 3)
  expect_equal(result$anova_table[["Df"]][3], 6)
})
test_that("anova_rbd stops with invalid columns", {
  data <- data.frame(x=1:10, y=rnorm(10))
  expect_error(anova_rbd(data, response="yield", treatment="trt", block="blk"))
})
test_that("anova_rbd_pooled produces correct output", {
  set.seed(789)
  data <- expand.grid(trt = paste0("T",1:4), env = paste0("E",1:3), block = 1:2)
  data$yield <- rnorm(nrow(data), 1500, 200)
  result <- anova_rbd_pooled(data, response="yield", treatment="trt", environment="env", block="block", verbose=FALSE)
  expect_type(result, "list")
  expect_true("pooled_anova" %in% names(result))
})
test_that("anova_latin produces correct output", {
  set.seed(404)
  data <- data.frame(
    trt = c("A","B","C","D","E","B","C","D","E","A","C","D","E","A","B","D","E","A","B","C","E","A","B","C","D"),
    row = rep(1:5, each=5), col = rep(1:5, times=5), yield = rnorm(25,500,50))
  result <- anova_latin(data, response="yield", treatment="trt", row="row", column="col", verbose=FALSE)
  expect_type(result, "list")
  expect_true("anova_table" %in% names(result))
})
test_that("anova_factorial produces correct output", {
  set.seed(303)
  data <- expand.grid(fert = c("F1","F2","F3"), variety = c("V1","V2"), rep = 1:4)
  data$yield <- rnorm(nrow(data), 1000, 100)
  result <- anova_factorial(data, response="yield", factor1="fert", factor2="variety", replication="rep", verbose=FALSE)
  expect_type(result, "list")
  expect_true("anova_table" %in% names(result))
})
test_that("anova_factorial_3way produces correct output", {
  set.seed(505)
  data <- expand.grid(f1 = c("A","B"), f2 = c("X","Y"), f3 = c("P","Q"), rep = 1:3)
  data$yield <- rnorm(nrow(data), 100, 15)
  result <- anova_factorial_3way(data, response="yield", factor1="f1", factor2="f2", factor3="f3", replication="rep", verbose=FALSE)
  expect_type(result, "list")
  expect_true("anova_table" %in% names(result))
})

