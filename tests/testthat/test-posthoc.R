
test_that("anova_crd with posthoc lsd works", {
  set.seed(1)
  data <- data.frame(trt = rep(c("A","B","C","D"), each=5), yield = c(rnorm(5,20,2), rnorm(5,25,2), rnorm(5,22,2), rnorm(5,30,2)))
  result <- anova_crd(data, response="yield", treatment="trt", posthoc="lsd", verbose=FALSE)
  expect_type(result, "list")
  expect_true("posthoc" %in% names(result))
})

test_that("anova_rbd with posthoc duncan works", {
  set.seed(2)
  data <- data.frame(trt = rep(c("A","B","C"), times=4), block = rep(1:4, each=3), yield = rnorm(12,50,5))
  result <- anova_rbd(data, response="yield", treatment="trt", block="block", posthoc="duncan", verbose=FALSE)
  expect_type(result, "list")
  expect_true("posthoc" %in% names(result))
})

test_that("anova_rbd with posthoc tukey works", {
  set.seed(3)
  data <- data.frame(trt = rep(c("A","B","C"), times=4), block = rep(1:4, each=3), yield = rnorm(12,50,5))
  result <- anova_rbd(data, response="yield", treatment="trt", block="block", posthoc="tukey", verbose=FALSE)
  expect_type(result, "list")
})

test_that("anova_crd with posthoc snk works", {
  set.seed(4)
  data <- data.frame(trt = rep(c("A","B","C"), each=5), yield = rnorm(15,100,10))
  result <- anova_crd(data, response="yield", treatment="trt", posthoc="snk", verbose=FALSE)
  expect_type(result, "list")
})

test_that("anova_crd with posthoc scheffe works", {
  set.seed(5)
  data <- data.frame(trt = rep(c("A","B","C"), each=5), yield = rnorm(15,100,10))
  result <- anova_crd(data, response="yield", treatment="trt", posthoc="scheffe", verbose=FALSE)
  expect_type(result, "list")
})

test_that("anova_crd with posthoc bonferroni works", {
  set.seed(6)
  data <- data.frame(trt = rep(c("A","B","C"), each=5), yield = rnorm(15,100,10))
  result <- anova_crd(data, response="yield", treatment="trt", posthoc="bonferroni", verbose=FALSE)
  expect_type(result, "list")
})

test_that("anova_crd with posthoc dunnett works", {
  set.seed(7)
  data <- data.frame(trt = rep(c("A","B","C"), each=5), yield = rnorm(15,100,10))
  result <- anova_crd(data, response="yield", treatment="trt", posthoc="dunnett", verbose=FALSE)
  expect_type(result, "list")
})

