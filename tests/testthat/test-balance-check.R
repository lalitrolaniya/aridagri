
test_that("anova_crd warns on unbalanced data", {
  data <- data.frame(trt = c(rep("A",5), rep("B",3), rep("C",4)), yield = rnorm(12,50,5))
  expect_warning(anova_crd(data, response="yield", treatment="trt", verbose=FALSE), "Unbalanced design detected")
})

test_that("anova_crd no warning on balanced data", {
  data <- data.frame(trt = rep(c("A","B","C"), each=4), yield = rnorm(12,50,5))
  expect_no_warning(anova_crd(data, response="yield", treatment="trt", verbose=FALSE))
})

test_that("anova_rbd warns on unbalanced data", {
  set.seed(1)
  data <- data.frame(trt = rep(c("A","B","C"), times=4), block = rep(1:4, each=3), yield = rnorm(12,50,5))
  data_unbal <- data[-5, ]
  expect_warning(anova_rbd(data_unbal, response="yield", treatment="trt", block="block", verbose=FALSE), "Unbalanced design detected")
})

test_that("anova_rbd no warning on balanced data", {
  set.seed(1)
  data <- data.frame(trt = rep(c("A","B","C"), times=4), block = rep(1:4, each=3), yield = rnorm(12,50,5))
  expect_no_warning(anova_rbd(data, response="yield", treatment="trt", block="block", verbose=FALSE))
})

test_that("anova_spd warns on unbalanced data", {
  set.seed(1)
  data <- expand.grid(main = c("M1","M2"), sub = c("S1","S2","S3"), rep = 1:3)
  data$yield <- rnorm(nrow(data), 100, 10)
  data_unbal <- data[-3, ]
  expect_warning(anova_spd(data_unbal, response="yield", main_plot="main", sub_plot="sub", replication="rep", verbose=FALSE), "Unbalanced design detected")
})

test_that("anova_spd no warning on balanced data", {
  set.seed(1)
  data <- expand.grid(main = c("M1","M2"), sub = c("S1","S2","S3"), rep = 1:3)
  data$yield <- rnorm(nrow(data), 100, 10)
  expect_no_warning(anova_spd(data, response="yield", main_plot="main", sub_plot="sub", replication="rep", verbose=FALSE))
})

