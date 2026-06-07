
test_that("anova_spd produces correct output", {
  set.seed(789)
  data <- expand.grid(irrigation = c("I1","I2","I3"), variety = c("V1","V2","V3","V4"), rep = 1:3)
  data$yield <- rnorm(nrow(data), 1500, 200)
  result <- anova_spd(data, response="yield", main_plot="irrigation", sub_plot="variety", replication="rep", verbose=FALSE)
  expect_type(result, "list")
  expect_true("anova_table" %in% names(result))
  expect_true("cv_a" %in% names(result))
  expect_true("cv_b" %in% names(result))
})

test_that("anova_spd returns correct error df", {
  set.seed(1)
  data <- expand.grid(main = c("A","B"), sub = c("X","Y","Z"), rep = 1:4)
  data$yield <- rnorm(nrow(data), 100, 15)
  result <- anova_spd(data, response="yield", main_plot="main", sub_plot="sub", replication="rep", verbose=FALSE)
  expect_equal(result$df_error_a, 3)
})

test_that("anova_spd_ab_main produces correct output", {
  set.seed(111)
  data <- expand.grid(f1 = c("A","B"), f2 = c("X","Y"), sub = c("S1","S2","S3"), rep = 1:3)
  data$yield <- rnorm(nrow(data), 200, 25)
  result <- anova_spd_ab_main(data, response="yield", main_factor1="f1", main_factor2="f2", sub_plot="sub", replication="rep", verbose=FALSE)
  expect_type(result, "list")
  expect_true("anova_table" %in% names(result))
})

test_that("anova_spd_c_main_ab_sub produces correct output", {
  set.seed(222)
  data <- expand.grid(main = c("M1","M2","M3"), sf1 = c("A","B"), sf2 = c("X","Y"), rep = 1:3)
  data$yield <- rnorm(nrow(data), 300, 30)
  result <- anova_spd_c_main_ab_sub(data, response="yield", main_plot="main", sub_factor1="sf1", sub_factor2="sf2", replication="rep", verbose=FALSE)
  expect_type(result, "list")
  expect_true("anova_table" %in% names(result))
})

test_that("anova_spd_ab_cd produces correct output", {
  set.seed(333)
  data <- expand.grid(mf1 = c("A","B"), mf2 = c("X","Y"), sf1 = c("P","Q"), sf2 = c("R","S"), rep = 1:3)
  data$yield <- rnorm(nrow(data), 400, 40)
  result <- anova_spd_ab_cd(data, response="yield", main_factor1="mf1", main_factor2="mf2", sub_factor1="sf1", sub_factor2="sf2", replication="rep", verbose=FALSE)
  expect_type(result, "list")
})

test_that("anova_spd_pooled produces correct output", {
  set.seed(444)
  data <- expand.grid(main = c("A","B"), sub = c("X","Y","Z"), env = paste0("E",1:3), rep = 1:2)
  data$yield <- rnorm(nrow(data), 500, 50)
  result <- anova_spd_pooled(data, response="yield", main_plot="main", sub_plot="sub", environment="env", replication="rep", verbose=FALSE)
  expect_type(result, "list")
  expect_true("anova_table" %in% names(result))
})

test_that("anova_sspd produces correct output", {
  set.seed(101)
  data <- expand.grid(main = c("M1","M2"), sub = c("S1","S2","S3"), subsub = c("SS1","SS2"), rep = 1:3)
  data$yield <- rnorm(nrow(data), 200, 30)
  result <- anova_sspd(data, response="yield", main_plot="main", sub_plot="sub", sub_sub_plot="subsub", replication="rep", verbose=FALSE)
  expect_type(result, "list")
  expect_true("anova_table" %in% names(result))
})

test_that("anova_sspd_pooled produces correct output", {
  set.seed(555)
  data <- expand.grid(main = c("A","B"), sub = c("X","Y"), subsub = c("P","Q"), env = paste0("E",1:2), rep = 1:2)
  data$yield <- rnorm(nrow(data), 600, 60)
  result <- anova_sspd_pooled(data, response="yield", main_plot="main", sub_plot="sub", sub_sub_plot="subsub", environment="env", replication="rep", verbose=FALSE)
  expect_type(result, "list")
})

test_that("anova_strip produces correct output", {
  set.seed(202)
  data <- expand.grid(tillage = c("ZT","CT","MT"), nitrogen = c("N0","N60","N120"), rep = 1:3)
  data$yield <- rnorm(nrow(data), 2000, 250)
  result <- anova_strip(data, response="yield", horizontal_factor="tillage", vertical_factor="nitrogen", replication="rep", verbose=FALSE)
  expect_type(result, "list")
  expect_true("anova_table" %in% names(result))
})

