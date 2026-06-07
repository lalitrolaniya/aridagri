
test_that("anova_augmented produces correct output", {
  set.seed(666)
  checks <- expand.grid(genotype = c("C1","C2","C3"), block = 1:4)
  tests <- data.frame(genotype = paste0("G",1:12), block = rep(1:4, each=3))
  data <- rbind(checks, tests)
  data$yield <- rnorm(nrow(data), 1000, 120)
  result <- anova_augmented(data, response="yield", genotype="genotype", block="block", check_names=c("C1","C2","C3"), verbose=FALSE)
  expect_type(result, "list")
})

test_that("anova_alpha_lattice produces correct output", {
  set.seed(777)
  data <- expand.grid(genotype = paste0("G",1:9), rep = 1:3)
  data$block <- rep(rep(1:3, each=3), 3)
  data$yield <- rnorm(nrow(data), 800, 80)
  result <- anova_alpha_lattice(data, response="yield", genotype="genotype", replication="rep", block="block", verbose=FALSE)
  expect_type(result, "list")
})

