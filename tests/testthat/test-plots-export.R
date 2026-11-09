test_that("ggplot2 figures are built for the main result types", {
  skip_if_not_installed("ggplot2")
  data(mothbean_rbd, package = "aridagri")
  data(clusterbean_spd, package = "aridagri")
  data(pearlmillet_met, package = "aridagri")
  data(mungbean_repeated, package = "aridagri")
  r1 <- anova_rbd(mothbean_rbd, "seed_yield", "treatment", "block",
                  verbose = FALSE)
  expect_s3_class(plot_means(r1), "ggplot")
  expect_s3_class(plot(r1), "ggplot")
  r2 <- anova_spd(clusterbean_spd, "seed_yield", "irrigation", "nitrogen",
                  "rep", verbose = FALSE)
  expect_s3_class(plot_interaction(r2), "ggplot")
  expect_s3_class(plot_means(r2, which = "nitrogen"), "ggplot")
  st <- stability_analysis(pearlmillet_met, "genotype", "environment", "rep",
                           "grain_yield", verbose = FALSE)
  for (tp in c("eberhart", "ammi1", "ammi2")) {
    expect_s3_class(plot_stability(st, tp), "ggplot")
  }
  gg <- gge_biplot(pearlmillet_met, "genotype", "environment", "grain_yield",
                   verbose = FALSE)
  expect_s3_class(plot(gg), "ggplot")
  expect_s3_class(plot(gg, type = "mean_stability"), "ggplot")
  rm <- anova_repeated(mungbean_repeated, "spad", "treatment", "das",
                       block = "block", covariance = FALSE, verbose = FALSE)
  expect_s3_class(plot(rm), "ggplot")
  nr <- nutrient_response(c(0, 20, 40, 60, 80), c(900, 1100, 1250, 1320, 1330),
                          model = "all", verbose = FALSE)
  expect_s3_class(plot_response(nr, c(0, 20, 40, 60, 80),
                                c(900, 1100, 1250, 1320, 1330)), "ggplot")
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  for (p in list(plot(r1), plot(r2), plot_stability(st, "ammi2"), plot(gg))) {
    expect_silent(print(p))
  }
})

test_that("base-graphics arid_plot draws without error", {
  data(mothbean_rbd, package = "aridagri")
  r1 <- anova_rbd(mothbean_rbd, "seed_yield", "treatment", "block",
                  verbose = FALSE)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(arid_plot(r1))
  r2 <- anova_spd(aridagri::clusterbean_spd, "seed_yield", "irrigation",
                  "nitrogen", "rep", verbose = FALSE)
  expect_silent(arid_plot(r2))
})

test_that("export_results writes CSV and Excel files", {
  data(mothbean_rbd, package = "aridagri")
  r1 <- anova_rbd(mothbean_rbd, "seed_yield", "treatment", "block",
                  verbose = FALSE)
  f <- tempfile(fileext = ".csv")
  export_results(r1, f, format = "csv", verbose = FALSE)
  expect_true(file.exists(f))
  tb <- anova_traits(mothbean_rbd, c("seed_yield", "plant_height"),
                     treatment = "treatment", block = "block",
                     verbose = FALSE)
  f2 <- tempfile(fileext = ".csv")
  export_results(tb, f2, format = "csv", verbose = FALSE)
  expect_equal(nrow(utils::read.csv(f2)), nrow(tb$table))
  skip_if_not_installed("writexl")
  f3 <- tempfile(fileext = ".xlsx")
  export_results(r1, f3, verbose = FALSE)
  expect_true(file.exists(f3))
})
