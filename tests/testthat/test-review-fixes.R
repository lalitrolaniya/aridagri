test_that("anova_repeated runs on incomplete data and with posthoc none", {
  data(mungbean_repeated, package = "aridagri")
  d <- mungbean_repeated
  d$spad[5] <- NA
  r <- suppressWarnings(anova_repeated(d, "spad", "treatment", "das",
                                       subject = "plot", block = "block",
                                       verbose = FALSE))
  expect_false(r$complete)
  expect_null(r$anova_table)
  expect_true(r$best_structure %in% c("CS", "AR1", "UN"))
  r2 <- suppressWarnings(anova_repeated(d, "spad", "treatment", "das",
                                        subject = "plot", posthoc = "none",
                                        verbose = FALSE))
  expect_s3_class(r2, "aridagri_repeated")
  r3 <- anova_repeated(mungbean_repeated[mungbean_repeated$plot != "R4-4", ],
                       "spad", "treatment", "das", subject = "plot",
                       verbose = FALSE)
  expect_false(is.null(r3$anova_table))
})

test_that("anova_repeated handles non-syntactic names and tibbles", {
  data(mungbean_repeated, package = "aridagri")
  d <- mungbean_repeated
  names(d)[3] <- "Days after sowing"
  r <- anova_repeated(d, "spad", "treatment", "Days after sowing",
                      block = "block", covariance = FALSE, verbose = FALSE)
  expect_equal(nrow(r$corrected), 2)
  skip_if_not_installed("tibble")
  r2 <- anova_repeated(tibble::as_tibble(mungbean_repeated), "spad",
                       "treatment", "das", block = "block",
                       covariance = FALSE, verbose = FALSE)
  expect_s3_class(r2, "aridagri_repeated")
})

test_that("Dunnett averages over interactions and uses the given error", {
  skip_if_not_installed("multcomp")
  data(clusterbean_spd, package = "aridagri")
  fit <- aov(seed_yield ~ rep + irrigation * nitrogen, clusterbean_spd)
  ph <- perform_posthoc(fit, clusterbean_spd, "seed_yield", "irrigation",
                        mse = 3200.82, df_error = 4,
                        posthoc = c("lsd", "dunnett"), verbose = FALSE)
  m <- tapply(clusterbean_spd$seed_yield, clusterbean_spd$irrigation, mean)
  expect_equal(ph$dunnett$estimate, as.numeric(m[2:3] - m[1]),
               tolerance = 1e-8)
  expect_equal(ph$dunnett$SE, rep(sqrt(2 * 3200.82 / 12), 2),
               tolerance = 1e-6)
  expect_equal(ph$dunnett$SE[1], ph$lsd$statistics$SE_diff, tolerance = 1e-8)
})

test_that("response models are compared with the right parameter count", {
  x <- c(0, 20, 40, 60, 80, 100)
  y <- c(1100, 1350, 1520, 1610, 1650, 1660)
  r <- nutrient_response(x, y, "all", 60, 15, verbose = FALSE)
  aq <- r$comparison$AIC[r$comparison$Model == "quadratic"]
  am <- r$comparison$AIC[r$comparison$Model == "mitscherlich"]
  ref <- AIC(lm(y ~ x + I(x^2))) - AIC(nls(y ~ SSasymp(x, Asym, R0, lrc)))
  expect_equal(aq - am, ref, tolerance = 1e-6)
  mi <- r$models$mitscherlich
  expect_match(mi$note, "extrapolation")
})

test_that("unbalanced factorial means are least-squares means", {
  skip_if_not_installed("emmeans")
  fa <- expand.grid(rep = factor(1:4), N = c("N0", "N40", "N80"),
                    P = c("P0", "P30"))
  set.seed(1)
  fa$y <- 900 + 3 * as.numeric(sub("N", "", fa$N)) +
    2 * as.numeric(sub("P", "", fa$P)) + rnorm(nrow(fa), 0, 40)
  fa <- fa[-c(2, 9, 17), ]
  f <- suppressWarnings(anova_factorial(fa, "y", "N", "P", "rep",
                                        verbose = FALSE))
  em <- suppressMessages(summary(emmeans::emmeans(lm(y ~ rep + N * P, fa),
                                                  "N")))
  expect_equal(f$mean_a$Mean, em$emmean, tolerance = 1e-8)
  expect_equal(f$mean_a$SE, em$SE, tolerance = 1e-8)
})

test_that("non-estimable least-squares means are dropped with a warning", {
  fe <- expand.grid(rep = factor(1:3), N = c("N0", "N40", "N80"),
                    P = c("P0", "P30"))
  set.seed(2)
  fe$y <- rnorm(nrow(fe), 100, 5)
  fe <- fe[!(fe$N == "N80" & fe$P == "P30"), ]
  expect_warning(ph <- perform_posthoc(aov(y ~ rep + N * P, fe), fe, "y",
                                       "N", verbose = FALSE),
                 "not estimable")
  expect_equal(nrow(ph$lsd$table), 2)
})

test_that("anova_traits works with non-syntactic factor names", {
  fa <- expand.grid(rep = factor(1:3), `N dose` = c("N0", "N40"),
                    `P dose` = c("P0", "P30"))
  set.seed(4)
  fa$y <- rnorm(nrow(fa), 100, 5)
  t <- anova_traits(fa, "y", design = "factorial", factor1 = "N dose",
                    factor2 = "P dose", block = "rep", verbose = FALSE)
  expect_equal(nrow(t$stats), 3)
  expect_false(anyNA(t$stats$F_value))
})

test_that("SNK and Duncan p-values agree with the letters", {
  cm <- compare_means(c(A = 3.3, B = 3, C = 0), 4, 1000, 4, "snk",
                      verbose = FALSE)
  expect_equal(cm$comparisons$significant, cm$comparisons$p_value < 0.05)
})

test_that("augmented genetic parameters use the block-adjusted test MS", {
  data(chickpea_augmented, package = "aridagri")
  au <- anova_augmented(chickpea_augmented, "seed_yield", "entry", "block",
                        c("C1", "C2", "C3", "C4"), verbose = FALSE)
  at <- au$anova_adjusted
  expect_equal(sum(at[c("  Tests vs Checks", "  Checks", "  Tests"),
                      "Sum Sq"]),
               at["Treatments (eliminating blocks)", "Sum Sq"])
  gp <- au$genetic_parameters
  vp <- gp$Value[gp$Parameter == "Phenotypic variance"]
  expect_equal(vp, at["  Tests", "Mean Sq"])
  expect_equal(gp$Value[gp$Parameter == "Genotypic variance"],
               vp - au$mse)
})
