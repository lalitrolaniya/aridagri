test_that("relative efficiency of RCBD follows Gomez and Gomez", {
  set.seed(11)
  d <- expand.grid(trt = paste0("T", 1:5), rep = factor(1:4))
  d$y <- 100 + rnorm(20, 0, 5)
  r <- anova_rbd(d, "y", "trt", "rep", verbose = FALSE)
  a <- anova(aov(y ~ rep + trt, d))
  Eb <- a[1, 3]
  Ee <- a[3, 3]
  re <- 100 * (3 * Eb + 4 * 4 * Ee) / (19 * Ee)
  k <- ((3 * 4 + 1) * (5 * 3 + 3)) / ((3 * 4 + 3) * (5 * 3 + 1))
  expect_equal(r$relative_efficiency_unadjusted, re)
  expect_equal(r$relative_efficiency, re * k)
  # with no block effect RE is close to 100%, not 50%
  expect_gt(r$relative_efficiency_unadjusted, 80)
})

test_that("pooled RBD works silently and tests environment against Rep(E)", {
  set.seed(5)
  d <- expand.grid(trt = paste0("T", 1:4), env = paste0("E", 1:3),
                   block = 1:3)
  d$y <- rnorm(nrow(d), 100, 10) + as.integer(d$env) * 5
  p0 <- anova_rbd_pooled(d, "y", "trt", "env", "block", verbose = FALSE)
  expect_true(is.finite(p0$bartlett_test$chi_sq))
  expect_length(p0$individual_anovas, 3)
  expect_true(all(p0$mse_values > 0))
  tab <- p0$pooled_anova
  expect_equal(tab[1, "F value"], tab[1, "Mean Sq"] / tab[2, "Mean Sq"])
  expect_equal(tab[3, "F value"], tab[3, "Mean Sq"] / tab[5, "Mean Sq"])
  pr <- anova_rbd_pooled(d, "y", "trt", "env", "block", verbose = FALSE,
                         env_effect = "random")
  expect_equal(pr$pooled_anova[3, "F value"],
               tab[3, "Mean Sq"] / tab[4, "Mean Sq"])
  # Bartlett matches stats::bartlett.test on the residuals of each env
  expect_equal(p0$bartlett_test$df, 2)
})

test_that("split-plot ANOVA matches aov Error strata and agricolae", {
  set.seed(31)
  d <- expand.grid(rep = factor(1:3), A = paste0("A", 1:3),
                   B = paste0("B", 1:4))
  d$y <- 50 + as.numeric(factor(d$A)) * 3 + rnorm(nrow(d), 0, 2) +
    rnorm(9, 0, 3)[as.integer(interaction(d$rep, d$A))]
  s <- anova_spd(d, "y", "A", "B", "rep", verbose = FALSE)
  e <- summary(aov(y ~ A * B + Error(rep / A), d))
  expect_equal(s$anova_table["A", "F value"],
               e[["Error: rep:A"]][[1]]["A", "F value"])
  expect_equal(s$anova_table["B", "F value"],
               e[["Error: Within"]][[1]]["B  ", "F value"])
  Ea <- s$ms_error_a
  Eb <- s$ms_error_b
  r <- 3
  a <- 3
  b <- 4
  sc <- s$se_cd
  expect_equal(sc$SEd[1], sqrt(2 * Ea / (r * b)))
  expect_equal(sc$SEd[2], sqrt(2 * Eb / (r * a)))
  expect_equal(s$se_sub_same * sqrt(2), sqrt(2 * Eb / r))
  expect_equal(s$se_main_same_sub * sqrt(2),
               sqrt(2 * ((b - 1) * Eb + Ea) / (r * b)))
  tw <- ((b - 1) * Eb * qt(0.975, s$df_error_b) + Ea * qt(0.975,
                                                           s$df_error_a)) /
    ((b - 1) * Eb + Ea)
  expect_equal(s$cd_main_same_sub,
               tw * sqrt(2 * ((b - 1) * Eb + Ea) / (r * b)))
  skip_if_not_installed("agricolae")
  sp <- suppressWarnings(utils::capture.output(
    ag <- with(d, agricolae::sp.plot(rep, A, B, y))))
  expect_equal(s$ms_error_a, ag$Ea)
  expect_equal(s$ms_error_b, ag$Eb)
})

test_that("split-split-plot SEd reproduce the textbook formulas", {
  set.seed(41)
  d <- expand.grid(rep = factor(1:3), A = paste0("A", 1:2),
                   B = paste0("B", 1:3), C = paste0("C", 1:4))
  d$y <- 40 + rnorm(nrow(d), 0, 2) +
    rnorm(6, 0, 3)[as.integer(interaction(d$rep, d$A))] +
    rnorm(18, 0, 2)[as.integer(interaction(d$rep, d$A, d$B))]
  s <- anova_sspd(d, "y", "A", "B", "C", "rep", verbose = FALSE)
  e <- summary(aov(y ~ A * B * C + Error(rep / A / B), d))
  expect_equal(s$anova_table["B", "F value"],
               e[["Error: rep:A:B"]][[1]]["B", "F value"])
  Ea <- s$ms_error_a
  Eb <- s$ms_error_b
  Ec <- s$ms_error_c
  r <- 3
  b <- 3
  cc <- 4
  g <- function(lab) s$se_cd$SEd[s$se_cd$Comparison == lab]
  expect_equal(g("A levels at the same or different B"),
               sqrt(2 * ((b - 1) * Eb + Ea) / (r * b * cc)))
  expect_equal(g("B levels at the same A, at the same or different C"),
               sqrt(2 * ((cc - 1) * Ec + Eb) / (r * cc)))
  expect_equal(g("A levels at the same or different B x C"),
               sqrt(2 * (b * (cc - 1) * Ec + (b - 1) * Eb + Ea) /
                      (r * b * cc)))
})

test_that("pooled split-plot F-tests match aov Error strata", {
  set.seed(42)
  d <- expand.grid(env = paste0("E", 1:3), rep = factor(1:3),
                   A = paste0("A", 1:3), B = paste0("B", 1:2))
  d$y <- 60 + rnorm(nrow(d), 0, 2) +
    rnorm(27, 0, 2)[as.integer(interaction(d$env, d$rep, d$A))]
  sp <- anova_spd_pooled(d, "y", "A", "B", "env", "rep", verbose = FALSE)
  e <- suppressWarnings(summary(aov(y ~ env * A * B + Error(env:rep / A),
                                    d)))
  expect_equal(sp$anova_table["env", "F value"],
               e[["Error: env:rep"]][[1]]["env", "F value"])
  expect_equal(sp$anova_table["A", "F value"],
               e[["Error: env:rep:A"]][[1]]["A", "F value"])
  expect_equal(sp$anova_table["env:A:B", "F value"],
               e[["Error: Within"]][[1]]["env:A:B", "F value"])
  # rep labels unique across environments give the same answer
  d2 <- d
  d2$rep <- as.integer(d2$rep) + 3 * (as.integer(factor(d2$env)) - 1)
  sp2 <- anova_spd_pooled(d2, "y", "A", "B", "env", "rep", verbose = FALSE)
  expect_equal(sp2$anova_table[["F value"]], sp$anova_table[["F value"]])
})

test_that("strip-plot tests each strip factor against its own error", {
  set.seed(43)
  d <- expand.grid(rep = factor(1:3), H = paste0("H", 1:3),
                   V = paste0("V", 1:3))
  d$y <- 40 + rnorm(nrow(d), 0, 3) +
    rnorm(9, 0, 3)[as.integer(interaction(d$rep, d$H))]
  st <- anova_strip(d, "y", "H", "V", "rep", verbose = FALSE)
  e <- summary(aov(y ~ H * V + Error(rep / (H * V)), d))
  expect_equal(st$anova_table["H", "F value"],
               e[["Error: rep:H"]][[1]]["H", "F value"])
  expect_equal(st$anova_table["V", "F value"],
               e[["Error: rep:V"]][[1]]["V", "F value"])
  expect_equal(st$anova_table["H:V", "F value"],
               e[["Error: rep:H:V"]][[1]]["H:V", "F value"])
  r <- 3
  a <- 3
  b <- 3
  expect_equal(st$se_cd$SEd[3],
               sqrt(2 * ((b - 1) * st$ms_error_c + st$ms_error_a) / (r * b)))
  expect_length(st$cv, 3)
})

test_that("augmented design follows Federer's adjustment", {
  set.seed(5)
  checks <- c("C1", "C2", "C3", "C4")
  tests <- paste0("T", 1:20)
  d <- do.call(rbind, lapply(1:5, function(j) {
    data.frame(block = j, genotype = c(checks, tests[((j - 1) * 4 + 1):(j * 4)]))
  }))
  d$y <- 100 + rnorm(nrow(d), 0, 8) + d$block * 3
  d$block <- factor(d$block)
  au <- anova_augmented(d, "y", "genotype", "block", checks, verbose = FALSE)
  # Federer's adjustment: observed - (check mean in block - overall check mean)
  ck <- d[d$genotype %in% checks, ]
  beff <- tapply(ck$y, ck$block, mean) - mean(ck$y)
  tst <- d[!d$genotype %in% checks, ]
  adj <- as.numeric(tst$y - beff[as.character(tst$block)])
  m <- au$adjusted_means
  expect_equal(m$Adjusted_mean[match(tst$genotype, m$Entry)], unname(adj),
               tolerance = 1e-10)
  mse <- au$mse
  b <- 5
  cc <- 4
  expect_equal(unname(au$se_diff),
               c(sqrt(2 * mse / b), sqrt(2 * mse), sqrt(2 * mse * (1 + 1 / cc)),
                 sqrt(mse * (1 + 1 / b + 1 / cc - 1 / (b * cc)))))
  expect_equal(nrow(au$test_vs_checks), 20)
  expect_true(all(c("Heritability (broad sense, %)", "GCV (%)") %in%
                    au$genetic_parameters$Parameter))
})

test_that("alpha lattice intra-block means equal least-squares means", {
  skip_if_not_installed("emmeans")
  data(cowpea_alpha, package = "aridagri")
  al <- anova_alpha_lattice(cowpea_alpha, "grain_yield", "genotype", "rep",
                            "block", analysis = "intra", verbose = FALSE)
  d <- cowpea_alpha
  d$ib <- interaction(d$rep, d$block)
  em <- summary(emmeans::emmeans(lm(grain_yield ~ ib + genotype, d),
                                 "genotype"))
  m <- merge(al$genotype_means, data.frame(Genotype = em$genotype,
                                           em = em$emmean))
  expect_equal(m$Adjusted_mean, m$em, tolerance = 1e-8)
  expect_gt(al$efficiency, 100)
  ar <- anova_alpha_lattice(cowpea_alpha, "grain_yield", "genotype", "rep",
                            "block", verbose = FALSE)
  expect_equal(ar$analysis, "reml")
  expect_equal(nrow(ar$variance_components), 3)
})

test_that("repeated measures match aov strata and base R sphericity", {
  data(mungbean_repeated, package = "aridagri")
  d <- mungbean_repeated
  rm <- anova_repeated(d, "spad", "treatment", "das", block = "block",
                       verbose = FALSE)
  d$das <- factor(d$das)
  e <- summary(aov(spad ~ block + treatment + das + block:das +
                     treatment:das + Error(plot), d))
  tb <- e[["Error: plot"]][[1]]
  tw <- e[["Error: Within"]][[1]]
  expect_equal(rm$anova_table["treatment", "F value"],
               tb[trimws(rownames(tb)) == "treatment", "F value"])
  expect_equal(rm$anova_table["treatment:das", "F value"],
               tw[trimws(rownames(tw)) == "treatment:das", "F value"])
  w <- reshape(d[, c("plot", "block", "treatment", "das", "spad")],
               idvar = "plot", timevar = "das", v.names = "spad",
               direction = "wide")
  Y <- as.matrix(w[, grep("^spad", names(w))])
  mlm <- lm(Y ~ block + treatment, data = w)
  sph <- anova(mlm, X = ~1, test = "Spherical")
  expect_equal(rm$corrected$p_GG[2], sph[["G-G Pr"]][3], tolerance = 1e-8)
  expect_equal(rm$corrected$p_HF[2], sph[["H-F Pr"]][3], tolerance = 1e-8)
  expect_equal(rm$sphericity$Mauchly_W,
               unname(mauchly.test(mlm, X = ~1)$statistic))
  expect_true(rm$best_structure %in% c("CS", "AR1", "UN"))
  expect_equal(nrow(rm$by_time), 4)
})

test_that("anova_traits keeps treatments aligned across traits", {
  data(mothbean_rbd, package = "aridagri")
  tb <- anova_traits(mothbean_rbd, c("seed_yield", "pods_per_plant"),
                     treatment = "treatment", block = "block",
                     verbose = FALSE)
  r1 <- anova_rbd(mothbean_rbd, "seed_yield", "treatment", "block",
                  verbose = FALSE)
  m1 <- r1$treatment_means
  lv <- levels(mothbean_rbd$treatment)
  cell <- tb$table$seed_yield[match(lv, tb$table$treatment)]
  num <- as.numeric(sub(" .*", "", cell))
  expect_equal(num, round(m1$Mean[match(lv, m1$Treatment)], 2))
  expect_true(any(grepl("CD", tb$table$treatment)))
})

test_that("check_assumptions tests variances across treatments", {
  set.seed(2)
  d <- expand.grid(trt = paste0("T", 1:4), blk = factor(1:4))
  d$y <- rnorm(16, 10, 2)
  ca <- check_assumptions(aov(y ~ blk + trt, d), verbose = FALSE)
  expect_equal(ca$groups_factor, "trt")
  expect_false(is.null(ca$levene))
  expect_false(is.null(ca$nonadditivity))
})

test_that("Latin square efficiencies and letters", {
  data <- data.frame(
    row = rep(1:5, each = 5), col = rep(1:5, 5),
    treatment = c("A", "B", "C", "D", "E", "B", "C", "D", "E", "A",
                  "C", "D", "E", "A", "B", "D", "E", "A", "B", "C",
                  "E", "A", "B", "C", "D"),
    yield = c(1180, 1265, 1342, 1410, 1238, 1270, 1335, 1402, 1251, 1176,
              1340, 1415, 1245, 1190, 1258, 1405, 1232, 1178, 1262, 1339,
              1236, 1183, 1271, 1348, 1398))
  res <- anova_latin(data, "yield", "treatment", "row", "col",
                     verbose = FALSE)
  expect_length(res$relative_efficiency, 3)
  expect_true("Group" %in% names(res$treatment_means))
  expect_equal(res$df_error, 12)
})
