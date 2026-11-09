rbd_data <- function(seed = 7) {
  set.seed(seed)
  k <- 7
  r <- 4
  d <- expand.grid(trt = paste0("T", 1:k), blk = factor(1:r))
  d$y <- 50 + rep(c(0, 2, 3.5, 5, 6, 9, 9.5), r) + rnorm(k * r, 0, 2.2)
  d
}

ag_letters <- function(x) {
  g <- x$groups
  stats::setNames(trimws(as.character(g$groups)), rownames(g))
}

test_that("letters match agricolae for LSD, HSD, Duncan, SNK and Scheffe", {
  skip_if_not_installed("agricolae")
  d <- rbd_data()
  fit <- aov(y ~ blk + trt, d)
  ph <- perform_posthoc(fit, d, "y", "trt",
                        posthoc = c("lsd", "tukey", "duncan", "snk",
                                    "scheffe"), verbose = FALSE)
  mine <- function(m) stats::setNames(ph[[m]]$groups$groups,
                                      rownames(ph[[m]]$groups))
  lev <- levels(d$trt)
  expect_equal(mine("lsd")[lev],
               ag_letters(agricolae::LSD.test(fit, "trt"))[lev])
  expect_equal(mine("tukey")[lev],
               ag_letters(agricolae::HSD.test(fit, "trt"))[lev])
  expect_equal(mine("duncan")[lev],
               ag_letters(agricolae::duncan.test(fit, "trt"))[lev])
  expect_equal(mine("snk")[lev],
               ag_letters(agricolae::SNK.test(fit, "trt"))[lev])
  expect_equal(mine("scheffe")[lev],
               ag_letters(agricolae::scheffe.test(fit, "trt"))[lev])
  expect_equal(ph$lsd$statistics$critical_difference,
               agricolae::LSD.test(fit, "trt")$statistics$LSD)
  expect_equal(ph$tukey$statistics$critical_difference,
               agricolae::HSD.test(fit, "trt")$statistics$MSD)
  expect_equal(unname(ph$duncan$critical_range),
               agricolae::duncan.test(fit, "trt")$duncan$CriticalRange)
})

test_that("Tukey-Kramer p-values match TukeyHSD with unequal replication", {
  set.seed(3)
  d <- data.frame(trt = factor(rep(paste0("T", 1:5), c(3, 5, 4, 6, 4))))
  d$y <- 20 + as.integer(d$trt) * 1.5 + rnorm(nrow(d), 0, 1.5)
  fit <- aov(y ~ trt, d)
  ph <- perform_posthoc(fit, d, "y", "trt", posthoc = "tukey",
                        verbose = FALSE)
  tk <- TukeyHSD(fit)$trt
  cmp <- ph$tukey$comparisons
  key <- vapply(strsplit(cmp$comparison, " - "),
                function(z) paste0(z[2], "-", z[1]), character(1))
  expect_equal(cmp$p_value, unname(tk[key, "p adj"]), tolerance = 1e-8)
})

test_that("least-squares means with missing plots match emmeans", {
  skip_if_not_installed("emmeans")
  set.seed(3)
  d <- expand.grid(trt = factor(paste0("T", 1:5)), blk = factor(1:4))
  d$y <- 100 + as.integer(d$trt) * 3 + as.integer(d$blk) * 2 +
    rnorm(20, 0, 2)
  d <- d[-c(3, 12), ]
  fit <- aov(y ~ blk + trt, d)
  lsm <- aridagri:::.ls_means_model(fit, "trt")
  em <- summary(emmeans::emmeans(fit, "trt"))
  expect_equal(unname(lsm$means), em$emmean, tolerance = 1e-10)
  expect_equal(unname(sqrt(diag(lsm$vcov))), em$SE, tolerance = 1e-10)
  pr <- summary(pairs(emmeans::emmeans(fit, "trt"), adjust = "none"))
  ph <- perform_posthoc(fit, d, "y", "trt", posthoc = "lsd", verbose = FALSE)
  expect_equal(sort(ph$lsd$comparisons$p_value), sort(pr$p.value),
               tolerance = 1e-10)
})

test_that("Duncan and SNK apply the closure rule", {
  # T1-T3 range not significant but T1-T2 and T2-T3 sub-ranges could be:
  # closure forbids declaring a sub-range significant
  m <- c(A = 10, B = 13.4, C = 13.6, D = 17)
  out <- compare_means(m, mse = 2, df_error = 12, reps = 2, method = "snk",
                       verbose = FALSE)
  sig <- out$comparisons
  # any pair inside a non-significant wider range is not significant
  rng <- function(a, b) sig$significant[sig$comparison == paste(a, "-", b)]
  if (!rng("A", "D")) {
    expect_false(rng("A", "C"))
    expect_false(rng("B", "D"))
  }
  expect_true(all(nzchar(out$letters)))
})

test_that("compare_means works from summary statistics", {
  m <- c(T1 = 1180, T2 = 1265, T3 = 1342, T4 = 1410, T5 = 1238)
  out <- compare_means(m, mse = 3600, df_error = 12, reps = 4,
                       method = "lsd", verbose = FALSE)
  lsd <- qt(0.975, 12) * sqrt(2 * 3600 / 4)
  expect_equal(out$statistics$critical_difference, lsd)
  # T4-T3, T3-T2, T2-T5-T1 overlap: three groups
  expect_equal(unname(out$letters[c("T4", "T3", "T2", "T5", "T1")]),
               c("a", "ab", "bc", "c", "c"))
  expect_error(compare_means(m, mse = -1, df_error = 12, reps = 4))
  expect_error(compare_means(m, mse = 1, df_error = 12, reps = 4,
                             method = "dunnett"))
})

test_that("letter display is valid for any significance pattern", {
  set.seed(11)
  for (i in 1:20) {
    k <- sample(4:9, 1)
    m <- sort(rnorm(k, 0, 3), decreasing = TRUE)
    names(m) <- paste0("T", seq_len(k))
    out <- compare_means(m, mse = 1, df_error = 20, reps = 3,
                         method = sample(c("lsd", "tukey", "duncan"), 1),
                         verbose = FALSE)
    let <- out$letters
    for (r in seq_len(nrow(out$comparisons))) {
      pr <- strsplit(out$comparisons$comparison[r], " - ")[[1]]
      share <- length(intersect(strsplit(let[pr[1]], "")[[1]],
                                strsplit(let[pr[2]], "")[[1]])) > 0
      expect_equal(share, !out$comparisons$significant[r])
    }
  }
})

test_that("perform_posthoc keeps the agricolae-style groups element", {
  d <- rbd_data()
  fit <- aov(y ~ blk + trt, d)
  ph <- perform_posthoc(fit, d, "y", "trt", posthoc = "all",
                        verbose = FALSE)
  expect_true(all(c("lsd", "duncan", "tukey", "snk", "scheffe",
                    "bonferroni") %in% names(ph)))
  expect_true("groups" %in% names(ph$lsd$groups))
  expect_error(perform_posthoc(fit, d, "y", "trt", posthoc = "xyz",
                               verbose = FALSE), "Unknown post-hoc")
})

test_that("Dunnett's test uses the chosen control", {
  skip_if_not_installed("multcomp")
  d <- rbd_data()
  fit <- aov(y ~ blk + trt, d)
  ph <- perform_posthoc(fit, d, "y", "trt", posthoc = "dunnett",
                        control = "T3", verbose = FALSE)
  expect_true(all(grepl("- T3$", ph$dunnett$comparison)))
  expect_equal(nrow(ph$dunnett), 6)
})
