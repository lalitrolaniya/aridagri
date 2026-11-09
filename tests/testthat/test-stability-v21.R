met_data <- function(seed = 21) {
  set.seed(seed)
  g <- 8
  e <- 5
  r <- 3
  d <- expand.grid(gen = paste0("G", 1:g), env = paste0("E", 1:e), rep = 1:r)
  ge <- matrix(rnorm(g * e, 0, 60), g, e)
  d$y <- 1200 + rep(rnorm(g, 0, 80), e * r) +
    rep(rep(rnorm(e, 0, 150), each = g), r) +
    ge[cbind(as.integer(d$gen), as.integer(d$env))] + rnorm(nrow(d), 0, 50)
  d
}

test_that("Shukla stability variance matches agricolae", {
  skip_if_not_installed("agricolae")
  d <- met_data()
  s <- stability_analysis(d, "gen", "env", "rep", "y", verbose = FALSE)
  m <- tapply(d$y, list(d$gen, d$env), mean)
  sp <- agricolae::stability.par(as.data.frame(m), rep = 3,
                                 MSerror = s$mse, console = FALSE)
  ag <- stats::setNames(sp$statistics[["Sigma-square"]], rownames(m)) / 3
  mine <- stats::setNames(s$shukla$Sigma2_i, s$shukla$Genotype)[rownames(m)]
  expect_equal(unname(mine), unname(ag), tolerance = 1e-3)
})

test_that("AMMI stability value matches agricolae index.AMMI", {
  skip_if_not_installed("agricolae")
  d <- met_data()
  s <- stability_analysis(d, "gen", "env", "rep", "y", method = "ammi",
                          verbose = FALSE)
  am <- with(d, agricolae::AMMI(env, gen, rep, y, console = FALSE))
  ix <- agricolae::index.AMMI(am)
  mine <- stats::setNames(s$ammi$stability$ASV,
                          s$ammi$stability$Genotype)[rownames(ix)]
  expect_equal(unname(mine), ix$ASV, tolerance = 1e-4)
  # IPCA sums of squares add up to the G x E sum of squares
  expect_equal(sum(s$ammi$components$Sum_Sq) +
                 0, s$anova["Genotype x Environment", "Sum Sq"],
               tolerance = 1e-6)
})

test_that("Eberhart-Russell parameters follow the published formulas", {
  d <- met_data()
  s <- stability_analysis(d, "gen", "env", "rep", "y", method = "eberhart",
                          verbose = FALSE)
  m <- tapply(d$y, list(d$gen, d$env), mean)
  I <- colMeans(m) - mean(m)
  b <- as.numeric(m %*% I) / sum(I^2)
  dev <- m - rowMeans(m) - outer(b, I)
  s2d <- rowSums(dev^2) / (ncol(m) - 2) - s$mse / 3
  er <- s$eberhart[match(rownames(m), s$eberhart$Genotype), ]
  expect_equal(er$bi, b)
  expect_equal(er$S2di, unname(s2d))
  era <- s$eberhart_anova
  expect_equal(era["Env + (G x E)", "Sum Sq"],
               era["Environment (linear)", "Sum Sq"] +
                 era["G x E (linear)", "Sum Sq"] +
                 era["Pooled deviation", "Sum Sq"])
})

test_that("stability analysis works with genotype x environment means", {
  d <- met_data()
  m <- aggregate(y ~ gen + env, d, mean)
  s <- stability_analysis(m, "gen", "env", NULL, "y", verbose = FALSE)
  expect_null(s$anova)
  expect_true(all(is.na(s$shukla$p_value)))
  expect_equal(nrow(s$integrated), 8)
})

test_that("stability analysis rejects incomplete tables", {
  d <- met_data()
  d <- d[!(d$gen == "G1" & d$env == "E2"), ]
  expect_error(stability_analysis(d, "gen", "env", "rep", "y",
                                  verbose = FALSE), "every environment")
})

test_that("GGE biplot finds the winner in each environment", {
  d <- met_data()
  gg <- gge_biplot(d, "gen", "env", "y", verbose = FALSE)
  expect_equal(nrow(gg$which_won_where), 5)
  expect_true(all(gg$which_won_where$Winner %in% levels(d$gen)))
  expect_equal(sum(gg$variance_percent), 100)
  expect_equal(nrow(gg$aec), 8)
  g2 <- gge_biplot(d, "gen", "env", "y", scaling = "sd", svp = "genotype",
                   verbose = FALSE)
  expect_equal(g2$svp, "genotype")
})
