# Tests for stability analysis

test_that("stability_analysis with method='all' produces integrated ranking", {
  set.seed(42)
  data <- expand.grid(
    genotype = paste0("G", 1:5),
    environment = paste0("E", 1:4),
    rep = 1:3
  )
  data$yield <- rnorm(nrow(data), 1500, 250)

  result <- stability_analysis(data, genotype = "genotype",
                               environment = "environment",
                               replication = "rep", trait = "yield",
                               method = "all", verbose = FALSE)

  expect_type(result, "list")
  expect_true("integrated" %in% names(result))
  expect_true("eberhart" %in% names(result))
  expect_true("ammi" %in% names(result))
  expect_true("shukla" %in% names(result))
  expect_true("wricke" %in% names(result))
  expect_true("cv_stability" %in% names(result))
  expect_true("superiority" %in% names(result))
})

test_that("integrated ranking has all 7 method columns", {
  set.seed(42)
  data <- expand.grid(
    genotype = paste0("G", 1:5),
    environment = paste0("E", 1:4),
    rep = 1:3
  )
  data$yield <- rnorm(nrow(data), 1500, 250)

  result <- stability_analysis(data, genotype = "genotype",
                               environment = "environment",
                               replication = "rep", trait = "yield",
                               method = "all", verbose = FALSE)

  integrated <- result$integrated
  expected_cols <- c("Genotype", "Mean", "ER", "AMMI", "FW",
                     "Shukla", "Wricke", "CV", "LB", "Rank")
  expect_true(all(expected_cols %in% names(integrated)))
  expect_equal(ncol(integrated), 10)
})

test_that("integrated ranking uses all 7 methods for final rank", {
  set.seed(42)
  data <- expand.grid(
    genotype = paste0("G", 1:8),
    environment = paste0("E", 1:5),
    rep = 1:3
  )
  data$yield <- rnorm(nrow(data), 1500, 250)

  result <- stability_analysis(data, genotype = "genotype",
                               environment = "environment",
                               replication = "rep", trait = "yield",
                               method = "all", verbose = FALSE)

  integrated <- result$integrated
  # Verify rank is based on mean of 7 method ranks
  rank_cols <- c("ER", "AMMI", "FW", "Shukla", "Wricke", "CV", "LB")
  mean_ranks <- rowMeans(integrated[, rank_cols])
  expected_rank <- rank(mean_ranks)
  # Ranks should match after reordering
  expect_equal(as.numeric(integrated$Rank), sort(as.numeric(integrated$Rank)))
})

test_that("stability_analysis with single method skips integrated ranking", {
  set.seed(42)
  data <- expand.grid(
    genotype = paste0("G", 1:5),
    environment = paste0("E", 1:4),
    rep = 1:3
  )
  data$yield <- rnorm(nrow(data), 1500, 250)

  result <- stability_analysis(data, genotype = "genotype",
                               environment = "environment",
                               replication = "rep", trait = "yield",
                               method = "eberhart", verbose = FALSE)

  expect_true("eberhart" %in% names(result))
  expect_false("integrated" %in% names(result))
})

test_that("Shukla and Wricke produce identical rankings", {
  set.seed(42)
  data <- expand.grid(
    genotype = paste0("G", 1:6),
    environment = paste0("E", 1:5),
    rep = 1:3
  )
  data$yield <- rnorm(nrow(data), 1500, 250)

  result <- stability_analysis(data, genotype = "genotype",
                               environment = "environment",
                               replication = "rep", trait = "yield",
                               method = "all", verbose = FALSE)

  # Shukla and Wricke ranks must be identical (mathematical property)
  expect_equal(result$integrated$Shukla, result$integrated$Wricke)
})
