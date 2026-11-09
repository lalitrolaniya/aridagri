# Simulated example datasets for aridagri.
# All values are generated from the statistical models below with means and
# variances typical of arid western Rajasthan. They are for teaching and for
# the examples only and are not results of real experiments.
# Run from the package root: Rscript data-raw/make_datasets.R

set.seed(2026)

## 1. mothbean_rbd: RCBD, 8 nutrient treatments x 4 blocks, 8 traits
trt <- c("T1 Control", "T2 RDF", "T3 RDF+Rhizobium", "T4 RDF+PSB",
         "T5 RDF+Rhizobium+PSB", "T6 75% RDF+FYM", "T7 RDF+S20",
         "T8 RDF+Zn5")
yield_eff <- c(0, 95, 140, 128, 205, 168, 150, 118)
d1 <- expand.grid(treatment = trt, block = paste0("R", 1:4),
                  stringsAsFactors = FALSE)
d1 <- d1[, c("block", "treatment")]
it <- match(d1$treatment, trt)
ib <- match(d1$block, paste0("R", 1:4))
blk <- c(-18, 6, 14, -2)[ib]
d1$seed_yield <- round(505 + yield_eff[it] + blk + rnorm(nrow(d1), 0, 34))
d1$plant_height <- round(31 + 0.022 * yield_eff[it] + rnorm(nrow(d1), 0, 1.6), 1)
d1$branches_per_plant <- round(3.1 + 0.004 * yield_eff[it] +
                                 rnorm(nrow(d1), 0, 0.25), 1)
d1$pods_per_plant <- round(18.5 + 0.055 * yield_eff[it] + 0.02 * blk +
                             rnorm(nrow(d1), 0, 1.4), 1)
d1$seeds_per_pod <- round(5.6 + 0.0035 * yield_eff[it] +
                            rnorm(nrow(d1), 0, 0.22), 2)
d1$test_weight <- round(28.4 + 0.009 * yield_eff[it] +
                          rnorm(nrow(d1), 0, 0.55), 2)
d1$stover_yield <- round(1180 + 1.55 * yield_eff[it] + 1.8 * blk +
                           rnorm(nrow(d1), 0, 70))
d1$harvest_index <- round(100 * d1$seed_yield /
                            (d1$seed_yield + d1$stover_yield), 2)
d1$block <- factor(d1$block)
d1$treatment <- factor(d1$treatment, levels = trt)
mothbean_rbd <- d1

## 2. clusterbean_spd: split plot, irrigation (main) x nitrogen (sub), 3 reps
irr <- c("I1 Rainfed", "I2 One irrigation", "I3 Two irrigations")
nit <- c("N0", "N10", "N20", "N20+Foliar")
d2 <- expand.grid(nitrogen = nit, irrigation = irr, rep = paste0("R", 1:3),
                  stringsAsFactors = FALSE)
d2 <- d2[, c("rep", "irrigation", "nitrogen")]
ii <- match(d2$irrigation, irr)
inn <- match(d2$nitrogen, nit)
ir <- match(d2$rep, paste0("R", 1:3))
main_err <- rnorm(9, 0, 38)[(ir - 1) * 3 + ii]
eff_i <- c(0, 210, 330)[ii]
eff_n <- c(0, 75, 120, 165)[inn]
inter <- c(0, 0, 0, 0, 0, 20, 35, 55, 0, 30, 55, 80)[(ii - 1) * 4 + inn]
d2$seed_yield <- round(820 + eff_i + eff_n + inter + c(-15, 10, 5)[ir] +
                         main_err + rnorm(nrow(d2), 0, 30))
water <- c(212, 272, 332)[ii]
d2$water_use <- water + round(rnorm(nrow(d2), 0, 4))
d2$wue <- round(d2$seed_yield / d2$water_use, 3)
d2$stover_yield <- round(2.05 * d2$seed_yield + rnorm(nrow(d2), 0, 90))
d2$rep <- factor(d2$rep)
d2$irrigation <- factor(d2$irrigation, levels = irr)
d2$nitrogen <- factor(d2$nitrogen, levels = nit)
clusterbean_spd <- d2

## 3. pearlmillet_met: 10 genotypes x 6 environments x 3 reps
gen <- sprintf("PM-%02d", 1:10)
env <- c("BKN-2025", "BKN-2026", "JDH-2025", "JDH-2026", "JSM-2025",
         "JSM-2026")
g_eff <- c(120, 60, -40, 200, -110, 30, 150, -60, -20, -10)
e_eff <- c(-180, 60, 320, 420, -520, -380)
# interaction with structure: some genotypes respond to good environments
sens <- c(1.25, 0.95, 0.80, 1.10, 0.70, 1.00, 1.30, 0.85, 1.05, 1.00)
d3 <- expand.grid(rep = 1:3, environment = env, genotype = gen,
                  stringsAsFactors = FALSE)
d3 <- d3[, c("environment", "rep", "genotype")]
ig <- match(d3$genotype, gen)
ie <- match(d3$environment, env)
ge <- matrix(rnorm(60, 0, 55), 10, 6)
rep_eff <- matrix(rnorm(18, 0, 40), 6, 3)
d3$grain_yield <- round(1850 + g_eff[ig] + sens[ig] * e_eff[ie] +
                          ge[cbind(ig, ie)] + rep_eff[cbind(ie, d3$rep)] +
                          rnorm(nrow(d3), 0, 95))
d3$days_to_maturity <- round(78 + c(2, 0, -3, 4, -5, 1, 3, -2, 0, -1)[ig] +
                               c(1, 2, -1, 0, -3, -2)[ie] +
                               rnorm(nrow(d3), 0, 1.5))
d3$environment <- factor(d3$environment, levels = env)
d3$genotype <- factor(d3$genotype, levels = gen)
d3$rep <- factor(d3$rep)
pearlmillet_met <- d3

## 4. chickpea_augmented: augmented RCBD, 6 blocks, 4 checks, 48 tests
checks <- c("C1", "C2", "C3", "C4")
tests <- sprintf("IC-%03d", 1:48)
chk_eff <- c(0, 85, 140, 40)
d4 <- do.call(rbind, lapply(1:6, function(b) {
  ent <- c(checks, tests[((b - 1) * 8 + 1):(b * 8)])
  data.frame(block = paste0("B", b), entry = sample(ent),
             stringsAsFactors = FALSE)
}))
bb <- c(-60, -20, 15, 40, 70, -45)[match(d4$block, paste0("B", 1:6))]
test_eff <- stats::setNames(rnorm(48, 0, 120), tests)
eff <- ifelse(d4$entry %in% checks, chk_eff[match(d4$entry, checks)],
              test_eff[d4$entry])
d4$seed_yield <- round(1420 + eff + bb + rnorm(nrow(d4), 0, 55))
d4$plant_height <- round(46 + 0.012 * eff + 0.02 * bb +
                           rnorm(nrow(d4), 0, 2.2), 1)
d4$block <- factor(d4$block)
d4$entry <- factor(d4$entry, levels = c(checks, tests))
chickpea_augmented <- d4

## 5. cowpea_alpha: alpha lattice, 24 genotypes, 3 reps, 6 blocks of 4
## Layout generated with agricolae::design.alpha(trt, k = 4, r = 3, seed = 7)
genos <- sprintf("CP-%02d", 1:24)
lay <- agricolae::design.alpha(genos, k = 4, r = 3, seed = 7)$book
names(lay)[names(lay) == "genos"] <- "genotype"
d5 <- data.frame(rep = factor(lay$replication),
                 block = factor(as.integer(as.character(lay$block)) -
                                  (as.integer(as.character(lay$replication)) -
                                     1) * 6),
                 plot = lay$plots, genotype = factor(lay$genotype,
                                                     levels = genos))
gef <- stats::setNames(rnorm(24, 0, 70), genos)
bef <- rnorm(18, 0, 60)
bid <- (as.integer(d5$rep) - 1) * 6 + as.integer(as.character(d5$block))
d5$grain_yield <- round(1080 + gef[as.character(d5$genotype)] + bef[bid] +
                          c(-30, 10, 20)[as.integer(d5$rep)] +
                          rnorm(nrow(d5), 0, 45))
cowpea_alpha <- d5

## 6. mungbean_repeated: RCBD, 4 foliar treatments x 4 blocks, SPAD at 4 dates
trt6 <- c("Control", "Water spray", "Urea 2%", "ZnSO4 0.5%")
das <- c(25, 35, 45, 55)
d6 <- expand.grid(das = das, treatment = trt6, block = paste0("R", 1:4),
                  stringsAsFactors = FALSE)
d6 <- d6[, c("block", "treatment", "das")]
d6$plot <- paste(d6$block, match(d6$treatment, trt6), sep = "-")
plots <- unique(d6$plot)
ar <- 0.6
e6 <- unlist(lapply(plots, function(p) {
  z <- numeric(4)
  z[1] <- rnorm(1, 0, 1.3)
  for (k in 2:4) z[k] <- ar * z[k - 1] + rnorm(1, 0, 1.3 * sqrt(1 - ar^2))
  z
}))
d6 <- d6[order(d6$plot, d6$das), ]
it6 <- match(d6$treatment, trt6)
id6 <- match(d6$das, das)
d6$spad <- round(31.5 + c(0, 6.5, 10.2, 7.1)[id6] +
                   matrix(c(0, 0, 0, 0, 0, 0.4, 0.6, 0.3, 0, 1.8, 3.6, 3.0,
                            0, 1.1, 2.4, 2.2), 4, 4,
                          byrow = TRUE)[cbind(it6, id6)] +
                   c(-0.6, 0.2, 0.5, -0.1)[match(d6$block, paste0("R", 1:4))] +
                   e6, 1)
d6$block <- factor(d6$block)
d6$treatment <- factor(d6$treatment, levels = trt6)
d6$plot <- factor(d6$plot)
rownames(d6) <- NULL
mungbean_repeated <- d6

for (nm in c("mothbean_rbd", "clusterbean_spd", "pearlmillet_met",
             "chickpea_augmented", "cowpea_alpha", "mungbean_repeated")) {
  obj <- get(nm)
  rownames(obj) <- NULL
  assign(nm, obj)
  save(list = nm, file = file.path("data", paste0(nm, ".rda")),
       compress = "xz")
}
