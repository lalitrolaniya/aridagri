# aridagri 2.1.0

This release is dedicated to the memory of co-authors Dr. Monika Punia
and Dr. Hemant Poonia.

Version 2.1.0 corrects several calculations, adds ten functions and six
example datasets, and removes the dependence on 'agricolae' for letter
groupings. Results from earlier versions should be rechecked for the
functions listed under "Corrected results".

## Corrected results

* `anova_rbd()`: the relative efficiency of the RCBD over a CRD used the
  coefficient (t + r - 2) instead of r(t - 1) on the error mean square and
  replaced a small block mean square by the error mean square. With five
  treatments, four blocks and no block effect it reported 52.6% instead of
  about 100%. It now follows Gomez and Gomez (1984), with the precision
  adjustment for fewer than 20 error df.
* `anova_rbd_pooled()`: with `verbose = FALSE` the analyses of the
  individual environments were skipped, so Bartlett's test returned NaN and
  `individual_anovas` was empty. The environment effect was tested against
  the pooled residual instead of blocks within environments. Both are fixed,
  and a new `env_effect = "random"` option tests treatments against the
  treatment x environment interaction (McIntosh, 1983).
* `anova_spd()`: the returned `anova_table` tested the main-plot factor
  against Error(b); in one test case F changed from 2.91 (p = 0.080) to the
  correct 0.60 (p = 0.59). The printed table was already correct.
* `anova_strip()`: both strip factors were tested against Error(c). The
  horizontal factor is now tested against Error(a) = Rep x A and the
  vertical factor against Error(b) = Rep x B; in one test case a false
  p = 0.007 became p = 0.30. CV(a), CV(b), CV(c) and the four SE/CD values of
  Gomez and Gomez are now reported.
* `stability_analysis()`:
  * Shukla's stability variance was W/(e - 1). It is now Shukla's unbiased
    estimator with an F-test (rankings were not affected).
  * The AMMI stability value put the SS ratio inside the square and used
    IPCA scores scaled by the singular value instead of its square root.
    ASV now follows Purchase et al. (2000) and agrees with 'agricolae'.
  * Eberhart-Russell S2di now subtracts the pooled error (MSE/r), and the
    tests of bi = 1 and S2di = 0 are reported.
* `check_assumptions()`: the homogeneity test used the first factor of the
  model, which for an RCBD was the block. It now uses the treatment (or a
  factor you choose) and the residuals. It no longer claims to run an
  Anderson-Darling test for more than 5000 residuals.
* `nue_calculate()`: the "utilization efficiency" column was AE x ARE, which
  is not a defined index. UE is now the biomass-based physiological
  efficiency times ARE (Fageria and Baligar, 2005), given when biomass is
  supplied. The function also failed for vector input.
* `crop_growth_analysis()`: LAI was computed but never returned.
* `nutrient_response()`: no check that the quadratic curve has a maximum.
* `anova_augmented()` did not compute adjusted means of the test entries;
  `anova_alpha_lattice()` reported unadjusted genotype means.
* Many functions ignored `alpha` (CD was always at 5%) or `posthoc`, and
  counted rows with a missing response.

## New functions

* `anova_repeated()`: repeated measures (split-plot in time) with Mauchly's
  test, Greenhouse-Geisser and Huynh-Feldt corrections, a comparison of
  CS, AR(1) and unstructured covariance models, SE/CD for every comparison
  and an analysis at each time.
* `anova_traits()`: the treatment-means table for many traits that journals
  ask for (means with letters, SEm, CD or NS, CV) for CRD, RCBD, factorial
  and split-plot designs.
* `compare_means()`: letters from means, error MS and df alone.
* `gge_biplot()`: GGE biplot with which-won-where, mega-environments and the
  mean vs stability view.
* `day_length()`: day length from latitude and date (FAO-56).
* `plot_means()`, `plot_interaction()`, `plot_stability()`,
  `plot_response()`, `theme_aridagri()` and `plot()` methods for
  publication figures with 'ggplot2'.

## Improvements

* Letter groupings for LSD, Duncan, Tukey, SNK, Scheffe and Bonferroni are
  computed inside the package (no 'agricolae' needed), with exact pairwise
  tests for unequal replication (Tukey-Kramer) and missing plots, the
  closure rule for Duncan and SNK (the reported p-values follow it), and
  the letter algorithm of Piepho (2004). Results agree with 'agricolae' for
  balanced data. Least-squares means that cannot be estimated (empty cells)
  are left out with a warning.
* `perform_posthoc()`: Dunnett's test now averages over interacting factors
  and uses `mse` and `df_error` when given; a new `control` argument sets
  the control level.
* `anova_rbd()`, `anova_crd()` and `anova_latin()` handle missing plots and
  unequal replication by least squares; means are least-squares means.
* All split-plot functions report SEm, SEd and CD for every type of
  comparison, including those that mix error terms (weighted t), and give
  letters for each factor and for each factor within the levels of another.
  The pooled split-plot functions now report SE and CD too.
* `anova_augmented()`: adjusted means, the four standard errors of
  difference, ANOVA with tests vs checks, checks and tests (adjusted and
  unadjusted for blocks), comparison of each test entry with the check mean
  and the best check, and genetic parameters (PCV, GCV, heritability,
  genetic advance) from the block-adjusted mean square of tests.
* `anova_alpha_lattice()`: adjusted means with exact SEs and letters,
  efficiency relative to an RCBD, and REML analysis with recovery of
  inter-block information (default; `analysis = "intra"` for the
  intra-block analysis).
* `stability_analysis()`: combined ANOVA, Eberhart-Russell ANOVA, AMMI ANOVA
  with Gollob's df and F-tests, Kendall's W among the seven rankings, and
  input of genotype x environment means (`replication = NULL`).
* `thermal_indices()`: upper threshold temperature, McMaster-Wilhelm
  method 2, helio-thermal use efficiency, input checks.
* `crop_growth_analysis()`: ground area for CGR and LAI, leaf area duration.
* `yield_gap_analysis()`: technology gap, extension gap and technology index
  for front-line demonstrations.
* `nutrient_response()`: linear-plateau, quadratic-plateau and Mitscherlich
  models and `model = "all"` for a comparison by AIC; optimum doses beyond
  the tested range are flagged.
* `anova_factorial()` and `anova_factorial_3way()`: SE/CD for every effect
  and letters for main effects and treatment combinations; with unbalanced
  data the means are least-squares means.
* `anova_latin()`: efficiency relative to CRD and RCBD.
* `check_assumptions()`: Levene (Brown-Forsythe) test, Tukey's test for
  non-additivity, studentized residuals and optional residual plots.
* Print methods for all results; `export_results()` writes every table to a
  separate Excel sheet.
* All functions accept tibbles (for example from 'readxl') and column names
  with spaces.

## Behaviour changes

* Split-plot type functions stop with an informative error when data are
  incomplete, instead of returning unreliable results with a warning.
* `stability_analysis()` needs at least three genotypes and environments.
* `nlme` (shipped with R) is now imported. Unused suggested packages were
  removed.

## Other changes

* Six simulated example datasets: `mothbean_rbd`, `clusterbean_spd`,
  `pearlmillet_met`, `chickpea_augmented`, `cowpea_alpha`,
  `mungbean_repeated`.
* A vignette, a pkgdown site and 607 test expectations (from 82), many
  checking results against 'agricolae', 'emmeans', `aov()` error strata and
  base R's sphericity tests.
* Hemant Poonia, co-author of the article describing the package
  (Software Impacts, 2026), is added as an author. He and Monika Punia
  are marked as deceased in the author list.

# aridagri 2.0.4

## Bug fixes

* `anova_spd()` split-plot ANOVA: the main-plot error term (Error a) and the
  sub-plot sum of squares are now extracted from the ANOVA table by term name
  rather than by fixed row position. Because R orders main-effect terms ahead of
  interaction terms, the previous positional extraction swapped the sub-plot
  effect with Error(a), which produced incorrect main-plot and sub-plot F-tests
  and an incorrect Error(a) mean square, main-plot standard error, critical
  difference and CV(a). Sub-plot standard error, critical difference and CV(b),
  and the returned ANOVA table itself, were not affected.
* `anova_spd_ab_main()`, `anova_spd_c_main_ab_sub()`, `anova_spd_ab_cd()` and
  `anova_sspd()` now test each effect against its correct error term. Previously
  these functions reported the raw combined ANOVA, in which the main-plot (and,
  for the split-split design, sub-plot) factors were tested against the residual
  rather than against the properly pooled main-plot / sub-plot error. The
  returned `anova_table` now gives the correct F-values and p-values for every
  stratum, and the objects carry the corresponding error mean squares, degrees
  of freedom and CVs. Verified against independent `aov()` `Error()`-strata
  models.
* `anova_spd_pooled()` and `anova_sspd_pooled()` (combined analysis over
  environments) now test each source against the correct pooled error. Both
  functions previously reported a single-stratum `anova()` table in which every
  effect was tested against one residual, so the environment main effect and the
  main-plot (and, for the split-split design, sub-plot) factors were tested
  against the sub-sub-plot residual rather than against their own pooled errors.
  Treating environments as fixed and replications as random within environment,
  the environment effect is now tested against the environment/replication error
  (E:Rep), the main-plot factors against the pooled main-plot error, and the
  sub-plot and sub-sub-plot factors against their respective pooled errors. The
  returned `anova_table` gives the correct F- and p-values for every stratum, the
  objects now carry per-stratum coefficients of variation, and `anova_sspd_pooled()`
  now returns factor means for plotting. Verified against independent `aov()`
  `Error()`-strata models.
* `perform_posthoc()` can now be called directly without pre-computing the error
  mean square and error degrees of freedom. When `mse` or `df_error` are not
  supplied they are derived from the fitted model. This fixes errors in the
  Bonferroni method and in the manual (no `agricolae`) fallbacks for LSD and
  Tukey.
* Fixed Dunnett's test: the comparison is now built against the actual treatment
  factor supplied by the user instead of a hard-coded factor name, so it works
  correctly when `multcomp` is installed.
* `verbose = FALSE` is now fully respected by `anova_crd()`, `anova_rbd()` and
  `anova_latin()`; these no longer print post-hoc and assumption-check output
  when silent operation is requested. Result objects are unchanged and still
  contain the full post-hoc and diagnostic components regardless of `verbose`.
* `arid_plot()` now draws factor or treatment means for every ANOVA design
  (factorial, three-way factorial, Latin square, all split-plot and
  split-split-plot variants, strip-plot, alpha lattice and the pooled designs),
  not only CRD and RBD. Designs with several factors are shown as a multi-panel
  bar chart, one panel per factor, with standard-deviation whiskers where
  available. Objects that expose no factor-means table return quietly with an
  informative message instead of an error.
* `arid_plot()` produces plots using base graphics for correlation objects (a
  heatmap of the correlation matrix) and stability objects (a bar chart of the
  integrated stability ranking); it previously returned without plotting.
* `export_results()` now supports `format = "csv"` in addition to `"xlsx"` and
  raises a clear error for unsupported formats.

## Other changes

* Added `grDevices`, `graphics`, `stats` and `utils` to `Imports`.

# aridagri 2.0.3

* Initial CRAN release (24 February 2026)
* 33 exported functions across 6 modules
* 16 ANOVA designs with proper error terms
* 7 post-hoc comparison tests with letter groupings
* 7 stability analysis methods with integrated ranking
* Thermal indices (GDD, HTU, PTU, HUE)
* Crop growth analysis (CGR, RGR, NAR, LAI)
* Nutrient use efficiency indices
* Path analysis and SEM via lavaan
* Publication-ready output with SE, CD, CV
* Unified verbose parameter for all functions
* Zero external dependencies (base R only)
