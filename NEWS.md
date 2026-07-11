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
