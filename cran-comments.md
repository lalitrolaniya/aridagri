## Submission

This is an update from 2.0.4 to 2.1.0.

Authors: Hemant Poonia, a co-author of the article describing the package
<doi:10.1016/j.simpa.2026.100851>, is added as an author. He and Monika
Punia have passed away and are marked as deceased in Authors@R.

Corrected: the relative efficiency of the RCBD over a CRD, the error terms
and silent mode of the combined RCBD analysis, the F-tests of the strip-plot
design, the ANOVA table returned by the split-plot function, Shukla's
stability variance, the AMMI stability value, the Eberhart-Russell
deviation variance, the factor used for the variance test in
check_assumptions(), and an undefined efficiency index in nue_calculate().
Details are in NEWS.md.

New: repeated measures ANOVA with sphericity corrections, a multi-trait
publication table, GGE biplots, letter groupings computed in the package,
adjusted means for augmented and alpha lattice designs, response curve
models, ggplot2 figures, six simulated datasets and a vignette.

'nlme' (a recommended package) is now imported. Unused suggested packages
were removed.

Words that may be flagged as possibly misspelled in DESCRIPTION are author
and method names (Federer, Piepho, Eberhart, Shukla, Wricke, Binns, Yan,
Geisser, Huynh, Feldt, Keuls, Scheffe, Tukey, Dunnett, Bonferroni,
Mitscherlich), abbreviations defined in the text (CRD, RCBD, SEm, SEd,
REML, AMMI, GGE) and technical terms (biplots, genotype, helio-thermal,
intra-block, et al.).

## Test environments

* Local: Ubuntu 24.04, R 4.6.2 (R CMD check --as-cran)
* Local: Ubuntu 24.04, R 4.6.2 with none of the suggested packages
* (to run before submission) win-builder R-devel, and the GitHub Actions
  workflow for macOS, Windows and Ubuntu (release, devel, oldrel-1)

## R CMD check results

0 errors | 0 warnings | 1 note

* The note lists the package website https://lalitrolaniya.github.io/aridagri/
  as not found when checked before the site was first published. The site
  is published by GitHub Actions and is live at submission time.

## Reverse dependencies

There are no reverse dependencies.
