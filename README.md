
<!-- README.md is generated from README.Rmd. Please edit that file -->

# success <img src="hexagon/success_hex.png" align="right" width="120" />

<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/d-gomon/success/workflows/R-CMD-check/badge.svg)](https://github.com/d-gomon/success/actions/)
[![CRAN
status](https://www.r-pkg.org/badges/version/success)](https://CRAN.R-project.org/package=success)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/success)](https://cran.r-project.org/package=success)
[![GitHub Repo
stars](https://img.shields.io/github/stars/d-gomon/success?style=social)](https://github.com/d-gomon/success)
[![Biostatistics](https://img.shields.io/badge/Biostatistics-kxac041-%23003365)](https://doi.org/10.1093/biostatistics/kxac041)
<!-- badges: end -->

<!-- [![arXiv](https://img.shields.io/badge/stat.AP-arXiv%3A2205.07618-B31B1B)](https://doi.org/10.48550/arXiv.2205.07618) -->

# SUrvival Control Chart EStimation Software

The goal of the package is to allow easy applications of continuous time
CUSUM procedures on survival data. Specifically, the Biswas &
Kalbfleisch CUSUM (2008) and the CGR-CUSUM (Gomon et al. 2022).

Besides continuous time procedures, it is also possible to construct the
Bernoulli (binary) CUSUM and funnel plot (Spiegelhalter 2005) on
survival data.

## Installation

You can install the released version of success from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("success")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("d-gomon/success")
```

## CGR-CUSUM Example

This is a basic example which shows you how to construct a CGR-CUSUM
chart on a hospital from the attached data set “surgerydat”:

``` r
dat <- subset(surgerydat, unit == 1)
exprfit <- as.formula("Surv(survtime, censorid) ~ age + sex + BMI")
tcoxmod <- coxph(exprfit, data = surgerydat)

cgr <- cgr_cusum(data = dat, coxphmod = tcoxmod, stoptime = 200)
plot(cgr)
```

<img src="man/figures/README-success-1.png" width="100%" />

You can plot the figure with control limit `h = 10` by using:

``` r
plot(cgr, h = 10)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

And determine the runlength of the chart when using control limit
`h = 10`:

``` r
runlength(cgr, h = 10)
#> [1] 151
```

Using a control limit of `h = 10` Hospital 1 would be detected by a
CGR-CUSUM 151 days after the first patient entered the study.

## References

Gomon D., Putter H., Nelissen R.G.H.H., van der Pas S (2022):
[CGR-CUSUM: A Continuous time Generalized Rapid Response Cumulative Sum
chart](https://doi.org/10.1093/biostatistics/kxac041), *Biostatistics*

Biswas P. and Kalbfleisch J.D. (2008): [A risk-adjusted CUSUM in
continuous time based on the Cox
model](https://doi.org/10.1002/sim.3216), *Statistics in Medicine*

Spiegelhalter D.J. (2005): [Funnel plots for comparing institutional
performance](https://doi.org/10.1002/sim.1970), *Statistics in Medicine*
