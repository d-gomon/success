
<!-- README.md is generated from README.Rmd. Please edit that file -->

# success

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
[![arXiv](https://img.shields.io/badge/stat.AP-arXiv%3A2205.07618-B31B1B)](https://doi.org/10.48550/arXiv.2205.07618)
<!-- badges: end -->

# SUrvival Control Chart EStimation Software

The goal of the package is to allow easy applications of continuous time
CUSUM procedures on survival data. Specifically, the Biswas &
Kalbfleisch CUSUM (2008) and the CGR-CUSUM (2021).

Besides this, it allows for the construction of the Binary CUSUM chart
and funnel plot on survival data as well.

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

You can plot the figure with control limit
![h = 10](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;h%20%3D%2010 "h = 10")
by using:

``` r
plot(cgr, h = 10)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

And determine the runlength of the chart when using control limit
![h = 10](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;h%20%3D%2010 "h = 10"):

``` r
runlength(cgr, h = 10)
#> [1] 151
```

Hospital 1 would be detected by a CGR-CUSUM with control limit
![h = 10](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;h%20%3D%2010 "h = 10")
after
![151](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;151 "151")
days.

Alternatively, you can construct the CGR-CUSUM only until it crosses
control limit
![h = 10](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;h%20%3D%2010 "h = 10")
by:

``` r
cgr <- cgr_cusum(data = dat, coxphmod = tcoxmod, h = 10)
plot(cgr)
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

## References

The theory behind the methods in this package can be found in:

Gomon D., Putter H., Nelissen R.G.H.H., van der Pas S (2022):
[CGR-CUSUM: A Continuous time Generalized Rapid Response Cumulative Sum
chart](https://doi.org/10.48550/arXiv.2205.07618), *arXiv: a preprint*
