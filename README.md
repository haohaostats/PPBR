# ppbr

`ppbr` implements posterior predictive bracketing and resolution for phase I
dose-finding trials with ordered doses and binary toxicity outcomes.

[Launch the PPBR Shiny web application](https://01a02ea1-beac-e827-0c00-bcca8d89e653.share.connect.posit.cloud/)

## Installation

```r
remotes::install_github("haohaostats/PPBR")
```

## Shiny web application

Use PPBR through the public
[Shiny web application](https://01a02ea1-beac-e827-0c00-bcca8d89e653.share.connect.posit.cloud/).
No local R installation is required.

## Core workflow

```r
library(ppbr)

design <- ppbr_design(
  dose = c(1, 3, 6, 8, 10),
  target = 0.10,
  cohort_size = 3,
  max_sample = 33,
  overdose_threshold = 0.20
)

fit <- ppbr(
  design,
  n = c(3, 10, 12, 8, 0),
  dlt = c(0, 0, 2, 0, 0),
  current_dose = 8
)

summary(fit)
ppbr_next(fit)
plot(fit, type = "all")
```

## Evaluate a proposed design

```r
oc <- ppbr_simulate(
  design,
  scenarios = list(
    c(0.03, 0.08, 0.15, 0.25, 0.40),
    c(0.01, 0.04, 0.10, 0.20, 0.35)
  ),
  nsim = 1000,
  seed = 2026
)

summary(oc)
plot(oc)
```

The V1.0 public interface is intentionally small: `ppbr_design()`, `ppbr()`,
`ppbr_next()`, and `ppbr_simulate()`, plus standard `summary()` and `plot()`
methods.
