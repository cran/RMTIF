# RMTIF

Estimate restricted mean time in favor of treatment from ordered state paths.
Lower numeric states denote better outcomes. Estimation is available for
complete paths and paths subject to independent censoring; inference includes
projection standard errors, normal intervals, two-sided tests, and stratified
bootstrap methods that re-estimate the censoring distribution.

Install the source archive with base R:

```r
install.packages("RMTIF_0.3.1.tar.gz", repos = NULL, type = "source")
library(RMTIF)
```

For example, obtain projection inference from a matrix of pairwise net times
whose rows identify control subjects and columns identify treatment subjects:

```r
net <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
rmtif_projection_inference(net, conf_level = 0.95)
```

The R source is copied byte for byte from the accompanying
`functions_v1.2.R`. The separate `examples_v1.2.R` contains manuscript design
settings and writes the manuscript tables. Long manuscript Monte Carlo runs
are excluded from package examples and tests. Function help includes concise
examples using independently chosen illustration inputs.

## Method references

The restricted-time estimand follows Mao (2023), "On restricted mean time in
favor of treatment", Biometrics 79, 61-72, doi:10.1111/biom.13570. The
pairwise regression is motivated by Mao and Wang (2021), "A class of
proportional win-fractions regression models for composite outcomes",
Biometrics 77, 1265-1275, doi:10.1111/biom.13382.

The manuscript describing this package's reversible-path extensions is
unpublished and currently has no public URL or DOI. The references above
identify the published foundations; package help documents the implemented
estimators and functions.

Study functions return summary, replicate, design, and runtime data frames.
They run quietly in both serial and parallel execution. Users can explicitly
print or save the returned information they need.
