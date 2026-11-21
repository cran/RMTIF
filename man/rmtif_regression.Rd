\name{rmtif_regression}
\alias{rmtif_regression}
\title{Fit covariate-adjusted pairwise rmt-IF regression}
\description{
The model uses `treatment + (X_treatment - X_control)` as its pairwise design, so the first coefficient is the treatment effect for equal baseline covariates.  The score is exactly the manuscript score `(win - loss * exp(beta^T difference)) * difference`.  With arm-specific IPCW, a stratified bootstrap refits both censoring curves and is the primary inference; the fixed-weight sandwich is retained as a diagnostic. A singular or nearly singular score-sensitivity matrix makes coefficients and intervals unidentified: the function raises an error instead of reporting a pseudoinverse-based, spuriously narrow interval. Near singular means the smallest singular value is at most the largest times square-root machine precision; remove redundant covariates and refit.
}
\usage{
rmtif_regression(subjects,
    tau,
    censoring = c("none",
    "arm_km"),
    bootstrap,
    conf_level,
    
    bootstrap_success_fraction,
    seed = NULL,
    minimum_survival = 0,
    beta_initial = NULL)
}
\arguments{
\item{subjects}{List of `rmtif_subject` objects.}
\item{tau}{Restriction horizon.}
\item{censoring}{Either `"none"` or terminal-aware `"arm_km"`.}
\item{bootstrap}{Stratified bootstrap replicates. `0` returns fixed-weight   diagnostic inference only.}
\item{conf_level}{Confidence level.}
\item{bootstrap_success_fraction}{Minimum fraction of finite bootstrap   fits needed to report bootstrap inference.}
\item{seed}{Optional reproducibility seed.}
\item{minimum_survival}{Positivity threshold for each censoring curve.}
\item{beta_initial}{Optional starting coefficient vector.}
}
\value{
An `rmtif_regression` fit with treatment and covariate coefficients.
}
