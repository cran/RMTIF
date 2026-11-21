\name{rmtif_ipcw}
\alias{rmtif_ipcw}
\title{Estimate terminal-aware arm-specific IPCW rmt-IF}
\description{
The aim is unbiased marginal rmt-IF estimation under unconditional independent loss to follow-up within treatment arm.  The estimator uses two arm-specific censoring Kaplan-Meier curves, carries observed death forward, and integrates over the union of state, death, and censoring times.  Primary inference is a treatment-stratified bootstrap that refits both censoring curves in every replicate.  The fixed-weight projection is returned only as a diagnostic because it omits censoring-curve estimation uncertainty.
}
\usage{
rmtif_ipcw(subjects,
    tau,
    bootstrap,
    conf_level,
    bootstrap_success_fraction,
    
    seed = NULL,
    minimum_survival = 0,
    projection_diagnostic = TRUE,
    allow_insufficient_bootstrap = FALSE)
}
\arguments{
\item{subjects}{List of observed `rmtif_subject` objects.}
\item{tau}{Restriction horizon.}
\item{bootstrap}{Number of stratified bootstrap replicates.}
\item{conf_level}{Confidence level.}
\item{bootstrap_success_fraction}{Minimum fraction of finite bootstrap   estimates needed to report bootstrap inference.}
\item{seed}{Optional reproducibility seed.}
\item{minimum_survival}{Positivity threshold for each arm-specific curve.}
\item{projection_diagnostic}{Whether to form the dense original-sample pair   matrix for the fixed-weight projection diagnostic.  Set to `FALSE` in a   Monte Carlo study when only the IPCW point estimate and refitted-bootstrap   inference are required.}
\item{allow_insufficient_bootstrap}{If `TRUE`, return a diagnostic object   when too few bootstrap resamples retain support.  Its bootstrap   standard error and intervals are left undefined; this option is intended   only to audit a boundary horizon, not to report primary inference.}
}
\value{
An `rmtif_fit` list with bootstrap inference when requested.
}
