\name{rmtif_gumbel_recovery_study}
\alias{rmtif_gumbel_recovery_study}
\title{Summarize Gumbel regression recovery over simulated trials}
\description{
The summary contains the mean estimate, empirical SD, mean estimated SE, interval coverage, and convergence counts for every parameter and sample size. Trials with singular or nearly singular sensitivity are excluded from summary estimates and counted as failed fits.
}
\usage{
rmtif_gumbel_recovery_study(n_per_arm_values,
    replications,
    beta_true,
    tau,
    sigma_x,
    conf_level,
    
    seed)
}
\arguments{
\item{n_per_arm_values}{Per-arm sample sizes.}
\item{replications}{Monte Carlo trials per sample size.}
\item{beta_true}{True treatment and covariate coefficients.}
\item{tau}{Restriction horizon.}
\item{sigma_x}{Baseline covariate covariance matrix.}
\item{conf_level}{Confidence level.}
\item{seed}{Reproducibility seed.}
}
\value{
List containing replicate-level and summary data frames.
}
