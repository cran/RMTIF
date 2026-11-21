\name{rmtif_simulate_gumbel_regression}
\alias{rmtif_simulate_gumbel_regression}
\title{Simulate the constant-outcome Gumbel regression DGM}
\description{
This descriptive entry point generates the same Gaussian-covariate, Gumbel-error pairwise-regression trial as `rmtif_simulate_gumbel_trial`. All treatment, slope, covariance, horizon, and sample-size inputs come from the caller.
}
\usage{
rmtif_simulate_gumbel_regression(n_per_arm,
    beta_true,
    tau,
    sigma_x,
    seed = NULL)
}
\arguments{
\item{n_per_arm}{Subjects per treatment arm.}
\item{beta_true}{Named treatment-intercept and covariate-slope vector.}
\item{tau}{Restricted time horizon.}
\item{sigma_x}{Positive-definite baseline-covariate covariance matrix.}
\item{seed}{Optional reproducibility seed.}
}
\value{
List with complete ordered-state subjects, true parameter vector, and restriction horizon.
}
