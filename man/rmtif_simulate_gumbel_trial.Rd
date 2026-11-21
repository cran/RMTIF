\name{rmtif_simulate_gumbel_trial}
\alias{rmtif_simulate_gumbel_trial}
\title{Simulate a Gumbel pairwise-regression data-generating mechanism}
\description{
Generate Gaussian baseline covariates with the supplied covariance matrix, treatment allocation, and independent standard Gumbel errors. The resulting continuous latent utility is constant over follow-up and is represented by its negative as an ordered state, so smaller states mean better outcomes. This is a regression-equation sanity check rather than a multistate DGM.
}
\usage{
rmtif_simulate_gumbel_trial(n_per_arm,
    beta_true,
    tau,
    sigma_x,
    seed = NULL)
}
\arguments{
\item{n_per_arm}{Subjects per arm.}
\item{beta_true}{Named vector: treatment intercept followed by slopes.}
\item{tau}{Restriction horizon.}
\item{sigma_x}{Positive-definite covariance matrix for baseline covariates.}
\item{seed}{Optional seed.}
}
\value{
Complete `rmtif_subject` list and true parameter.
}
\examples{
trial <- rmtif_simulate_gumbel_trial(
  n_per_arm = 12, beta_true = c(0.2, 0.3), tau = 3.7,
  sigma_x = matrix(1), seed = 397)
rmtif_complete(trial$subjects, tau = 3.7, conf_level = 0.95)
}
