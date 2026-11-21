\name{rmtif_simulate_oncology_trial}
\alias{rmtif_simulate_oncology_trial}
\title{Simulate one shared Appendix C oncology trial}
\description{
All six methods are derived from this same patient-level realization. Nominal visits are generated far enough to cover the jitter window and then visits outside `(0,tau]`, after loss to follow-up, or at/after death are dropped.  The RMT-IF representation is compressed to actual observed state changes; the complete visit history is retained for OR and BOR confirmation.
}
\usage{
rmtif_simulate_oncology_trial(n_per_arm,
    design,
    tau,
    assessment_interval,
    assessment_jitter,
    
    censoring_rate,
    covariate_mean,
    covariate_sd,
    seed = NULL)
}
\arguments{
\item{n_per_arm}{Exact subjects in each randomized arm.}
\item{design}{Validated output of `rmtif_oncology_design()`.}
\item{tau}{Administrative restriction horizon.}
\item{assessment_interval}{Nominal imaging interval.}
\item{assessment_jitter}{Half-width of independent uniform visit jitter.}
\item{censoring_rate}{Exponential LTFU rate; zero specifies complete follow-up.}
\item{covariate_mean}{Mean of the Gaussian baseline covariate.}
\item{covariate_sd}{Positive standard deviation of the Gaussian covariate.}
\item{seed}{Reproducibility seed.}
}
\value{
Shared endpoint data, RMT-IF subjects, visits, and latent paths.
}
\examples{
transition <- data.frame(
  from = c(2, 2, 1, 0, 3), to = c(1, 3, 0, 4, 4),
  baseline_rate = c(0.17, 0.23, 0.19, 0.13, 0.27),
  covariate_log_rate = rep(0, 5),
  treatment_log_rate = rep(0, 5))
design <- rmtif_oncology_design(transition,
  states = c(CR = 0L, PR = 1L, SD = 2L, PD = 3L, death = 4L),
  initial_state = 2L, death_state = 4L, scenario = "illustration")
trial <- rmtif_simulate_oncology_trial(
  n_per_arm = 8, design = design, tau = 3.7,
  assessment_interval = 0.8, assessment_jitter = 0.07,
  censoring_rate = 0, covariate_mean = 0, covariate_sd = 1,
  seed = 397)
names(trial)
}
