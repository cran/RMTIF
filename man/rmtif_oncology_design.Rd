\name{rmtif_oncology_design}
\alias{rmtif_oncology_design}
\alias{rmtif_oncology_parameters}
\title{Validate an ordered-state oncology simulation design}
\description{
All baseline intensities, covariate and treatment log-rate effects, state labels, and the initial and terminal states are supplied by the caller. The named clinical states define the response-confirmation endpoints. No manuscript simulation settings are inserted by this constructor.
}
\usage{
rmtif_oncology_design(transition,
    states,
    initial_state,
    death_state,
    scenario)
}
\arguments{
\item{transition}{Data frame with `from`, `to`, `baseline_rate`,   `covariate_log_rate`, and `treatment_log_rate` columns, one row per allowed   transition. Rates are nonnegative; log-rate effects are finite.}
\item{states}{Named, increasing, zero-based integer codes with names `CR`,   `PR`, `SD`, `PD`, and `death`, in clinical order from better to worse.}
\item{initial_state}{Initial state code.}
\item{death_state}{Absorbing terminal state code, equal to `states["death"]`.}
\item{scenario}{Nonempty scenario label used only in returned metadata.}
}
\value{
A validated list with `scenario`, `transition`, `states`, `initial_state`, and `death_state`.
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
design$transition
}
