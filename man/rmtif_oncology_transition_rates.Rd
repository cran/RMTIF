\name{rmtif_oncology_transition_rates}
\alias{rmtif_oncology_transition_rates}
\title{Calculate conditional oncology transition intensities and jump probabilities}
\description{
Apply the supplied Cox-form transition model at a specified treatment arm and baseline covariate value. The next-state probabilities divide each transition intensity by the total exit rate of its origin state. They are undefined (`NA`) when all allowed exit rates from that state are zero.
}
\usage{
rmtif_oncology_transition_rates(design,
    arm,
    covariate)
}
\arguments{
\item{design}{Validated output of `rmtif_oncology_design()`.}
\item{arm}{Treatment-arm code, either control `0` or treatment `1`.}
\item{covariate}{Finite scalar baseline covariate value.}
}
\value{
The design's transition data frame with `intensity`, `total_exit_rate`, and `next_transition_probability` columns added.
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
rmtif_oncology_transition_rates(design, arm = 1, covariate = 0)
}
