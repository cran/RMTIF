\name{rmtif_subject}
\alias{rmtif_subject}
\title{Construct one observed ordered-state trajectory}
\description{
The function creates the common subject representation used by every rmt-IF estimator in this file.  `times` and `states` define a right-continuous, piecewise-constant observed path.  An observed terminal event is inserted at `death_time` and carried forward; it never truncates the outcome path. `censor_time` is loss to follow-up or administrative end, and `censor_event` distinguishes loss to follow-up from administrative end.
}
\usage{
rmtif_subject(arm,
    times,
    states,
    tau,
    covariates = numeric(),
    censor_time = tau,
    
    censor_event = FALSE,
    death_time = Inf,
    death_state = NULL,
    id = NULL)
}
\arguments{
\item{arm}{Treatment indicator, 0 for control or 1 for treatment.}
\item{times}{Numeric state-observation/change times and must include baseline.}
\item{states}{Ordered states corresponding to `times`; smaller is better.}
\item{tau}{Planned restriction horizon.}
\item{covariates}{Optional numeric baseline covariate vector.}
\item{censor_time}{Loss-to-follow-up or administrative time.}
\item{censor_event}{Whether `censor_time` is a loss-to-follow-up event.}
\item{death_time}{Observed terminal-event time, or `Inf` if unobserved/absent.}
\item{death_state}{Ordered absorbing terminal state, or `NULL` if not used.}
\item{id}{Optional subject identifier.}
}
\value{
An object of class `rmtif_subject`.
}
\examples{
path <- rmtif_subject(arm = 1, times = c(0, 1.5),
                      states = c(2, 1), tau = 3.7)
rmtif_state_at(path, c(0, 1.5, 3))
}
