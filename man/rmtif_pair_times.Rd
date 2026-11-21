\name{rmtif_pair_times}
\alias{rmtif_pair_times}
\title{Compute pairwise restricted win, loss, and net times}
\description{
This is the exact interval-by-interval calculation for two piecewise-constant paths. With `censoring = "none"`, both paths must be complete through `tau`. With `censoring = "arm_km"`, the grid also contains every censoring time from `reference_subjects`, and the integrand uses the terminal-aware, arm-specific inverse censoring weight described in the revised manuscript.
}
\usage{
rmtif_pair_times(control,
    treatment,
    tau,
    censoring = c("none",
    "arm_km"),
    reference_subjects = NULL,
    
    minimum_survival = 0)
}
\arguments{
\item{control}{One arm-0 `rmtif_subject`.}
\item{treatment}{One arm-1 `rmtif_subject`.}
\item{tau}{Restriction horizon.}
\item{censoring}{Either `"none"` or `"arm_km"`.}
\item{reference_subjects}{Full sample used to estimate censoring curves.}
\item{minimum_survival}{Positivity threshold for each censoring survival.}
}
\value{
Named vector with `win`, `loss`, and `net` in time units.
}
