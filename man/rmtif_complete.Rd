\name{rmtif_complete}
\alias{rmtif_complete}
\title{Estimate complete-data rmt-IF with projection inference}
\description{
Average pairwise treatment-minus-control net times over every control and treatment subject pair. Subject-level row and column projections yield the two-sample U-statistic standard error. The function reports a normal confidence interval and the z statistic and two-sided p-value for a zero contrast; no independent-increment assumption across time is used.
}
\usage{
rmtif_complete(subjects,
    tau,
    conf_level)
}
\arguments{
\item{subjects}{List of `rmtif_subject` objects with complete paths to `tau`.}
\item{tau}{Restriction horizon.}
\item{conf_level}{Confidence level.}
}
\value{
An `rmtif_fit` list containing `estimate`, projection `se`, `conf_int`, zero-contrast `statistic` and `p_value`, mean pairwise `win` and `loss` times, and row and column projections.
}
\examples{
subjects <- list(
  rmtif_subject(0, 0, 5, tau = 3.7),
  rmtif_subject(0, 0, 5, tau = 3.7),
  rmtif_subject(1, 0, 4, tau = 3.7),
  rmtif_subject(1, 0, 4, tau = 3.7)
)
rmtif_complete(subjects, tau = 3.7, conf_level = 0.95)
}
