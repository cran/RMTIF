\name{rmtif_ipcw_equivalence_check}
\alias{rmtif_ipcw_equivalence_check}
\title{Check the exact IPCW state-mass implementation against dense pair averaging}
\description{
This audit helper evaluates one observed data set by both algebraically identical routes.  It is useful for regression testing the faster bootstrap calculation; it does not fit a bootstrap or change the estimator.
}
\usage{
rmtif_ipcw_equivalence_check(subjects,
    tau,
    minimum_survival = 0,
    tolerance = 1e-10)
}
\arguments{
\item{subjects}{List of observed `rmtif_subject` objects.}
\item{tau}{Restriction horizon.}
\item{minimum_survival}{Positivity threshold for each arm-specific curve.}
\item{tolerance}{Absolute numerical comparison tolerance.}
}
\value{
One-row data frame with dense, direct state-mass, and cached state-mass IPCW means, including a deterministic duplicate-subject bootstrap resample.
}
