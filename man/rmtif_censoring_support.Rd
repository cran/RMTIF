\name{rmtif_censoring_support}
\alias{rmtif_censoring_support}
\title{Audit arm-specific censoring support for IPCW}
\description{
The aim is to check the empirical positivity condition before estimation. Survival is evaluated as the predictable left limit at `tau`; therefore a Kaplan--Meier curve that first reaches zero exactly at `tau` still supports integration through that endpoint, whereas a zero before `tau` does not.
}
\usage{
rmtif_censoring_support(subjects,
    tau)
}
\arguments{
\item{subjects}{List of `rmtif_subject` objects.}
\item{tau}{Candidate restriction horizon.}
}
\value{
One row per arm with the left-limit survival at `tau`, the first zero observed by `tau`, maximum empirically verified horizon, and a support flag. Absence of a zero by `tau` is not extrapolated beyond `tau`.
}
