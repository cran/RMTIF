\name{rmtif_state_at}
\alias{rmtif_state_at}
\title{Evaluate a subject's right-continuous state}
\description{
This helper is intended for transparent checking of the path construction. It returns `NA` after unobserved loss to follow-up, but continues to return the absorbing terminal state after an observed death.
}
\usage{
rmtif_state_at(subject,
    time)
}
\arguments{
\item{subject}{An `rmtif_subject`.}
\item{time}{One or more nonnegative evaluation times.}
}
\value{
A numeric vector of states, with `NA` where the path is unavailable.
}
