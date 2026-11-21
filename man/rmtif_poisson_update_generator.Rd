\name{rmtif_poisson_update_generator}
\alias{rmtif_poisson_update_generator}
\title{Convert a Poisson attempted-event update model to a CTMC generator}
\description{
Attempted updates occur at the state-specific Poisson rates. Each attempt draws a destination from the corresponding row of `transition`. A diagonal draw leaves the state unchanged; it is not a state transition. The resulting generator has off-diagonal entries `rates[i] * transition[i, j]` and rows summing to zero. This construction accommodates reversible movement and positive diagonal update probabilities.
}
\usage{
rmtif_poisson_update_generator(transition,
    rates)
}
\arguments{
\item{transition}{Row-stochastic attempted-update probability matrix.}
\item{rates}{Positive attempted-update rates, one per state.}
}
\value{
CTMC generator matrix with the same row and column names.
}
\examples{
transition <- matrix(c(0.37, 0.63, 0.44, 0.56),
                     nrow = 2, byrow = TRUE)
rmtif_poisson_update_generator(transition, rates = c(1.3, 2.1))
}
