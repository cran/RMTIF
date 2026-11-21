\name{rmtif_ctmc_generator}
\alias{rmtif_ctmc_generator}
\title{Build a continuous-time Markov generator from transition intensities}
\description{
The caller supplies nonnegative off-diagonal transition intensities in a square matrix. Diagonal values are ignored and replaced by the negative off-diagonal row sums. This permits reversible movement among nonterminal states and an absorbing terminal state when its off-diagonal row is zero.
}
\usage{
rmtif_ctmc_generator(off_diagonal_rates)
}
\arguments{
\item{off_diagonal_rates}{Square matrix of transition intensities.}
}
\value{
Valid CTMC generator with zero row sums.
}
\examples{
rates <- matrix(c(0, 0.17, 0.29, 0), nrow = 2, byrow = TRUE)
rmtif_ctmc_generator(rates)
}
