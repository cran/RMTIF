\name{rmtif_simulate_ctmc_path}
\alias{rmtif_simulate_ctmc_path}
\title{Simulate a complete continuous-time Markov state path}
\description{
Wait at each state according to its generator exit intensity, then choose the next state in proportion to the off-diagonal intensities. An absorbing state ends transitions while retaining its state through the horizon.
}
\usage{
rmtif_simulate_ctmc_path(generator,
    initial_state,
    tau)
}
\arguments{
\item{generator}{Valid CTMC generator matrix.}
\item{initial_state}{Initial state index starting at zero.}
\item{tau}{Follow-up horizon.}
}
\value{
List with right-continuous `times` and `states` vectors.
}
\examples{
rates <- matrix(c(0, 0.17, 0.29, 0), nrow = 2, byrow = TRUE)
generator <- rmtif_ctmc_generator(rates)
set.seed(397)
rmtif_simulate_ctmc_path(generator, initial_state = 0, tau = 3.7)
}
