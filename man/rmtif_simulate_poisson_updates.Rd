\name{rmtif_simulate_poisson_updates}
\alias{rmtif_simulate_poisson_updates}
\title{Simulate finite-state Poisson attempted-event update trajectories}
\description{
At each state an exponential waiting time with the specified attempted-event rate precedes a destination draw from `transition`. A self-update advances time while retaining the same state. The output is one complete observed ordered-state trajectory per subject; the supplied state order determines favorable versus unfavorable states in later rmt-IF calculations.
}
\usage{
rmtif_simulate_poisson_updates(n,
    transition,
    rates,
    tau,
    states,
    initial_state,
    arm,
    absorbing_states,
    
    seed = NULL)
}
\arguments{
\item{n}{Number of subjects.}
\item{transition}{Row-stochastic attempted-update probability matrix.}
\item{rates}{Positive attempted-update rates, one per state.}
\item{tau}{Follow-up horizon.}
\item{states}{Ordered numeric state labels.}
\item{initial_state}{Initial state label.}
\item{arm}{Treatment arm, 0 or 1.}
\item{absorbing_states}{Optional explicit absorbing state labels.}
\item{seed}{Optional seed.}
}
\value{
List of complete `rmtif_subject` trajectories.
}
\examples{
transition <- matrix(c(0.37, 0.63, 0.44, 0.56),
                     nrow = 2, byrow = TRUE)
paths <- rmtif_simulate_poisson_updates(
  n = 4, transition = transition, rates = c(1.3, 2.1),
  tau = 3.7, states = c(0, 1), initial_state = 0,
  arm = 0, absorbing_states = numeric(), seed = 397)
length(paths)
}
