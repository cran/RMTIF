\name{rmtif_observe_ctmc_path}
\alias{rmtif_observe_ctmc_path}
\title{Apply loss to follow-up and terminal-event observation to a CTMC path}
\description{
A terminal state reached before censoring is observed exactly and carried forward. Earlier loss to follow-up truncates the observed nonterminal path.
}
\usage{
rmtif_observe_ctmc_path(path,
    arm,
    tau,
    censor_time,
    death_state,
    covariates = numeric(),
    
    id = NULL)
}
\arguments{
\item{path}{Complete CTMC path from `rmtif_simulate_ctmc_path`.}
\item{arm}{Treatment indicator.}
\item{tau}{Restriction horizon.}
\item{censor_time}{Independent loss-to-follow-up time.}
\item{death_state}{Absorbing terminal-state label.}
\item{covariates}{Optional numeric baseline covariate vector.}
\item{id}{Optional subject identifier.}
}
\value{
One observed `rmtif_subject`.
}
