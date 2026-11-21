\name{rmtif_simulate_ipcw_trial}
\alias{rmtif_simulate_ipcw_trial}
\title{Simulate one dedicated IPCW trial}
\description{
The aim is to evaluate the marginal arm-specific IPCW estimator under a reversible nonterminal process, absorbing death, unequal allocation, and differential but within-arm independent loss to follow-up.
}
\usage{
rmtif_simulate_ipcw_trial(n0,
    n1,
    scenario,
    tau,
    censoring_rates,
    parameters,
    seed = NULL)
}
\arguments{
\item{n0}{Control sample size.}
\item{n1}{Treatment sample size.}
\item{scenario}{`"Null"` or `"Treatment benefit"`.}
\item{tau}{Restriction horizon.}
\item{censoring_rates}{Arm-specific exponential loss-to-follow-up rates.}
\item{parameters}{Markov generators and state labels.}
\item{seed}{Optional seed.}
}
\value{
List of observed subjects and the two population generators.
}
