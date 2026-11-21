\name{rmtif_ipcw_simulation_study}
\alias{rmtif_ipcw_simulation_study}
\title{Run the dedicated differential-censoring IPCW study}
\description{
The returned manuscript-ready summary contains truth, mean, bias, empirical SD, mean bootstrap SE, normal-interval coverage, rejection probability, and fit counts for every requested scenario/allocation combination.
}
\usage{
rmtif_ipcw_simulation_study(allocations,
    scenarios,
    trials,
    bootstrap,
    tau,
    seed,
    minimum_survival,
    
    censoring_rates,
    parameters,
    conf_level,
    bootstrap_success_fraction,
    
    test_level,
    truth_grid_points)
}
\arguments{
\item{allocations}{List of named `c(n0, n1)` vectors.}
\item{scenarios}{Null and/or treatment-benefit scenarios.}
\item{trials}{Number of independently simulated trials.}
\item{bootstrap}{Stratified bootstrap replicates per trial.}
\item{tau}{Restriction horizon.}
\item{seed}{Reproducibility seed.}
\item{minimum_survival}{Positivity threshold.}
\item{censoring_rates}{Named or ordered arm-specific exponential   loss-to-follow-up rates used in every generated trial.}
\item{parameters}{Markov generators and state labels.}
\item{conf_level}{Confidence level.}
\item{bootstrap_success_fraction}{Required successful bootstrap fraction.}
\item{test_level}{Rejection threshold.}
\item{truth_grid_points}{Number of numerical integration points for the   population CTMC contrast.}
}
\value{
List containing replicate-level and summary data frames.
}
