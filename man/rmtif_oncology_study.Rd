\name{rmtif_oncology_study}
\alias{rmtif_oncology_study}
\title{Monte Carlo comparison of oncology endpoints on shared trials}
\description{
Simulate each scenario using a separate seed for every trial and bootstrap. All six methods analyse the same subjects within a replicate. Failed fits remain in the replicate file; rejection proportions and their Monte Carlo standard errors use the number of valid fits as the denominator. Optional scenario checkpoints permit an interrupted run to resume without changing the random streams. Parallel workers are allocated by scenario, not method.
}
\usage{
rmtif_oncology_study(scenarios,
    designs,
    n_per_arm,
    replications,
    tau,
    assessment_interval,
    
    assessment_jitter,
    censoring_rate,
    covariate_mean,
    covariate_sd,
    rmt_inference,
    
    regression_bootstrap,
    seed,
    conf_level,
    bootstrap_success_fraction,
    
    minimum_survival,
    test_level,
    workers = 1L,
    checkpoint_dir = NULL)
}
\arguments{
\item{scenarios}{Distinct scenario labels to simulate.}
\item{designs}{Named list of validated `rmtif_oncology_design()` results,   with one entry per requested scenario. List names and each design's   `scenario` label must agree.}
\item{n_per_arm}{Number of subjects in each arm.}
\item{replications}{Number of attempted trials per scenario.}
\item{tau}{Administrative horizon, in months.}
\item{assessment_interval}{Nominal interval between assessments.}
\item{assessment_jitter}{Half-width of uniform assessment-time jitter.}
\item{censoring_rate}{Exponential loss-to-follow-up rate; zero gives complete   follow-up through the administrative horizon.}
\item{covariate_mean}{Mean of the Gaussian baseline covariate.}
\item{covariate_sd}{Positive standard deviation of the Gaussian covariate.}
\item{rmt_inference}{Complete-path projection or censored IPCW bootstrap.}
\item{regression_bootstrap}{Bootstrap resamples for censored RMT-IF inference.}
\item{seed}{Master seed; replicate seeds are retained in the output.}
\item{conf_level}{Confidence level for fitted models.}
\item{bootstrap_success_fraction}{Required successful bootstrap fraction.}
\item{minimum_survival}{IPCW positivity threshold.}
\item{test_level}{Two-sided rejection threshold.}
\item{workers}{Number of scenario-level PSOCK workers.}
\item{checkpoint_dir}{Optional directory for resumable scenario checkpoints.}
}
\value{
List containing summary, replicate, transition-design, and runtime data frames. The example driver saves these outside the package.
}
