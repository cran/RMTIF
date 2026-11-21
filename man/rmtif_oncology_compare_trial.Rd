\name{rmtif_oncology_compare_trial}
\alias{rmtif_oncology_compare_trial}
\title{Compare six covariate-adjusted tests in one oncology trial}
\description{
OS/PFS use Cox Wald tests, OR and CDC use logistic Wald tests, study-specific ordinal BOR uses a proportional-odds Wald test, and RMT-IF tests its pairwise treatment coefficient.  Every method uses the same generated subjects.
}
\usage{
rmtif_oncology_compare_trial(simulation,
    rmt_inference,
    regression_bootstrap,
    conf_level,
    
    bootstrap_success_fraction,
    minimum_survival,
    test_level,
    bootstrap_seed = NULL)
}
\arguments{
\item{simulation}{Output of `rmtif_simulate_oncology_trial()`.}
\item{rmt_inference}{Complete-path projection or censored IPCW bootstrap.}
\item{regression_bootstrap}{IPCW bootstrap count (ignored for complete data).}
\item{conf_level}{Confidence level.}
\item{bootstrap_success_fraction}{Required successful bootstrap fraction.}
\item{minimum_survival}{IPCW positivity threshold.}
\item{test_level}{Two-sided rejection threshold.}
\item{bootstrap_seed}{Separate seed for the stratified bootstrap.}
}
\value{
Six-row long-format method result.
}
