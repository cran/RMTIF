\name{rmtif_hfmock_analysis}
\alias{rmtif_hfmock_analysis}
\title{Analyze the public `hfmock` scheduled-visit process}
\description{
If `data` is omitted, the function loads `hfmock` only when the optional `Wcompo` package is installed.  No synthetic replacement is silently used. The analysis reports arm counts and terminal-aware IPCW bootstrap inference.
}
\usage{
rmtif_hfmock_analysis(data = NULL,
    tau,
    visit_interval,
    bootstrap,
    seed,
    raw_status,
    
    ordered_state,
    quantile_probs,
    conf_level,
    bootstrap_success_fraction,
    
    minimum_survival = 0,
    projection_diagnostic = FALSE)
}
\arguments{
\item{data}{Optional `hfmock` data frame.}
\item{tau}{Restriction horizon.}
\item{visit_interval}{Scheduled visit spacing.}
\item{bootstrap}{Stratified bootstrap replicates.}
\item{seed}{Reproducibility seed.}
\item{raw_status}{Named raw event-code vector; see `rmtif_hfmock_subjects`.}
\item{ordered_state}{Named ordered-state vector; see `rmtif_hfmock_subjects`.}
\item{quantile_probs}{Censor-time audit probabilities; see   `rmtif_hfmock_subjects`.}
\item{conf_level}{Normal confidence level.}
\item{bootstrap_success_fraction}{Required fraction of successful bootstrap   resamples before reporting inference.}
\item{minimum_survival}{Positivity threshold.}
\item{projection_diagnostic}{Whether to retain the dense fixed-weight   projection diagnostic.}
}
\value{
List with construction audit and an `rmtif_ipcw` fit.
}
