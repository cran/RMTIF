\name{rmtif_hfmock_subjects}
\alias{rmtif_hfmock_subjects}
\title{Convert `Wcompo::hfmock` event records to scheduled-visit state paths}
\description{
At caller-specified visits, the better nonterminal state means no hospitalization since the prior visit, the worse nonterminal state means at least one hospitalization, and death is an exact absorbing state. An end-of-record observation between visits censors the scheduled-visit process at the most recent completed visit. The raw data do not distinguish every such observation as loss to follow-up versus administrative closure. A record extending to or beyond the requested horizon supports carrying the most recent scheduled state to that horizon, including a horizon between scheduled visits.
}
\usage{
rmtif_hfmock_subjects(data,
    tau,
    visit_interval,
    raw_status,
    ordered_state,
    quantile_probs)
}
\arguments{
\item{data}{A data frame with `id`, `time`, `status`, `Training`, and   `HF.etiology`, such as `Wcompo::hfmock`.}
\item{tau}{Restriction horizon.}
\item{visit_interval}{Scheduled-visit spacing.}
\item{raw_status}{Named vector with `censor`, `death`, and   `hospitalization` raw event codes.}
\item{ordered_state}{Named vector with `no_hospitalization`,   `hospitalization`, and `death` ordered-state codes.}
\item{quantile_probs}{Three requested probabilities for the arm audit's   lower, median, and upper censor-time summaries.}
}
\value{
List with `rmtif_subject` objects and raw-data audit summaries.
}
