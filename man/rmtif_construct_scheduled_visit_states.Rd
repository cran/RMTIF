\name{rmtif_construct_scheduled_visit_states}
\alias{rmtif_construct_scheduled_visit_states}
\title{Construct scheduled-visit ordered states from recurrent HF event records}
\description{
This entry point accepts the event-record schema used by the public HF example: patient identifier, event time, raw event code, treatment arm, and baseline etiology. Hospitalization since the previous scheduled visit defines the nonterminal state; observed death is an absorbing exact event.
}
\usage{
rmtif_construct_scheduled_visit_states(data,
    tau,
    visit_interval,
    raw_status,
    ordered_state,
    quantile_probs)
}
\arguments{
\item{data}{Event records with `id`, `time`, `status`, `Training`, and   `HF.etiology` columns.}
\item{tau}{Restriction horizon.}
\item{visit_interval}{Scheduled-visit spacing.}
\item{raw_status}{Named `censor`, `death`, `hospitalization` code vector.}
\item{ordered_state}{Named ordered-state code vector.}
\item{quantile_probs}{Three probabilities for the censor-time arm audit.}
}
\value{
List with constructed `rmtif_subject` paths and subject and arm audit data frames.
}
