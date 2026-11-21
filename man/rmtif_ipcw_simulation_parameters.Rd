\name{rmtif_ipcw_simulation_parameters}
\alias{rmtif_ipcw_simulation_parameters}
\title{Prepare control and treatment CTMC generators for a censoring study}
\description{
The design supplies every transition intensity and treatment multiplier. Additional trial design fields are passed through for convenience.
}
\usage{
rmtif_ipcw_simulation_parameters(design)
}
\arguments{
\item{design}{List containing `rates_control`, `benefit_multiplier`,   `initial_state`, `death_state`, `tau`, `primary_censoring_rates`,   `allocations`, `trials`, and `bootstrap`.}
}
\value{
Named list of control, null-treatment and benefit-treatment CTMC generators, state labels, and supplied trial settings.
}
