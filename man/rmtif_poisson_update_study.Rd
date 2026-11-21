\name{rmtif_poisson_update_study}
\alias{rmtif_poisson_update_study}
\title{Calibrate complete-path projection inference under Poisson updates}
\description{
Repeatedly simulate two arms of a reversible ordered-state process and compute the complete-path rmt-IF curve at each requested horizon. The population truth follows from the corresponding CTMC generators. The returned summary reports bias, empirical SD, projection SE, coverage and a two-sided zero-contrast rejection proportion for every horizon.
}
\usage{
rmtif_poisson_update_study(scenario,
    n_per_arm,
    replications,
    tau,
    increment,
    seed,
    design,
    
    conf_level,
    test_level,
    truth_grid_points)
}
\arguments{
\item{scenario}{Either `"Null"` or an alternative design.}
\item{n_per_arm}{Subjects per arm.}
\item{replications}{Monte Carlo replications.}
\item{tau}{Restricted horizon.}
\item{increment}{Evaluation-grid spacing.}
\item{seed}{Reproducibility seed.}
\item{design}{List with `transition_control`, `transition_treatment`,   `rates_control`, `rates_treatment`, `states`, `initial_state`, and   `absorbing_states`; all study values are supplied by the caller.}
\item{conf_level}{Confidence level.}
\item{test_level}{Rejection threshold.}
\item{truth_grid_points}{Number of numerical integration points for the   population CTMC contrast.}
}
\value{
List with replicate-level curves and a pointwise summary containing population truth, bias, empirical SD, mean projection SE, normal-interval coverage, and zero-null rejection proportion.
}
