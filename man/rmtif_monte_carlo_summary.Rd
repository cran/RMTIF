\name{rmtif_monte_carlo_summary}
\alias{rmtif_monte_carlo_summary}
\title{Summarize Monte Carlo estimator and interval performance}
\description{
For each requested design cell, calculate the mean estimate, bias, empirical standard deviation, mean reported standard error, interval coverage of the supplied truth, two-sided zero-test rejection fraction, and successful-fit count. Additional trial-level columns can be averaged over all generated trials, as with observed loss-to-follow-up fractions.
}
\usage{
rmtif_monte_carlo_summary(replicates,
    group_columns,
    estimate_column,
    se_column,
    truth_column,
    
    lower_column,
    upper_column,
    p_value_column,
    test_level,
    extra_mean_columns = character())
}
\arguments{
\item{replicates}{Data frame with one row per simulated trial.}
\item{group_columns}{Character names of the columns defining design cells.}
\item{estimate_column}{Name of the numeric estimate column.}
\item{se_column}{Name of the numeric estimated-standard-error column.}
\item{truth_column}{Name of the numeric population-truth column.}
\item{lower_column}{Name of the lower confidence-limit column.}
\item{upper_column}{Name of the upper confidence-limit column.}
\item{p_value_column}{Name of the two-sided zero-test p-value column.}
\item{test_level}{Rejection threshold for the p-value.}
\item{extra_mean_columns}{Optional trial-level numeric columns to average.}
}
\value{
One row per design cell with truth, mean estimate, bias, empirical SD, mean estimated SE, coverage, rejection, successful-fit count, total trial count, and means of requested extra columns.
}
