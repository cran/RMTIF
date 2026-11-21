\name{rmtif_complete_state_count_curve}
\alias{rmtif_complete_state_count_curve}
\title{Evaluate the exact finite-state projection curve by state counts}
\description{
This is the fast implementation of the complete-path pairwise estimator. For each interval it obtains row and column projections from strict empirical state counts in the opposite arm, algebraically matching the dense pair matrix without assuming independent increments across time.
}
\usage{
rmtif_complete_state_count_curve(subjects,
    times,
    conf_level)
}
\arguments{
\item{subjects}{Complete observed ordered-state paths.}
\item{times}{Increasing positive restriction horizons.}
\item{conf_level}{Pointwise normal confidence level.}
}
\value{
Data frame with horizon, rmt-IF estimate, projection SE, interval endpoints, and two-sided zero-contrast p-value.
}
