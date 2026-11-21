\name{rmtif_complete_curve}
\alias{rmtif_complete_curve}
\title{Estimate a complete-data rmt-IF curve}
\description{
The function generates pointwise estimates and projection intervals for complete ordered-state paths.
}
\usage{
rmtif_complete_curve(subjects,
    times,
    conf_level)
}
\arguments{
\item{subjects}{Complete observed paths.}
\item{times}{Increasing restriction horizons.}
\item{conf_level}{Pointwise confidence level.}
}
\value{
Data frame with estimate, standard error, interval, and p-value.
}
