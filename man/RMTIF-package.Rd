\name{RMTIF-package}
\alias{RMTIF-package}
\alias{RMTIF}
\docType{package}
\title{Restricted Mean Time in Favor of Treatment}
\description{Estimation, inference, and parameterized simulation for
ordered multistate trajectories. Lower numeric state values mean better
outcomes; net time equals treatment win time minus treatment loss time.}
\details{Use \code{rmtif_subject} to represent observed trajectories,
\code{rmtif_complete} for complete-path inference, \code{rmtif_ipcw}
for censoring-adjusted inference with a refitted stratified bootstrap,
and \code{rmtif_regression} for pairwise regression. Design settings and
Monte Carlo sizes are supplied explicitly by the caller.}
\references{
Mao, L. (2023). "On restricted mean time in favor of treatment."
Biometrics, 79(1), 61-72. \doi{10.1111/biom.13570}.

Mao, L. and Wang, T. (2021). "A class of proportional win-fractions
regression models for composite outcomes." Biometrics, 77(4), 1265-1275.
\doi{10.1111/biom.13382}.
}
\seealso{\code{\link{rmtif_subject}}, \code{\link{rmtif_complete}},
\code{\link{rmtif_ipcw}}, \code{\link{rmtif_regression}}}
