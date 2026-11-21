\name{rmtif_projection_inference}
\alias{rmtif_projection_inference}
\title{Compute two-sample projection inference from pairwise net times}
\description{
The estimate is the mean of a control-by-treatment net-time matrix. Its standard error combines the empirical variance of row means divided by the control sample size with the variance of column means divided by the treatment sample size. A normal confidence interval and two-sided test of zero contrast use the supplied confidence level.
}
\usage{
rmtif_projection_inference(pairwise_net_times,
    conf_level)
}
\arguments{
\item{pairwise_net_times}{Finite numeric matrix with control subjects in   rows and treatment subjects in columns.}
\item{conf_level}{Confidence level for the normal interval.}
}
\value{
List with `estimate`, `se`, `conf_int`, `statistic`, `p_value`, `row_projection`, and `column_projection`.
}
\examples{
net <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
rmtif_projection_inference(net, conf_level = 0.95)
}
