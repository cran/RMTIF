\name{rmtif_true_ctmc_contrast}
\alias{rmtif_true_ctmc_contrast}
\title{Compute a population CTMC rmt-IF contrast}
\description{
The aim is to obtain simulation truth without Monte Carlo error.  Transition probabilities are computed by uniformization and the pairwise state contrast is integrated numerically.
}
\usage{
rmtif_true_ctmc_contrast(generator_control,
    generator_treatment,
    initial_state,
    tau,
    grid_points)
}
\arguments{
\item{generator_control}{Control CTMC generator.}
\item{generator_treatment}{Treatment CTMC generator.}
\item{initial_state}{Common initial state, numbered from zero.}
\item{tau}{Restriction horizon.}
\item{grid_points}{Number of equally spaced integration points, supplied   by the caller to control numerical accuracy.}
}
\value{
Population treatment-minus-control rmt-IF in time units.
}
\examples{
rates <- matrix(c(0, 0.17, 0.29, 0), nrow = 2, byrow = TRUE)
generator <- rmtif_ctmc_generator(rates)
rmtif_true_ctmc_contrast(generator, generator, initial_state = 0,
                         tau = 3.7, grid_points = 51)
}
