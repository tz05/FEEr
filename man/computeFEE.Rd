% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fee_ind.R
\name{computeFEE}
\alias{computeFEE}
\title{Calculate Functional Extension and Evenness (FEE) index}
\usage{
computeFEE(
  pool_traits,
  abund,
  dis_metric = c("euclidean", "manhattan"),
  abundWeighted = TRUE,
  poolBased = FALSE,
  doFEE0 = FALSE,
  user_ecdf = NULL,
  doParallel = FALSE
)
}
\arguments{
\item{pool_traits}{a matrix (or data frame) of numeric functional traits for
species pool. Columns are for different traits and rows are for different
species in the species pool.}

\item{abund}{matrix (or data frame) of species abundances in the communities
of interest. Columns are for different species in the species pool (in the
same species order as in \code{pool_traits}. Thus its number of columns
needs to be equal to the number of rows in \code{pool_traits}). Each row
is for a community. For species not presented in a community, their
abundance values should be assigned 0 or \code{NA}.}

\item{dis_metric}{string specifying the scheme of quantifying species
dissimilarity in trait space. The currently available options are
\code{"euclidean"} (the default) and \code{"manhattan"}. Euclidean
distances are root sum-of-squares of differences, and manhattan distances
are the sum of absolute differences.}

\item{abundWeighted}{logical; indicating whether the abundance values in
\code{abund} are considered as weights in FEE calculation. The default is
\code{TRUE}.}

\item{poolBased}{logical; indicating whether the calculation is based on
species pool. The default is \code{FALSE}, indicating no prior knowledge
about trait distribution (i.e., adopting uniform distribution for trait
values to build null communities). If \code{TRUE}, \code{pool_traits}
is used to build null communities.}

\item{doFEE0}{logical; indicating whether the output includes FEE0 (raw
functional extension and evenness metric). The default is \code{FALSE}.}

\item{user_ecdf}{string specifying the path and name of the file that saves
user's results of empirical cumulative distribution function (eCDF) of
FEE0. Under the same setting (\code{dis_metric} and \code{pool_traits}),
calculated eCDFs could be reused to save time. If the argument is
\code{NULL} (default), no user's eCDF will be adopted or saved, and the
species-pool-based eCDFs will be calculated from scratch. If the file path
is provided, the eCDFs contained in this file will be loaded and used as
needed, and when the eCDF of a specific combination of species number and
trait number is not available in the file, it will be calculated from
scratch and then saved to the file.}

\item{doParallel}{logical; indicating whether the calculation of FEE is
conducted in parallel. The default value is \code{FALSE}, i.e.,
sequentially calculating the FEE index for the communities specified in
\code{abund}.}
}
\value{
The FEE index values of the communities. If \code{doFEE0} is
  \code{TRUE}, the output includes two lists: the raw FEE0 metrics and the
  FEE indices. Specifically, FEE values for monocultures are \code{NA}.
}
\description{
Computes the FEE index for the communities specified in \code{abund}, based
on the trait data specified in \code{pool_traits}.
}
\details{
This function is the core function of the package. The \code{pool_traits}
and \code{abund} are key arguments. The \code{dis_metric},
\code{abundWeighted}, and \code{poolBased} are essential arguments, which
together define the FEE calculation scheme. The other three arguments are
auxiliary arguments, which might affect calculation efficiency or alter
output format, but do not affect FEE calculation results.

Preparing eCDFs usually takes the most time of FEE calculation. To increase
calculation speed, some pre-calculated eCDFs are enclosed in the package.
For the cases with 1 to 6 traits, 2 to 90 species, and under the assumption
of no prior knowledge about trait distribution (i.e., \code{poolBased} is
\code{FALSE}, thus species pool does not affect eCDFs), these pre-calculated
eCDFs are used in the FEE calculation. For other scenarios, especially
species-pool-based calculation (\code{poolBased = TRUE}), eCDFs are not
prepared and need to be calculated from scratch. In this case, specifying
\code{user_ecdf} allows the newly calculated eCDFs to be saved and reused
as needed. As a result, the saved eCDFs can considerably reduce calculation
time for communities that are under the same settings of \code{pool_traits},
\code{dis_metric}, and \code{poolBased}.
}
\examples{
d_trait <- cbind(tr1 = rnorm(10, 0, 1), tr2 = rnorm(10, 3, 2))
d_comm <- rbind(comm1 = c(0:9),
                comm2 = rep(1, 10),
                comm3 = runif(10, 0, 5),
                comm4 = rchisq(10, 2))
computeFEE(d_trait, d_comm, abundWeighted = TRUE)
}
