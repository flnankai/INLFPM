\name{inmax_test}
\alias{inmax_test}
\title{Joint intercept tests for a panel with observed factors}
\description{Compute sum, maximum and Cauchy-combination tests of the joint null that all intercepts are zero.}
\usage{
inmax_test(Y, F, bootstrap_reps = 100L, seed = NULL, multipliers = NULL)
}
\arguments{
  \item{Y}{Numeric matrix of returns, with times in rows and assets in columns. At least two assets are required.}
  \item{F}{Numeric factor matrix with the same rows as \code{Y}. Omit the intercept. The factor design must have full rank and leave more than four residual degrees of freedom after adding the intercept.}
  \item{bootstrap_reps}{Positive integer number of entrywise wild-bootstrap draws for SS-SUM.}
  \item{seed}{Optional random seed for the SS-SUM bootstrap.}
  \item{multipliers}{Optional array of fixed signs with dimensions \code{c(bootstrap_reps, nrow(Y), ncol(Y))}. Supplying the array makes the calculation deterministic.}
}
\details{
LS-SUM uses thresholded residual correlations. LS-MAX uses the centered squared OLS maximum and a Gumbel reference law.
SS-SUM uses a bias-corrected spatial-sign statistic with the ordered-pair split trace. IN-SUM uses a separate sign scale for each deleted pair and the design-normalized variance.
SS-MAX and IN-MAX use their fitted finite-sample reference tails and the Sidak maximum mapping. Each CC combines the two p-values from its own family and panel with equal weights. SS-CC truncates component Cauchy contributions at zero. LS-CC and IN-CC retain the signs of both contributions.
}
\value{A list with named vectors \code{p.value} and \code{statistic} in LS, SS and IN order (SUM, MAX, CC within each family); \code{failures}, a list of messages keyed by method; and \code{diagnostics}, containing the successful robust MAX fits. The robust MAX statistics are the raw, uncentered squared maxima. Method-specific failures produce missing values for the affected tests.}
\seealso{\code{\link{ac_ss_analytic}}, \code{\link{ac_in_analytic}}}
