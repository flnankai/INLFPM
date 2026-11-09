\name{ac_ss_analytic}
\alias{ac_ss_analytic}
\alias{ac_in_analytic}
\title{Finite-sample analytic reference tails for robust maximum tests}
\description{Fit a single robust maximum test with the finite-sample working calibration. Both computations are deterministic.}
\usage{
ac_ss_analytic(Y, F, calibration = "sidak", fit_tol = 1e-7,
               max_iter = 500L)
ac_in_analytic(Y, F, calibration = "sidak",
               pqr_profile = "heterogeneous")
}
\arguments{
  \item{Y}{Finite numeric time-by-asset matrix of returns.}
  \item{F}{Full-rank matrix of observed factors, with the same rows as \code{Y}, excluding the intercept.}
  \item{calibration}{Mapping of the marginal two-sided reference tail to the maximum. The default \code{"sidak"} uses one minus the Nth power of the marginal nonrejection probability. The diagnostic alternative \code{"gumbel"} maps the calibrated marginal tail through a normal quantile and a Gumbel approximation.}
  \item{fit_tol}{Convergence tolerance for the spatial-sign estimating equations.}
  \item{max_iter}{Maximum number of spatial-sign fit iterations.}
  \item{pqr_profile}{Use \code{"heterogeneous"} for the fitted time-variance profile or \code{"gaussian"} for the constant-profile PQR moment.}
}
\details{The SS procedure uses the formal published radial normalization, fits common time variances from projected residuals, and evaluates a frozen-weight Gaussian quadratic-ratio tail with an angular correction. IN uses PQR coordinate scales, an angular correction and a moment-matched Student reference. Both procedures estimate nuisance parameters from data. The resulting maximum p-values are working approximations. A fit or reference computation outside its domain raises an error.}
\value{A list with \code{p} (selected calibrated maximum p-value), \code{original_p} (original Gumbel p-value of the raw test), \code{raw} (uncentered squared maximum), \code{marginal_p}, estimated \code{lambda}, angular quantities \code{gamma} and \code{cJ}, and fit-specific diagnostics. For SS the latter include the radial moments and integration diagnostics. For IN they include the inverse-norm root, PQR moment, effective degrees of freedom and squared reference scale.}
\seealso{\code{\link{inmax_test}}}
