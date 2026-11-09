\name{inmax_run}
\alias{inmax_run}
\title{Run the calibrated nine-method experiment}
\description{Generate the null and signal designs, combining component p-values on each panel and accumulating rejection and failure counts in memory.}
\usage{
inmax_run(replications = 1000L, workers = 25L, bootstrap_reps = 100L,
          master_seed = 20261109L, level = 0.05, output = NULL)
}
\arguments{
  \item{replications}{Number of independent base panels per null design cell.}
  \item{workers}{Number of PSOCK worker processes. Tasks are split across design cells and replication batches on all platforms.}
  \item{bootstrap_reps}{Wild-bootstrap draws per base panel for SS-SUM. The bias estimate is shared across intercept shifts of the same panel.}
  \item{master_seed}{Positive integer used to derive a reproducible seed for each design cell and replication.}
  \item{level}{Rejection threshold strictly between zero and one.}
  \item{output}{Optional path for a single final summary CSV, in an existing directory.}
}
\details{Six time-length and panel-width combinations are crossed with four error laws. Short panels of width 100 and 200 also enter the varying-sparsity design. Short panels of width 200 enter the signal-strength design. R uses its own random stream; independent implementations of the same design can have different realized rejection counts. Computation starts only on an explicit call.}
\value{A data frame by design and method, including number attempted, valid fits, failures, rejections, rejection rate per attempted panel, its Monte Carlo standard error, the conditional-on-valid rate and standard error, and bounds treating all failures as nonrejections or rejections.}
\seealso{\code{\link{inmax_test}}}
