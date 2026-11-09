# Explicit opt-in to the full experiment. See ?INLFPM::inmax_run.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L || args[1L] != "--run") {
  stop("Usage: Rscript run.R --run output.csv [workers]", call. = FALSE)
}
workers <- if (length(args) == 3L) as.numeric(args[3L]) else 25L
invisible(INLFPM::inmax_run(workers = workers, output = args[2L]))
