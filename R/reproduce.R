inmax_counts <- function(cell, bootstrap_reps, master_seed, level) {
  T <- cell$T
  N <- cell$N
  distribution <- cell$distribution
  settings <- inmax_settings(T, N)
  valid <- rejected <- matrix(0, nrow(settings), length(INMAX_METHODS))

  for (r in seq.int(cell$first, cell$last)) {
    set.seed(ac_simulation_seed(T, N, distribution, r, master_seed))
    panel <- ac_draw_data(T, N, distribution)
    prepared <- inmax_prepare(panel$Y, panel$F, bootstrap_reps, NULL)
    for (i in seq_len(nrow(settings))) {
      Y <- panel$Y
      s <- settings$s[i]
      if (s > 0L) {
        active <- seq_len(s)
        Y[, active] <- Y[, active, drop = FALSE] + sqrt(settings$delta[i] / s)
      }
      p <- inmax_evaluate(Y, panel$F, prepared)$p.value
      usable <- is.finite(p)
      valid[i, ] <- valid[i, ] + usable
      rejected[i, ] <- rejected[i, ] + (usable & p < level)
    }
  }

  cbind(valid = as.vector(t(valid)), rejected = as.vector(t(rejected)))
}

inmax_summary <- function(cell, counts, replications) {
  settings <- inmax_settings(cell$T, cell$N)
  cells <- settings[rep(seq_len(nrow(settings)), each = length(INMAX_METHODS)), , drop = FALSE]
  n <- counts[, "valid"]
  k <- counts[, "rejected"]
  rate <- k / replications
  rate_valid <- ifelse(n > 0, k / n, NA_real_)
  data.frame(T = cell$T, N = cell$N, distribution = cell$distribution,
             cells, method = rep(INMAX_METHODS, nrow(settings)),
             replications, valid = n, failures = replications - n, rejections = k, rate,
             MCSE = sqrt(rate * (1 - rate) / replications), rate_among_valid = rate_valid,
             MCSE_among_valid = sqrt(rate_valid * (1 - rate_valid) / n),
             rejection_lower = rate, rejection_upper = (k + replications - n) / replications,
             row.names = NULL)
}

inmax_run <- function(replications = 1000L, workers = 25L, bootstrap_reps = 100L,
                      master_seed = 20261109L, level = 0.05, output = NULL) {
  replications <- ac_driver_positive_integer(replications, "replications")
  workers <- ac_driver_positive_integer(workers, "workers")
  bootstrap_reps <- ac_driver_positive_integer(bootstrap_reps, "bootstrap_reps")
  master_seed <- ac_driver_positive_integer(master_seed, "master_seed")
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("level must be between zero and one.")
  if (!is.null(output) && (!is.character(output) || length(output) != 1L ||
                           is.na(output) || !nzchar(output) || !dir.exists(dirname(output))))
    stop("output must be a file path in an existing directory.")
  grid <- ac_simulation_grid()
  # Keep several independent batches per worker for heterogeneous panel costs.
  per_batch <- max(1L, ceiling(replications / workers))
  work <- expand.grid(cell = seq_len(nrow(grid)), first = seq.int(1, replications, by = per_batch))
  batches <- grid[work$cell, , drop = FALSE]
  batches$first <- work$first
  batches$last <- pmin(work$first + per_batch - 1, replications)
  tasks <- lapply(seq_len(nrow(batches)), function(i) batches[i, , drop = FALSE])
  if (workers == 1L) {
    counts <- lapply(tasks, inmax_counts, bootstrap_reps, master_seed, level)
  } else {
    cluster <- parallel::makePSOCKcluster(min(workers, length(tasks)))
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    # Load from the same installation even when it is in a private library.
    parallel::clusterCall(cluster, function(paths) {
      .libPaths(paths)
      loadNamespace("INLFPM")
      NULL
    }, .libPaths())
    counts <- parallel::parLapplyLB(cluster, tasks, inmax_counts,
                                    bootstrap_reps, master_seed, level, chunk.size = 1L)
  }
  cells <- lapply(seq_len(nrow(grid)), function(i) {
    total <- Reduce(`+`, counts[work$cell == i])
    inmax_summary(grid[i, , drop = FALSE], total, replications)
  })
  summary <- do.call(rbind, cells)
  rownames(summary) <- NULL
  if (!is.null(output)) utils::write.csv(summary, output, row.names = FALSE)
  summary
}
