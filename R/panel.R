# Common factor projection, unrestricted residuals and nuisance estimates.
# All three are invariant under additive shifts of the asset intercepts.
inmax_prepare <- function(Y, F, bootstrap_reps, multipliers) {
  F <- as.matrix(F)
  Y <- as.matrix(Y)
  if (!is.numeric(Y) || !is.numeric(F) || any(!is.finite(Y)) ||
      any(!is.finite(F)) || nrow(Y) != nrow(F) || ncol(Y) < 2L ||
      ncol(F) < 1L || nrow(F) <= ncol(F) + 5L) {
    stop("Supply finite numeric Y and F with matching rows, N >= 2, at least one factor, and residual df > 4.")
  }
  bootstrap_reps <- ac_driver_positive_integer(bootstrap_reps, "bootstrap_reps")
  if (!is.null(multipliers) &&
      (!identical(dim(multipliers), c(bootstrap_reps, dim(Y))) ||
       anyNA(multipliers) || any(abs(multipliers) != 1))) {
    stop("multipliers must contain signs (-1 or 1) in a B by T by N array.")
  }
  pr <- inmax_projection(F)
  Z <- inmax_project(Y, pr)
  E <- Z - tcrossprod(pr$h, as.numeric(crossprod(pr$h, Z)) / pr$h2)
  bias <- tryCatch(inmax_bootstrap_bias(E, pr, bootstrap_reps, multipliers)$bias,
                   error = identity)
  pairs <- tryCatch(inmax_stat_pair_scales(E), error = identity)
  list(projection = pr, bias = bias, pairs = pairs)
}

inmax_evaluate <- function(Y, F, prepared) {
  pr <- prepared$projection
  Z <- inmax_project(Y, pr)
  p <- statistics <- setNames(rep(NA_real_, length(INMAX_METHODS)), INMAX_METHODS)
  failures <- list()
  diagnostics <- list()

  ols <- tryCatch(inmax_ols(Z, pr), error = identity)
  if (inherits(ols, "error")) {
    failures[INMAX_METHODS[1:3]] <- rep(list(conditionMessage(ols)), 3L)
  } else {
    p[1:3] <- ols$p
    statistics[1:3] <- c(ols$S, ols$M, ols$C)
  }

  ss <- tryCatch({
    if (inherits(prepared$bias, "error")) stop(prepared$bias)
    fit <- inmax_diag_sign(Z)
    trace <- inmax_split_trace(Y, pr, fit$D)
    (inmax_sign_q(fit$U, pr$h) - prepared$bias) / sqrt(2 * trace)
  }, error = identity)
  if (inherits(ss, "error")) {
    failures[["SS-SUM"]] <- conditionMessage(ss)
  } else {
    statistics["SS-SUM"] <- ss
    p["SS-SUM"] <- pnorm(ss, lower.tail = FALSE)
  }
  ss_max <- tryCatch(ac_ss_analytic(Y, F), error = identity)
  if (inherits(ss_max, "error")) {
    failures[["SS-MAX"]] <- conditionMessage(ss_max)
  } else {
    statistics["SS-MAX"] <- ss_max$raw
    p["SS-MAX"] <- ss_max$p
    diagnostics[["SS-MAX"]] <- ss_max
  }

  ins <- tryCatch({
    if (inherits(prepared$pairs, "error")) stop(prepared$pairs)
    inmax_stat_sum(Z, pr, cache = prepared$pairs)
  }, error = identity)
  if (inherits(ins, "error")) {
    failures[["IN-SUM"]] <- conditionMessage(ins)
  } else {
    statistics["IN-SUM"] <- ins$S
    p["IN-SUM"] <- ins$p
  }
  in_max <- tryCatch(ac_in_analytic(Y, F), error = identity)
  if (inherits(in_max, "error")) {
    failures[["IN-MAX"]] <- conditionMessage(in_max)
  } else {
    statistics["IN-MAX"] <- in_max$raw
    p["IN-MAX"] <- in_max$p
    diagnostics[["IN-MAX"]] <- in_max
  }

  for (family in c("SS", "IN")) {
    components <- p[paste0(family, c("-SUM", "-MAX"))]
    name <- paste0(family, "-CC")
    if (anyNA(components)) {
      failures[[name]] <- paste("Unavailable component:", paste(names(components)[is.na(components)], collapse = ", "))
    } else {
      cc <- inmax_cauchy(components, truncated = family == "SS")
      statistics[name] <- cc["statistic"]
      p[name] <- cc["p"]
    }
  }
  list(p.value = p, statistic = statistics, failures = failures, diagnostics = diagnostics)
}

inmax_test <- function(Y, F, bootstrap_reps = 100L, seed = NULL, multipliers = NULL) {
  Y <- as.matrix(Y)
  F <- as.matrix(F)
  if (!is.null(seed)) set.seed(seed)
  prepared <- inmax_prepare(Y, F, bootstrap_reps, multipliers)
  inmax_evaluate(Y, F, prepared)
}
