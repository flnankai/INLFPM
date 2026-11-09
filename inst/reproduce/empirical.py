"""Rolling and full-sample pricing tests on locally supplied monthly data."""
import os
for name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"

import argparse
import csv
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtr, stdtr

from calibration_math import Reference, sidak
from panel_math import (prepare_projection, ols_projected, sign_bootstrap_bias,
                        diagonal_sign_fit, split_trace, sign_q, cauchy_p, FitError)
from pair_scales import pair_diagonals, signal_pvalues, StatFitError

METHODS = ("LS-SUM", "LS-MAX", "LS-CC", "SS-SUM", "SS-MAX", "SS-CC",
           "IN-SUM", "IN-MAX", "IN-CC")
SEED = 20261109


def read_panel(prices, factors):
    """Read GBK prices (two metadata columns) and dated percentage-point factors."""
    raw = pd.read_csv(prices, encoding="gbk")
    px = raw.iloc[:, 2:].apply(pd.to_numeric, errors="raise")
    months = pd.to_datetime(px.columns, format="%Y/%m/%d").to_period("M")
    if (not months.is_monotonic_increasing or months.has_duplicates or
            not months.equals(pd.period_range(months[0], months[-1], freq="M"))):
        raise ValueError("Price columns must be unique, consecutive and chronological.")
    values = px.to_numpy(dtype=float).T
    if (not np.isfinite(values).all() or np.any(values <= 0)):
        raise ValueError("Supply a complete, positive monthly price panel.")
    fac = pd.read_csv(factors)
    fac.index = pd.to_datetime(fac["month"].astype(str), format="%Y%m").dt.to_period("M")
    if fac.index.has_duplicates:
        raise ValueError("Factor months must be unique.")
    months = months[1:]
    fac = fac.loc[months]
    rf = fac["RF"].to_numpy(dtype=float) / 100
    y = values[1:] / values[:-1] - 1 - rf[:, None]
    # Follow the market-column construction of the published application.
    f = np.column_stack([fac["Mkt-RF"].to_numpy(dtype=float) / 100 - rf, fac["SMB"], fac["HML"]])
    if not np.isfinite(f).all() or not np.isfinite(y).all():
        raise ValueError("Returns or factors contain nonfinite values.")
    return np.ascontiguousarray(y), np.ascontiguousarray(f), months


def panel_pvalues(y, f, start):
    """Nine tests on one window; start is its zero-based position in the panel."""
    pr = prepare_projection(f)
    ref = Reference(y, pr)
    ols = ols_projected(ref.z, pr)
    p = np.full(len(METHODS), np.nan)
    p[:3] = ols["p"]
    # Fixed by calendar position and window length, independent of scheduling.
    rng = np.random.default_rng(np.random.SeedSequence([SEED, len(y), start]))
    try:
        bias, _ = sign_bootstrap_bias(ref.e, pr, 100, rng)
        sign = diagonal_sign_fit(ref.z)
        trace = split_trace(y, pr, sign["D"])
        p[3] = ndtr(-(sign_q(sign["U"], pr.h) - bias) / np.sqrt(2 * trace))
    except FitError:
        pass
    try:
        p[4], _ = ref.spatial(ref.z)
    except FitError:
        pass
    p[5] = cauchy_p(p[3:5], truncated=True)[0]
    try:
        scales = pair_diagonals(ref.e)
        p[6] = signal_pvalues(ref.z, pr, [("null", 0, 0.0)], scales)[0]
    except StatFitError:
        pass
    try:
        alpha, sc, adj, m = ref.inverse()
        p[7] = sidak(2 * stdtr(m, -np.max(np.abs(alpha * sc)) * adj), y.shape[1])
    except FitError:
        pass
    p[8] = cauchy_p(p[6:8])[0]
    return p


def initialize(y, f):
    global RETURNS, FACTORS
    RETURNS, FACTORS = y, f


def window(task):
    length, start = task
    end = start + length
    return length, panel_pvalues(RETURNS[start:end], FACTORS[start:end], start)


def summarize(y, f, workers):
    windows = (48, 60, 72)
    n = len(y)
    if n <= max(windows):
        raise ValueError("The rolling design requires more than 72 return months.")
    tasks = [(length, start) for length in windows for start in range(n - length)]
    values = {length: [] for length in (*windows, n)}
    with Pool(workers, initializer=initialize, initargs=(y, f)) as pool:
        for length, p in pool.imap_unordered(window, tasks, chunksize=1):
            values[length].append(p)
    # The full-sample fit uses larger pairwise covariance arrays. Compute it
    # after releasing worker memory.
    values[n] = [panel_pvalues(y, f, 0)]
    for length, panel in values.items():
        p = np.asarray(panel)
        valid = np.isfinite(p).sum(axis=0)
        for level in (.01, .05):
            counts = np.sum(p < level, axis=0)
            for j, method in enumerate(METHODS):
                yield {"sample": "full" if length == n else "rolling",
                       "assets": y.shape[1], "months": length, "level": level,
                       "method": method, "attempted": len(panel), "valid": valid[j],
                       "rejections": counts[j], "rate": counts[j] / len(panel),
                       "full_pvalue": p[0, j] if length == n else ""}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prices", required=True, type=Path)
    parser.add_argument("--factors", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=25)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    y, f, dates = read_panel(args.prices, args.factors)
    records = list(summarize(y, f, args.workers))
    with args.output.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    print(f"Completed {dates[0]} to {dates[-1]}, {y.shape[1]} stocks, {len(records)} summary rows.")


if __name__ == "__main__":
    main()
