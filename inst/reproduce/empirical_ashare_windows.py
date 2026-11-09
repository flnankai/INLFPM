"""Fit CH3/CH4 on all stock series complete within each rolling window."""
import argparse
from pathlib import Path

from empirical import window_values, summary_records
from empirical_models import in_overlap, write_table
from empirical_ashares import mingshi_models
from ashare_inputs import load_ashare_returns
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", required=True, type=Path)
    parser.add_argument("--ch3", required=True, type=Path)
    parser.add_argument("--ch4", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--overlap", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=25)
    args = parser.parse_args()
    if args.workers < 1 or args.summary.resolve() == args.overlap.resolve():
        parser.error("Supply positive workers and two separate output paths.")
    returns = load_ashare_returns(args.source_dir)
    months = returns.index
    observed = np.isfinite(returns.to_numpy(dtype=float))
    summaries, overlaps = [], []
    for model, rf, factors in mingshi_models(args.ch3, args.ch4, months):
        y = np.ascontiguousarray(returns.subtract(rf, axis="index"), dtype=float)
        print(f"Fitting {model} on stocks observed throughout each window", flush=True)
        windows = window_values(y, np.ascontiguousarray(factors), args.workers,
                                include_last_window=True, full_sample=True,
                                windowwise_complete=True)
        ranges = {}
        for length, values in windows.items():
            counts = [int(observed[i:i + length].all(axis=0).sum()) for i in range(len(values))]
            ranges[length] = {"assets_min": min(counts), "assets_max": max(counts)}
        for row in summary_records(windows, len(months), assets=""):
            row.update(ranges[row["months"]])
            if row["assets_min"] == row["assets_max"]:
                row["assets"] = row["assets_min"]
            row.update(factor_model=model, first_month=str(months[0]), last_month=str(months[-1]))
            summaries.append(row)
        for row in in_overlap(windows, months, model):
            row.update(ranges[row["months"]])
            overlaps.append(row)
    write_table(args.summary, summaries)
    write_table(args.overlap, overlaps)
    print(f"Completed two models with {len(summaries)} estimates "
          f"and {len(overlaps)} paired comparisons.", flush=True)


if __name__ == "__main__":
    main()
