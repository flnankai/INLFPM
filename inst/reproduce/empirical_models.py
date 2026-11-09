"""Paired nine-test comparisons across five common U.S. factor models."""
import argparse
import csv
from pathlib import Path

from empirical import window_values, summary_records
from empirical_updated import monthly_factors, read_updated_panel
import numpy as np


def designs(three, rf, months, ff5, momentum):
    five = monthly_factors(ff5, months, ["Mkt-RF", "SMB", "HML", "RMW", "CMA", "RF"])
    mom = monthly_factors(momentum, months, ["Mom"])
    if not np.array_equal(rf, five[:, 5]):
        raise ValueError("Monthly risk-free rates in FF3 and FF5 differ.")
    # Use the five-factor file's definitions, including the FF5 size factor.
    return {
        "CAPM": three[:, :1],
        "FF3": three,
        "Carhart4": np.column_stack([three, mom]),
        "FF5": five[:, :5],
        "FF5+UMD": np.column_stack([five[:, :5], mom]),
    }


def first_window(mask, p, months, length, label):
    indices = np.flatnonzero(mask)
    example = dict(start="", end="", sum_p="", max_p="", cc_p="")
    if len(indices):
        i = indices[0]
        example.update(start=str(months[i]), end=str(months[i + length - 1]),
                       sum_p=p[i, 0], max_p=p[i, 1], cc_p=p[i, 2])
    return {f"first_{label}_{key}": value for key, value in example.items()}


def in_overlap(windows, months, model):
    # All comparisons use complete p-value triples on the same windows.
    for length, samples in windows.items():
        if length == len(months):
            continue
        p = np.asarray(samples)[:, 6:9]
        valid = np.isfinite(p).all(axis=1)
        for level in (.01, .05):
            reject_sum, reject_max, reject_cc = (p < level).T
            sum_only = valid & reject_sum & ~reject_max
            max_only = valid & reject_max & ~reject_sum
            both = valid & reject_sum & reject_max
            neither = valid & ~reject_sum & ~reject_max
            cc_only = neither & reject_cc
            row = {
                "factor_model": model, "months": length, "level": level,
                "attempted": len(p), "valid": int(valid.sum()),
                "sum_rejections": int(np.sum(valid & reject_sum)),
                "max_rejections": int(np.sum(valid & reject_max)),
                "cc_rejections": int(np.sum(valid & reject_cc)),
                "sum_only": int(sum_only.sum()), "max_only": int(max_only.sum()),
                "both_components": int(both.sum()), "neither_component": int(neither.sum()),
                "cc_only": int(cc_only.sum()),
                "cc_adds_vs_sum": int(np.sum(valid & reject_cc & ~reject_sum)),
                "cc_loses_vs_sum": int(np.sum(valid & ~reject_cc & reject_sum)),
                "cc_adds_vs_max": int(np.sum(valid & reject_cc & ~reject_max)),
                "cc_loses_vs_max": int(np.sum(valid & ~reject_cc & reject_max)),
                "cc_keeps_sum_only": int(np.sum(sum_only & reject_cc)),
                "cc_keeps_max_only": int(np.sum(max_only & reject_cc)),
            }
            if (row["cc_only"] or np.any(both & ~reject_cc) or
                row["cc_rejections"] - row["sum_rejections"] !=
                row["cc_keeps_max_only"] - (row["sum_only"] - row["cc_keeps_sum_only"])):
                raise ArithmeticError("Inconsistent paired Cauchy rejection accounting.")
            row.update(first_window(max_only, p, months, length, "max_only"))
            row.update(first_window(cc_only, p, months, length, "cc_only"))
            yield row


def write_table(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fit_models(y, models, months, workers):
    """Apply all nine tests to identical full and rolling samples for every model."""
    summaries, overlaps = [], []
    for model, f in models.items():
        print(f"Fitting {model} on {len(months)} months and {y.shape[1]} stocks", flush=True)
        windows = window_values(y, np.ascontiguousarray(f), workers,
                                include_last_window=True, full_sample=True)
        for row in summary_records(windows, len(months), y.shape[1]):
            row.update(factor_model=model, first_month=str(months[0]), last_month=str(months[-1]))
            summaries.append(row)
        overlaps.extend(in_overlap(windows, months, model))
    return summaries, overlaps


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prices", required=True, type=Path)
    parser.add_argument("--ff3", required=True, type=Path)
    parser.add_argument("--ff5", required=True, type=Path)
    parser.add_argument("--momentum", required=True, type=Path)
    parser.add_argument("--calendar", type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--overlap", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=25)
    args = parser.parse_args()
    if args.workers < 1 or args.summary.resolve() == args.overlap.resolve():
        parser.error("Supply positive workers and two separate output paths.")
    y, three, rf, months = read_updated_panel(args.prices, args.ff3, "2005-02",
                                               "2026-09", args.calendar)
    models = designs(three, rf, months, args.ff5, args.momentum)
    summaries, overlaps = fit_models(y, models, months, args.workers)
    write_table(args.summary, summaries)
    write_table(args.overlap, overlaps)
    print(f"Completed {len(models)} models with {len(summaries)} rates and {len(overlaps)} paired comparisons.")


if __name__ == "__main__":
    main()
