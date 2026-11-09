"""Run paired A-share pricing tests on the common monthly panel."""
import argparse
from pathlib import Path

from empirical_models import fit_models, write_table
from ashare_inputs import load_ashare_inputs
import numpy as np
import pandas as pd


def mingshi_design(path, months, field_names):
    table = pd.read_excel(path, sheet_name="Returnseries",
                          usecols=["mnthdt", "rf_mon", *field_names])
    dates = pd.to_datetime(table.pop("mnthdt").astype(str), format="%Y%m%d")
    if not dates.dt.is_month_end.all():
        raise ValueError(f"{path.name} contains a non-month-end date.")
    table.index = dates.dt.to_period("M")
    if table.index.has_duplicates:
        raise ValueError(f"Duplicate monthly factor row in {path.name}.")
    table = table.reindex(months).astype(float)
    if not np.isfinite(table.to_numpy()).all():
        raise ValueError(f"Incomplete or invalid factor rows in {path.name}.")
    return table.pop("rf_mon"), table[list(field_names)]


def mingshi_models(ch3, ch4, months):
    for name, path, keys in [("CH3", ch3, ("mktrf", "SMB", "VMG")),
                             ("CH4", ch4, ("mktrf", "SMB", "VMG", "PMO"))]:
        yield name, *mingshi_design(path, months, keys)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", required=True, type=Path,
                        help="Directory with the documented return and factor inputs")
    parser.add_argument("--summary", required=True, type=Path,
                        help="Output CSV for full-sample p-values and rolling rejection rates")
    parser.add_argument("--overlap", required=True, type=Path,
                        help="Output CSV for paired IN-SUM, IN-MAX and IN-CC rejection counts")
    parser.add_argument("--workers", type=int, default=25)
    parser.add_argument("--ch3", type=Path, help="Mingshi CH3 monthly factor workbook")
    parser.add_argument("--ch4", type=Path, help="Mingshi CH4 monthly factor workbook")
    args = parser.parse_args()
    if args.workers < 1 or args.summary.resolve() == args.overlap.resolve():
        parser.error("Supply positive workers and two separate output paths.")
    if (args.ch3 is None) != (args.ch4 is None):
        parser.error("Provide both monthly CH3 and CH4 workbooks for the supplementary models.")
    returns, rf, models = load_ashare_inputs(args.source_dir)
    months = returns.index
    if not rf.index.equals(months) or any(not f.index.equals(months) for f in models.values()):
        raise ValueError("Every model must align with the same stock return months.")
    y = np.ascontiguousarray(returns.subtract(rf, axis="index"), dtype=float)
    summaries, overlaps = fit_models(y, models, months, args.workers)
    if args.ch3 is not None:
        for name, model_rf, factors in mingshi_models(args.ch3, args.ch4, months):
            y = np.ascontiguousarray(returns.subtract(model_rf, axis="index"), dtype=float)
            rows, counts = fit_models(y, {name: factors}, months, args.workers)
            summaries.extend(rows)
            overlaps.extend(counts)
    write_table(args.summary, summaries)
    write_table(args.overlap, overlaps)
    model_count = len({row["factor_model"] for row in summaries})
    print(f"Completed {model_count} models with {len(summaries)} estimates "
          f"and {len(overlaps)} paired comparisons.", flush=True)


if __name__ == "__main__":
    main()
