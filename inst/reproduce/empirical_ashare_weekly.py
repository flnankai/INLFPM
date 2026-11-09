"""Compare CSI 300 and CSI 500 on common 52-trading-week windows."""
import argparse
from pathlib import Path

from empirical import window_values, summary_records
from empirical_models import in_overlap, write_table
from ashare_weekly_inputs import load_weekly_inputs
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csi300", required=True, type=Path)
    parser.add_argument("--csi500", required=True, type=Path)
    parser.add_argument("--factors", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--overlap", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=25)
    args = parser.parse_args()
    if args.workers < 1 or args.summary.resolve() == args.overlap.resolve():
        parser.error("Supply positive workers and two separate output paths.")

    panels = {market: load_weekly_inputs(path, args.factors)
              for market, path in [("CSI300", args.csi300), ("CSI500", args.csi500)]}
    weeks = panels["CSI300"][0].index
    summaries, overlaps = [], []
    for market, (returns, rf, models) in panels.items():
        if not returns.index.equals(weeks) or not rf.index.equals(weeks):
            raise ValueError("Both stock universes must share the same dates and risk-free series.")
        if not rf.equals(panels["CSI300"][1]):
            raise ValueError("The weekly risk-free returns differ between stock universes.")
        y = np.ascontiguousarray(returns.subtract(rf, axis="index"), dtype=float)
        for model in ("FF3", "FF5"):
            f = models[model]
            if not f.index.equals(weeks) or not f.equals(panels["CSI300"][2][model]):
                raise ValueError("Factor observations differ between the two universes.")
            print(f"Fitting {model} on {len(weeks)} weeks of {market}, N={y.shape[1]}", flush=True)
            windows = window_values(y, np.ascontiguousarray(f), args.workers, lengths=(52,),
                                    include_last_window=True, full_sample=False)
            for row in summary_records(windows, len(weeks), y.shape[1]):
                row["weeks"] = row.pop("months")
                del row["full_pvalue"]
                row.update(universe=market, factor_model=model,
                           first_week=str(weeks[0]), last_week=str(weeks[-1]))
                summaries.append(row)
            for row in in_overlap(windows, weeks, model):
                row["weeks"] = row.pop("months")
                row.update(universe=market, assets=y.shape[1])
                overlaps.append(row)
    write_table(args.summary, summaries)
    write_table(args.overlap, overlaps)
    print(f"Completed {len(summaries)} rates and {len(overlaps)} paired comparisons.", flush=True)


if __name__ == "__main__":
    main()
