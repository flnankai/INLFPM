"""Pricing tests for a fixed stock panel with adjusted closes and standard FF3."""
import argparse
import csv
import io
from pathlib import Path
from zipfile import ZipFile

from empirical import summarize
import numpy as np
import pandas as pd


def csv_text(path):
    if path.suffix.lower() == ".zip":
        with ZipFile(path) as archive:
            csv_names = [n for n in archive.namelist() if n.lower().endswith(".csv")]
            if len(csv_names) != 1:
                raise ValueError("A factor ZIP must contain one CSV file.")
            return archive.read(csv_names[0]).decode("utf-8-sig")
    return path.read_text(encoding="utf-8-sig")


def monthly_ff3(path, months):
    # French distributes monthly and annual observations together. Only
    # six-digit dates belong to the monthly table.
    reader = csv.reader(io.StringIO(csv_text(path)))
    fields = ["Mkt-RF", "SMB", "HML", "RF"]
    header = None
    rows = {}
    for raw in reader:
        row = [s.strip() for s in raw]
        if header is None:
            if all(k in row for k in fields):
                header = [row.index(k) for k in fields]
            continue
        if not row:
            continue
        key = row[0].replace("-", "")
        if len(key) != 6 or not key.isdigit():
            continue
        if key in rows:
            raise ValueError(f"Duplicate monthly factor row: {key}")
        rows[key] = [float(row[i]) for i in header]
    if header is None:
        raise ValueError("Expected columns Mkt-RF, SMB, HML and RF.")
    missing = [str(m) for m in months if m.strftime("%Y%m") not in rows]
    if missing:
        raise ValueError(f"Factors missing for {', '.join(missing)}.")
    values = np.array([rows[m.strftime("%Y%m")] for m in months])
    if not np.isfinite(values).all() or np.any(values <= -99.99):
        raise ValueError("Missing or invalid monthly factor values.")
    return values / 100


def calendar_ends(path, month_range):
    # Accept ISO session dates or the eight-digit dates in French daily CSVs.
    sessions = []
    for row in csv.reader(io.StringIO(csv_text(path))):
        if row:
            key = row[0].strip().replace("-", "")
            if len(key) == 8 and key.isdigit():
                sessions.append(key)
    dates = pd.to_datetime(sessions, format="%Y%m%d")
    if not len(dates) or dates.has_duplicates:
        raise ValueError("Expected a trading calendar of unique dated sessions.")
    return pd.Series(dates).groupby(dates.to_period("M")).max().reindex(month_range)


def read_updated_panel(prices, factors, first, last, calendar=None):
    with prices.open(encoding="utf-8-sig", newline="") as stream:
        header = next(csv.reader(stream))
    if len(header) < 3 or len(set(header)) != len(header):
        raise ValueError("Provide a date column followed by distinct stock columns.")
    data = pd.read_csv(prices)
    dates = pd.to_datetime(data.iloc[:, 0], format="%Y-%m-%d", errors="raise")
    if not dates.is_monotonic_increasing or dates.duplicated().any():
        raise ValueError("Price dates must be unique and chronological.")
    months = pd.period_range(first, last, freq="M")
    if len(months) < 73:
        raise ValueError("The application requires at least 73 months of returns.")
    month_range = pd.period_range(months[0] - 1, months[-1], freq="M")
    # Pick the same last observed date for every stock. Missing closes at the
    # selected date stay missing, even when a prior-day close is present.
    periods = dates.dt.to_period("M")
    ends = dates.groupby(periods).max().reindex(month_range)
    if ends.isna().any() or (ends.dt.weekday >= 5).any():
        raise ValueError("Missing month or weekend session among price month ends.")
    if calendar is not None:
        expected = calendar_ends(calendar, month_range)
        mismatch = ends != expected
        if mismatch.any():
            raise ValueError(f"Calendar/price month-end mismatch in {list(month_range[mismatch])}")
    px = data.iloc[:, 1:].set_axis(dates).loc[ends].to_numpy(dtype=float)
    if not np.isfinite(px).all() or (px <= 0).any():
        raise ValueError("The fixed panel must have positive finite prices at every month end.")
    fac = monthly_ff3(factors, months)
    y = px[1:] / px[:-1] - 1 - fac[:, 3, None]
    return np.ascontiguousarray(y), np.ascontiguousarray(fac[:, :3]), fac[:, 3], months


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prices", required=True, type=Path, help="Dated adjusted-close wide CSV")
    parser.add_argument("--factors", required=True, type=Path, help="Monthly FF3 CSV/ZIP in percent")
    parser.add_argument("--calendar", type=Path, help="Dated trading-session CSV/ZIP")
    parser.add_argument("--first", default="2005-02", help="First return month (YYYY-MM)")
    parser.add_argument("--last", default="2026-09", help="Last return month (YYYY-MM)")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=25)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    y, standard, rf, months = read_updated_panel(args.prices, args.factors, args.first,
                                                 args.last, args.calendar)
    published = standard.copy()
    published[:, 0] -= rf
    records = []
    for spec, design in [("standard_FF3", standard), ("source_market_minus_RF", published)]:
        for row in summarize(y, design, args.workers, include_last_window=True,
                             full_sample=spec == "standard_FF3"):
            row.update(market_spec=spec, first_month=str(months[0]), last_month=str(months[-1]))
            records.append(row)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    print(f"Processed {y.shape[1]} stocks over {len(months)} months ending in {months[-1]}.")


if __name__ == "__main__":
    main()
