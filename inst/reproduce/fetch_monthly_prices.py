"""Build an endpoint-matched stock panel from Yahoo Finance daily adjusted closes.

Dependencies: pandas, numpy and requests. The archive is Wind's
SP500_2005_now_month.csv (ticker, name, then dated monthly prices).
The input constituent list has one column named ticker. Only public ticker
symbols, the date range and chart parameters are sent to Yahoo.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import re
import time
from urllib.parse import quote

import numpy as np
import pandas as pd
import requests


# The 2018 issuers used these tickers before different securities began trading.
FOX_SOURCE = "https://investor.foxcorporation.com/establishment-of-fox/"
IDENTITY_BREAKS = {
    "FOX": ("New Fox Corporation Class B shares began trading in March 2019.", FOX_SOURCE),
    "FOXA": ("New Fox Corporation Class A shares began trading in March 2019.", FOX_SOURCE),
    "IR": ("Gardner Denver changed its symbol from GDI to IR in March 2020.",
           "https://www.sec.gov/Archives/edgar/data/1699150/000114036120004816/form8k.htm"),
}


def screen_identities(universe, prices):
    universe = universe.copy()
    universe["series_status"] = universe["status"]
    universe["identity_note"] = ""
    universe["identity_source"] = ""
    for ticker, (note, source) in IDENTITY_BREAKS.items():
        matched = (universe.ticker == ticker) & universe.in_endpoint_roster
        universe.loc[matched, ["status", "identity_note", "identity_source"]] = [
            "documented_identity_break", note, source]
    included = universe.loc[universe.status == "included", "ticker"]
    return universe, prices.loc[:, included]


def symbol(code, wind=False):
    code = str(code).strip().upper()
    if wind:
        code = re.sub(r"\.(N|O|A)$", "", code)
    code = code.replace(".", "-").replace("_", "-")
    if not re.fullmatch("[A-Z0-9]+(?:-[A-Z0-9]+)*", code):
        raise ValueError(f"Unrecognized ticker {code!r}")
    return code


def get_daily(code, start, stop):
    url = "https://query1.finance.yahoo.com/v8/finance/chart/" + quote(code)
    params = {"period1": int(start.timestamp()), "period2": int(stop.timestamp()),
              "interval": "1d", "events": "history", "includeAdjustedClose": "true"}
    for attempt in range(5):
        try:
            r = requests.get(url, params=params, timeout=(12, 40),
                             headers={"User-Agent": "Mozilla/5.0", "Accept": "application/json"})
            r.raise_for_status()
            chart = r.json()["chart"]
            if chart.get("error"):
                raise ValueError(str(chart["error"]))
            result = chart["result"][0]
            meta = result["meta"]
            if meta["symbol"].upper() != code:
                raise ValueError("Yahoo returned a different symbol: " + meta["symbol"])
            if meta["currency"] != "USD":
                raise ValueError("Yahoo returned a non-USD series")
            dates = pd.to_datetime(result["timestamp"], unit="s", utc=True)
            dates = dates.tz_convert(meta["exchangeTimezoneName"]).tz_localize(None).normalize()
            prices = pd.Series(result["indicators"]["adjclose"][0]["adjclose"],
                               index=dates, dtype=float)
            if prices.index.has_duplicates:
                raise ValueError("duplicate exchange dates")
            return prices.sort_index(), meta
        except requests.RequestException as exc:
            status = exc.response.status_code if exc.response is not None else None
            if attempt == 4 or status in (400, 401, 403, 404, 410, 422):
                raise
            time.sleep(min(20, 2 ** attempt))
    raise RuntimeError("Failed to fetch chart")


def get_archive(path):
    archive = pd.read_csv(path.resolve(), encoding="gb18030")
    codes = archive.iloc[:, 0].map(lambda value: symbol(value, wind=True))
    if codes.duplicated().any():
        raise ValueError("Nonunique normalized archive tickers")
    dated = archive.iloc[:, 2:].apply(pd.to_numeric, errors="raise")
    dates = pd.to_datetime(dated.columns, format="%Y/%m/%d").to_period("M")
    if not dates.is_unique or not dates.is_monotonic_increasing:
        raise ValueError("Duplicate or unordered source months")
    dated.index, dated.columns = codes, dates
    return archive, codes, dated


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--roster", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=5)
    parser.add_argument("--start-month", default="2005-01")
    parser.add_argument("--end-month", default="2026-09")
    parser.add_argument("--min-spearman", type=float, default=.80)
    args = parser.parse_args()
    if not 1 <= args.workers <= 10 or not 0 <= args.min_spearman <= 1:
        parser.error("Use 1 to 10 workers and a rank-correlation threshold in [0, 1].")
    archive, codes, history = get_archive(args.archive)
    roster = pd.read_csv(args.roster)["ticker"].map(symbol)
    if roster.duplicated().any():
        raise ValueError("Nonunique normalized roster tickers")
    months = pd.period_range(args.start_month, args.end_month, freq="M")
    start = months[0].start_time.tz_localize("UTC")
    stop = (months[-1]+1).start_time.tz_localize("UTC") + pd.Timedelta(days=7)
    spy, _ = get_daily("SPY", start, stop)
    spy = spy[spy.notna() & (spy > 0)]
    if not any(spy.index.to_period("M") > months[-1]):
        raise ValueError("Trading calendar ends before the final month can be verified.")
    calendar = spy.index.to_series().groupby(spy.index.to_period("M")).last().reindex(months)
    if calendar.isna().any():
        raise ValueError("Incomplete SPY trading calendar")
    dates = pd.DatetimeIndex(calendar, name="date")
    hist_months = history.columns.intersection(months)
    expected = len(hist_months) - 1
    if expected < 24:
        raise ValueError("Insufficient overlapping history")
    roster_set = set(roster)
    args.out.mkdir(parents=True, exist_ok=True)
    print(f"{len(codes)} original and {len(roster_set)} endpoint symbols, "
          f"{sum(codes.isin(roster_set))} matched symbols and {len(dates)} monthly endpoints",
          flush=True)

    def one(i):
        code = codes.iloc[i]
        info = {"archive_code": archive.iloc[i,0], "archive_name": archive.iloc[i,1],
                "ticker": code, "in_endpoint_roster": code in roster_set,
                "status": "not_in_endpoint_roster", "source": "Yahoo Finance chart v8",
                "yahoo_long_name": None, "yahoo_exchange": None,
                "matched_months": None, "missing_months": None,
                "overlap_returns": None, "pearson": None, "spearman": None,
                "first_missing_month": None, "error": None}
        if code not in roster_set:
            return info, None
        try:
            daily, meta = get_daily(code, start, stop)
        except requests.RequestException as exc:
            info.update(status="http_failure" if exc.response is not None else "network_failure",
                        error=(f"HTTP {exc.response.status_code}" if exc.response is not None
                               else type(exc).__name__))
            return info, None
        except (KeyError, IndexError, ValueError, TypeError) as exc:
            info.update(status="invalid_chart", error=str(exc)[:200])
            return info, None
        monthly = daily.reindex(dates)
        monthly = monthly.where(np.isfinite(monthly) & (monthly > 0))
        monthly.index = months
        miss = monthly.isna()
        calc = pd.concat([history.loc[code, hist_months].pct_change(fill_method=None),
                          monthly.reindex(hist_months).pct_change(fill_method=None)], axis=1)
        calc = calc.replace([np.inf,-np.inf], np.nan).dropna()
        pearson = calc.iloc[:,0].corr(calc.iloc[:,1]) if len(calc) >= 2 else np.nan
        ranks = calc.rank()
        spearman = ranks.iloc[:,0].corr(ranks.iloc[:,1]) if len(calc) >= 2 else np.nan
        info.update(yahoo_long_name=meta.get("longName") or meta.get("shortName"),
                    yahoo_exchange=meta.get("fullExchangeName"), matched_months=int((~miss).sum()),
                    missing_months=int(miss.sum()), overlap_returns=len(calc),
                    pearson=pearson, spearman=spearman,
                    first_missing_month=str(months[miss][0]) if miss.any() else None)
        if miss.any():
            info["status"] = "incomplete_month_ends"
        elif len(calc) != expected or not np.isfinite(spearman):
            info["status"] = "incomplete_overlap"
        elif spearman < args.min_spearman:
            info["status"] = "low_return_correlation"
        else:
            info["status"] = "included"
        return info, monthly.to_numpy() if info["status"] == "included" else None

    records, columns = [], {}
    with ThreadPoolExecutor(args.workers) as pool:
        for i, (item, values) in enumerate(pool.map(one, range(len(codes))), 1):
            records.append(item)
            if values is not None:
                columns[item["ticker"]] = values
            if i % 50 == 0:
                print(f"Reviewed {i}/{len(codes)} archive securities", flush=True)
    universe = pd.DataFrame(records)
    prices = pd.DataFrame(columns, index=dates)
    universe, prices = screen_identities(universe, prices)
    universe.to_csv(args.out / "universe.csv", index=False)
    if prices.empty:
        raise RuntimeError("No qualified securities. Inspect the security audit.")
    prices.to_csv(args.out / "monthly_adjusted_prices.csv", float_format="%.12g")
    print(universe["status"].value_counts().to_string(), flush=True)
    print(f"Price matrix: {prices.shape}, endpoints {dates[0].date()} to {dates[-1].date()}", flush=True)


if __name__ == "__main__":
    main()
