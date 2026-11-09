"""Read weekly CSI constituent returns and match FF factors by trading week."""
import re

import numpy as np
import pandas as pd


def load_weekly_inputs(stock_path, factor_path):
    """Return raw decimal stock returns, weekly RF and two matched factor designs."""
    stocks = pd.read_excel(stock_path, sheet_name="Sheet1", index_col="Date").T
    stocks.index = pd.to_datetime(stocks.index)
    dates = stocks.index
    weeks = dates.to_period("W-FRI")
    if (len(dates) != 270 or dates[0] != pd.Timestamp("2020-01-03") or
            dates[-1] != pd.Timestamp("2025-04-11") or not dates.is_monotonic_increasing
            or not weeks.is_unique or (dates.dayofweek > 4).any()):
        raise ValueError("Expected one date per trading week in the archived stock workbook.")
    if (not stocks.columns.is_unique or len(stocks.columns) not in (300, 500) or
            any(not re.fullmatch(r"[0-9]{6}\.(SH|SZ)", str(code)) for code in stocks)):
        raise ValueError("Expected one unique exchange-qualified code per stock.")
    pct = stocks.apply(pd.to_numeric, errors="raise")
    if np.isinf(pct.to_numpy()).any() or (pct < -100).any().any():
        raise ValueError("Invalid weekly percentage change in the stock workbook.")
    # Pre-listing periods can be coded zero. Require observed trading before the sample.
    run_in = pct.iloc[:4]
    traded = (run_in.notna() & run_in.ne(0)).any()
    sample = pct.iloc[4:].loc[:, traded].dropna(axis="columns")
    flat_window = sample.rolling(52).max().eq(sample.rolling(52).min()).any()
    returns = sample.loc[:, ~flat_window] / 100

    factors = pd.read_csv(factor_path, index_col="trdwk", float_precision="round_trip")
    factors.index = pd.to_datetime(factors.index, format="%Y-%m-%d").to_period("W-FRI")
    if not factors.index.is_unique or not factors.index.is_monotonic_increasing:
        raise ValueError("Expected one dated factor row per trading week, in time order.")
    covered = factors.loc[weeks[0]:weeks[-1]]
    if not covered.index.equals(weeks):
        raise ValueError("Stock and factor trading-week calendars differ.")
    # Holiday weeks have early stock closes. Factor labels denote the same trading week.
    columns = ["mkt_rf", "smb", "hml", "rmw", "cma", "rf"]
    matched = covered.iloc[4:][columns].copy()
    if not np.isfinite(matched.to_numpy(dtype=float)).all():
        raise ValueError("Incomplete factor or RF series for a stock-return week.")
    matched.index = returns.index
    return returns, matched["rf"], {"FF3": matched[columns[:3]], "FF5": matched[columns[:5]]}
