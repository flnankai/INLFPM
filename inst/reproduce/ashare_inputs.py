"""Read the two licensed/local CSV inputs used in the CSI 300 application.

The archive's monthly_returns file contains simple decimal stock returns
before subtraction of the risk-free rate. Its Tencent monthly price returns
are supplemented by returns from an archived price panel where unavailable.
The external monthly factor file supplies decimal factor returns and RF.
No local source data are bundled with INLFPM.
"""
from pathlib import Path

import numpy as np
import pandas as pd


RETURNS_FILE = "hs300_ff4_external_rf_hybrid_2007_2025_monthly_returns.csv"
FACTORS_FILE = (
    Path("factor_data_20260607") / "three_four_five_factor_monthly"
    / "fivefactor_monthly.csv"
)


def load_ashare_returns(source_dir):
    """Read all 299 raw monthly series, retaining their original missing values."""
    source_dir = Path(source_dir)
    months = pd.period_range("2007-02", "2025-12", freq="M", name="month")
    raw = pd.read_csv(source_dir / RETURNS_FILE, index_col="Date",
                      float_precision="round_trip")
    raw.index = pd.to_datetime(raw.index, format="%Y-%m-%d").to_period("M")
    if not raw.index.equals(months):
        raise ValueError("Expected one stock-return row per month, Feb 2007 to Dec 2025.")
    if not all(len(code) == 6 and code.isdigit() for code in raw.columns):
        raise ValueError("Expected six-digit A-share stock codes.")
    returns = raw.apply(pd.to_numeric, errors="raise")
    if np.isinf(returns.to_numpy()).any() or (returns < -1).any().any():
        raise ValueError("Invalid monthly stock returns in the source panel.")
    if returns.shape[1] != 299:
        raise ValueError("Expected the 299 stock series in the archived return panel.")
    returns.index = months
    return returns


def load_ashare_inputs(source_dir):
    """Return raw monthly stock returns, monthly RF and five matched factor designs."""
    returns = load_ashare_returns(source_dir)
    months = returns.index
    # Freeze the stock universe before fitting any models or rolling windows.
    complete = returns.dropna(axis="columns")
    if complete.shape[1] != 103:
        raise ValueError("The article used 103 complete histories among 299 stock series.")
    factors = pd.read_csv(Path(source_dir) / FACTORS_FILE, index_col="trdmn",
                          float_precision="round_trip")
    factors.index = pd.to_datetime(factors.index.astype(str), format="%Y%m").to_period("M")
    if not factors.index.is_unique:
        raise ValueError("Duplicate monthly observations in the factor archive.")
    factors = factors.rename(columns={"mkt_rf": "MKT", "smb": "SMB", "hml": "HML",
                                      "umd": "UMD", "rmw": "RMW", "cma": "CMA", "rf": "RF"})
    factors = factors.loc[months, ["MKT", "SMB", "HML", "RMW", "CMA", "UMD", "RF"]]
    if not np.isfinite(factors.to_numpy(dtype=float)).all():
        raise ValueError("Nonfinite factor return or RF on an analysis month.")
    # The market factor is already an excess return. Callers subtract RF only from stocks.
    complete.index = months
    factors.index = months
    specifications = {
        "CAPM": ["MKT"],
        "FF3": ["MKT", "SMB", "HML"],
        "Carhart4": ["MKT", "SMB", "HML", "UMD"],
        "FF5": ["MKT", "SMB", "HML", "RMW", "CMA"],
        "FF5+UMD": ["MKT", "SMB", "HML", "RMW", "CMA", "UMD"],
    }
    designs = {name: factors[columns] for name, columns in specifications.items()}
    return complete, factors["RF"], designs
