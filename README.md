# INLFPM

INLFPM tests whether every intercept in a high-dimensional factor model is zero. It implements the LS, SS and IN families, each with a sum, maximum and Cauchy combination test.

## Installation

Install from the GitHub repository with `remotes::install_github("flnankai/INLFPM")`, or run `R CMD INSTALL INLFPM` on a local source directory. The package uses standard R libraries.

## Test a panel

Pass a matrix of returns with one row per date and one column per asset. The factor matrix uses the same dates and excludes the intercept.

```r
library(INLFPM)
fit <- inmax_test(Y, F, seed = 123)
fit$p.value
fit$failures

ss <- ac_ss_analytic(Y, F)
inv <- ac_in_analytic(Y, F)
c(calibrated = inv$p, original = inv$original_p)
```

SS-MAX evaluates a Gaussian quadratic-ratio tail at the angular-adjusted formal statistic. IN-MAX uses a fitted time-variance profile, a PQR moment correction and a matched Student tail. Both default to the Sidak maximum mapping. SS-CC uses the truncated Cauchy statistic and IN-CC uses the ordinary Cauchy statistic, each with the calibrated MAX and SUM p-values from the same panel. The LS maximum uses its original Gumbel reference.

The fitted finite-sample reference is approximate under cross-sectional dependence. The `diagnostics` element provides fitted nuisance parameters and the marginal tail calculation. Failed component fits are identified in `failures` and return missing p-values for the affected methods.

## Manuscript results

Read the reported results directly from the installed package. Each row identifies the design, method, number of attempted and valid fits, number of rejections and rejection rate.

```r
paper <- read.csv(system.file("extdata", "paper-summary.csv", package = "INLFPM"))
subset(paper, design == "null" & T == 60 & N == 200)
```

To regenerate the manuscript counts, run the Python script below from the repository root. It uses the manuscript's NumPy random streams and bootstrap seeds to compute all nine tests, including the leave-two-out IN-SUM statistic and both calibrated maxima. It runs 1,000 replications per cell with 25 processes by default and writes a single summary CSV.

```sh
python inst/reproduce/manuscript.py --run --workers 25 --output calibrated_summary.csv
```

Python 3.14.4, NumPy 2.4.4 and SciPy 1.17.1 were used for the calibrated results. The independent-component design divides each Student t draw with three degrees of freedom by `sqrt(3)` before multiplying by the symmetric correlation-matrix square root. All four designs therefore have the same error covariance.

If the original `INMAX_Research_Package_20261108.zip` is available locally, add `--archive path/to/INMAX_Research_Package_20261108.zip` to reuse archived component p-values for the three elliptical designs. This mode also uses pandas 3.0.2, verifies the archive's panel seeds and fitted maxima for those designs, and recomputes all nine procedures for the standardized independent-component design. Both modes use a strict five-percent rejection threshold.

## Independent implementation in R

`inmax_run()` implements all 24 null cells and the manuscript's sparsity and signal-strength designs with R's random stream. This produces an independent Monte Carlo realization with its own rejection counts. Each base panel has a seed determined by its design and replication index, independent of worker scheduling. The driver retains only cell counts in memory and can write a final summary CSV.

```r
# Run explicitly when the full experiment is needed.
results <- inmax_run(replications = 1000L, workers = 25L,
                     output = "calibrated_summary.csv")
```

The final summary gives rejection rates per attempted panel, rejected and valid counts, failures, Monte Carlo standard errors and bounds for the effect of failed fits. An installed command-line entry point is available through `system.file("reproduce", "run.R", package = "INLFPM")`.

## Monthly stock application

The application uses 375 stocks with complete prices from January 2005 through October 2018. Supply local copies of the monthly stock-price and factor CSVs to the Python entry point. Prices are GBK-encoded, with two metadata columns followed by consecutive dates. The factor file has columns `month` (YYYYMM), `Mkt-RF`, `SMB`, `HML` and `RF`, in percentage points.

```sh
python inst/reproduce/empirical.py --prices SP500_2005_now_month.csv --factors F-F_Research_Data_Factors.CSV --workers 25 --output empirical_summary.csv
```

Simple stock returns and factor rows are paired by their ending month. The source application's market regressor is constructed as `(Mkt-RF column - RF column)/100`; the response subtracts `RF/100` from each simple return. SMB and HML enter in their source units. For 165 return months and window length T, the application uses the first 165 - T start positions, yielding 117, 105 and 93 windows of length 48, 60 and 72. Each window's 100 bootstrap draws use a NumPy SeedSequence with entropy `[20261109, T, zero_based_start]`. One additional row per method and level describes the full 165-month sample. The `full_pvalue` column records those full-sample p-values.

```r
application <- read.csv(system.file("extdata", "empirical-summary.csv", package = "INLFPM"))
subset(application, sample == "rolling" & level == 0.05)
```

## Updated monthly stock application

The extended application starts with 375 identifiers in the earlier Wind panel. Of these, 275 ticker symbols match the September 30, 2026 S&P 500 roster. Three matches (FOX, FOXA and IR) have documented changes of security identity. The remaining 272 have positive month-end adjusted closes throughout and a Spearman correlation of at least 0.80 between the Wind and Yahoo monthly returns during the overlap. This fixed panel uses information from both endpoints. Membership can change between those endpoints.

The 503 distinct symbols in `inst/extdata/constituents_20260930.csv` come from the September 30, 2026 observation in Chinobing's historical S&P 500 constituent archive. [Frozen constituent history](https://github.com/chinobing/historical_sp500_constituents/blob/edf70ac0cb1237ef2cc585dcbb79f294a8fb282a/sp_500_historical_components.csv)

Given the local Wind file and a single-column endpoint roster headed `ticker`, the price-preparation command requests Yahoo daily adjusted closes, obtains month ends from SPY sessions and writes a security audit and the month-end matrix. It uses Python with NumPy, pandas and requests. The data used for the published update have adjusted-price CSV SHA-256 `01c321e50f2b0e9dfd59cbe0fbade255f54e433abb18617eff612d397d0224e8`. Historical adjusted prices can change after vendor revisions.

```sh
python inst/reproduce/fetch_monthly_prices.py --archive SP500_2005_now_month.csv --roster inst/extdata/constituents_20260930.csv --out panel
```

Supply an adjusted-close CSV with one row per trading date or month end. The first column contains ISO dates and the remaining columns identify individual stocks. For each calendar month the script selects one common closing date for every stock. Supply the monthly Fama–French U.S. research factors as the official CSV or ZIP in percentage units. An optional calendar file accepts the French daily CSV/ZIP or a CSV with trading dates in its first column.

```sh
python inst/reproduce/empirical_updated.py --prices panel/monthly_adjusted_prices.csv --factors F-F_Research_Data_Factors_CSV.zip --calendar F-F_Research_Data_Factors_daily_CSV.zip --workers 25 --output empirical_updated_summary.csv
```

January 2005 supplies the starting adjusted price. The 260 simple monthly returns cover February 2005 through September 2026. For each stock, the script subtracts RF from the return. Under `standard_FF3` the regressors use Mkt-RF, SMB and HML as published, all in decimal units. Under `source_market_minus_RF` the market column is Mkt-RF minus RF, reproducing the earlier application's convention. The 48-, 60- and 72-month analyses include all 213, 201 and 189 available windows under both conventions with common bootstrap draws. Full-sample p-values use standard FF3.

```r
updated <- read.csv(system.file("extdata", "empirical-updated-summary.csv", package = "INLFPM"))
subset(updated, sample == "rolling" & level == 0.05 & market_spec == "standard_FF3")
```

## Five pricing models on matched windows

The model comparison keeps the same 272 stocks, 260 return months, 603 rolling windows, risk-free series and bootstrap draws. CAPM uses Mkt-RF. FF3 adds SMB and HML from the three-factor file. Carhart4 adds the monthly momentum factor (Mom, also called UMD) to FF3. FF5 uses Mkt-RF, SMB, HML, RMW and CMA from the five-factor file, including that file's five-factor SMB definition. FF5+UMD adds Mom to FF5. All factors use the September 2026 U.S. data release. [Kenneth French Data Library](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html)

```sh
python inst/reproduce/empirical_models.py --prices panel/monthly_adjusted_prices.csv --ff3 F-F_Research_Data_Factors_CSV.zip --ff5 F-F_Research_Data_5_Factors_2x3_CSV.zip --momentum F-F_Momentum_Factor_CSV.zip --calendar F-F_Research_Data_Factors_daily_CSV.zip --workers 25 --summary factor-model-summary.csv --overlap factor-model-IN-overlap.csv
```

The model summary provides valid fits, attempted windows and rejection counts at each nominal level, plus full-sample p-values. The overlap table uses the common valid windows of IN-SUM, IN-MAX and IN-CC. Its `max_only` and `sum_only` columns count component disagreements. The `cc_keeps_max_only` and `cc_keeps_sum_only` columns count the disagreements for which IN-CC rejects. Each available `first_max_only` example reports the earliest start, end and three p-values. Nested rolling windows overlap in calendar time.

```r
models <- read.csv(system.file("extdata", "factor-model-summary.csv", package = "INLFPM"))
overlap <- read.csv(system.file("extdata", "factor-model-IN-overlap.csv", package = "INLFPM"))
subset(overlap, months == 60 & level == 0.05)
```

## Monthly CSI 300 A-share application

The fixed panel contains 103 stocks with uninterrupted histories from February 2007 through December 2025. The June 2026 Eastmoney snapshot lists 300 CSI 300 constituents. The archived return file covers 299 of them. Selecting complete histories from this later roster conditions the sample on both membership and data availability. Historical index membership changed over the sample.

Of the 23,381 monthly returns in the fixed panel, 23,053 come from ratios of consecutive Tencent monthly closing prices with its `hfq` adjustment (后复权). The archive obtains the other 328 from consecutive prices in a local adjusted-price file and fills the missing returns by stock and month. Both prices in each ratio come from the same source. The two sources use their own adjustment conventions.

The seven models use the same raw returns, 180 four-year windows, 168 five-year windows and 156 six-year windows. CAPM, FF3, Carhart4, FF5 and FF5+UMD use a common archive of Chinese factors constructed with CSMAR data. The market factor is the excess return on the float-cap-weighted A-share market. SMB, HML, UMD, RMW and CMA are the archive's value-weighted factor series. Each stock's return is reduced by the supplied monthly RF once. For CH3 and CH4, each model uses its own workbook's market, factor and RF series. CH3 uses size and earnings-to-price value. CH4 also uses abnormal turnover. [Liu, Stambaugh and Yuan (2019), *Size and Value in China*](https://en.mingshiim.com/static/sizevaluechina.pdf)

To reproduce all seven specifications, install Python with NumPy, SciPy, pandas and openpyxl. Supply the two local archive files under the following paths relative to `--source-dir`. The stock-return panel and CSMAR-based archive must be supplied by the user.

| Input | Data layout | SHA-256 |
| --- | --- | --- |
| `hs300_ff4_external_rf_hybrid_2007_2025_monthly_returns.csv` | `Date` (ISO month end), followed by 299 six-digit stock codes. Values are simple decimal monthly returns before RF. | `7a1a4ecd932800c16a18428a5503ef0b34d0fddb5228665d5d6f7ef2c1206864` |
| `factor_data_20260607/three_four_five_factor_monthly/fivefactor_monthly.csv` | `trdmn` (YYYYMM), `mkt_rf`, `smb`, `hml`, `umd`, `rmw`, `cma`, `rf` in decimal monthly units. | `60ac677b3692dd29028ccee55ba9468c08b3be651d66dd8fd7de03568a6daeb1` |
| `CH3_factors_monthly_202608.xlsx` | `Returnseries` sheet with `mnthdt` (YYYYMMDD), `rf_mon`, `mktrf`, `SMB`, `VMG`. | `8d956043a1a39d10bc1b777a865fca8cf40b1b628c1d94f900ba2b23c22819a4` |
| `CH4_factors_monthly_202608.xlsx` | The CH4 `Returnseries` sheet also contains `PMO`. | `9ebb2d93fb8730f100251a2fe511ea4749c332d33183662eee75558f0b632555` |

The two monthly workbooks are available from the authors' data library. The article uses the August 2026 release, restricted to the stock sample dates. [Mingshi monthly CH3 and CH4 downloads](https://en.mingshiim.com/wq/)

```sh
python inst/reproduce/empirical_ashares.py --source-dir path/to/local/archive --ch3 CH3_factors_monthly_202608.xlsx --ch4 CH4_factors_monthly_202608.xlsx --workers 25 --summary inst/extdata/ashare-summary.csv --overlap inst/extdata/ashare-IN-overlap.csv
```

Omit both workbook flags to run the five specifications supplied by the local archive. The summary includes full-sample p-values and all nine tests' rolling rejection counts. The overlap table records counts for the paired IN-SUM, IN-MAX and IN-CC decisions. All windows include both endpoints and move one month at a time.

```r
cn <- read.csv(system.file("extdata", "ashare-summary.csv", package = "INLFPM"))
cn_pairs <- read.csv(system.file("extdata", "ashare-IN-overlap.csv", package = "INLFPM"))
subset(cn_pairs, months == 60 & level == 0.05)
```

For the CH3 and CH4 coverage check, start from all 299 series and retain each stock with returns in every month of the current window. The resulting panels contain 104–280 stocks over four-year windows, 104–270 over five-year windows and 104–256 over six-year windows. Within a given window, both models and all nine tests use the same stocks. The full-sample check uses the original 103 complete histories.

```sh
python inst/reproduce/empirical_ashare_windows.py --source-dir path/to/local/archive --ch3 CH3_factors_monthly_202608.xlsx --ch4 CH4_factors_monthly_202608.xlsx --workers 25 --summary inst/extdata/ashare-window-summary.csv --overlap inst/extdata/ashare-window-IN-overlap.csv
```

## Weekly CSI 300 and CSI 500 comparison

Two archived constituent workbooks each report 270 trading weeks of percentage changes, from January 3, 2020 through April 11, 2025. The stock lists are disjoint and fixed at the rosters in the respective exports. Each workbook's `Sheet1` has a first column headed `Date` containing exchange-qualified stock codes, with week-ending trading dates across the remaining columns. Divide each percentage change by 100. The third file supplies decimal weekly returns for `mkt_rf`, `smb`, `hml`, `rmw`, `cma` and `rf`, indexed by `trdwk`. The source places each week's factors on Friday. Stock dates can end earlier in holiday weeks. The reader requires a one-to-one match of all trading weeks in chronological order.

The constituent exports contain zeros before some stocks first trade. The first four observations establish eligibility. A stock must have at least one recorded nonzero return during this January 2020 run-in, finite returns on every subsequent date and some variation in every 52-week analysis window. One CSI 300 security, 000792.SZ, had a 66-week stretch of zero returns during a listing suspension. Shorter zero runs remain in the analysis as reported. The resulting fixed panels contain 265 CSI 300 and 417 CSI 500 stocks. Each has 266 weekly observations from February 7, 2020, giving 215 overlapping 52-trading-week windows. [Shenzhen Stock Exchange filing on 000792.SZ](https://disc.static.szse.cn/disc/disk02/finalpage/2021-08-03/226f764b-539a-4c13-80b8-a2453ccf49ff.PDF)

Each index is tested with the three-factor and five-factor specifications. The response subtracts weekly RF from each stock once. Supply the archived stock workbooks and CSMAR-based weekly factor file with the following SHA-256 fingerprints.

| Filename | SHA-256 |
| --- | --- |
| `沪深300成分股近五年周数据PctChange.xlsx` | `0a5192f1b8bb1161a843402e2763c4f8807270c1319e5651c22948a4d00af9c4` |
| `中证500成分股近五年周数据PctChange.xlsx` | `15807d1c981883d97e25a2de26860d868397e2cd1af390c6d5b3124d788b2f4f` |
| `fivefactor_weekly.csv` | `cc6262cc4ed31d5aa7b09111703037777efdcd2590b7fb8cb4c5c436d3f470c7` |

```sh
python inst/reproduce/empirical_ashare_weekly.py --csi300 "沪深300成分股近五年周数据PctChange.xlsx" --csi500 "中证500成分股近五年周数据PctChange.xlsx" --factors fivefactor_weekly.csv --workers 25 --summary inst/extdata/ashare-weekly-summary.csv --overlap inst/extdata/ashare-weekly-IN-overlap.csv
```

## License

GPL-3.
