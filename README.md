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

To regenerate the manuscript counts, run the Python script below from the repository root. It generates the original NumPy panels, computes the five unchanged p-values using the original bootstrap seed and corrected leave-two-out IN-SUM kernel, then evaluates both calibrated maxima and their combinations on each panel. It runs 1,000 replications per cell with 25 processes by default. The package contains the required numerical routines. The computation writes a single summary CSV.

```sh
python inst/reproduce/manuscript.py --run --workers 25 --output calibrated_summary.csv
```

Python 3.14.4, NumPy 2.4.4 and SciPy 1.17.1 were used for the calibrated results. If the original `INMAX_Research_Package_20261108.zip` is available locally, add `--archive path/to/INMAX_Research_Package_20261108.zip` to reuse the five unchanged p-values stored for each panel. This mode also uses pandas 3.0.2, regenerates the original data, verifies panel seeds and fitted maxima, and recomputes the four changed procedures. Both modes use a strict five-percent rejection threshold.

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

## License

GPL-3.
