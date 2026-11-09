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

## License

GPL-3.
