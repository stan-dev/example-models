# Foreign Exchange Rate Data

Data and example from Fernando Tussell, "Kalman Filtering in R", *Journal of
Statistical Software*, http://www.jstatsoft.org/v39/i02/paper.

The data are three time series of exchange rates of BEF (Belgian
Franc), CHF (Swiss Franc), and DEM (German Deutchmark)

Two models are estimated.

- `fx_equicorr.stan`: local level models for each series with a common
  measurement correlation.
- `fx_factor.stan`: Single factor model with local level factor.

Both models use the data in `fx.data.R` as input.
