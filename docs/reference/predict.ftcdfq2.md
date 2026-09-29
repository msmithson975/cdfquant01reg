# Predicting quantiles Obtains predicted values (quantiles) from a fitted CDF quantile regression model object.

Predicting quantiles Obtains predicted values (quantiles) from a fitted
CDF quantile regression model object.

## Usage

``` r
# S3 method for class 'ftcdfq2'
predict(object, newdata = NULL, q = 0.5, ...)
```

## Arguments

- object:

  An object of class \\`"ftcdfq2"`.\\

- newdata:

  An optional data frame or vector of new values for out-of-sample
  prediction. If omitted, fitted values for the original dataset are
  returned.

- q:

  A numeric value or vector indicating the target quantile(s). Defaults
  to `0.5` (median).

- ...:

  Additional arguments (not used).

## Value

A numeric vector of predicted values.
