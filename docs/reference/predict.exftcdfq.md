# Predicting quantiles Obtains predicted values (quantiles) from a fitted extended CDF quantile regression model object.

Predict Method for exftcdfq Model Fits

## Usage

``` r
# S3 method for class 'exftcdfq'
predict(object, newdata = NULL, q = 0.5, ...)
```

## Arguments

- object:

  An object of class \\`"exftcdfq"`.\\

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

## Details

Obtains predicted values (quantiles) from a fitted extended CDF quantile
regression model object.
