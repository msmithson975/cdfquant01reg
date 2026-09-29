# Maximum Likelihood Estimation for 3-Parameter CDF Quantile Regression

Fits a parametric frequentist regression model for bounded outcomes on
the unit interval using maximum likelihood estimation across three
distribution parameters: mu, sigma, and theta.

## Usage

``` r
ftcdfqMLE3(formula, data, start = NULL)

# S3 method for class 'ftcdfq3'
print(x, ...)

# S3 method for class 'ftcdfq3'
summary(object, ...)

# S3 method for class 'summary.ftcdfq3'
print(x, ...)

# S3 method for class 'ftcdfq3'
coef(object, ...)

# S3 method for class 'ftcdfq3'
logLik(object, ...)
```

## Arguments

- formula:

  A three-part formula separated by pipes (`|`), where the first part
  specifies predictors for mu (location), the second for sigma
  (dispersion), and the third for theta (skew).

- data:

  A data frame containing the variables specified in the formula.

- start:

  An optional numeric vector of starting values for the parameters. If
  `NULL`, defaults to `0.1` for all parameters.

- x:

  An object of class `"ftcdfq3"` or `"summary.ftcdfq3"` (for print
  methods).

- ...:

  Additional arguments (not used).

- object:

  An object of class `"ftcdfq3"`.

## Value

An object of class (`"ftcdfq3"`.) This is a list containing:

An object of class `"summary.ftcdfq3"`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Select a distribution first
select_ftcdfq3()

# Fit an MLE model
m1 <- ftcdfqMLE3(ptotpalest1 ~ palestpk | 1 | 1 | 1, data = blame)
summary(m1)
} # }
# Fit a model
```
