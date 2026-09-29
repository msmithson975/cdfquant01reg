# Maximum Likelihood Estimation for 2-Parameter CDF Quantile Regression

Fits a parametric frequentist regression model for bounded outcomes on
the unit interval using maximum likelihood estimation across two
distribution parameters: mu, sigma.

## Usage

``` r
ftcdfqMLE2(formula, data, start = NULL)

# S3 method for class 'ftcdfq2'
print(x, ...)

# S3 method for class 'ftcdfq2'
summary(object, ...)

# S3 method for class 'summary.ftcdfq2'
print(x, ...)

# S3 method for class 'ftcdfq2'
coef(object, ...)

# S3 method for class 'ftcdfq2'
logLik(object, ...)
```

## Arguments

- formula:

  A two-part formula separated by pipes (`|`), where the first part
  specifies predictors for mu (location) and the second for sigma
  (dispersion).

- data:

  A data frame containing the variables specified in the formula.

- start:

  An optional numeric vector of starting values for the parameters. If
  `NULL`, defaults to `0.1` for all parameters.

- x:

  An object of class `"ftcdfq2"` or `"summary.ftcdfq2"` (for print
  methods).

- ...:

  Additional arguments (not used).

- object:

  An object of class `"ftcdfq2"`.

## Value

An object of class (`"ftcdfq2"`.) This is a list containing:

An object of class `"summary.ftcdfq2"`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Select a distribution first
select_ftcdfq2()

# Fit an MLE model
m1 <- ftcdfqMLE2(ptotpalest1 ~ palestpk | 1 | 1 , data = blame)
summary(m1)
} # }
# Fit a model
```
