# Maximum Likelihood Estimation for Extended CDF Quantile Regression

Fits a parametric frequentist regression model for bounded outcomes on
the unit interval using maximum likelihood estimation across three
distribution parameters: mu, sigma, and u.

## Usage

``` r
exftcdfqMLE(formula, data, start = NULL)

# S3 method for class 'exftcdfq'
print(x, ...)

# S3 method for class 'exftcdfq'
summary(object, ...)

# S3 method for class 'summary.exftcdfq'
print(x, ...)

# S3 method for class 'exftcdfq'
coef(object, ...)

# S3 method for class 'exftcdfq'
logLik(object, ...)
```

## Arguments

- formula:

  A three-part formula separated by pipes (`|`), where the first part
  specifies predictors for mu (skew and location), the second for sigma
  (dispersion), and the third for u (exceedance).

- data:

  A data frame containing the variables specified in the formula.

- start:

  An optional numeric vector of starting values for the parameters. If
  `NULL`, defaults to `0.1` for all parameters.

- x:

  An object of class `"exftcdfq"` or `"summary.exftcdfq"` (for print
  methods).

- ...:

  Additional arguments (not used).

- object:

  An object of class `"exftcdfq"`.

## Value

An object of class (`"exftcdfq"`.) This is a list containing:

An object of class `"summary.exftcdfq"`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Select a distribution first
select_exftcdfq()

# Fit an MLE model
m1 <- exftcdfqMLE(ptotpalest1 ~ palestpk | 1 | 1, data = blame)
summary(m1)
} # }
```
