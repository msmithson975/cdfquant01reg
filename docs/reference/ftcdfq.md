# Fit an FTCDFQ Model using brms

Fit an FTCDFQ Model using brms

## Usage

``` r
ftcdfq(formula, data, backend = "cmdstanr", silent = 2, ...)
```

## Arguments

- formula:

  A brms formula object.

- data:

  A data frame containing the variables in the model.

- backend:

  The backend engine to use for compilation (defaults to "cmdstanr").

- silent:

  A number controlling the warnings output from stan.

- ...:

  Additional arguments passed to brms::brm.
