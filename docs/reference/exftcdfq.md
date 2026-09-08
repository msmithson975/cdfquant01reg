# Fit an EXFTCDFQ Model using brms

Fit an EXFTCDFQ Model using brms

## Usage

``` r
exftcdfq(formula, data, backend = "cmdstanr", ...)
```

## Arguments

- formula:

  A brms formula object.

- data:

  A data frame containing the variables in the model.

- backend:

  The backend engine to use for compilation (defaults to "cmdstanr").

- ...:

  Additional arguments passed to brms::brm.
