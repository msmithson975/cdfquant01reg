#' Fit an FTCDFQ Model using brms
#'
#' @param formula A brms formula object.
#' @param data A data frame containing the variables in the model.
#' @param backend The backend engine to use for compilation (defaults to "cmdstanr").
#' @param ... Additional arguments passed to brms::brm.
#' @export
ftcdfq <- function(formula, data, backend = "cmdstanr", ...) {
  # Look for the loaded family/stanvars objects in the user's workspace
  if (!exists("my_family", envir = .GlobalEnv) || !exists("stanvars", envir = .GlobalEnv)) {
    stop("No model has been loaded yet. Please run select_ftcdfq() first.")
  }

  if (backend == "cmdstanr" && !requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("The 'cmdstanr' package is required for this backend. Please install it.")
  }

  # Fetch objects safely from the global space
  fam <- get("my_family", envir = .GlobalEnv)
  sv  <- get("stanvars", envir = .GlobalEnv)

  brms::brm(
    formula = formula,
    family = fam,
    stanvars = sv,
    data = data,
    backend = backend,
    ...
  )
}

#' Fit an EXFTCDFQ Model using brms
#'
#' @param formula A brms formula object.
#' @param data A data frame containing the variables in the model.
#' @param backend The backend engine to use for compilation (defaults to "cmdstanr").
#' @param ... Additional arguments passed to brms::brm.
#' @export
exftcdfq <- function(formula, data, backend = "cmdstanr", ...) {
  if (!exists("my_family", envir = .GlobalEnv) || !exists("stanvars", envir = .GlobalEnv)) {
    stop("No model has been loaded yet. Please run select_exftcdfq() first.")
  }

  if (backend == "cmdstanr" && !requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("The 'cmdstanr' package is required for this backend. Please install it.")
  }

  fam <- get("my_family", envir = .GlobalEnv)
  sv  <- get("stanvars", envir = .GlobalEnv)

  brms::brm(
    formula = formula,
    family = fam,
    stanvars = sv,
    data = data,
    backend = backend,
    ...
  )
}
#
