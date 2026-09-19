#' Fit an FTCDFQ Model using brms
#'
#' @param formula A brms formula object.
#' @param data A data frame containing the variables in the model.
#' @param backend The backend engine to use for compilation (defaults to "cmdstanr").
#' @param silent A number controlling the warnings output from stan.
#' @param ... Additional arguments passed to brms::brm.
#' @export
ftcdfq <- function(formula, data, backend = "cmdstanr",  silent = 2, ...) {
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
    silent = silent,
    ...
  )
}

#' Fit an EXFTCDFQ Model using brms
#'
#' @param formula A brms formula object.
#' @param data A data frame containing the variables in the model.
#' @param backend The backend engine to use for compilation (defaults to "cmdstanr").
#' @param silent A number controlling the warnings output from stan.
#' @param ... Additional arguments passed to brms::brm.
#' @export
exftcdfq <- function(formula, data, backend = "cmdstanr", silent = 2, ...) {
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
    silent = silent,
    ...
  )
}

 # Fit an FTCDFQ Model using optim

 # This log likelihood function uses the dnsity2 function in the selected distribution.
 # It assumes that a design matrix has been constructed, containing y, x, z, and w.
 # Note that mu is really theta.
cdfloglik <- function(h, y, x, z, w, density_fun) {
  k_x <- ncol(x)
  k_z <- ncol(z)
  k_w <- ncol(w)
  hx <- x %*% h[1:k_x]
  mu <- hx
  gz <- z %*% h[(k_x + 1):(k_x + k_z)]
  sigma <- exp(gz)
  gw <- w %*% h[(k_x + k_z + 1):(k_x + k_z + k_w)]
  u <- exp(gw)
  # Evaluate the version of dnsity2 that belongs to this specific model run
  loglik <- log(density_fun(y, mu, sigma, u))
  return(-sum(loglik, na.rm = TRUE))
}

# Here is the optimizer function:
# I use Nelder-Mead first and then BFGS:
optimfunc <- function(start, y, x, z, w, density_fun) {
  # Pass density_fun downstream through the extra arguments (...) of optim()
  precdf <- optim(start, cdfloglik, hessian = TRUE, y = y, x = x, z = z, w = w,
                  density_fun = density_fun, method = "Nelder-Mead")
  cdfmle <- optim(precdf$par, cdfloglik, hessian = TRUE, y = y, x = x, z = z, w = w,
                  density_fun = density_fun, method = "BFGS")
  return(cdfmle)
}

#' Maximum Likelihood Estimation for Extended CDF Quantile Regression
#'
#' Fits a parametric frequentist regression model for bounded outcomes on the unit interval
#' using maximum likelihood estimation across three distribution parameters: mu, sigma, and u.
#'
#' @param formula A three-part formula separated by pipes (`|`),
#'   where the first part specifies predictors
#'   for mu (median), the second for sigma (dispersion), and the third for u (shape).
#' @param data A data frame containing the variables specified in the formula.
#' @param start An optional numeric vector of starting values for the parameters.
#'   If `NULL`, defaults to `0.1` for all parameters.
#'
#' @param x An object of class \code{"exftcdfq"} or \code{"summary.exftcdfq"} (for print methods).
#' @return An object of class (\code{"exftcdfq"}.) This is a list containing:
#' \describe{
#('   \item{coefficients}{A named vector of parameter estimates.})
#('   \item{vcov}{The variance-covariance matrix derived from the inverted Hessian.})
#('   \item{summary_table}{A matrix of estimates, standard errors, z-values, and p-values.})
#('   \item{loglik}{The log-likelihood value at convergence.})
#('   \item{optim_raw}){(The raw output list returned by the \code{optim}) function.}
#('   \item{formula}{The formula used in the model.})
#('   \item{nobs}{The number of observations used in the model.})
#('   \item{y}{The observed response vector.})
#('   \item{data_matrices}{A list containing the design matrices for mu, sigma, and u.})
#('   \item{dnsity2}{The specific density function active at estimation time.})
#('   \item{inverse_CDF}{The specific inverse CDF function active at estimation time.})
#' }
#'
#' @importFrom Formula Formula
#' @importFrom stats optim model.matrix model.response pnorm
#' @importFrom stats optim model.matrix model.response pnorm model.frame predict printCoefmat
#' @export
#'
#' @examples
#' \dontrun{
#' # Select a distribution first
#' select_exftcdfq()
#'
#' # Fit an MLE model
#' m1 <- exftcdfqMLE(ptotpalest1 ~ palestpk | 1 | 1, data = blame)
#' summary(m1)
#' }
exftcdfqMLE <- function(formula, data, start = NULL) {
 # 1. Force the formula into a multi-part Formula object
  f <- Formula(formula)

 # Ensure the user provided exactly 3 parts on the right-hand side (mu, sigma, u)
  if (length(f)[2] != 3) {
    stop("Formula must have exactly 3 right-hand parts separated by '|'. Example: y ~ x | z | w")
  }

  # 2. Extract the response variable (y)
  # drop = TRUE extracts it as a clean vector
  y <- model.response(model.frame(f, data = data, rhs = 1))

  # 3. Dynamically build the design matrices for mu, sigma, and u
  X_mu    <- model.matrix(f, data = data, rhs = 1)
  X_sigma <- model.matrix(f, data = data, rhs = 2)
  X_u     <- model.matrix(f, data = data, rhs = 3)

  # 4. Generate clean parameter names for tracking and final output
  names_mu    <- paste0("mu_", colnames(X_mu))
  names_sigma <- paste0("sigma_", colnames(X_sigma))
  names_u     <- paste0("u_", colnames(X_u))
  all_param_names <- c(names_mu, names_sigma, names_u)

  # 5. Handle starting values automatically if not provided
  total_params <- ncol(X_mu) + ncol(X_sigma) + ncol(X_u)
  if (is.null(start)) {
    # Default to 0.1 for all parameters if user doesn't specify them
    start <- rep(0.1, total_params)
  } else if (length(start) != total_params) {
    stop(paste("Length of start values must equal", total_params))
  }

  # 5b. Capture the EXACT versions of the functions active in the session right now
  if (!exists("dnsity2") || !exists("inverse_CDF")) {
    stop("No distribution is currently active. Please run select_exftcdfq() first.")
  }
  active_density <- get("dnsity2")
  active_inverse <- get("inverse_CDF")

  # 6. Run the optimization pipeline
  # (Passing the matrices directly into your robust optimfunc with the captured density)
  fit <- optimfunc(start = start, y = y, x = X_mu, z = X_sigma, w = X_u, density_fun = active_density)

  # 7. Polish the final output structure
  names(fit$par) <- all_param_names
  rownames(fit$hessian) <- all_param_names
  colnames(fit$hessian) <- all_param_names

  # Calculate standard errors and summary statistics safely
  # (Wrapped in tryCatch in case Hessian is non-invertible)
  summary_table <- tryCatch({
    vcov_mat <- solve(fit$hessian)
    serr     <- sqrt(diag(vcov_mat))
    zstat    <- fit$par / serr
    prob     <- 2 * (1 - pnorm(abs(zstat)))

    cbind(Estimate = fit$par, `Std. Error` = serr, `z value` = zstat, `Pr(>|z|)` = prob)
  }, error = function(e) {
    warning("Hessian could not be inverted to calculate Standard Errors.")
    cbind(Estimate = fit$par)
  })

  # Return a structured list that acts like a model object
  output <- list(
    coefficients  = fit$par,
    vcov          = if(exists("vcov_mat")) vcov_mat else NULL,
    summary_table = summary_table,
    loglik        = -fit$value, # Convert back to positive log-likelihood
    optim_raw     = fit,
    formula       = formula
  )
  # Return a structured list that acts like a model object
  output <- list(
    coefficients  = fit$par,
    vcov          = if(exists("vcov_mat")) vcov_mat else NULL,
    summary_table = summary_table,
    loglik        = -fit$value,
    optim_raw     = fit,
    formula       = formula,
    nobs          = length(y),
    y             = as.numeric(y),
    data_matrices = list(mu = X_mu, sigma = X_sigma, u = X_u),
    # Save the functions themselves inside the model object!
    dnsity2       = active_density,
    inverse_CDF   = active_inverse
  )
  # Assign the custom S3 class
  class(output) <- "exftcdfq"
  return(output)
}

#' S3 Methods

# Print function
#' @rdname exftcdfqMLE
#' @export
print.exftcdfq <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$formula), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Log-Likelihood:", round(x$loglik, 4), "\n\n")
  cat("Coefficients:\n")
  print(x$coefficients)
  cat("\n")
  invisible(x)
}

#' Summary function
#' Summarizing Extended CDF Quantile Regression Fits
#'
#' @param object An object of class \code{"exftcdfq"}.
#' @param ... Additional arguments (not used).
#'
#' @return An object of class \code{"summary.exftcdfq"}.
#' @rdname exftcdfqMLE
#' @export
summary.exftcdfq <- function(object, ...) {
  # Calculate information criteria
  k <- length(object$coefficients)
  aic <- 2 * k - 2 * object$loglik
  bic <- k * log(object$nobs) - 2 * object$loglik

  ans <- list(
    formula       = object$formula,
    summary_table = object$summary_table,
    loglik        = object$loglik,
    aic           = aic,
    bic           = bic
  )
  class(ans) <- "summary.exftcdfq"
  return(ans)
}

#' The companion print method for the summary object
#' @rdname exftcdfqMLE
#' @export
print.summary.exftcdfq <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$formula), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("Fit Statistics:\n")
  cat("  Log-Likelihood:", round(x$loglik, 4), "\n")
  cat("  AIC:           ", round(x$aic, 4), "\n")
  cat("  BIC:           ", round(x$bic, 4), "\n\n")

  cat("Coefficients:\n")
  printCoefmat(x$summary_table, P.values = TRUE, has.Pvalue = TRUE)
  cat("\n")
  invisible(x)
}

#' coef, vcov extraction
#' @importFrom stats coef
#' @rdname exftcdfqMLE
#' @export
coef.exftcdfq <- function(object, ...) {
  return(object$coefficients)
}
#' @importFrom stats vcov
#' @export
vcov.exftcdfq <- function(object, ...) {
  return(object$vcov)
}

#' Extract Parameter Correlation Matrix for exftcdfq Fits
#'
#' Calculates the correlation matrix of the estimated parameters from a fitted model.
#'
#' @param object An object of class \code{"exftcdfq"}.
#' @param ... Additional arguments (not used).
#'
#' @return A numeric correlation matrix.
#' @importFrom stats cov2cor
#' @export
parameter_cor <- function(object, ...) {
  if (is.null(object$vcov)) {
    stop("Variance-covariance matrix is missing; cannot compute correlations.")
  }
  return(stats::cov2cor(object$vcov))
}

#' Fit criteria extraction helper function
#' Extract Log-Likelihood from an exftcdfq Model
#'
#' @param object An object of class \code{"exftcdfq"}.
#' @param ... Additional arguments (not used).
#'
#' @importFrom stats logLik
#' @rdname exftcdfqMLE
#' @export
logLik.exftcdfq <- function(object, ...) {
  # Calculate degrees of freedom (number of estimated parameters)
  df_count <- length(object$coefficients)

     # Structure it as an official logLik object class
  val <- object$loglik
  attr(val, "df") <- df_count
  attr(val, "nobs") <- object$nobs
  class(val) <- "logLik"

  return(val)
}

#' Predicting quantiles
#' Obtains predicted values (quantiles) from a fitted extended CDF quantile regression model object.
#'
#' Predict Method for exftcdfq Model Fits
#'
#' Obtains predicted values (quantiles) from a fitted extended CDF quantile regression model object.
#'
#' @param object An object of class \(\code{"exftcdfq"}.\)
#' @param newdata An optional data frame or vector of new values for out-of-sample prediction.
#'   If omitted, fitted values for the original dataset are returned.
#' @param q A numeric value or vector indicating the target quantile(s). Defaults to `0.5` (median).
#' @param ... Additional arguments (not used).
#'
#' @return A numeric vector of predicted values.
#' @importFrom Formula Formula
#' @importFrom stats model.matrix formula
#' @export
predict.exftcdfq <- function(object, newdata = NULL, q = 0.5, ...) {
  # 1. Parse the multi-part formula structure cleanly from the model object
  f <- Formula::Formula(object$formula)

  # 2. Extract or reconstruct the design matrices safely
  if (is.null(newdata)) {
    X_mu    <- object$data_matrices$mu
    X_sigma <- object$data_matrices$sigma
    X_u     <- object$data_matrices$u
  } else {
   # Defensive check: if the user passed a raw vector, convert it to a data.frame
    if (!is.data.frame(newdata)) {
      # Use base R to get all variables on the right-hand side of the formula
      vars_needed <- all.vars(formula(f, lhs = 0, rhs = 1))

      if (length(vars_needed) == 1) {
        newdata <- data.frame(newdata)
        names(newdata) <- vars_needed
      } else {
        stop("`newdata` must be a data.frame with column names matching your formula predictors.")
      }
    }

    # Safely build matrices from the structured data.frame
    X_mu    <- model.matrix(f, data = newdata, rhs = 1)
    X_sigma <- model.matrix(f, data = newdata, rhs = 2)
    X_u     <- model.matrix(f, data = newdata, rhs = 3)
  }

  # 3. Extract counts using NCOL
  k_mu    <- NCOL(X_mu)
  k_sigma <- NCOL(X_sigma)
  k_u     <- NCOL(X_u)

  h <- object$coefficients

  # 4. Calculate linear predictors and apply parameter link functions
  mu    <- X_mu %*% h[1:k_mu]
  sigma <- exp(X_sigma %*% h[(k_mu + 1):(k_mu + k_sigma)])
  u     <- exp(X_u %*% h[(k_mu + k_sigma + 1):(k_mu + k_sigma + k_u)])

  # 5. Pass parameters to your inverse_CDF function by position
  predictions <- object$inverse_CDF(q, as.vector(mu), as.vector(sigma), as.vector(u))

  return(predictions)
}

#' Extract Model Residuals
#'
#' Calculates the raw residuals (Observed minus Predicted Medians) for an exftcdfq object.
#'
#' @param object An object of class \(\code{"exftcdfq"}.\)
#' @param ... Additional arguments (not used).
#'
#' @return A numeric vector of residuals.
#' @importFrom stats residuals
#' @export
residuals.exftcdfq <- function(object, ...) {
   # Standard residuals are: Observed minus Predicted Median
   # Coerce observed values to a clean, flat numeric vector
  observed <- as.vector(object$y)
  predicted <- predict(object, q = 0.5)
  res <- observed - predicted
  return(res)
}
