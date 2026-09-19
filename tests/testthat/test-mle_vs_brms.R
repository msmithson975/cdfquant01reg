# Testing unit
library(testthat)
library(Formula)
library(brms)
test_that("MLE and brms estimates are statistically consistent (Z-score < 2)", {
  # 1. Generate a tiny synthetic dataset to keep the test exceptionally fast
  set.seed(123)
  n_test <- 34
  test_data <- data.frame(
    x_var = rbinom(n_test, 1, 0.5),
    # Generate bounded values strictly on the closed unit interval [0,1]
    y_var = c(rep(0,2),rep(1,2),runif(n_test-4, 0.01, 0.99))
  )

  # 2. Loop through all 16 distributions to test the entire family
  # Replace these placeholder strings with your actual package distribution names
  distributions_to_test <- c(
       "Ex_Arcsinh_Arcsinh_outer_W",
       "Ex_Arcsinh_Cauchy_outer_W",
       "Ex_Cauchit_Arcsinh_outer_W",
       "Ex_Cauchit_Cauchy_outer_W",
       "Ex_t2_t2_outer_W",
       "Ex_Cauchit_Arcsinh_outer_V",
       "Ex_Cauchit_Cauchy_outer_V",
       "Ex_t2_t2_outer_V",
       "Ex_Arcsinh_Arcsinh_inner_W",
       "Ex_Arcsinh_Cauchy_inner_W",
       "Ex_Cauchit_Cauchy_inner_W",
       "Ex_t2_t2_inner_W",
       "Ex_Arcsinh_Cauchy_inner_V",
       "Ex_Cauchit_Arcsinh_inner_V",
       "Ex_Cauchit_Cauchy_inner_V",
       "Ex_t2_t2_inner_V"
  )

  for (dist_name in distributions_to_test) {

    # 3. Use your package's selector mechanism to load the active distribution functions
    # (Assuming select_exftcdfq or an internal variant can accept a string programmatically)
    # If select_exftcdfq doesn't take strings, we explicitly source the file:
    rel_path <- paste0("FTCDFQ_SourceFiles/", dist_name, ".R")

    if (requireNamespace("pkgload", quietly = TRUE) && pkgload::is_dev_package("cdfquant01reg")) {
      chosen_file <- pkgload::package_file("inst", rel_path)
    } else {
      chosen_file <- system.file(rel_path, package = "cdfquant01reg")
    }
    source(chosen_file, local = .GlobalEnv)

    # 4. Fit the Frequentist MLE model
    fit_mle <- exftcdfqMLE(y_var ~ x_var | 1 | 1, data = test_data)

    # 5. Fit the Bayesian brms model
    # Using 1 chain and 500 iterations keeps Stan blazingly fast for a unit test
    fit_brms <- exftcdfq(
      brms::bf(y_var ~ x_var, sigma ~ 1, u ~ 1),
      data = test_data,
      backend = "cmdstanr",
      chains = 1,
      iter = 500,
      silent = 2
    )

    # 6. Extract estimates and standard errors
    mle_estimates  <- coef(fit_mle)
    mle_se         <- fit_mle$summary_table[, "Std. Error"]
    brms_estimates <- brms::fixef(fit_brms)[, "Estimate"]

    # 7. Map the brms coefficients to match MLE naming conventions
    brms_names <- names(brms_estimates)
    brms_names[brms_names == "Intercept"]       <- "mu_(Intercept)"
    brms_names[brms_names == "sigma_Intercept"] <- "sigma_(Intercept)"
    brms_names[brms_names == "u_Intercept"]     <- "u_(Intercept)"

    is_bare_predictor <- !grepl("^mu_|^sigma_|^u_", brms_names)
    brms_names[is_bare_predictor] <- paste0("mu_", brms_names[is_bare_predictor])
    names(brms_estimates) <- brms_names

    # Re-order the brms vector to exactly match the MLE vector layout
    brms_estimates_aligned <- brms_estimates[names(mle_estimates)]

    # 8. Calculate the absolute Z-score difference
    z_diff <- abs(mle_estimates - brms_estimates_aligned) / mle_se

    # 9. Execute the self-scaling test criteria
    # If any parameter diverges by more than 2 SEs, testthat halts and tells you which distribution failed.
    expect_true(
      all(z_diff < 2),
      info = paste("Paradigm mismatch in distribution:", dist_name, "- Max SE distance:", max(z_diff))
    )
  }
})
