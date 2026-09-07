##### Cauchit_Arcsinh_inner_V_3 #####
#  real Cauchit_Arcsinh_inner_V_rng(real x, real theta, real sigma) {
#    return Vtransinv(Arcsinhfunc(Utransinv(InvCauchy(uniform_rng(0, 1)),sigma)),theta); }
# Helper functions
Cauchyf <- function(x){
  return ((atan(x)/pi) + 0.5) }

InvCauchyf <- function(x){
  return (tan(pi*x - pi/2)) }

signum <- function(x){
  if(x < 0){
    return (-1)}
  if(x > 0){
    return (1)}
  else{ 
    return (0)}
}

Arcsinhf <- function(x){
  return (1/(exp(-asinh(x))+1)) }

InvArcsinhf <- function(x){
  return ((1-2*x)/(2*x*(x-1)))}

Vtransf <- function(x, theta){
  return ((exp(theta)*x)/(1-x+exp(theta)*x))}

Vinv <- function(x, theta){
  return (x/(exp(theta) + x - x*exp(theta)))}

Utransf <- function(x, mu, sigma){
  return ((x-mu)/sigma)}

Uinv <- function(x,mu, sigma){
  return (mu+x*sigma)}

# CDF function
CDF <- function(x, mu, sigma, theta){
  if(x==0){
    return ((2*sigma*exp(theta))/pi)}
  if(x==1){
    return ((2*sigma*exp(-theta))/pi)}
  else{
    return ((Cauchyf(Utransf(InvArcsinhf(Vtransf(x,theta)),mu,sigma))))}
}

# density function
dnsity <- function(x, mu, sigma, theta){
  if(x==0){
    return ((2*sigma*exp(theta))/pi)}
  if(x==1){
    return ((2*sigma*exp(-theta))/pi)}
  else{
    return ((Cauchyf(Utransf(InvArcsinhf(Vtransf(x + 0.000001,theta)),mu,sigma)) - Cauchyf(Utransf(InvArcsinhf(Vtransf(x,theta)),mu,sigma)))/0.000001)}
}

inverse_CDF <- function(x, mu, sigma, theta){
    return (Vinv(Arcsinhf(Uinv(InvCauchyf(x),mu,sigma)),theta))}
  
  # random draws function
  rsamp <- function(n, mu, sigma, theta){
    return (inverse_CDF(runif(n),mu, sigma, theta)) }
  
  # log_lik function
  log_lik_Cauchit_Arcsinh_inner_V_3 <- function(i, prep) {
    mu <- brms::get_dpar(prep, "mu", i = i)
    sigma <- brms::get_dpar(prep, "sigma", i = i)
    theta <- brms::get_dpar(prep, "theta", i = i)
    y <- prep$data$Y[i]
    log(dnsity(y, mu, sigma, theta))
  }
  
  # posterior_predict function
  posterior_predict_Cauchit_Arcsinh_inner_V_3 <- function(i, prep, ...) {
    mu <- brms::get_dpar(prep, "mu", i = i)
    sigma <- brms::get_dpar(prep, "sigma", i = i)
    theta <- brms::get_dpar(prep, "theta", i = i)
    ndraws <- prep$ndraws
    rsamp(ndraws, mu, sigma, theta)
  }
  
  # posterior_epred function
  posterior_epred_Cauchit_Arcsinh_inner_V_3 <- function(prep) {
    mu <- brms::get_dpar(prep, "mu")
    sigma <- brms::get_dpar(prep, "sigma")
    theta <- brms::get_dpar(prep, "theta")
    inverse_CDF(0.5,mu, sigma, theta)
  }
  
  # This defines the custom family.
  Cauchit_Arcsinh_inner_V_3 <- custom_family(
    "Cauchit_Arcsinh_inner_V_3", dpars = c("mu", "sigma", "theta"),
    links = c("identity", "log", "identity"),
    lb = c(NA, 0, NA), ub = c(NA, NA, NA),
    type = "real"
  )
  # This is the piecewise log-likelihood. 
  stan_funs <- '
  real Cauchy(real x){
    return (atan(x)/pi()) + 0.5;}

  real InvCauchy(real x){
    return tan(pi()*x - pi()/2); }
    
  real signum(real x) {
    if(x < 0){
      return -1;}
    if(x > 0){
      return 1;}
    else{ 
      return 0;}
    }
  
  real Arcsinhfunc(real x){
    return 1/(exp(-asinh(x))+1); }

  real InvArcsinh(real x){
    return (1-2*x)/(2*x*(x-1));}
    
  real VTransform(real x, real theta){
    return (exp(theta)*x)/(1-x+exp(theta)*x);}

  real Vtransinv(real x, real theta){
    return x/(exp(theta) + x - x*exp(theta));}
    
  real U(real x, real mu, real sigma){
    return (x-mu)/sigma;}

  real Utransinv(real x, real mu, real sigma){
    return mu+x*sigma;}
  
  real Cauchit_Arcsinh_inner_V_3_lpdf(real x, real mu, real sigma, real theta){
    if(x==0){
      return log((2*sigma*exp(theta))/pi());
    }
    if(x==1){
      return log((2*sigma*exp(-theta))/pi());
    }
    else{
       return log(Cauchy(U(InvArcsinh(VTransform(x+0.00001,theta)),mu,sigma)) - Cauchy(U(InvArcsinh(VTransform(x,theta)),mu,sigma)));}
}
'
# We make an object that tells brms that the object stan_funs contains 
# STAN-code from the functions block of the STAN program. 
my_family <- Cauchit_Arcsinh_inner_V_3
stanvars <- stanvar(scode = stan_funs, block = "functions")