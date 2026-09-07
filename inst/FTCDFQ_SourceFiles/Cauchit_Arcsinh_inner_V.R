##### Cauchit_Arcsinh_inner_V #####
#  real Cauchit_Arcsinh_inner_V_rng(real x, real mu, real sigma) {
#    return Vtransinv(Arcsinhfunc(Utransinv(InvCauchy(uniform_rng(0, 1)),sigma)),mu); }
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

Vtransf <- function(x, mu){
  return ((exp(mu)*x)/(1-x+exp(mu)*x))}

Vinv <- function(x, mu){
  return (x/(exp(mu) + x - x*exp(mu)))}

Utransf <- function(x, sigma){
  return (x/sigma)}

Uinv <- function(x, sigma){
  return (x*sigma)}

# CDF function
CDF <- function(x, mu, sigma){
  if(x==0){
    return ((2*sigma*exp(mu))/pi)}
  if(x==1){
    return ((2*sigma*exp(-mu))/pi)}
  else{
    return ((Cauchyf(Utransf(InvArcsinhf(Vtransf(x,mu)),sigma))))}
}

# density function
dnsity <- function(x, mu, sigma){
  if(x==0){
    return ((2*sigma*exp(mu))/pi)}
  if(x==1){
    return ((2*sigma*exp(-mu))/pi)}
  else{
    return ((Cauchyf(Utransf(InvArcsinhf(Vtransf(x + 0.000001,mu)),sigma)) - Cauchyf(Utransf(InvArcsinhf(Vtransf(x,mu)),sigma)))/0.000001)}
}

inverse_CDF <- function(x, mu, sigma){
    return (Vinv(Arcsinhf(Uinv(InvCauchyf(x),sigma)),mu))}
  
  # random draws function
  rsamp <- function(n, mu, sigma){
    return (inverse_CDF(runif(n),mu,sigma)) }
  
  # log_lik function
  log_lik_Cauchit_Arcsinh_inner_V <- function(i, prep) {
    mu <- brms::get_dpar(prep, "mu", i = i)
    sigma <- brms::get_dpar(prep, "sigma", i = i)
    y <- prep$data$Y[i]
    log(dnsity(y, mu, sigma))
  }
  
  # posterior_predict function
  posterior_predict_Cauchit_Arcsinh_inner_V <- function(i, prep, ...) {
    mu <- brms::get_dpar(prep, "mu", i = i)
    sigma <- brms::get_dpar(prep, "sigma", i = i)
    ndraws <- prep$ndraws
    rsamp(ndraws, mu, sigma)
  }
  
  # posterior_epred function
  posterior_epred_Cauchit_Arcsinh_inner_V <- function(prep) {
    mu <- brms::get_dpar(prep, "mu")
    sigma <- brms::get_dpar(prep, "sigma")
    inverse_CDF(0.5,mu,sigma)
  }
  
  # This defines the custom family.
  Cauchit_Arcsinh_inner_V <- custom_family(
    "Cauchit_Arcsinh_inner_V", dpars = c("mu", "sigma"),
    links = c("identity", "log"),
    lb = c(NA, 0), ub = c(NA, NA),
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
    
  real VTransform(real x, real mu){
    return (exp(mu)*x)/(1-x+exp(mu)*x);}

  real Vtransinv(real x, real mu){
    return x/(exp(mu) + x - x*exp(mu));}
    
  real U(real x, real sigma){
    return x/sigma;}

  real Utransinv(real x, real sigma){
    return x*sigma;}
  
  real Cauchit_Arcsinh_inner_V_lpdf(real x, real mu, real sigma){
    if(x==0){
      return log((2*sigma*exp(mu))/pi());
    }
    if(x==1){
      return log((2*sigma*exp(-mu))/pi());
    }
    else{
       return log(Cauchy(U(InvArcsinh(VTransform(x+0.00001,mu)),sigma)) - Cauchy(U(InvArcsinh(VTransform(x,mu)),sigma)));}
}
'
# We make an object that tells brms that the object stan_funs contains 
# STAN-code from the functions block of the STAN program. 
my_family <- Cauchit_Arcsinh_inner_V
stanvars <- stanvar(scode = stan_funs, block = "functions")