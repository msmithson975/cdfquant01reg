##### Arcsinh_Cauchy_inner_V #####
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

Wtransf <- function(x, mu){
  return (sinh(asinh(x) + mu))}

Winv <- function(x, mu){
  return (sinh(asinh(x) - mu))}

Utransf <- function(x, sigma){
  return (x/sigma)}

Uinv <- function(x, sigma){
  return (x*sigma)}

# CDF function
CDF <- function(x, mu, sigma){
  if(x==0){
    return ((pi*sigma*exp(mu))/2)}
  if(x==1){
    return ((pi*sigma*exp(-mu))/2)}
  else{
    return ((Arcsinhf(Utransf(InvCauchyf(Vtransf(x,mu)),sigma))))}
}

# density function
dnsity <- function(x, mu, sigma){
  if(x==0){
    return ((pi*sigma*exp(mu))/2)}
  if(x==1){
    return ((pi*sigma*exp(-mu))/2)}
  else{
    return ((Arcsinhf(Utransf(InvCauchyf(Vtransf(x + 0.000001,mu)),sigma)) - Arcsinhf(Utransf(InvCauchyf(Vtransf(x,mu)),sigma)))/0.000001)}
}

inverse_CDF <- function(x, mu, sigma){
  return(Vinv(Cauchyf(Uinv(InvArcsinhf(x),sigma)),mu))
}

# random draws function
rsamp <- function(n, mu, sigma){
  return (inverse_CDF(runif(n, min = 0, max = 1),mu,sigma)) }

# log_lik function
log_lik_Arcsinh_Cauchy_inner_V <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  y <- prep$data$Y[i]
  log(dnsity(y, mu, sigma))
}

# posterior_predict function
posterior_predict_Arcsinh_Cauchy_inner_V <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  ndraws <- prep$ndraws
  rsamp(ndraws, mu, sigma)
}

# posterior_epred function
posterior_epred_Arcsinh_Cauchy_inner_V <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  sigma <- brms::get_dpar(prep, "sigma")
  inverse_CDF(0.5,mu,sigma)
}

# This defines the custom family.
Arcsinh_Cauchy_inner_V <- custom_family(
  "Arcsinh_Cauchy_inner_V", dpars = c("mu", "sigma"),
  links = c("identity", "log"),
  lb = c(NA, 0), ub = c(NA, NA),
  type = "real"
)
# This is the piecewise log-likelihood. 
# STAN lacks native cotangent or cosecant functions so these must be defined. 
stan_funs <- '
  real Arcsinh(real x){
    return 1 / (exp(-asinh(x))+1);}
    
  real signum(real x) {
    if(x < 0){
      return -1;}
    if(x > 0){
      return 1;}
    else{ 
      return 0;}
    }
  
  real InvCauchy(real x){
    return tan(pi()*x - pi()/2);}
    
  real VTransform(real x, real mu){
    return (exp(mu)*x)/(1-x+exp(mu)*x);}
    
  real U(real x, real sigma){
    return x/sigma;}
  
  real Arcsinh_Cauchy_inner_V_lpdf(real x, real mu, real sigma){
    if(x==0){
      return log((pi()*sigma*exp(mu))/2);
    }
    if(x==1){
      return log((pi()*sigma*exp(-mu))/2);
    }
    else{
       return log(Arcsinh(U(InvCauchy(VTransform(x+0.00001,mu)),sigma)) - Arcsinh(U(InvCauchy(VTransform(x,mu)),sigma)));}
}'
# We make an object that tells brms that the object stan_funs contains 
# STAN-code from the functions block of the STAN program. 
my_family <- Arcsinh_Cauchy_inner_V
stanvars <- stanvar(scode = stan_funs, block = "functions")
#
